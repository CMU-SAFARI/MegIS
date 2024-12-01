#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <omp.h>
#include <cstdint>
#include <tuple>
#include <cstring>
#include <memory>
#include <mutex>

#include "progress_bar.hpp"

using uint64 = uint64_t;
using uint32 = uint32_t;
using int32 = int32_t;
using uchar = unsigned char;
static std::mutex mu;

#define MAX_K 256
#define SIZE ((MAX_K + 31) / 32)
#define _bswap_uint64(X) __builtin_bswap64(X)

class KMCDB {
  public:
    KMCDB(const std::string &prefix)
          : pre_file_name(prefix + ".kmc_pre"), suf_file_name(prefix + ".kmc_suf") {
        FILE *file_pre = fopen(pre_file_name.c_str(), "rb");

        fseek(file_pre, -12, SEEK_END);

        uint32 kmc_version;
        fread(&kmc_version, sizeof(uint32), 1, file_pre);
        if (kmc_version != 0) {
            std::cerr << "Error: Only KMC1 supported";
            exit(1);
        }

        fseek(file_pre, -8, SEEK_END);

        header_offset = fgetc(file_pre);

        fseek(file_pre, (0LL - (header_offset + 8)), SEEK_END);
        fread(&kmer_length, 1, sizeof(uint32), file_pre);
        fread(&mode, 1, sizeof(uint32), file_pre);
        if (mode != 0) {
            std::cerr << "Error: Quake quake compatible counters are not supported anymore\n";
            exit(1);
        }

        fread(&counter_size, 1, sizeof(uint32), file_pre);
        fread(&lut_prefix_len, 1, sizeof(uint32), file_pre);

        fread(&min_count, 1, sizeof(uint32), file_pre);

        fread(&max_count_lo, 1, sizeof(uint32), file_pre);

        fread(&total_kmers, 1, sizeof(uint64), file_pre);

        fread(&both_strands, 1, 1, file_pre);
        both_strands = !both_strands;

        fread(&max_count_hi, 1, sizeof(uint32), file_pre);

        fclose(file_pre);

        suffix_bytes = (kmer_length - lut_prefix_len) / 4;
        record_size = suffix_bytes + counter_size;

        std::cout << "Opened " << prefix << " containing "
                  << total_kmers << " " << kmer_length << "-mers\tprefix length " << lut_prefix_len
                  << "\trecord size " << record_size << "\tcounter size " << counter_size << "\n";
    }

    using Record = std::unique_ptr<uchar[]>;

    class Generator {
      public:
        Generator(const std::string &pre_file_name,
                  const std::string &suf_file_name,
                  uint64 min_prefix,
                  uint64 max_prefix,
                  uint64 total_kmers,
                  uint64 record_size,
                  uint32 lut_prefix_len)
              : file_pre(pre_file_name, std::ios::binary),
                file_suf(suf_file_name, std::ios::binary),
                min_prefix(min_prefix),
                max_prefix(max_prefix),
                record_size(record_size),
                //iobuf_pre(std::make_unique<char[]>(1000000)),
                //iobuf_suf(std::make_unique<char[]>(1000000)),
                cur_record(std::make_unique<uchar[]>(record_size)) {
            //file_pre.rdbuf()->pubsetbuf(iobuf_pre.get(), sizeof iobuf_pre.get());
            //file_suf.rdbuf()->pubsetbuf(iobuf_suf.get(), sizeof iobuf_suf.get());
            if (!file_pre.good()) {
                std::cerr << "Error: failed to open pre: " << pre_file_name << "\n";
                exit(1);
            }
            if (!file_suf.good()) {
                std::cerr << "Error: failed to open suf: " << suf_file_name << "\n";
                exit(1);
            }

            char marker[4];
            file_suf.read(marker, 4);

            if (strncmp(marker, "KMCS", 4) != 0) {
                std::cerr << "Error: wrong start marker in file: " << suf_file_name << "\n";
                exit(1);
            }

            file_suf.seekg(-4, std::ios_base::end);
            file_suf.read(marker, 4);
            if (strncmp(marker, "KMCS", 4) != 0) {
                std::cerr << "Error: wrong end marker in file: " << suf_file_name << "\n";
                exit(1);
            }

            file_suf.seekg(4, std::ios_base::beg);

            file_pre.read(marker, 4);
            if (strncmp(marker, "KMCP", 4) != 0) {
                std::cerr << "Error: wrong start marker in file: " << pre_file_name << "\n";
                exit(1);
            }

            file_pre.seekg(-4, std::ios_base::end);
            file_pre.read(marker, 4);
            if (strncmp(marker, "KMCP", 4) != 0) {
                std::cerr << "Error: wrong end marker in file: " << pre_file_name << "\n";
                exit(1);
            }

            file_pre.seekg(min_prefix * sizeof(uint64) + 4, std::ios_base::beg);
            file_pre.read((char*)&num_to_skip_front, sizeof(num_to_skip_front));

            prefix_left_to_read = 1llu << (lut_prefix_len * 2);
            auto total_kmers_global = total_kmers;
            if (max_prefix < prefix_left_to_read) {
                file_pre.seekg(max_prefix * sizeof(uint64) + 4, std::ios_base::beg);
                uint64 total_kmers_check;
                file_pre.read((char*)&total_kmers_check, sizeof(total_kmers_check));
                if (total_kmers_check > total_kmers) {
                    std::cerr << "Error: failed to read total k-mers in middle\n";
                    exit(1);
                }
                total_kmers = total_kmers_check;
            }

            if (num_to_skip_front > total_kmers) {
                std::cerr << "Error: failed to read num to skip\n";
                exit(1);
            }

            total_kmers_left = total_kmers - num_to_skip_front;

            prefix_left_to_read = max_prefix - min_prefix;

            file_pre.seekg(min_prefix * sizeof(uint64) + 4, std::ios_base::beg);
            file_suf.seekg(num_to_skip_front * record_size + 4, std::ios_base::beg);

            file_pre.read((char*)&prefix_boundary, sizeof(prefix_boundary));
            if (prefix_boundary != num_to_skip_front) {
                std::cerr << "Error: failed to read number of k-mers to skip in front\n";
                exit(1);
            }

            file_pre.read((char*)&prefix_boundary, sizeof(prefix_boundary));
            if (num_to_skip_front > prefix_boundary) {
                std::cerr << "Error: invalid read of range from prefix\n";
                exit(1);
            }

            kmers_left_in_prefix = prefix_boundary - num_to_skip_front;
            ++min_prefix;

            if (kmers_left_in_prefix) {
                file_suf.read((char*)cur_record.get(), record_size);
            }

            {
                std::lock_guard<std::mutex> lock(mu);
                std::cout << "Reading " << total_kmers_left << " / " << total_kmers_global << " k-mers from prefix range ["
                          << min_prefix << ", " << max_prefix << ") in file "
                          << suf_file_name << "\n";
            }
        }

        std::ifstream file_pre;
        std::ifstream file_suf;
        uint64 total_kmers_left;
        uint64 min_prefix;
        uint64 max_prefix;
        uint64 prefix_left_to_read;
        uint64 prefix_boundary;
        uint64 kmers_left_in_prefix;
        uint64 record_size;
        uint64 counter_size;
        uint64 num_to_skip_front;
        //std::unique_ptr<char[]> iobuf_pre;
        //std::unique_ptr<char[]> iobuf_suf;

        Record cur_record;

        uint64 top_prefix() const { return min_prefix - 1; }

        void pop_prefix() {
            if (min_prefix == max_prefix) {
                std::cerr << "Error: read too much from prefix\n";
                exit(1);
            }

            ++min_prefix;
            if (kmers_left_in_prefix > total_kmers_left) {
                std::cerr << "Error: overreading\n";
                exit(1);
            }
            total_kmers_left -= kmers_left_in_prefix;
            if (total_kmers_left) {
#ifndef NDEBUG
                auto old_g = file_suf.tellg();
#endif
                file_suf.seekg(prefix_boundary * record_size + 4, std::ios_base::beg);
#ifndef NDEBUG
                if (kmers_left_in_prefix && file_suf.tellg() - old_g != static_cast<long long>((kmers_left_in_prefix - 1) * record_size)) {
                    std::cerr << "Error: bad jump\n";
                    exit(1);
                }
#endif
                if (min_prefix + 1 < max_prefix) {
                    uint64 next_prefix_boundary;
                    file_pre.read((char*)&next_prefix_boundary, sizeof(next_prefix_boundary));
#ifndef NDEBUG
                    if (prefix_boundary > next_prefix_boundary) {
                        std::cerr << "Error: failed to read next boundary after reading " << min_prefix << " / " << max_prefix << "\n";
                        exit(1);
                    }
#endif
                    kmers_left_in_prefix = next_prefix_boundary - prefix_boundary;
                    prefix_boundary = next_prefix_boundary;
#ifndef NDEBUG
                    if (kmers_left_in_prefix > total_kmers_left) {
                        std::cerr << "Error: about to read " << kmers_left_in_prefix << " when there are only " << total_kmers_left << " k-mers left\n";
                        exit(1);
                    }
#endif
                } else {
                    kmers_left_in_prefix = total_kmers_left;
                }
                if (kmers_left_in_prefix) {
                    file_suf.read((char*)cur_record.get(), record_size);
                }
            }
        }
        bool has_prefix() const {
            return min_prefix < max_prefix;
        }

        const Record& top_suffix() const { return cur_record; }
        void pop_suffix() {
#ifndef NDEBUG
            if (!kmers_left_in_prefix) {
                std::cerr << "Error: no kmers left in prefix\n";
                exit(1);
            }
            if (!total_kmers_left) {
                std::cerr << "Error: no kmers left, but there are still " << kmers_left_in_prefix << " kmers in the prefix\n";
                exit(1);
            }
#endif
            --total_kmers_left;
            if (--kmers_left_in_prefix) {
                file_suf.read((char*)cur_record.get(), record_size);
            }
        }
        bool has_suffix() const { return kmers_left_in_prefix; }

        uint64 kmers_left() const { return total_kmers_left; }
    };

    Generator get_generator(uint64 min_prefix, uint64 max_prefix) const {
        return Generator(pre_file_name, suf_file_name, min_prefix, max_prefix,
                         total_kmers, record_size, lut_prefix_len);
    }

    std::string pre_file_name;
    std::string suf_file_name;

    uint32 kmer_length;
    uint32 lut_prefix_len;
    uint64 total_kmers;
    uint64 suffix_bytes;
    uint64 record_size;
    uchar header_offset;
    uint32 mode;
    uint32 counter_size;
    uint32 min_count;
    uint32 max_count_lo;
    bool both_strands;
    uint32 max_count_hi;
};


int8_t suffix_cmp(const KMCDB::Record &a, const KMCDB::Record &b, uint32 n, bool little_endian) {
    /*
    auto a_buf = std::make_unique<unsigned long long[]>(SIZE);
    auto b_buf = std::make_unique<unsigned long long[]>(SIZE);

    auto copy_over = [&](unsigned long long *data, uchar *buffer) {
        uint32 p = ((n + 7) >> 3) - 1;
        for (uint32 i = SIZE - 1; i > p; --i)
            data[i] = 0;

        if (!(n & 7))
            ++p;
        else
        {
            memcpy(&data[p], buffer, sizeof(uint64));
            if (little_endian)
                data[p] = _bswap_uint64(data[p]);
            data[p] >>= (sizeof(uint64)-(n & 7)) << 3;
            buffer += n & 7;
        }
        for (int i = p - 1; i >= 0; --i)
        {
            memcpy(&data[i], buffer, sizeof(uint64));
            if (little_endian)
                data[i] = _bswap_uint64(data[i]);
            buffer += 8;
        }
    };

    copy_over(a_buf.get(), a.get());
    copy_over(b_buf.get(), b.get());

    for (int32 i = SIZE - 1; i >= 0; --i) {
        if (a_buf.get()[i] < b_buf.get()[i])
            return -1;

        if (a_buf.get()[i] > b_buf.get()[i])
            return 1;
    }

    return 0;
    */
    auto begin_a = a.get();
    auto end_a = a.get() + n;
    auto begin_b = b.get();
    auto end_b = b.get() + n;

    auto [a_it, b_it] = std::mismatch(begin_a, end_a, begin_b, end_b);
    if (a_it == end_a)
        return 0;

    if (*a_it > *b_it)
        return 1;

    return -1;
}


uint64 intersect_DBs(KMCDB &db1, KMCDB &db2, uint64 min_prefix, uint64 max_prefix, uint64 end,
                   const std::string &out_name_prefix, bool little_endian) {
    std::ignore = end;
    std::ofstream pre_fout(out_name_prefix + ".kmc_pre", std::ios::binary);
    std::ofstream suff_fout(out_name_prefix + ".kmc_suf", std::ios::binary);

    // std::string pre_marker = "KMCP";
    // std::string suf_marker = "KMCS";
    // pre_fout.write(pre_marker.c_str(), 4);
    // suff_fout.write(suf_marker.c_str(), 4);

    auto generator1 = db1.get_generator(min_prefix, max_prefix);
    auto generator2 = db2.get_generator(min_prefix, max_prefix);

    uint64 num_kmers_before = 0;

    //ProgressBar progress_bar(max_prefix - min_prefix, "Processing",
    //                         std::cerr, false);
    while (generator1.has_prefix() && generator2.has_prefix()) {
        pre_fout.write((char*)&num_kmers_before, sizeof(num_kmers_before));
        while (generator1.has_suffix() && generator2.has_suffix()) {
            int8_t cmp = suffix_cmp(generator1.top_suffix(), generator2.top_suffix(), db1.suffix_bytes, little_endian);
            if (cmp == 0) {
                suff_fout.write((const char*)generator1.top_suffix().get(), db1.record_size);
                generator1.pop_suffix();
                generator2.pop_suffix();
                ++num_kmers_before;
            } else if (cmp == -1) {
                generator1.pop_suffix();
            } else {
                generator2.pop_suffix();
            }
        }

        generator1.pop_prefix();
        generator2.pop_prefix();
        //++progress_bar;
    }

    // pre_fout.write((char*)&db1.kmer_length, sizeof(db1.kmer_length));
    // pre_fout.write((char*)&db1.mode, sizeof(db1.mode));
    // pre_fout.write((char*)&db1.counter_size, sizeof(db1.counter_size));
    // pre_fout.write((char*)&db1.lut_prefix_len, sizeof(db1.lut_prefix_len));
    // pre_fout.write((char*)&db1.min_count, sizeof(db1.min_count));
    // pre_fout.write((char*)&db1.max_count_lo, sizeof(db1.max_count_lo));
    // pre_fout.write((char*)&num_kmers_before, sizeof(num_kmers_before));

    // bool both_strands = false;
    // pre_fout.write((char*)&both_strands, sizeof(both_strands));

    // // 3 uchar padding
    // pre_fout.put('\0');
    // pre_fout.put('\0');
    // pre_fout.put('\0');

    // pre_fout.write((char*)&db1.max_count_hi, sizeof(db1.max_count_hi));

    // // 4 uint32s for padding
    // pre_fout.put('\0');
    // pre_fout.put('\0');
    // pre_fout.put('\0');
    // pre_fout.put('\0');
    // pre_fout.put('\0');
    // pre_fout.put('\0');
    // pre_fout.put('\0');
    // pre_fout.put('\0');
    // pre_fout.put('\0');
    // pre_fout.put('\0');
    // pre_fout.put('\0');
    // pre_fout.put('\0');
    // pre_fout.put('\0');
    // pre_fout.put('\0');
    // pre_fout.put('\0');
    // pre_fout.put('\0');
    // pre_fout.put('\0');
    // pre_fout.put('\0');
    // pre_fout.put('\0');
    // pre_fout.put('\0');

    // // version (uint32)
    // pre_fout.put('\0');
    // pre_fout.put('\0');
    // pre_fout.put('\0');
    // pre_fout.put('\0');

    // uint32 tmp = 64;
    // pre_fout.write((char*)&tmp, sizeof(tmp));

    // pre_fout.write(pre_marker.c_str(), 4);
    // suff_fout.write(suf_marker.c_str(), 4);

    return num_kmers_before;
}

int main(int argc, char **argv) {
    KMCDB db1(argv[1]);
    KMCDB db2(argv[2]);

    size_t num_threads = std::stol(argv[3]);

    if (db1.kmer_length != db2.kmer_length) {
        std::cerr << "K-mer lengths don't match: "
                  << db1.kmer_length << " != "
                  << db2.kmer_length << std::endl;
        exit(1);
    }

    if (db1.lut_prefix_len != db2.lut_prefix_len) {
        std::cerr << "Prefix lengths don't match: "
                  << db1.lut_prefix_len << " != "
                  << db2.lut_prefix_len << std::endl;
        exit(1);
    }


    uint64 prefix_len = db1.lut_prefix_len;
    uint64 max_prefix = (1llu << (prefix_len * 2));
    uint64 block_size = max_prefix / num_threads;
    uint64 total_kmer_count = 0;

    uint64 a = 0x0807060504030201;
    bool little_endian = ((uchar*)&a)[0] == 1 && ((uchar*)&a)[1] == 2 && ((uchar*)&a)[2] == 3 && ((uchar*)&a)[3] == 4 && ((uchar*)&a)[5] == 6 && ((uchar*)&a)[7] == 8;

    std::cout << "Num threads: " << num_threads << "\tlittle endian: " << little_endian << std::endl;

    #pragma omp parallel for num_threads(num_threads)
    for (uint64 min_prefix = 0; min_prefix < max_prefix; min_prefix += block_size) {
        uint64 local_max_prefix = std::min(max_prefix, min_prefix + block_size);
        std::string out_name_prefix = std::string(argv[4]) + std::to_string(min_prefix) + "-" + std::to_string(max_prefix);

        uint64 num_kmers_output = intersect_DBs(db1, db2, min_prefix, local_max_prefix, max_prefix, out_name_prefix, little_endian);

        std::lock_guard<std::mutex> lock(mu);
        std::cout << "Chunk " << out_name_prefix << " contains " << num_kmers_output << " kmers\n";
        total_kmer_count += num_kmers_output;
    }

    std::cout << "Intersection contains " << total_kmer_count << " kmers\n";
}
