# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.25

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.25.2/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.25.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/hmustafa/git/KMC/3rd_party/cloudflare

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/hmustafa/git/KMC/3rd_party/cloudflare

# Include any dependencies generated for this target.
include CMakeFiles/zlib.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/zlib.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/zlib.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/zlib.dir/flags.make

CMakeFiles/zlib.dir/adler32.c.o: CMakeFiles/zlib.dir/flags.make
CMakeFiles/zlib.dir/adler32.c.o: adler32.c
CMakeFiles/zlib.dir/adler32.c.o: CMakeFiles/zlib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/hmustafa/git/KMC/3rd_party/cloudflare/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/zlib.dir/adler32.c.o"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/zlib.dir/adler32.c.o -MF CMakeFiles/zlib.dir/adler32.c.o.d -o CMakeFiles/zlib.dir/adler32.c.o -c /Users/hmustafa/git/KMC/3rd_party/cloudflare/adler32.c

CMakeFiles/zlib.dir/adler32.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/adler32.c.i"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/hmustafa/git/KMC/3rd_party/cloudflare/adler32.c > CMakeFiles/zlib.dir/adler32.c.i

CMakeFiles/zlib.dir/adler32.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/adler32.c.s"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/hmustafa/git/KMC/3rd_party/cloudflare/adler32.c -o CMakeFiles/zlib.dir/adler32.c.s

CMakeFiles/zlib.dir/compress.c.o: CMakeFiles/zlib.dir/flags.make
CMakeFiles/zlib.dir/compress.c.o: compress.c
CMakeFiles/zlib.dir/compress.c.o: CMakeFiles/zlib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/hmustafa/git/KMC/3rd_party/cloudflare/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/zlib.dir/compress.c.o"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/zlib.dir/compress.c.o -MF CMakeFiles/zlib.dir/compress.c.o.d -o CMakeFiles/zlib.dir/compress.c.o -c /Users/hmustafa/git/KMC/3rd_party/cloudflare/compress.c

CMakeFiles/zlib.dir/compress.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/compress.c.i"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/hmustafa/git/KMC/3rd_party/cloudflare/compress.c > CMakeFiles/zlib.dir/compress.c.i

CMakeFiles/zlib.dir/compress.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/compress.c.s"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/hmustafa/git/KMC/3rd_party/cloudflare/compress.c -o CMakeFiles/zlib.dir/compress.c.s

CMakeFiles/zlib.dir/crc32.c.o: CMakeFiles/zlib.dir/flags.make
CMakeFiles/zlib.dir/crc32.c.o: crc32.c
CMakeFiles/zlib.dir/crc32.c.o: CMakeFiles/zlib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/hmustafa/git/KMC/3rd_party/cloudflare/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object CMakeFiles/zlib.dir/crc32.c.o"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/zlib.dir/crc32.c.o -MF CMakeFiles/zlib.dir/crc32.c.o.d -o CMakeFiles/zlib.dir/crc32.c.o -c /Users/hmustafa/git/KMC/3rd_party/cloudflare/crc32.c

CMakeFiles/zlib.dir/crc32.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/crc32.c.i"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/hmustafa/git/KMC/3rd_party/cloudflare/crc32.c > CMakeFiles/zlib.dir/crc32.c.i

CMakeFiles/zlib.dir/crc32.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/crc32.c.s"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/hmustafa/git/KMC/3rd_party/cloudflare/crc32.c -o CMakeFiles/zlib.dir/crc32.c.s

CMakeFiles/zlib.dir/deflate.c.o: CMakeFiles/zlib.dir/flags.make
CMakeFiles/zlib.dir/deflate.c.o: deflate.c
CMakeFiles/zlib.dir/deflate.c.o: CMakeFiles/zlib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/hmustafa/git/KMC/3rd_party/cloudflare/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object CMakeFiles/zlib.dir/deflate.c.o"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/zlib.dir/deflate.c.o -MF CMakeFiles/zlib.dir/deflate.c.o.d -o CMakeFiles/zlib.dir/deflate.c.o -c /Users/hmustafa/git/KMC/3rd_party/cloudflare/deflate.c

CMakeFiles/zlib.dir/deflate.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/deflate.c.i"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/hmustafa/git/KMC/3rd_party/cloudflare/deflate.c > CMakeFiles/zlib.dir/deflate.c.i

CMakeFiles/zlib.dir/deflate.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/deflate.c.s"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/hmustafa/git/KMC/3rd_party/cloudflare/deflate.c -o CMakeFiles/zlib.dir/deflate.c.s

CMakeFiles/zlib.dir/gzclose.c.o: CMakeFiles/zlib.dir/flags.make
CMakeFiles/zlib.dir/gzclose.c.o: gzclose.c
CMakeFiles/zlib.dir/gzclose.c.o: CMakeFiles/zlib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/hmustafa/git/KMC/3rd_party/cloudflare/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object CMakeFiles/zlib.dir/gzclose.c.o"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/zlib.dir/gzclose.c.o -MF CMakeFiles/zlib.dir/gzclose.c.o.d -o CMakeFiles/zlib.dir/gzclose.c.o -c /Users/hmustafa/git/KMC/3rd_party/cloudflare/gzclose.c

CMakeFiles/zlib.dir/gzclose.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/gzclose.c.i"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/hmustafa/git/KMC/3rd_party/cloudflare/gzclose.c > CMakeFiles/zlib.dir/gzclose.c.i

CMakeFiles/zlib.dir/gzclose.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/gzclose.c.s"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/hmustafa/git/KMC/3rd_party/cloudflare/gzclose.c -o CMakeFiles/zlib.dir/gzclose.c.s

CMakeFiles/zlib.dir/gzlib.c.o: CMakeFiles/zlib.dir/flags.make
CMakeFiles/zlib.dir/gzlib.c.o: gzlib.c
CMakeFiles/zlib.dir/gzlib.c.o: CMakeFiles/zlib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/hmustafa/git/KMC/3rd_party/cloudflare/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building C object CMakeFiles/zlib.dir/gzlib.c.o"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/zlib.dir/gzlib.c.o -MF CMakeFiles/zlib.dir/gzlib.c.o.d -o CMakeFiles/zlib.dir/gzlib.c.o -c /Users/hmustafa/git/KMC/3rd_party/cloudflare/gzlib.c

CMakeFiles/zlib.dir/gzlib.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/gzlib.c.i"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/hmustafa/git/KMC/3rd_party/cloudflare/gzlib.c > CMakeFiles/zlib.dir/gzlib.c.i

CMakeFiles/zlib.dir/gzlib.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/gzlib.c.s"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/hmustafa/git/KMC/3rd_party/cloudflare/gzlib.c -o CMakeFiles/zlib.dir/gzlib.c.s

CMakeFiles/zlib.dir/gzread.c.o: CMakeFiles/zlib.dir/flags.make
CMakeFiles/zlib.dir/gzread.c.o: gzread.c
CMakeFiles/zlib.dir/gzread.c.o: CMakeFiles/zlib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/hmustafa/git/KMC/3rd_party/cloudflare/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building C object CMakeFiles/zlib.dir/gzread.c.o"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/zlib.dir/gzread.c.o -MF CMakeFiles/zlib.dir/gzread.c.o.d -o CMakeFiles/zlib.dir/gzread.c.o -c /Users/hmustafa/git/KMC/3rd_party/cloudflare/gzread.c

CMakeFiles/zlib.dir/gzread.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/gzread.c.i"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/hmustafa/git/KMC/3rd_party/cloudflare/gzread.c > CMakeFiles/zlib.dir/gzread.c.i

CMakeFiles/zlib.dir/gzread.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/gzread.c.s"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/hmustafa/git/KMC/3rd_party/cloudflare/gzread.c -o CMakeFiles/zlib.dir/gzread.c.s

CMakeFiles/zlib.dir/gzwrite.c.o: CMakeFiles/zlib.dir/flags.make
CMakeFiles/zlib.dir/gzwrite.c.o: gzwrite.c
CMakeFiles/zlib.dir/gzwrite.c.o: CMakeFiles/zlib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/hmustafa/git/KMC/3rd_party/cloudflare/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building C object CMakeFiles/zlib.dir/gzwrite.c.o"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/zlib.dir/gzwrite.c.o -MF CMakeFiles/zlib.dir/gzwrite.c.o.d -o CMakeFiles/zlib.dir/gzwrite.c.o -c /Users/hmustafa/git/KMC/3rd_party/cloudflare/gzwrite.c

CMakeFiles/zlib.dir/gzwrite.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/gzwrite.c.i"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/hmustafa/git/KMC/3rd_party/cloudflare/gzwrite.c > CMakeFiles/zlib.dir/gzwrite.c.i

CMakeFiles/zlib.dir/gzwrite.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/gzwrite.c.s"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/hmustafa/git/KMC/3rd_party/cloudflare/gzwrite.c -o CMakeFiles/zlib.dir/gzwrite.c.s

CMakeFiles/zlib.dir/inflate.c.o: CMakeFiles/zlib.dir/flags.make
CMakeFiles/zlib.dir/inflate.c.o: inflate.c
CMakeFiles/zlib.dir/inflate.c.o: CMakeFiles/zlib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/hmustafa/git/KMC/3rd_party/cloudflare/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building C object CMakeFiles/zlib.dir/inflate.c.o"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/zlib.dir/inflate.c.o -MF CMakeFiles/zlib.dir/inflate.c.o.d -o CMakeFiles/zlib.dir/inflate.c.o -c /Users/hmustafa/git/KMC/3rd_party/cloudflare/inflate.c

CMakeFiles/zlib.dir/inflate.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/inflate.c.i"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/hmustafa/git/KMC/3rd_party/cloudflare/inflate.c > CMakeFiles/zlib.dir/inflate.c.i

CMakeFiles/zlib.dir/inflate.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/inflate.c.s"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/hmustafa/git/KMC/3rd_party/cloudflare/inflate.c -o CMakeFiles/zlib.dir/inflate.c.s

CMakeFiles/zlib.dir/infback.c.o: CMakeFiles/zlib.dir/flags.make
CMakeFiles/zlib.dir/infback.c.o: infback.c
CMakeFiles/zlib.dir/infback.c.o: CMakeFiles/zlib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/hmustafa/git/KMC/3rd_party/cloudflare/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building C object CMakeFiles/zlib.dir/infback.c.o"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/zlib.dir/infback.c.o -MF CMakeFiles/zlib.dir/infback.c.o.d -o CMakeFiles/zlib.dir/infback.c.o -c /Users/hmustafa/git/KMC/3rd_party/cloudflare/infback.c

CMakeFiles/zlib.dir/infback.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/infback.c.i"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/hmustafa/git/KMC/3rd_party/cloudflare/infback.c > CMakeFiles/zlib.dir/infback.c.i

CMakeFiles/zlib.dir/infback.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/infback.c.s"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/hmustafa/git/KMC/3rd_party/cloudflare/infback.c -o CMakeFiles/zlib.dir/infback.c.s

CMakeFiles/zlib.dir/inftrees.c.o: CMakeFiles/zlib.dir/flags.make
CMakeFiles/zlib.dir/inftrees.c.o: inftrees.c
CMakeFiles/zlib.dir/inftrees.c.o: CMakeFiles/zlib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/hmustafa/git/KMC/3rd_party/cloudflare/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building C object CMakeFiles/zlib.dir/inftrees.c.o"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/zlib.dir/inftrees.c.o -MF CMakeFiles/zlib.dir/inftrees.c.o.d -o CMakeFiles/zlib.dir/inftrees.c.o -c /Users/hmustafa/git/KMC/3rd_party/cloudflare/inftrees.c

CMakeFiles/zlib.dir/inftrees.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/inftrees.c.i"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/hmustafa/git/KMC/3rd_party/cloudflare/inftrees.c > CMakeFiles/zlib.dir/inftrees.c.i

CMakeFiles/zlib.dir/inftrees.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/inftrees.c.s"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/hmustafa/git/KMC/3rd_party/cloudflare/inftrees.c -o CMakeFiles/zlib.dir/inftrees.c.s

CMakeFiles/zlib.dir/inffast.c.o: CMakeFiles/zlib.dir/flags.make
CMakeFiles/zlib.dir/inffast.c.o: inffast.c
CMakeFiles/zlib.dir/inffast.c.o: CMakeFiles/zlib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/hmustafa/git/KMC/3rd_party/cloudflare/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building C object CMakeFiles/zlib.dir/inffast.c.o"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/zlib.dir/inffast.c.o -MF CMakeFiles/zlib.dir/inffast.c.o.d -o CMakeFiles/zlib.dir/inffast.c.o -c /Users/hmustafa/git/KMC/3rd_party/cloudflare/inffast.c

CMakeFiles/zlib.dir/inffast.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/inffast.c.i"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/hmustafa/git/KMC/3rd_party/cloudflare/inffast.c > CMakeFiles/zlib.dir/inffast.c.i

CMakeFiles/zlib.dir/inffast.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/inffast.c.s"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/hmustafa/git/KMC/3rd_party/cloudflare/inffast.c -o CMakeFiles/zlib.dir/inffast.c.s

CMakeFiles/zlib.dir/trees.c.o: CMakeFiles/zlib.dir/flags.make
CMakeFiles/zlib.dir/trees.c.o: trees.c
CMakeFiles/zlib.dir/trees.c.o: CMakeFiles/zlib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/hmustafa/git/KMC/3rd_party/cloudflare/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building C object CMakeFiles/zlib.dir/trees.c.o"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/zlib.dir/trees.c.o -MF CMakeFiles/zlib.dir/trees.c.o.d -o CMakeFiles/zlib.dir/trees.c.o -c /Users/hmustafa/git/KMC/3rd_party/cloudflare/trees.c

CMakeFiles/zlib.dir/trees.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/trees.c.i"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/hmustafa/git/KMC/3rd_party/cloudflare/trees.c > CMakeFiles/zlib.dir/trees.c.i

CMakeFiles/zlib.dir/trees.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/trees.c.s"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/hmustafa/git/KMC/3rd_party/cloudflare/trees.c -o CMakeFiles/zlib.dir/trees.c.s

CMakeFiles/zlib.dir/uncompr.c.o: CMakeFiles/zlib.dir/flags.make
CMakeFiles/zlib.dir/uncompr.c.o: uncompr.c
CMakeFiles/zlib.dir/uncompr.c.o: CMakeFiles/zlib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/hmustafa/git/KMC/3rd_party/cloudflare/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Building C object CMakeFiles/zlib.dir/uncompr.c.o"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/zlib.dir/uncompr.c.o -MF CMakeFiles/zlib.dir/uncompr.c.o.d -o CMakeFiles/zlib.dir/uncompr.c.o -c /Users/hmustafa/git/KMC/3rd_party/cloudflare/uncompr.c

CMakeFiles/zlib.dir/uncompr.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/uncompr.c.i"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/hmustafa/git/KMC/3rd_party/cloudflare/uncompr.c > CMakeFiles/zlib.dir/uncompr.c.i

CMakeFiles/zlib.dir/uncompr.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/uncompr.c.s"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/hmustafa/git/KMC/3rd_party/cloudflare/uncompr.c -o CMakeFiles/zlib.dir/uncompr.c.s

CMakeFiles/zlib.dir/zutil.c.o: CMakeFiles/zlib.dir/flags.make
CMakeFiles/zlib.dir/zutil.c.o: zutil.c
CMakeFiles/zlib.dir/zutil.c.o: CMakeFiles/zlib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/hmustafa/git/KMC/3rd_party/cloudflare/CMakeFiles --progress-num=$(CMAKE_PROGRESS_15) "Building C object CMakeFiles/zlib.dir/zutil.c.o"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/zlib.dir/zutil.c.o -MF CMakeFiles/zlib.dir/zutil.c.o.d -o CMakeFiles/zlib.dir/zutil.c.o -c /Users/hmustafa/git/KMC/3rd_party/cloudflare/zutil.c

CMakeFiles/zlib.dir/zutil.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/zutil.c.i"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/hmustafa/git/KMC/3rd_party/cloudflare/zutil.c > CMakeFiles/zlib.dir/zutil.c.i

CMakeFiles/zlib.dir/zutil.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/zutil.c.s"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/hmustafa/git/KMC/3rd_party/cloudflare/zutil.c -o CMakeFiles/zlib.dir/zutil.c.s

CMakeFiles/zlib.dir/inffast_chunk.c.o: CMakeFiles/zlib.dir/flags.make
CMakeFiles/zlib.dir/inffast_chunk.c.o: inffast_chunk.c
CMakeFiles/zlib.dir/inffast_chunk.c.o: CMakeFiles/zlib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/hmustafa/git/KMC/3rd_party/cloudflare/CMakeFiles --progress-num=$(CMAKE_PROGRESS_16) "Building C object CMakeFiles/zlib.dir/inffast_chunk.c.o"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/zlib.dir/inffast_chunk.c.o -MF CMakeFiles/zlib.dir/inffast_chunk.c.o.d -o CMakeFiles/zlib.dir/inffast_chunk.c.o -c /Users/hmustafa/git/KMC/3rd_party/cloudflare/inffast_chunk.c

CMakeFiles/zlib.dir/inffast_chunk.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/inffast_chunk.c.i"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/hmustafa/git/KMC/3rd_party/cloudflare/inffast_chunk.c > CMakeFiles/zlib.dir/inffast_chunk.c.i

CMakeFiles/zlib.dir/inffast_chunk.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/inffast_chunk.c.s"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/hmustafa/git/KMC/3rd_party/cloudflare/inffast_chunk.c -o CMakeFiles/zlib.dir/inffast_chunk.c.s

CMakeFiles/zlib.dir/adler32_simd.c.o: CMakeFiles/zlib.dir/flags.make
CMakeFiles/zlib.dir/adler32_simd.c.o: adler32_simd.c
CMakeFiles/zlib.dir/adler32_simd.c.o: CMakeFiles/zlib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/hmustafa/git/KMC/3rd_party/cloudflare/CMakeFiles --progress-num=$(CMAKE_PROGRESS_17) "Building C object CMakeFiles/zlib.dir/adler32_simd.c.o"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/zlib.dir/adler32_simd.c.o -MF CMakeFiles/zlib.dir/adler32_simd.c.o.d -o CMakeFiles/zlib.dir/adler32_simd.c.o -c /Users/hmustafa/git/KMC/3rd_party/cloudflare/adler32_simd.c

CMakeFiles/zlib.dir/adler32_simd.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/adler32_simd.c.i"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/hmustafa/git/KMC/3rd_party/cloudflare/adler32_simd.c > CMakeFiles/zlib.dir/adler32_simd.c.i

CMakeFiles/zlib.dir/adler32_simd.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/adler32_simd.c.s"
	/opt/homebrew/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/hmustafa/git/KMC/3rd_party/cloudflare/adler32_simd.c -o CMakeFiles/zlib.dir/adler32_simd.c.s

# Object files for target zlib
zlib_OBJECTS = \
"CMakeFiles/zlib.dir/adler32.c.o" \
"CMakeFiles/zlib.dir/compress.c.o" \
"CMakeFiles/zlib.dir/crc32.c.o" \
"CMakeFiles/zlib.dir/deflate.c.o" \
"CMakeFiles/zlib.dir/gzclose.c.o" \
"CMakeFiles/zlib.dir/gzlib.c.o" \
"CMakeFiles/zlib.dir/gzread.c.o" \
"CMakeFiles/zlib.dir/gzwrite.c.o" \
"CMakeFiles/zlib.dir/inflate.c.o" \
"CMakeFiles/zlib.dir/infback.c.o" \
"CMakeFiles/zlib.dir/inftrees.c.o" \
"CMakeFiles/zlib.dir/inffast.c.o" \
"CMakeFiles/zlib.dir/trees.c.o" \
"CMakeFiles/zlib.dir/uncompr.c.o" \
"CMakeFiles/zlib.dir/zutil.c.o" \
"CMakeFiles/zlib.dir/inffast_chunk.c.o" \
"CMakeFiles/zlib.dir/adler32_simd.c.o"

# External object files for target zlib
zlib_EXTERNAL_OBJECTS =

libz.a: CMakeFiles/zlib.dir/adler32.c.o
libz.a: CMakeFiles/zlib.dir/compress.c.o
libz.a: CMakeFiles/zlib.dir/crc32.c.o
libz.a: CMakeFiles/zlib.dir/deflate.c.o
libz.a: CMakeFiles/zlib.dir/gzclose.c.o
libz.a: CMakeFiles/zlib.dir/gzlib.c.o
libz.a: CMakeFiles/zlib.dir/gzread.c.o
libz.a: CMakeFiles/zlib.dir/gzwrite.c.o
libz.a: CMakeFiles/zlib.dir/inflate.c.o
libz.a: CMakeFiles/zlib.dir/infback.c.o
libz.a: CMakeFiles/zlib.dir/inftrees.c.o
libz.a: CMakeFiles/zlib.dir/inffast.c.o
libz.a: CMakeFiles/zlib.dir/trees.c.o
libz.a: CMakeFiles/zlib.dir/uncompr.c.o
libz.a: CMakeFiles/zlib.dir/zutil.c.o
libz.a: CMakeFiles/zlib.dir/inffast_chunk.c.o
libz.a: CMakeFiles/zlib.dir/adler32_simd.c.o
libz.a: CMakeFiles/zlib.dir/build.make
libz.a: CMakeFiles/zlib.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/hmustafa/git/KMC/3rd_party/cloudflare/CMakeFiles --progress-num=$(CMAKE_PROGRESS_18) "Linking C static library libz.a"
	$(CMAKE_COMMAND) -P CMakeFiles/zlib.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/zlib.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/zlib.dir/build: libz.a
.PHONY : CMakeFiles/zlib.dir/build

CMakeFiles/zlib.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/zlib.dir/cmake_clean.cmake
.PHONY : CMakeFiles/zlib.dir/clean

CMakeFiles/zlib.dir/depend:
	cd /Users/hmustafa/git/KMC/3rd_party/cloudflare && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/hmustafa/git/KMC/3rd_party/cloudflare /Users/hmustafa/git/KMC/3rd_party/cloudflare /Users/hmustafa/git/KMC/3rd_party/cloudflare /Users/hmustafa/git/KMC/3rd_party/cloudflare /Users/hmustafa/git/KMC/3rd_party/cloudflare/CMakeFiles/zlib.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/zlib.dir/depend
