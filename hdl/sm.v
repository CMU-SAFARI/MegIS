module sm
(
  input clk,
  input rst,

  // SM5
  input bckt_arrived,
  input endA,
  input endB,
  input AeqB,
  input AgtB,
  input AltB,
  input idx_end,
  
  output l2p_swap,
  output reg read_flash,
  output reg read_dram,
  output reg write_dram,
  output reg compare,
  output inc_switch_offset,
  output reg merger,
  output p_compare
);

  //SM 5
  reg [3:0] sm5_r, sm5_ns;
  reg [1:0] n_step_r, n_step_ns;

  assign inc_switch_offset = (sm5_r == 4'b1000) | (sm5_r == 4'b1001);
  assign l2p_swap = sm5_r == 4'b0001;
  assign p_compare = (sm5_r == 4'b0011) & (n_step_r == 2'b01);

  always @* begin
    sm5_ns = sm5_r;
    case (sm5_r)
      4'b0000: begin //idle
        read_flash = 0;
        read_dram = 0;
        write_dram = 0;
        compare = 0;
        merger = 0;

        if (bckt_arrived)
          sm5_ns = 4'b0001;
      end
      4'b0001: begin //prep
        sm5_ns = 4'b0010;
        n_step_ns = 2'b00;
      end
      4'b0010: begin //read_both
        read_dram = 1;
        read_flash = 1;
        write_dram = 0;

        if (endA | endB)
          sm5_ns = 4'b0111;
        else
          sm5_ns = 4'b0011;
      end
      4'b0011: begin //compare
        compare = 1;
        read_dram = 0;
        read_flash = 0;

        if (AeqB)
          sm5_ns = 4'b0100;
        else if(AgtB)
          sm5_ns = 4'b0101;
        else if(AltB)
          sm5_ns = 4'b0110;
      end
      4'b0100: begin //eq
        write_dram = 1;

        sm5_ns = 4'b0010;
      end
      4'b0101: begin //neq1
        read_dram = 1;
        read_flash = 0;

        if (endB)
          sm5_ns = 4'b0111;
        else
          sm5_ns = 4'b0011;
      end
      4'b0110: begin //neq2
        read_dram = 0;
        read_flash = 1;

        if (endA)
          sm5_ns = 4'b0111;
        else
          sm5_ns = 4'b0011;
      end
      4'b0111: begin //move_to_next_step
        n_step_ns = n_step_r+1;
        read_flash = 0;
        read_dram = 0;

        if (n_step_r == 2'b00)
          sm5_ns = 4'b1000;
        else if(n_step_r == 2'b01)
          sm5_ns = 4'b1001;
        else if(n_step_r == 2'b10)
          sm5_ns = 4'b1011;
      end
      4'b1000:
        sm5_ns = 4'b0010;
      4'b1001:
        sm5_ns = 4'b1010;
      4'b1010: begin //idx_merger
        merger = 1;
        write_dram = 1;

        if(idx_end)
          sm5_ns = 4'b0111;
        else
          sm5_ns = 4'b1010;
      end
      4'b1011:
        sm5_ns = 4'b0000;
    endcase
  end

  always @(posedge clk) begin
    if (rst) begin
      sm5_r <= 0;
      n_step_r <= 0;
    end
    else begin
      sm5_r <= sm5_ns;
      n_step_r <= n_step_ns;
    end
  end

endmodule