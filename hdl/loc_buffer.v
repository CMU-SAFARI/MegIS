module loc_buffer
(
	input clk,
	input rst,

	input [31:0] in,
	output [31:0] out
);

	reg [8:0] bigctr_r;
	reg smallctr_r;

	reg [63:0] buffer_r [511:0];

	reg [31:0] buffer_out;

	integer i;
	
	assign out = buffer_out;

	always @(posedge clk) begin
		if (rst) begin
			bigctr_r <= 0;
			smallctr_r <= 0;
		end
		else begin
			bigctr_r <= smallctr_r == 1'b1 ? bigctr_r + 1 : bigctr_r;
			smallctr_r <= smallctr_r == 1'b1 ? 1'b0 : smallctr_r + 1;
			
			for (i = 0 ; i < 512 ; i = i+1) begin
				if(bigctr_r == i) begin
					buffer_r[i][smallctr_r*32+:32] <= in;
				end
				else begin
					buffer_r[i] <= buffer_r[i];
				end
				if(bigctr_r + 1 == i) begin
					buffer_out <= buffer_r[i][smallctr_r*32+:32];
				end
				else begin
					buffer_out <= buffer_out;
				end
				
			end
		end
	end

endmodule
