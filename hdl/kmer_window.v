module kmer_window
(
	input clk,
	input rst,

	input [0:0] raddr,
	input [0:0] waddr,

	input [119:0] in,
	output [119:0] out
);

	reg [119:0] buffer_r [1:0], buffer_ns [1:0];

	reg [119:0] buffer_out;

	integer i;
	
	always @* begin
		for (i = 0 ; i < 2 ; i = i+1)
			buffer_ns[i]		=		buffer_r[i]	;
		
		buffer_ns[waddr]		=		in		;
		buffer_out			=		buffer_r[raddr] ;		
	end


	assign out = buffer_out;

	always @(posedge clk) begin
		if (rst) begin
			for (i = 0 ; i < 2 ; i = i+1) begin
				buffer_r[i] <= 120'b0;
			end
		end
		else begin
			
			for (i = 0 ; i < 2 ; i = i+1) begin
				buffer_r[i] <= buffer_ns[i];
			end
		end
	end

endmodule
