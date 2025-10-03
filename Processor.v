module testBench(); 
  reg CLK, CLR;

  // Instantiate the processor
  processor p (
    .CLK(CLK),
    .CLR(CLR)
  );

  initial begin
    // Initialize signals
    CLK = 0;
    CLR = 0;
	
	$display("++++++++++++++++++Starting+++++++++++++++++++++++++");
    // Apply reset (active low)
    #5 CLR = 0;
    #5 CLR = 1; // Release reset

    // Clock generation loop
    repeat (50) begin // Adjust iteration count as needed
      #5 CLK = ~CLK;  // Toggle clock every 5 time units
    end

    $finish;
  end
endmodule


// --------------------------------------------------------------------------------------------------------------
	   
	 module processor(
  input CLK,
  input CLR
);

  // Declare all intermediate wires
	FetchSegment fetch_segment (
	    .Instruction(Instruction_in),
	    .nextPC(NPC_in),
	    .BTA(BTA),
	    .JRA(JRA),
	    .JLA(JLA),
	    .PCsrc(PCsrc),
	    .INSTRUCTIONsrc(INSTRUCTIONsrc),
	    .NextDoubleInstruction(NextDoubleInstruction),
	    .Kill(Kill),
	    .CLK(CLK),
	    .CLR(CLR),
	    .LDW_SDW_Signal(ID_EX_LDW_SDW_Signal_in),
		.DisablePC (DisablePC)
	  );  
	  
	  IF_ID_Buffer if_id (
    .Instruction_out(Instruction_out),
    .NPC_out(NPC_out),
    .Instruction_in(Instruction_in),
    .NPC_in(NPC_in),
    .CLK(CLK),
    .CLR(CLR),
	.DisableIR(DisableIR)
  );   
  
  DecodeSegment decode_segment (
    .A(A_in),
    .B(B_in),
    .Imm32(Imm32_in),
    .NextDoubleInstruction(NextDoubleInstruction),
    .BTA(BTA),
    .JRA(JRA),
    .JLA(JLA),
    .Rd2_in(Rd2_in),
    .Branch_flag(Branch_flag),
    .INSTRUCTIONsrc(INSTRUCTIONsrc),
    .PCSrc(PCsrc),
    .ALUsrc(ALUsrc_in),
    .MemRd(ID_EX_MemRd_in),
    .RegWr(ID_EX_RegWr_in),
    .MemWr(ID_EX_MemWr_in),
    .WBdata(ID_EX_WBdata_in),
    .LDW_SDW_Signal(ID_EX_LDW_SDW_Signal_in),
    .Stall(Stall),
    .DisableIR(DisableIR),
    .DisablePC(DisablePC),
    .Instruction(Instruction_out),
    .nextPC(NPC_out),
    .ALU_Result(EX_MEM_ALU_result_in),
    .Memory_Result(MEM_WB_Data_wb_in),
    .Data_w(MEM_WB_Data_wb_out),
    .Rd2_f(Rd2_out),
    .Rd3_f(Rd3_out),
    .Rd4_f(Rd4_out),
    .EX_LDW_SDW_Signal(ID_EX_LDW_SDW_Signal_out),
    .EX_MemRd(ID_EX_MemRd_out),
    .EX_RegWr(ID_EX_RegWr_out),
    .MEM_RegWr(EX_MEM_RegWr_out),
    .WB_RegWr(MEM_WB_RegWr_out),
	.ALUop(ALUOp_in),
    .CLK(CLK),
    .CLR(CLR)
  );
  
  ID_EX_Buffer id_ex (
    .Imm32_out(Imm32_out),
    .A_out(A_out),
    .B_out(B_out),
    .Rd2_out(Rd2_out),
    .Imm32_in(Imm32_in),
    .A_in(A_in),
    .B_in(B_in),
    .Rd2_in(Rd2_in),
    .ALUOp_in(ALUOp_in),
    .ALUsrc_in(ALUsrc_in),
    .RegWr_in(ID_EX_RegWr_in),
    .MemRd_in(ID_EX_MemRd_in),
    .MemWr_in(ID_EX_MemWr_in),
    .WBdata_in(ID_EX_WBdata_in),
    .LDW_SDW_Signal_in(ID_EX_LDW_SDW_Signal_in),
    .ALUOp_out(ALUOp_out),
    .ALUsrc_out(ALUsrc_out),
    .RegWr_out(ID_EX_RegWr_out),
    .MemRd_out(ID_EX_MemRd_out),
    .MemWr_out(ID_EX_MemWr_out),
    .WBdata_out(ID_EX_WBdata_out),
    .LDW_SDW_Signal_out(ID_EX_LDW_SDW_Signal_out),
    .CLK(CLK),
    .CLR(CLR),
	.Stall(Stall)
  );
  
  ExecuteSegment execute_segment (
    .ALU_result(EX_MEM_ALU_result_in),
    .A(A_out),
    .B(B_out),
    .Imm32(Imm32_out),
    .zero_flag(zero_flag),
    .ALUop(ALUOp_out),
    .ALUsrc(ALUsrc_out),
	.CLK (CLK)
  ); 
  
  EX_MEM_Buffer ex_mem (
    .ALU_result_out(EX_MEM_ALU_result_out),
    .Data_in_mem_out(Data_in_mem_out),
    .Rd3_out(Rd3_out),
    .ALU_result_in(EX_MEM_ALU_result_in),
    .Data_in_mem_in(B_out),
    .Rd3_in(Rd2_out),
    .RegWr_in(ID_EX_RegWr_out),
    .MemRd_in(ID_EX_MemRd_out),
    .MemWr_in(ID_EX_MemWr_out),
    .WBdata_in(ID_EX_WBdata_out),
    .RegWr_out(EX_MEM_RegWr_out),
    .MemRd_out(EX_MEM_MemRd_out),
    .MemWr_out(EX_MEM_MemWr_out),
    .WBdata_out(EX_MEM_WBdata_out),
    .CLK(CLK),
    .CLR(CLR)
  ); 
  
  
  MemorySegment memory_segment (
    .Data_wb(MEM_WB_Data_wb_in),
    .ALU_result(EX_MEM_ALU_result_out),
    .Data_in_mem(Data_in_mem_out),
    .WBdata_signal(EX_MEM_WBdata_out),
    .MemRd(EX_MEM_MemRd_out),
    .MemWr(EX_MEM_MemWr_out),
    .CLR(CLR),
    .CLK(CLK)
  );
  
  
  MEM_WB_Buffer mem_wb (
    .Data_wb_out(MEM_WB_Data_wb_out),
    .Rd4_out(Rd4_out),
    .Data_wb_in(MEM_WB_Data_wb_in),
    .Rd4_in(Rd3_out),
    .RegWr_in(EX_MEM_RegWr_out),
    .RegWr_out(MEM_WB_RegWr_out),
    .CLK(CLK),
    .CLR(CLR)
  );


 
  
  // IF_ID buffer wires
  wire [31:0] Instruction_out, Instruction_in;
  wire [31:0] NPC_out, NPC_in; 

  // ID_EX buffer wires
  wire [31:0] A_out, A_in, B_out, B_in, Imm32_out, Imm32_in;
  wire [3:0]  Rd2_out, Rd2_in;
  wire [2:0]  ALUOp_out, ALUOp_in;
  wire  ALUsrc_out, ALUsrc_in,zero_flag, ID_EX_RegWr_out, ID_EX_RegWr_in, ID_EX_MemRd_out, 
  ID_EX_MemRd_in, ID_EX_MemWr_out, ID_EX_MemWr_in, ID_EX_WBdata_out, ID_EX_WBdata_in, ID_EX_LDW_SDW_Signal_out, ID_EX_LDW_SDW_Signal_in;

  // EX_MEM buffer wires
  wire [31:0] EX_MEM_ALU_result_out, EX_MEM_ALU_result_in, Data_in_mem_out, Data_in_mem_in;
  wire [3:0]  Rd3_out, Rd3_in;
  wire EX_MEM_RegWr_out, EX_MEM_RegWr_in,  EX_MEM_MemRd_out, EX_MEM_MemRd_in, EX_MEM_MemWr_out, EX_MEM_MemWr_in, EX_MEM_WBdata_out, EX_MEM_WBdata_in; 

  // MEM_WB buffer wires
  wire [31:0] MEM_WB_Data_wb_out, MEM_WB_Data_wb_in;
  wire [3:0]  Rd4_out, Rd4_in;
  wire  MEM_WB_RegWr_out, MEM_WB_RegWr_in; 


  // Decode segment signals
  wire [31:0] NextDoubleInstruction, BTA, JRA, JLA;
  wire Branch_flag;
  wire [1:0] INSTRUCTIONsrc, PCsrc;
  wire Stall, DisableIR, DisablePC;
  wire Kill;


  

endmodule




	  
// ************************ Data Buffers (without control signals)	   
		
module MEM_WB_Buffer (Data_wb_out, Rd4_out, Data_wb_in, Rd4_in, RegWr_in,RegWr_out,CLK, CLR);
	output reg [31:0] Data_wb_out;
	output reg [3:0] Rd4_out;
	output reg  RegWr_out ;
	input [31:0] Data_wb_in;
	input [3:0] Rd4_in;
	input CLK ,CLR;
	input  RegWr_in ;
	always @ (posedge CLK, negedge CLR)
		begin
			if (~CLR)
				begin
					Data_wb_out=0;
					Rd4_out=0;
					RegWr_out=0;
				end			  
			else
				begin
					Data_wb_out=Data_wb_in;
					Rd4_out=Rd4_in;
					RegWr_out = RegWr_in;
				end	
		end	 
endmodule 

// ------------			 



module EX_MEM_Buffer (ALU_result_out, Data_in_mem_out, Rd3_out, ALU_result_in, Data_in_mem_in, Rd3_in,RegWr_in, MemRd_in, MemWr_in, WBdata_in , 
	RegWr_out, MemRd_out, MemWr_out, WBdata_out , CLK ,CLR);
	output reg [31:0] ALU_result_out, Data_in_mem_out;
	output reg [3:0] Rd3_out;
	output reg RegWr_out, MemRd_out, MemWr_out, WBdata_out ;
	input [31:0] ALU_result_in, Data_in_mem_in;
	input [3:0] Rd3_in;
	input CLK ,CLR;
	input  RegWr_in, MemRd_in, MemWr_in, WBdata_in ;
	always @ (posedge CLK, negedge CLR)
		begin
			if (~CLR)
				begin
					ALU_result_out=0;
					Data_in_mem_out=0;
					Rd3_out=0;
					RegWr_out=0;
					MemRd_out=0;
					MemWr_out=0;
					WBdata_out=0;
				end			  
			else
				begin
					ALU_result_out=ALU_result_in;
					Data_in_mem_out=Data_in_mem_in;
					Rd3_out=Rd3_in;
					RegWr_out=RegWr_in;
					MemRd_out=MemRd_in;
					MemWr_out=MemWr_in;
					WBdata_out=WBdata_in;
		
				end	
		end	
		
		always @(posedge CLK)
		$display("[Execute]ALU_result= %h",ALU_result_out);
endmodule
// ------------	





  

module ID_EX_Buffer ( Imm32_out, A_out, B_out, Rd2_out, Imm32_in, A_in, B_in, Rd2_in, 
	ALUOp_in , ALUsrc_in, RegWr_in, MemRd_in, MemWr_in, WBdata_in , LDW_SDW_Signal_in , ALUOp_out,
	 ALUsrc_out, RegWr_out, MemRd_out, MemWr_out, WBdata_out , LDW_SDW_Signal_out ,CLK, CLR, Stall );	
	output reg [31:0]  Imm32_out, A_out, B_out; 
	output reg [3:0] Rd2_out;
	output reg ALUsrc_out, RegWr_out, MemRd_out, MemWr_out, WBdata_out , LDW_SDW_Signal_out;  
	output reg [2:0] ALUOp_out;
	input [31:0]  Imm32_in, A_in, B_in;	
	input [3:0] Rd2_in;
	input ALUsrc_in , RegWr_in, MemRd_in, MemWr_in, WBdata_in , LDW_SDW_Signal_in;
	input [2:0] ALUOp_in ;
	input CLK, CLR, Stall;
	
	
	always @ (posedge CLK, negedge CLR)
		begin
			if (~CLR)
				begin
					Imm32_out=0;
					A_out=0;
					B_out=0;
					Rd2_out=0;
					ALUOp_out = 0;
					ALUsrc_out = 0;
					RegWr_out=0;
					MemRd_out=0;
					MemWr_out=0;
					WBdata_out=0;
					LDW_SDW_Signal_out=0;
				end	
			else if (Stall)
				begin
					Imm32_out=0;
					A_out=0;
					B_out=0;
					Rd2_out=0;
					ALUOp_out = 0;
					ALUsrc_out = 0;
					RegWr_out=0;
					MemRd_out=0;
					MemWr_out=0;
					WBdata_out=0;
					LDW_SDW_Signal_out=0;
				end	
			else
				begin
					Imm32_out=Imm32_in;
					A_out=A_in;
					B_out=B_in;
					Rd2_out=Rd2_in;
					ALUOp_out = ALUOp_in;
					ALUsrc_out = ALUsrc_in;
					RegWr_out=RegWr_in;
					MemRd_out=MemRd_in;
					MemWr_out=MemWr_in;
					WBdata_out=WBdata_in;
					LDW_SDW_Signal_out=LDW_SDW_Signal_in;
				end	
		end	  
endmodule 	  













// ------------
module IF_ID_Buffer (Instruction_out, NPC_out, Instruction_in, NPC_in,DisableIR, CLK, CLR);
	output reg [31:0] Instruction_out, NPC_out;
	input [31:0] Instruction_in, NPC_in;
	input CLK, CLR, DisableIR;
	
	always @ (posedge CLK, negedge CLR)
		begin
			if (~CLR)
				begin
					Instruction_out=0;
					NPC_out=0;
				end	
			else if (DisableIR)
			begin
					Instruction_out= Instruction_out;
					NPC_out= NPC_out;
				end
			else
				begin
					Instruction_out= Instruction_in;
					NPC_out= NPC_in;
				end	
		end
endmodule  	  

// --------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------  

// ***** Pipeline stages (segments)

module MemorySegment (Data_wb, ALU_result, Data_in_mem, WBdata_signal, MemRd, MemWr, CLR, CLK);
	output [31:0] Data_wb; // Data to be passed to the pipeline reigster, either Data_out_mem or ALU_result
	input [31:0] ALU_result; // IN case of memory write, it is the address of the data memory
	input [31:0] Data_in_mem;
	input WBdata_signal, MemRd, MemWr; // pipelined control signals
	input CLK, CLR;
	
	wire [31:0] Data_mem_out;
	DataMemory dm (Data_mem_out, Data_in_mem, ALU_result, MemRd, MemWr, CLK, CLR); 
	
	assign Data_wb= WBdata_signal? Data_mem_out: ALU_result;  
	
	always @(posedge CLK)
		$display("[Memory] Data = %h, MemRead = %b, MemWrite = %b, Data_wb = %h", 
         Data_wb, MemRd, MemWr, Data_wb);

endmodule

// ---------

module ExecuteSegment (ALU_result, A, B, Imm32, zero_flag, ALUop, ALUsrc, CLK);
	output [31:0] ALU_result; 
	input [31:0] A, B, Imm32; //inputs from the pipeline registers
	output zero_flag; //this control signal is generated by the ALU
	input [2:0] ALUop; //input: pipelined control signal
	input  ALUsrc, CLK; //input: pipelined control signal
	
	reg [31:0] in2_alu;// determine the 2nd source of the ALU  
	
	assign in2_alu= ALUsrc? Imm32: B;
	
	/*
	always @ (*)
	begin
		case (ALUsrc)
			2'b00:in2_alu=B; // ALUsrc=0 (Register)
			2'b01: in2_alu=	Imm32; //ALUsrc=1 (Imm32)
			2'b10: in2_alu=	Imm32+1'b1; //ALUsrc=2 (Imm32 +1)	
		endcase
		
	end
	*/   
	

	ALU alu(ALU_result, zero_flag, A, in2_alu, ALUop);
endmodule

// --------
	
module DecodeSegment(
    A, B, Imm32, NextDoubleInstruction, BTA, JRA, JLA,Rd2_in, Branch_flag, INSTRUCTIONsrc, PCSrc,ALUsrc,MemRd,RegWr,  
	MemWr, WBdata , LDW_SDW_Signal, Stall, DisableIR, DisablePC, Instruction, nextPC, ALU_Result, 
	Memory_Result, Data_w, Rd2_f, Rd3_f, Rd4_f, EX_LDW_SDW_Signal,EX_MemRd, EX_RegWr,
	MEM_RegWr,WB_RegWr,ALUop, CLK, CLR
);
	input [3:0] Rd2_f, Rd3_f, Rd4_f;
	input EX_LDW_SDW_Signal, EX_MemRd, EX_RegWr,MEM_RegWr,WB_RegWr;
	output [2:0] ALUop;
	output [1:0] INSTRUCTIONsrc, PCSrc; 
	output ALUsrc,MemRd, MemWr, WBdata, LDW_SDW_Signal, Stall, DisableIR, DisablePC ;

	// Called control signals	  
	
	Main_Control_Signal mcs( OPcode , ALUop, RegDst,RBsrc,
	ExtOp, INSTRUCTIONsrc, Branch_flag, ALUsrc, RegWr, MemRd, MemWr, WBdata ,PCSrc , LDW_SDW_Signal , EX_LDW_SDW_Signal); 
	
	Stall_Unit su(OPcode,ForwardA, ForwardB, EX_MemRd, Branch_flag, LDW_SDW_Signal, Stall, DisableIR, DisablePC); 
	
	Forwarding_Unit fu(OPcode,Rd2_f, Rd3_f, Rd4_f, Rs, Rt, Rd,EX_RegWr,MEM_RegWr,WB_RegWr, ForwardA,ForwardB);
	
	
    // Outputs
    output [31:0] A, B, Imm32, BTA, JRA, JLA, NextDoubleInstruction;
    output [3:0] Rd2_in; // Desination register of the current instrcution, used when it is in WB stage
    output reg Branch_flag;

    // Inputs
    input [31:0] Instruction, nextPC;
    input [31:0] ALU_Result, Memory_Result, Data_w;	// Forwarded data from stages 3, 4 and 5
    output RegWr;  //control signal forwarded from stage 5 to determine if the instruction in stage 5 wishes to write on the register file
    wire RBsrc, ExtOp, RegDst;
    input CLK, CLR;
    wire [2:0] ForwardA, ForwardB;	//control signals from the control unit

    // Internal wires for instruction fields
    wire [5:0] OPcode;
    wire [3:0] Rd, Rs, Rt;
    wire [13:0] Imm14;

    // Register file read data wires
    wire [31:0] BusA, BusB;	//read registers in temporary vars, to support Forwarding to A, B pipeline-registers

    // Assign instruction fields
    assign OPcode = Instruction[31:26];
    assign Rd = Instruction[25:22];
    assign Rs = Instruction[21:18];
    assign Rt = Instruction[17:14];
    assign Imm14 = Instruction[13:0]; 
	
	wire [3:0] Rd_plus1;
	wire [13:0] Imm14_plus1;

	assign Rd_plus1 = Rd + 4'd1;
	assign Imm14_plus1 = Imm14 + 14'd1;

	assign NextDoubleInstruction = {OPcode, Rd_plus1, Rs, Rt, Imm14_plus1};

    // Destination register output 
	assign Rd2_in= RegDst? 4'b1110: Rd;
	
	/*
    always @ (*) 
	begin
		case (RegDst)
		2'b00: Rd2=Rd;
		2'b01: Rd2=Rd+1;
		2'b10: Rd2=4'b1110;	 // R14 for CLL
		endcase	
	end
	*/
	
	//Support the mux to choose RB's input
	reg [3:0] RB;
	assign RB= (RBsrc)? Rd2_in: Rt;

    // Instantiate RegisterFile
    RegisterFile rf (
        .BusA(BusA),
        .BusB(BusB),
        .BusW(Data_w),
        .RA(Rs),
        .RB(RB),
        .RW(Rd4_f),
        .RegWrite(WB_RegWr),
        .CLK(CLK),
        .CLR(CLR)
		);

    // Immediate extension
    assign Imm32 = ExtOp ? {{18{Imm14[13]}}, Imm14} : {18'b0, Imm14};

    // Branch Target Address calculation
    assign BTA = nextPC + Imm32; 
	
 	//jump lable address 
 	assign JRA = A; // JR Rs
    assign JLA = nextPC + Imm32; // for JL and CLL instructions
	
	

    // Forwarding muxes
    assign A = (ForwardA == 3'b000) ? BusA :
               (ForwardA == 3'b001) ? ALU_Result :
               (ForwardA == 3'b010) ? Memory_Result :
			   (ForwardA == 3'b011) ? Data_w: nextPC;  
			   
    assign B = (ForwardB == 3'b000) ? BusB :
               (ForwardB == 3'b001) ? ALU_Result :
               (ForwardB == 3'b010) ? Memory_Result:
               (ForwardB == 3'b011) ? Data_w: 0;

    // Branch flag and special register assignments
    always @(*) begin
        Branch_flag = 1'b0;
        case (OPcode)
            6'b001010: Branch_flag = (A == 0) ? 1'b1 : 1'b0; // BZ Rs, Label
            6'b001011: Branch_flag = (A > 0) ? 1'b1 : 1'b0;  // BGZ Rs, Label
            6'b001100: Branch_flag = (A < 0) ? 1'b1 : 1'b0;  // BLZ Rs, Label
        endcase

    end	  
	
	always @ (posedge CLK)
		$display("[Decode] Instruction = %h, Rs = %d, Rt = %d, Rd = %d, A = %h, B = %h, Imm32 = %h", 
         Instruction, Instruction[21:18], Instruction[17:14], Rd2_in, A, B, Imm32);

endmodule


// ------------

module FetchSegment (Instruction, nextPC, BTA, JRA, JLA, PCsrc, INSTRUCTIONsrc, NextDoubleInstruction, Kill, CLK ,CLR , LDW_SDW_Signal, DisablePC );
	//Fetch stage produces 
	output [31:0] nextPC; 
	output reg [31:0] Instruction;
	input [31:0] BTA, JRA, JLA, NextDoubleInstruction ; // Branch target address, Jump target address 
	wire [31:0] currentPC; 
	input [1:0] PCsrc, INSTRUCTIONsrc;
	input CLK, CLR, Kill;
	input LDW_SDW_Signal, DisablePC;
	
	handlingPC hpc (currentPC, nextPC, JLA, JRA, BTA, PCsrc,LDW_SDW_Signal, DisablePC, CLK, CLR); //get currentPC, nextPC
	
	//currentPC is the address to be passed to InstructionMemory to read the instruction  
	reg [31:0] InstructionRead; 
	InstructionMemory im(InstructionRead, currentPC);//here, the instruction is fetched	
	
	//Store either the instruction or a bubble (no Op) on the instruction pipeline bufffer
	always @(posedge CLK)
		$display("[Fetch] PC = %h, Instruction = %h", currentPC, Instruction);
	
	always @ (*)
	begin
		case (INSTRUCTIONsrc)
			2'b00: Instruction=InstructionRead; //  normal execution
			2'b01: Instruction= NextDoubleInstruction; // next instruction from a double-instruction gotten in decode stade
			2'b10: Instruction= 32'h00000000; //bubble	
		endcase	 

	end
endmodule



module handlingPC (currentPC, nextPC,JLA, JRA, BTA, PCsrc,LDW_SDW_Signal, DisablePC,  CLK, CLR);
	output reg [31:0] nextPC, currentPC;
	input [31:0] BTA, JRA, JLA;
	input [1:0] PCsrc;	
	input CLK, CLR, LDW_SDW_Signal, DisablePC;
	always @(posedge CLK, negedge CLR)
	begin
		if (~CLR)
		begin
			currentPC=0;
			nextPC=1; //increment next PC
		end	
		else if (LDW_SDW_Signal || DisablePC)
		begin
			currentPC=currentPC; //hold data (disable pc), if double instruction, or "Load Delay"
			nextPC=nextPC;
		end
		else 
		begin
			if (PCsrc==0)//Normal PC
			begin
				currentPC=nextPC;
				nextPC= currentPC+1;//increment next PC
			end	 
			else if (PCsrc==1)//JLA
			begin
				currentPC=JLA;
				nextPC= currentPC+1;//increment next PC
			end	
			else if (PCsrc==2)//JRA
			begin
				currentPC=JRA;
				nextPC= currentPC+1;//increment next PC
			end	
			else if (PCsrc==3)//BTA
			begin
				currentPC=BTA;
				nextPC= currentPC+1;//increment next PC
			end
		end	
		$display ("PCsrc= %d", PCsrc);
	end		
endmodule 	 


// --------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------
// ****** Components used in the datapath

module InstructionMemory (Instruction, Address);
    input  [31:0] Address;
    output [31:0] Instruction;
	
    reg [31:0] memory [1023:0]; // 4KB instruction memory

    assign Instruction = memory[Address]; // Word addressing

    initial begin
        $readmemh("Program_instructions.txt", memory); // Load instructions from file
    end
endmodule


// ------------

module DataMemory (Data_out, Data_in, Address, MemRead, MemWrite, CLK, CLR);
	output [31:0] Data_out;
	input [31:0] Data_in;
	input [31:0] Address;
	input MemRead, MemWrite, CLK, CLR;
	reg [31:0] memory [1023:0];   //can't allocate (2^32)*32 memory size
	
	// read doesn't need a CLK, but MemWrite should be 1
	assign Data_out = (MemRead) ? memory[Address] : 32'd0; 
	 
	//write operation needs a CLK
	always @(posedge CLK, negedge CLR)
	begin 
		if (~CLR) 
			begin
				integer i;
				for (i=0; i<1023; i=i+1)
					memory[i]= 0;	
			end	
			
		else if (MemWrite)	// &&CLK 
		begin
			memory[Address]= Data_in;
			$display ("(Write on address= %h) - (Data= %h)",Address,  Data_in);	
		end
	end	
endmodule

// ------------

module RegisterFile (BusA, BusB, BusW, RA, RB, RW, RegWrite, CLK, CLR);
	reg [31:0] Registers [14:0]; // R0-R14, 32-bit width (R15 is the PC)
	output  [31:0] BusA, BusB; 
	input [31:0] BusW;
	input [3:0] RA, RB, RW;
	input RegWrite, CLK, CLR;
	
	//read operation doesn't need a CLK
	assign BusA = Registers[RA];
	assign BusB = Registers[RB];
	
	integer i;
	integer j;
	//write operation needs a CLK
	always @(posedge CLK, negedge CLR)
	begin 
		if (~CLR) 
			begin
				
				for (i=0; i<15; i=i+1)
					Registers[i]= 0;									  
			end	
			
		else if (RegWrite && RW!=15)// && CLK	
			Registers[RW]= BusW; 
			
		// Display all register values
		
		for (j = 0; j < 15; j = j + 1)
			$display("Reg[%0d] = %h", j, Registers[j]);
	end	
endmodule


// ------------
							  // shamasneh add the carry_flag
module ALU (ALU_result, zero_flag, in1, in2, ALUop);
    output reg [31:0] ALU_result; 
    output reg zero_flag;
    input [31:0] in1, in2;           
    input [2:0] ALUop;
    reg signed [31:0] signed_in1, signed_in2;  // For signed operations
    
    always @(*) begin
        // signed version of inputs for arithmetic operations
        signed_in1 = in1;
        signed_in2 = in2;
        
        case (ALUop)
            3'b000: begin // ADD operation (signed)
                ALU_result = signed_in1 + signed_in2;
            end
            
            3'b001: begin // SUB operation (signed)
               ALU_result = signed_in1 - signed_in2;
            end    
            
            3'b010: begin // Bitwise OR operation (unsigned)
                ALU_result = in1 | in2; 
            end     
            
            3'b011: begin // Compare operation (signed)
                if (signed_in1 > signed_in2)
                    ALU_result = 32'd1;
                else if (signed_in1 < signed_in2)
                    ALU_result = -32'd1;
                else 
                    ALU_result = 32'd0;
            end
            
            default: begin
                ALU_result = 32'd0;
            end
        endcase
		
		zero_flag = (ALU_result == 0);	
		
		$display ("Testing ALU: result= %h, (in1=%h, in2=%h, ALUop=%d)",ALU_result, in1, in2, ALUop );

    end	 
	
endmodule	



// -------------------------------------------------------------------------------------------------------------- 
// --------------------------------------------------------------------------------------------------------------

 
module Stall_Unit(OPcode,ForwardA, ForwardB, EX_MemRd, Branch_flag, LDW_SDW_Signal, Stall, DisableIR, DisablePC);
  input [5:0] OPcode;
  input EX_MemRd, Branch_flag, LDW_SDW_Signal;	
  input [2:0] ForwardA, ForwardB;
  output reg Stall, DisableIR, DisablePC ;


always @(*) begin  
	 
    if ((EX_MemRd==1) && (ForwardA==1 || ForwardB==1))   // load delay	
	begin
        Stall = 1; 
		DisableIR=1;
		DisablePC=1;
	end	
	
    else if (((OPcode==8)||(OPcode==9)) && LDW_SDW_Signal)	// LDW, SDW
	begin
        Stall = 0; 
		DisableIR=0;
		DisablePC=1;
	end	
	else if (( (OPcode==10)||(OPcode==11)||(OPcode==12)) && Branch_flag)	// All branch instructions (if taken)
	begin
        Stall = 1; 
		DisableIR=0;
		DisablePC=0;
	end
	
	else if ((OPcode==13)||(OPcode==14))	// JR, J label
	begin
        Stall = 1; 
		DisableIR=0;
		DisablePC=0;
	end
  	else
    begin
        Stall = 0; 
		DisableIR=0;
		DisablePC=0;
	end
end
endmodule


module Main_Control_Signal( OPcode , ALUOp, RegDst,RBsrc, ExtOp, INSTRUCTIONsrc, Branch_flag, ALUsrc, RegWr, MemRd, MemWr, WBdata ,PCSrc , LDW_SDW_Signal , EX_LDW_SDW_Signal);
		   input [5:0] OPcode;
		   output reg  ExtOp, RBsrc, RegWr, MemRd, MemWr, WBdata , LDW_SDW_Signal;
		   output reg [2:0] ALUOp ;
		   output reg [1:0] PCSrc, RegDst , INSTRUCTIONsrc; 
		   output reg ALUsrc;
		   input EX_LDW_SDW_Signal, Branch_flag;
		   
		always @(*) begin
		   case (OPcode)
 			
        6'b000000: begin 	 // op =0 		OR Rd, Rs, Rt 
            RegDst = 0; RBsrc = 0; ExtOp = 0; PCSrc = 0;
			ALUsrc = 0; ALUOp = 2;  RegWr = 1;
			MemRd = 0; MemWr = 0; WBdata = 0;	
        end
        6'b000001: begin 	  //op = 1		   ADD Rd, Rs, Rt
            RegDst = 0; RBsrc = 0; ExtOp = 0; PCSrc = 0;
			ALUsrc = 0; ALUOp = 0; RegWr = 1;
			MemRd = 0; MemWr = 0; WBdata = 0;
        end
        6'b000010: begin 	  //op = 2			 SUB Rd, Rs, Rt
            RegDst = 0; RBsrc = 0; ExtOp = 0; PCSrc = 0;
			ALUsrc = 0; ALUOp = 1; RegWr = 1;
			MemRd = 0; MemWr = 0; WBdata = 0;
        end
        6'b000011: begin 			  //op = 3		  CMP Rd, Rs, Rt
            RegDst = 0; RBsrc = 0; ExtOp = 0; PCSrc = 0;
			ALUsrc = 0; ALUOp = 3; RegWr = 1;
			MemRd = 0; MemWr = 0; WBdata = 0;

        end
        6'b000100: begin 			   //op = 4		ORI Rd, Rs, Imm
            RegDst = 0; RBsrc = 0; ExtOp = 0; PCSrc = 0;
			ALUsrc = 1; ALUOp = 2; RegWr = 1;
			MemRd = 0; MemWr = 0; WBdata = 0;

        end
        6'b000101: begin		   //op = 5		  ADDI Rd, Rs, Imm 
           RegDst = 0; RBsrc = 0; ExtOp = 1; PCSrc = 0;
			ALUsrc = 1; ALUOp = 0; RegWr = 1;
			MemRd = 0; MemWr = 0; WBdata = 0;

        end
        6'b000110: begin 	   //op = 6			LW Rd, Imm(Rs) 
            RegDst = 0; RBsrc = 0; ExtOp = 1; PCSrc = 0;
			ALUsrc = 1; ALUOp = 0; RegWr = 1;
			MemRd = 1; MemWr = 0; WBdata = 1;

        end
        6'b000111: begin 	   //op =7			SW Rd, Imm(Rs)
            RegDst = 0; RBsrc = 0; ExtOp = 1; PCSrc = 0;
			ALUsrc = 1; ALUOp = 0; RegWr = 0;
			MemRd = 0; MemWr = 1; WBdata = 0;
			
        end		  
		
		//-------------------------------------------------------------------------------------------------------------------
		
		6'b001000: begin 	//op =8		  LDW Rd, Imm(Rs)
            RegDst = 0; RBsrc = 0; ExtOp = 1; PCSrc = 0;
			ALUsrc = 1; ALUOp = 0; RegWr = 1;
			MemRd = 1; MemWr = 0; WBdata = 1;

        end	
		6'b001001: begin 	  //op =9		SDW Rd, Imm(Rs)
            RegDst = 0; RBsrc = 1; ExtOp = 1;  PCSrc = 0;
			ALUsrc = 1; ALUOp = 0;  RegWr = 0;
			MemRd = 0; MemWr = 1; WBdata = 0;

        end	 
		
		//-------------------------------------------------------------------------------
		 6'b001010: begin		 //op =10		  BZ Rs, Label

            RegDst = 0; RBsrc = 0; ExtOp = 1; ; PCSrc = (Branch_flag == 1) ? 2'b11 : 2'b00;
			ALUsrc = 0; ALUOp = 0; RegWr = 0;
			MemRd = 0; MemWr = 0; WBdata = 0;

        end	 
	
														  
		6'b001011: begin 							 //op =11	   BGZ Rs, Labe
           							   
			RegDst = 0; RBsrc = 0; ExtOp = 1;  PCSrc = (Branch_flag == 1) ? 2'b11 : 2'b00;
			ALUsrc = 0; ALUOp = 0; RegWr = 0;
			MemRd = 0; MemWr = 0; WBdata = 0;

        end			  
		
		
		    6'b001100: begin 		//op =12		 BLZ Rs, Label

            RegDst = 0; RBsrc = 0; ExtOp = 1; PCSrc = (Branch_flag == 1) ? 2'b11 : 2'b00;
			ALUsrc = 0; ALUOp = 0; RegWr = 0;
			MemRd = 0; MemWr = 0; WBdata = 0;

        end	
		    6'b001101: begin 		  //op =  13		 JR Rs		 // done
            RegDst = 0; RBsrc = 0; ExtOp = 0; PCSrc = 2;
			ALUsrc = 0; ALUOp = 0; RegWr = 0;
			MemRd = 0; MemWr = 0; WBdata = 0;

        end	  
			6'b001110: begin 	  //op = 14		   J Label	   // done

            RegDst = 0; RBsrc = 0; ExtOp = 1; PCSrc = 1;
			ALUsrc = 0; ALUOp = 0; RegWr = 0;
			MemRd = 0; MemWr = 0; WBdata = 0;

        end	
			 6'b001111: begin 	   //op =  15		   CLL Label
            RegDst = 1; RBsrc = 0; ExtOp = 1; PCSrc = 1;
			ALUsrc = 0; ALUOp = 0;  RegWr = 1;
			MemRd = 0; MemWr = 0; WBdata = 0;
        end	
        endcase
		end	   
		
		// To find INSTRUCTIONsrc signal
		always @ (*)
		begin
			if (OPcode==13 ||OPcode==14 ||OPcode==15 ) // Jump instructions 
				INSTRUCTIONsrc=2;
			else if (OPcode==10 ||OPcode==11 ||OPcode==12)	 // Branch instructions
				INSTRUCTIONsrc= (Branch_flag == 1) ? 2'b10 : 2'b00;
			else if (OPcode==8 ||OPcode==9)	 // LDW, SDW instructions
				INSTRUCTIONsrc= (LDW_SDW_Signal== 1) ? 2'b01 : 2'b00;
			else INSTRUCTIONsrc=0;	
		end	
		
	// To find LDW_SDW_Signal signal
		always @ (*)
		begin
			if ((OPcode!=8) && (OPcode!=9))	
				LDW_SDW_Signal = 0;
			else if ((OPcode == 8) || (OPcode==9)) 
				LDW_SDW_Signal= (EX_LDW_SDW_Signal == 1)? 0: 1;	 
				// if LDW_SDW_Signal==1  -->> first double instuction is in decode stage
				// if LDW_SDW_Signal==0  -->> second double instuction is in decode stage	
		end	
		
endmodule
		   



	

module Forwarding_Unit (OPcode,Rd2, Rd3, Rd4, Rs, Rt, Rd,EX_RegWr,MEM_RegWr,WB_RegWr, ForwardA,ForwardB);
	
	input [5:0] OPcode;
    input [3:0] Rd2 ;// Rd2: Destination register from ALU stage    
    input [3:0] Rd3;  // Rd3: Destination register from MEMORY stage
    input [3:0] Rd4; // Rd4: Destination register from WRITE BACK stage
    input [3:0] Rs;  // Rs in the decode stage
    input [3:0] Rt;  // Rt in the decode stage 
	input [3:0] Rd;  // Rd in the decode stage
    input EX_RegWr;    // RegWrite signal from ID/EX stage
    input MEM_RegWr;    // RegWrite signal from EX/MEM stage
    input WB_RegWr;        // RegWrite signal from WB stage 
    output reg [2:0] ForwardA; // Forwarding control for Rs 
    output reg [2:0] ForwardB;  // Forwarding control for Rt


always @(*) begin
    
	// Forward A:
	
	// Case 1: CLL instruction
	if (OPcode == 6'b001111) 
		ForwardA = 4;
    // Case 2: Forward from EX/MEM (Rd2)
    else if ((Rs != 15) && (Rs == Rd2) && (EX_RegWr)) 
		 ForwardA = 1;

    // Case 3: Forward from MEM/WB (Rd3)
    else if ((Rs != 15) && (Rs == Rd3) && (MEM_RegWr)) 
		ForwardA = 2;

    // Case 4: Forward from WB (Rd4)
    else if ((Rs != 15) && (Rs == Rd4) && (WB_RegWr))
		ForwardA = 3 ;
	else 
		ForwardA = 0;// Default: No forwarding
		
		
	// Forward B:
	
	// Case 1: CLL instruction
	if (OPcode == 6'b001111) 
		ForwardB = 4;
    // Case 2: Forward from EX/MEM (Rd2)
    else if ( ((Rt != 15) && (Rt == Rd2) && (EX_RegWr)) || ((OPcode==7 ||OPcode==9) && (Rd==Rd2) &&  (EX_RegWr)) ) 
		 ForwardB = 1;

    // Case 3: Forward from MEM/WB (Rd3)
    else if ( ((Rt != 15) && (Rt == Rd3) && (MEM_RegWr)) || ((OPcode==7 ||OPcode==9) && (Rd==Rd3) &&  (MEM_RegWr)) )
		ForwardB = 2 ;

    // Case 4: Forward from WB (Rd4)
    else if ( ((Rt != 15) && (Rt == Rd4) && (WB_RegWr)) || ((OPcode==7 ||OPcode==9) && (Rd==Rd4) &&  (WB_RegWr)) )
		ForwardB = 3;
	else 
		ForwardB = 0;// Default: No forwarding
end
endmodule

												  	