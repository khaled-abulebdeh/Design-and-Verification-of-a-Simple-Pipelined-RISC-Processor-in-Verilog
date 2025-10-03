# Pipelined RISC Processor Design in Verilog 💻

## Project Overview

This project presents the **Register-Transfer Level (RTL) design** and verification of a simple **32-bit pipelined RISC processor** implemented using **Verilog HDL**.

Developed as a core requirement for a Computer Architecture course, the design implements a classic **5-stage pipeline** and features robust **hazard control mechanisms**, including **data hazard forwarding** and **control hazard flushing**, to ensure correct and efficient instruction execution based on a predefined Instruction Set Architecture (ISA).

## Key Features

* **Pipelined Architecture:** Implements the standard 5-stage pipeline: **Instruction Fetch (IF), Instruction Decode (ID), Execute (EX), Memory (MEM), and Write-Back (WB)**.
* **RISC Design:** Simple, 32-bit instruction size and word size, utilizing 16 general-purpose registers (R0-R15).
* **Data Hazard Control:** Features an aggressive **Forwarding Unit** to bypass results from the EX and MEM stages directly to the EX stage inputs, effectively minimizing stall cycles.
* **Control Hazard Handling:** Employs a **Kill signal** to flush (insert NOP) the instruction immediately following a branch or jump, preventing incorrect instruction execution.
* **Verilog Implementation:** Modular and verifiable RTL code structured for use with industry-standard simulators.

***

## Technologies and Prerequisites

To successfully simulate and verify this project, you will need the following software installed:

| Component | Type | Requirement |
| :--- | :--- | :--- |
| **Language** | **Verilog HDL** | Core implementation language. |
| **Toolchain** | **Verilog Simulator** | Active-HDL, ModelSim, or Icarus Verilog is required for compilation and simulation. |
| **Architecture** | **RISC Processor** | Custom ISA defined in the project specification. |

### Installation and Setup

1.  **Install a Verilog Simulator:** Ensure you have access to a tool capable of compiling and simulating Verilog HDL.
2.  **Load the Code:** Load the `Processor.v` file into your chosen simulator environment. This single file contains the entire processor design and the integrated testbench.

***

## Usage: Running the Simulation

The project's verification is handled by the internal **Testbench** within `Processor.v`, which loads instruction sequences and generates the necessary clock and reset signals.

### 1. The Core RTL (`Processor.v`)

This file contains the complete modular Verilog implementation of the 5-stage pipeline, including the datapath, the control unit, the **Forwarding Unit**, and the **Hazard Detection Unit**.

1.  Open `Processor.v` in your Verilog simulator.
2.  The Testbench module is at the top level and handles simulation setup.

### 2. The Testbench (`testBench` module)

The `testBench` module automatically sets up the environment and runs a pre-loaded sequence of instructions to test the processor's functionality, including hazard management.

1.  In your simulator, select the `testBench` module for execution.
2.  Run the simulation. The Testbench will automatically apply reset and start the clock.

### 3. Verification

Verification is done by analyzing the internal state of the processor during simulation.

1.  Use the simulator's waveform viewer to inspect signals after execution.
2.  Check the **Program Counter (PC)** flow and the final contents of the **Register File** to ensure all instructions, branches, and memory operations completed correctly.

***
