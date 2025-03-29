# SPICE Netlist Solver
## Overview

This repository contains two Python scripts for solving SPICE netlist files to determine node voltages in electrical circuits. The solvers use nodal analysis to construct a system of linear equations in the form $$G \cdot V = J$$, where:

-  G is the conductance matrix,
-  V  is the vector of node voltages, and
-  J  is the vector of current sources.

Depending on the script, the system is then solved using either dense or sparse matrix techniques.

- **`solver1.py`**: Solves simple circuit examples (`circuit1.sp` and `circuit2.sp`) using a dense matrix approach.
- **`solver2.py`**: Solves large-scale SPICE files from the [link](https://github.com/ASU-VDA-Lab/ML-for-IR-drop/tree/main/benchmarks/real-circuit-data) dataset using a sparse matrix approach to handle memory constraints.

The results for all processed SPICE files are saved in `results.txt`.

## Files in the Repository

- `circuit1.sp and circuit2.sp`: Two simple SPICE netlist files for test.
- `solver1.py`: Python script to solve simple circuits using a dense matrix approach.
- `solver2.py`: Python script to solve large-scale SPICE files using a sparse matrix approach.
- `results.txt`: Output file containing the node voltages for all processed SPICE files.

## Sanity Check
![image](https://github.com/user-attachments/assets/d36328be-d2e9-4caa-a757-ba24a1f11a7b)


## Algorithm
![image](https://github.com/user-attachments/assets/9e76e71e-c8b2-49f6-b2a4-f82ebadc0c1f)
