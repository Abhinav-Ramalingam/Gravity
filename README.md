# N-Body Simulation and Comparison Tool

This project implements a simulation for N-body dynamics and provides tools for comparing the simulation results with reference outputs.

## Project Overview

1. **galsim.c**: A C program that simulates the dynamics of particles under gravitational influence.
2. **Makefile**: The Makefile provides build and execution instructions for the program and comparison tool.
3. **Input Files**: Input data for the simulation (e.g., `ellipse_N_03000.gal`).
4. **Output Files**: Results of the simulation (e.g., `results.gal`) and comparison files (e.g., reference output files).
5. **cgf**: A tool used to compare the generated simulation results with reference data.

## Requirements

- GCC or any compatible C compiler.
- A Unix-based operating system (Linux/macOS) for running the Makefile and commands.
- Input files such as `ellipse_N_03000.gal` located in `./input_data/` directory.

## Directory Structure

```
├── Makefile                 # Makefile to build the program and run comparisons
├── galsim.c                 # C program for N-body simulation
├── input_data/              # Directory containing input files
│   └── ellipse_N_03000.gal  # Example input data file
├── ref_output_data/         # Directory containing reference output data
│   └── ellipse_N_03000_after100steps.gal  # Reference output file
├── results.gal              # Output file generated by the simulation
└── cgf                      # Comparison tool executable
```

## How to Use

### 1. Build the Program

To compile the simulation program (`galsim`), use the `make` command:

```bash
make
```

This will generate the `galsim` executable by compiling `galsim.c` with the necessary libraries.

### 2. Run the Simulation

To run the simulation program (`galsim`), use the following command:

```bash
make run
```

This will execute the simulation with the input file `ellipse_N_03000.gal`, running 100 timesteps with a timestep size of `1e-5`. The output will be saved in the `results.gal` file.

- **Command breakdown**:
  - `3000`: Number of particles in the simulation.
  - `./input_data/ellipse_N_03000.gal`: Input data file.
  - `100`: Number of timesteps.
  - `1e-5`: Timestep size.
  - `0`: Placeholder for graphics (not used in this case).

### 3. Compare the Results

To compare the results of the simulation (`results.gal`) with the reference data, use the `cgf` tool:

```bash
make compare
```

This will compare `results.gal` (generated by `galsim`) with the reference file `ellipse_N_03000_after100steps.gal` located in the `ref_output_data/` directory. The `cgf` tool will output the error (if any) between the two files.

### 4. Clean Up

To clean up the generated files (including the compiled executable), run:

```bash
make clean
```

This will remove the `galsim` executable.

## Notes

- Ensure the paths to input files and reference output files are correct.
- Modify the `Makefile` if different parameters or input files are required for different simulations.
- Make sure the `cgf` tool is properly compiled and in the correct location for comparison.
