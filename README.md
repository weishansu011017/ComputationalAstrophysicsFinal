# ComputationalAstrophysicsFinal

Final project for the **Computational Astrophysics** course in NTU (2025 spring)

------

## Installation

1. Clone this repository:

   ~~~bash
   git clone https://github.com/YOURUSERNAME/ComputationalAstrophysicsFinal.git
   ~~~

2. Enter the project directory:

   ~~~bash
   cd ComputationalAstrophysicsFinal
   ~~~

3. Download the `tomlplusplus` header-only library:

   ~~~bash
   git clone https://github.com/marzer/tomlplusplus.git util/tomlplusplus 
   ~~~

------

## Usage

### Makefile Support

All compilation processes are handled through `Makefile_templete`.
 For **Windows** users, use the Python-based alternative `Makefile.py`.

To configure your environment:

- Set the correct path for the `CXX` compiler and update `LIB_DIRS` accordingly.
- If using CUDA, ensure `CUDA_DIR` points to your CUDA installation.

------

3. - ### Initial Condition Sampler

     1. **Compile the setup tool**

        On **Unix/Linux/macOS**:
   
        ~~~bash
        make setup
        ~~~
        
        On **Windows**:
        
        ~~~powershell
        python .\Makefile.py setup
        ~~~
   
        This will generate an executable named `setup` (`setup.exe` on Windows).
        
     2. **Run the setup tool with mode and tag**
     
        On Unix/macOS:
     
        ~~~bash
        ./setup setup_mode SIMULATIONTAG
        ~~~
        
        On Windows (PowerShell or CMD):
        
        ~~~powershell
        .\setup setup_mode SIMULATIONTAG
        ~~~
     
        This creates an initial condition file named:
        
        ~~~bash
        SIMULATIONTAG_00000.h5
        ~~~
        
     3. **Automatic setup file generation**
   
        If `SIMULATIONTAG.setup` does not exist, a default version will be generated. You can then edit this file and rerun the command to regenerate the initial condition file.
     
        For example, to create a uniform box named `testuniform`:
     
        ~~~powershell
        .\setup uniform testuniform
        ~~~
        
        If `testuniform.setup` does not exist, it will be created.
         Modify it as needed, then re-run the same command to produce:
        
        - `testuniform_00000.h5` — initial condition
        - `testuniform.in` — simulation parameter file

------

### Main Simulation

To compile the main simulation:

~~~
make
~~~

To enable CUDA compilation:

~~~bash
make USE_CUDA=1
~~~

On **Windows**, use:

~~~powershell
python .\Makefile.py
~~~

To enable CUDA on Windows:

~~~powershell
python .\Makefile.py --gpu true
~~~

To run the simulation:

~~~bash
./simulation YOURPARAMS.in
~~~

Each time a dump file is created, the parameter file is automatically updated.

> **Note:**
>  GPU acceleration is only enabled if the `use_GPU` flag in `YOURPARAMS.in` is set to `true`.
>  **Do not** set `use_GPU = true` if the program was not compiled with CUDA support.

------

### Test Suite

To compile a single test file `TestHello.cpp` placed in the `./test` directory:

~~~bash
make TestHello
~~~

This creates:

~~~bash
./build/test/TestHello
~~~

To compile **all** test files in `./test`:

~~~bash
make testall
~~~

To clean up all compiled files:

~~~bash
make clean
~~~

The Python-based replacement (`Makefile.py`) on Windows supports the same build targets and functionalities as the standard `Makefile` used on Unix-based systems.

------

### Julia-based 2D simulation visualization tool

You can visualize the 2D simulation dumpfiles interactively using Julia.

- **Unix/macOS users**:

  ~~~
  julia ./util/visualization_2Dparticles_interactive.jl Sim_00*.h5
  ~~~

  This command will load all matching dumpfiles in the current folder.

- **Windows users**:

  ~~~powershell
  julia ./util/visualization_2Dparticles_interactive.jl PATH\TO\YOUR\RESULT\FOLDER
  ~~~
  
  Replace `PATH\TO\YOUR\RESULT\FOLDER` with the folder path that contains the simulation dumpfiles. The script will automatically search for `.h5` files inside that folder.

#### Julia Setup

This script requires the following Julia packages:

~~~julia
import Pkg
Pkg.add(["Makie", "GLMakie", "Observables", "GeometryBasics", "Statistics", "Glob", "Dates", "Printf", "HDF5
~~~

------

### Read dumpfile

The dumpfile reader is implemented in both **Python** and **Julia**, located in:

- `util/read_dumpfile.py` (Python version, wrapped as a class)
- `util/read_dumpfile.jl` (Julia version, implemented as a `struct`)

Both versions provide consistent access to the following data fields:

- `Table`: a table containing particle data with `Pandas`(Python)/`DataFrame`(Julia) (e.g., `x`, `y`, `vx`, `vy`, `m`, etc.)
- `params`: a dictionary (or `Dict`) of global simulation parameters (e.g., `t`, `utime`, `SimulationTag`, `Utot`, etc.)

The julia-based readers are used in interactive visualizations.
