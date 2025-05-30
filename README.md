# ComputationalAstrophysicsFinal
The final project of Computational Astrophysics course

## Installation
1. Downloding from this git

2. Enter the `./ComputationalAstrophysicsFinal/` folder. Downloading the `tomlplusplus` library via

   ```bash
   git clone https://github.com/marzer/tomlplusplus.git util/tomlplusplus 
   ```

   

## Usage

### Initial Condition sampler

1. Compiling the setup file via

   ```bash
   make setup
   ```

   which will generate a `./setup` binary file in your current folder.

2. Executing `./setup` with specific parameter

   ```bash
   ./setup setup_mode SIMULATIONTAG
   ```

   The output initial condition file would be labeled in `SIMULATIONTAG_00000.h5`

   If `SIMULATIONTAG.setup` does not exist in current folder, a sample setup file correspond to given `setup_mode` would be generated automatically.

   For example, sampling a uniform box with `SIMULATIONTAG` = `testuniform`

   ```bash
   ./setup uniform testuniform
   ```

   Since `testuniform.setup` does not exist, it will generate a`testuniform.setup`.  After modifing the setup file, typing `./setup uniform testuniform` again, the initial condition file `testuniform_00000.h5` would be generated. Also, a parameter file `testuniform.in` will also been generated

### Main simulation suite

The main loop of simulation is written inside `simulation.cpp`. To compile it, type in

```bash
make
```

in your terminal. To start the simulation, type in

```bash
./simulation YOURPARAMS.in
```

Whenever a dumpfile is generated, the parameter file will update automatically.

### Test code compiling suite

If you want to test your code, please put those source code into `./test` folder. For example, if you have a code named `TestHello.cpp` that would be tested. Move the `TestHello.cpp` into `./test`. Then type

``` bash
make TestHello
```

The compiled binary would be generated with path:  `./build/test/TestHello`. 

Or, if you wish to compile all the test files inside `./test`, just simply typing

```bash
make testall
```

Then all of the source code inside `./test` would be compiled inside `./build/test/`. 

If you want to clean all of the current compiled binaries. Entering

```bash
make clean
```

