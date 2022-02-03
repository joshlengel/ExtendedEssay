# Extended Essay Code Repository

## Setup

To properly execute this code and verify the results of my essay (found in the 'analysis' folder), one must satisfy the following requirements:

- Have a linux machine (Only tested on Ubuntu 21.04)
- Have a C++ compiler (Only tested with GCC 11.2.0)
- Have a Python interpreter (Only tested with Python 3.9.7)
- Have the packages 'matplotlib', 'numpy' and 'spiceypy' installed
- Have CMake accessible on the path (Only tested with CMake 3.18.4)
- Have a csh interpreter (Only tested with csh 20110502-6ubuntu1)
- Have Make accessible on the path (Only tested with GNU Make 4.3)

## Compiling

To compile the code, one must

1. cd into the 'deps/cspice' folder and execute `chmod +x makeall.csh && ./makeall.csh` (may require root privilege)
2. cd into the 'build' folder and execute `cmake ..`
3. in the same folder, execute `make`

Now the executable 'EE' should have been generated in the folder 'bin'.

## Execution

Note that this executable must be executed from the 'bin' folder for the relative paths in the program to work. This executable uses a considerable amount of RAM and CPU processing power so make sure to have around 2 GB of RAM and up to 10 minutes available before running. The results of the simulation will be available in the 'results.json' file under the 'out' directory (Note that this file can take up over 1 GB of disk space!).

## Analysis

To analyse the results in 'results.json', open the JupyterLab notebook 'EE.ipynb' in the 'analysis' folder in an editor of your choice (I used a Visual Studio Code extension which allowed me to switch between the notebook and my C++ code quickly and easily). Now click 'Run all' (or an equivalent button/command) and the appropriate data and graphs should be displayed. The graphs will also be saved to images in the 'out' directory.