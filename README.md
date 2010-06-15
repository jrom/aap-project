# Awareness of Architecture in Programming
## Final Project: Protein Dock

The main code lives in sources/3D_Dock/progs/ but to run it we will need the **fftw** library installed.
All the steps are automated in bash scripts:

To run the scripts, you must go into the bin/ directory and run them from there.

The first one compiles and installs the **fftw** library in the sources/fftw-2.1.3/installation/ directory.

    cd bin/
    ./setup.sh

Once the previous step is done, the main program must be compiled, you can do it by:

    cd bin/
    ./compile.sh

This will build the code with the default optimization flags defined in sources/3D_Dock/progs/Makefile.

The last step is to run the program. For that you can just use the runner script:

    cd bin/
    ./run.sh

This will by default run the three testcases and output the elapsed time from the timed section of the code.
