Executive summary
=================

This is the software I developed for Task 2.2 of the "TASTE Maintenance and
Evolutions" project of the European Space Agency. It is a port of the IDL
implementation of the StrayLight algorithm to C++, and it utilizes OpenMP,
Eigen and CUDA to achieve much better speeds than the original IDL code
(execution time in gdl: 18 seconds - execution time of my C++ code: 169 ms)

    Results for single precision speed - time per frame:

    $ ./configure && make && ./src/strayLight -b
    Execution time per image : 169 ms
    Average and std deviation: 171.23 +/- 1.6 ms
    
    Results for double precision speed - time per frame:

    $ ./configure --enable-double && make && ./src/strayLight -b
    Execution time per image : 222 ms
    Average and std deviation: 224.46 +/- 1.9 ms

Documentation
=============

Extensive documentation of the optimizations I performed to achieve these speedups
exists [inside this deliverable](doc/finalReport.pdf).

Usage:
======

    $ ./src/strayLight -h
    strayLight 1.1c (Thu May 30 14:57:08 2013)
    
    Usage: strayLight [OPTIONS]
    
      -h            this help
      -v            increase verbosity
      -V            show version
      -b            run benchmark (50 images)
      -i filename   instead of the ESA test image, process this file
      -c channel    process this channel
      -d N          dump log files from computation stages >= N
                    (default: 13, i.e. the final result is saved,
                     and never used during benchmarking (no output).

Validation of results
=====================
To create a log of the output image generated by the StrayLight algorithm:

    $ mkdir -p output
    $ ./src/strayLight

This will use the default value of option -d, and generate the stage13 output
file. If you wish to generate outputs from all stages, starting with stage 10:

    $ ./src/strayLight -d 10

This will generate log files for stages 10, 11, 12 and 13.

To compare the final output (stage13) with the outputs generated by IDL, run
this:

    $ ./contrib/elementsDiff.py \
        ./output/stage13_1 \
        ./output.from.IDL/stage13_1 | ./contrib/stats.py | grep Overall
    Overall: 6.61607071392e-06 +/- 375.7%

...which will verify that the difference between the values in the two outputs
(from IDL and from C++) is indeed minimal.

Contact point
=============
Thanassis Tsiodras, Dr.-Ing.
ttsiodras@gmail.com / ttsiodras@semantix.gr / a\_tsiodras@neuropublic.gr
