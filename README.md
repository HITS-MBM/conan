Requirements
-------------

In order to successfully run CONAN you need to have working installations of:

 - Python 3, with numpy and scipy installed. (Anaconda 3+ recommended). Using an older version of scipy can cause an error in testing. The error will read something like: ValueError: Valid methods when the raw observations are omitted are 'single', 'complete', 'weighted', and 'average'. If you are NOT interestred in cluster analysis, feel free to ignore this error in testing.
 - GROMACS molecular dynamics engine, version 5.0 or later. GROMACS can be downloaded from the official website. (NB: GROMACS is only used for the first step of analysis, the MD trajectory can be obtained by other programs too.)
 - Gnuplot plotting tool. It can be downloaded from the official website.
 - Mencoder is obtained through MPlayer, which can be found on the official website.

In fact, even gnuplot and mencoder are optional (but recommended!).

Usage
-------------

CONAN reads in a single input file and runs. For example, if the input file is `input.inp`, write:

    ./conan.py input.inp
And you are good to go. Note that, for a first run, CONAN will run `gmx mdmat`, which will require you to choose a group of atoms.

Memory requirements
-------------

Since the involved datasets scale, in principle, as N_(frames) x N_(res)^2, you should be mindful of disk space and memory. The main data file on the disk is `dmf.xpm`, which scales roughly as ~N_(frames) x N_(res)^2 x 2B. For example, a 100-frame simulation on a 100-residue protein will take up about 2MB of disk space.

The _output_ of CONAN will be many `.png` files and also many `.dat` files, which contain the raw data for the former. The `.dat` files can take up a lot of space; about ~N_(frames) x N_(res)^2 x 20B, if all the raw matrices are to be stored. If this is too much, you can try the keyword `CLEAN_MATRICES` (creating `.png` files only) or `MATRICES no` and `MAKE_MOVIE no` (creating only aggregate plots).

The memory requirements of CONAN are usually not too extreme, since data storage is done in a sparse way. Still, if you find that this does not help,  consider using the "economy mode" (keyword `ECONOMY`), which does only basic statistics on the data.

In general, it is a good idea to run CONAN on a subset of your data (setting `DT` high, or restricting `BEGIN` or `END` to a part of the simulation) before running it on the entire set.

Testing CONAN
-------------

You can test your CONAN distribution by:

    cd test/
    ./test.py

And read the results on the screen. The test script cleans up after itself, unless there are errors. Please note the previous remark -- older versions of Scipy are known to return a ValueError on cluster analysis.

CONAN tutorials
-------------

There are some sample applications (input files and trajectories) in the folder `applications/`. 
Check out our tutorials on https://contactmaps.blogspot.de/ for step-by-step instructions.

Also check out other example input files that you can use on your own trajectories at https://github.com/HITS-MBM/conan/blob/master/docs/Example_inputs.md

Using parts of CONAN
-------------

You can use CONAN as a library by `import conan` in a Python script. The variable names are relatively easy to understand, although not everything is 100% user-friendly at the moment.

Where to reach us?
-------------

You can contact me (Csaba) at csaba.daday at h-its.org. Check out the blog at https://contactmaps.blogspot.de/ and the YouTube channel at https://tinyurl.com/Con4nMD .
