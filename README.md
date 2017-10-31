Welcome to StackEdit!
===================


Hey! I'm your first Markdown document in **StackEdit**[^stackedit]. Don't delete me, I'm very helpful! I can be recovered anyway in the **Utils** tab of the <i class="icon-cog"></i> **Settings** dialog.

----------


Requirements
-------------

In order to successfully run CONAN you need to have working installations of:

 - Python 3, with numpy and scipy installed.
 - GROMACS molecular dynamics engine version 5.0 or later. GROMACS can be downloaded from the official website. (NB: GROMACS is only used for the first step of analysis, the MD trajectory can be obtained by other programs too.)
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

Using parts of CONAN
-------------

You can use CONAN as a library by `import conan` in a Python script. The variable names are relatively easy to understand, although not everything is 100% user-friendly at the moment.