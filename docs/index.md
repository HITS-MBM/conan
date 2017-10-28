# CONAN – Analysis of contacts in molecular dynamics trajectories

## Documentation

Contact analysis can be extremely useful to understand structure and dynamics of molecules and molecular complexes as they are atomistic-resolved descriptions of the interactions occurring within and between the investigated molecules.

CONAN is a tool developed for the statistical and dynamical analysis of contacts along molecular dynamics trajectories. It uses a combination of open-source tools for the computation of contacts along trajectories and for the creation of images and videos.

In particular:

- The computation of contacts over time is calculated using the tool mdmat,available in the molecular dynamics engine [GROMACS](http://www.gromacs.org/) (from version 5.0).
- The creation of graphs is achieved through the use of [gnuplot](http://www.gnuplot.info/).
- The creation of videos is achieved through the use of mencoder, which is part of the software [MPlayer](https://mplayerhq.hu/).
- The visualization of the cluster centers or the fitting of a series of structure output by CONAN can be achieved using python scripts interpreted by the molecular visualization software [UCSF Chimera](https://www.cgl.ucsf.edu/chimera/).

## What's in the distribution?

When you download CONAN, you will find two executables: `conan.py`, (the main program) and `conan_comp.py` (comparative CONAN). Additionally, you will find three directories called `example_input`, `gnuplot_inputscripts` and `chimera_visualization_scripts`. 

While the `gnuplot_inputscripts` contains gnuplot scripts. Feel free to change them.

In the `example_input` folder, you will find the following files:

- `mdmat.sample` – An example of the CONAN parameter file containing all the keywords that CONAN uses to perform several analyses. You can use it as a template. **Remember to set the PATH for each file that CONAN will use**.
- `compare.sample` – An example of the CONAN parameter file for the "comparative" CONAN script.
- `domains.txt` – An example of the domains file that contains the domains defined definitions for a molecule of ubiquitin.
- `stages.txt` – An example of the file used for calculating contact maps along defined intervals of the simulated trajectory.
- `observables_for_correlation.txt` - An example of the file containing a series of observables used for the calculation of the correlation coefficients between inter-residue distances and the chosen observables.
- `ubi.pdb` – An example coordinate file for ubiquitin. A coordinate file can be provided so that CONAN will write as B-factors the average interaction lifetime of contacts between residues. The structure can then be easily coloured by the B-factor in any software for molecular visualization (VMD, Pymol, Chimera...) so that the longer or shorter interaction lifetime can be readily visualized.
- `zoom_list.txt` – An example file with the list of residues to "zoom in" on.

In the folder named `chimera_visualization_scripts`, you will find the following scripts:

- `visualize_clusters.py`: this script allows you to visualize the clusters identified by CONAN if you require a cluster analysis of your trajectory. It will take two arguments: (i) the kmedoid\_summary.{metric}.txt fileand a coordinate file (in pdf format) matching your trajectory. The script can be launched by typing the following command:

```
chimera visualize_clusters.py summary.txt xxx.pdb
```

Launching the script will open chimera (that will prompt a message window asking which format are you wishing to open – no need for that, just click cancel) and will visualize, on the structure that you have specified as xxx.pdb, the cluster that have been identified by CONAN. Molecular structure (in cartoon representation) will be coloured according to each cluster and the atoms of the residues composing the center of the clusters will be shown using vdW radii. From the chimera window the user can then adjust the orientation of the molecule for the best visualisation and save a picture.

- `fit_structures.py`: this script opens all `*.gro` files that are in the same directory where the script is executed and returns the best structure alignment. The script can be launched by typing the following command:

```
chimera fit_structures.py
```

No arguments need to be specified as all the .gro files in the current working directory will be considered for the alignment.

## A few hints about what CONAN can do for you

CONAN has been developed in order to provide information about how inter-atomic contacts evolve over time during molecular dynamics simulations (MD). Beyond providing videos that show how contacts are formed or broken during MD, it can also read in a series of input files describing any observable against which the software calculates the correlation of the contacts made or broken. In this way it is easy for the user to identify events of interest along the simulated trajectory and relate them to the evolution of contact formation or ruptures. As an example, the formation of contacts can be related to increases of radius of gyration or end-to-end distance in a molecule, or to the formation, disruption or change of secondary structure (in this case, for example, the provided observables can be the φ or ψ angles of peptide chains). Additionally, contact maps can also be calculated in defined parts of the trajectory. As another example, in pulling simulations where the user identifies major rupture events of the tertiary and/or secondary structure of a molecule, CONAN can calculate contact maps for a defined amount of time, along the trajectory, defined by the user. The program in this case reads a file (through the keyword `STAGES` – see below) in which the user specifies beginning and end of the interval for which contact maps need to be calculated. In addition, a title for the set interval will be used to set the title of the contact map created.

Further functionalities (cluster analysis of residues or of frames in a trajectoryies, PCA, blocking analysis, ...) are also described below.

### Minimal requirements

In order to successfully run CONAN you need to have working installations of:

- Python 3, with [numpy](http://www.numpy.org/) and [scipy](https://scipy.org/) installed.
- GROMACS molecular dynamics engine version 5.0 or later. GROMACS can be downloaded from [the official website](http://www.gromacs.org/Downloads). (NB: GROMACS is only used for the first step of analysis, the MD trajectory can be obtained by other programs too.)
- Gnuplot plotting tool. It can be downloaded from [the official website](http://www.gnuplot.info/download.html).
- Mencoder is obtained through MPlayer, which can be found on [the official website](https://mplayerhq.hu/).

### Optional requirements

If you want to use the scripts provided for the visualization of the clusters onto a molecular structure or the fitting of the cluster conformations (see description of the scripts above)you will need:

- UCSF chimera, a molecular visualization software which can be found on [the official website at UCSF](https://www.cgl.ucsf.edu/chimera/download.html).

All of these softwares naturally run on Linux.

### Requirements on disk and memory

Given that CONAN is working in N^2 dimensional space, you should be aware that it **can get expensive** in terms of disk space and memory. As a rough guide, the following are the requirements of CONAN:

- Disk space: about N\_frames\*N\_res^2\*20B. For example, if your protein has 100 residues and you are including 100 frames in your analysis, expect about 20MB of storage to be taken. **Workaround**: the most of this disk space is taken up by the plain-text matrices. Consider the keyword `CLEAN_MATRICES` to get the animation and frames but no plain-text matrices (only one at a time will be produced, gnuplot will plot it, and CONAN will remove it, "cleaning up").
- Memory for general analysis: about N\_frames\*N\_pairs\*4B. N\_pairs is the total number of pairs that ever come within the main cutoff radius, `TRUNC`, of each other. In our experience, as a very rough guideline, N\_pairs is about 5-20% of the possible pairs (N\_pairs^2/2). **Workaround**: if this is too much, you can try the keyword `ECONOMY`, which will only store a tiny amount of data (perhaps a few megabytes), but no statistics are available.
- Memory for cluster analysis: about (N\_frames^2/2)\*4B. Unless you are reading in an extraordinary number of frames, this will not be the bottleneck.

## Functions

CONAN can perform the following types of analysis:

- Frame-wise contact maps and movies
  - _Related input keywords:_ `MATRICES`, `CLEAN_MATRICES`, `MAKE_MOVIE`
  - _Output:_ `frames/ matrices/ movies/`
- Differential contact maps and movies
  - _Related input keywords:_ `DR_MODE`, `MATRICES`, `CLEAN_MATRICES`, `MAKE_MOVIE`
  - _Output:_ `frames/ matrices/ movies/`
- "Aggregate" maps (average maps, time-encoded plots, ...)
  - _Related input keywords:_ `TRUNC`, `TRUNC_INTER`, `TRUNC_INTER_HIGH`,
  - _Output:_ `aggregate/`
- Pearson correlations with time, external observables, or inter-residue cross-correlation
  - _Related input keywords:_ `PEARSON_TIME`, `PEARSON_OBS`, `PEARSON_INTER`
  - _Output:_ `aggregate/`
- Cluster analysis of residues
  - _Related input keywords:_ `K_RES_CLUSTERS`
  - _Output:_ `cluster_res/`
- Cluster analysis of frames
  - _Related input keywords:_ `K_TRAJ_CLUSTERS`
  - _Output:_ `cluster_trj/`
- Principal component analysis
  - _Related input keywords:_ `PCA_MAKE`, `PCA_READ`
  - _Output:_ `pca/`
- "Zoomed-in" analysis on particular pairs of interest
  - _Related input keywords:_ `ZOOM_LIST`
  - _Output:_ `zoom/`
- "Blocking" analysis to estimate the correlation time of the contact map
  - _Related input keywords:_ `BLOCKING`
  - _Output:_ `blocking/`
- Asymmetric analysis, i.e., interactions between domains or between chains
  - _Related input keywords:_ `START_X`, `NRES_X`, `NTERM_X`, `DOMAINS_X`, (same for `Y`), `DIMER`

There is one more tools for contact analysis, distributed along with CONAN:

- Comparative CONAN: a small tool that compares two completed CONAN output directories. It creates asymmetric plots and plots of differences between various quantities (averages, standard deviations, interaction lifetimes, and contact formation/breaking). The runs should be done on the same number of residues and the same time interval.

## Input keywords

CONAN needs to have an input file that defines the parameters of the analysis.

An example input file containing most of the keywords and some explanation is contained in the distribution under: `example_input/mdmat.sample`

The keywords are **case-insensitive** (but are given in all-caps in this guide), but **the keys are case-sensitive** (since some of them are file names). If not specified, CONAN will either pick up the default value of that parameter and pass it to the routine used for the calculations, or will skip the kind of analysis defined by the omitted keyword. At the bottom of the list, we will give specific keywords applying only to the "asymmetric" runs or specific for the "comparative CONAN" tool.

### GNUS_PATH <path>

Path to the gnu files.

- _Output:_ -
- _Default:_ -
- _Remarks:_ If not given, CONAN will not produce any plots, just plain-text output.

### TRAJ <file name>

The trajectory file (`-f` option. Supported formats: `.pdb`, `.gro`, `.g96`, `.cpt`, `.xtc`, `.trr`, `.tng`).

- _Output:_ -
- _Default:_ -
- _Remarks:_ **compulsory** unless you set `RUN_MDMAT no` or `REREAD (path)`.

### COORD <file name>

The structure/coordinate file (`-s` option. Supported formats: `.pdb`, `.gro`, `.g96`, `.tpr`, `.brk`, `.ent`).

- _Output:_ -
- _Default:_ -
- _Remarks:_ **compulsory** unless you set `RUN_MDMAT no` or `REREAD (path)`.

### TRUNC <float (nm)>

The truncation value, i.e., any distance larger than this one will be set to this value.

- _Output:_ -
- _Default:_ 1.5
- _Remarks:_ -

### NLEVEL <int>

The number of discretization levels between 0 and `TRUNC` (maximum value: 7744). For example, if you want to use a truncation of 1.0 and a resolution of 0.001, set this to 1001. Setting `NLEVEL` to a high level will not affect the performance of mdmat or CONAN.

- _Output:_ -
- _Default:_ 40
- _Remarks:_ 40 is enough for visualization but probably not for statistical analysis!

### PATCH_TIME <yes|no>

Ignore time information from the `.xpm` file and just take the first frame and DT and keep increasing the time.

- _Output:_ -
- _Default:_ no
- _Remarks:_ This can be useful in case we analyze several concatenated trajectories or have duplicate time for some other reason. Can also be useful if the time delay between frames is less than 1 ps, since the .xpm file will contain values rounded to the nearest picosecond.

### MEAN <yes|no>

Should `gmx mdmat` create an average XPM contact map?

- _Output:_ `./dm.xpm`
- _Default:_ no
- _Remarks:_ ConAn will generate one anyway (`aggregate/average_mdmat.png`).

### INDEX <file name>

An index file.

- _Output:_ -
- _Default:_ (default Gromacs groups).
- _Remarks:_ In case you want to use a special group.

### BEGIN <time (ps)>

The first time in picoseconds (`-b`).

- _Output:_ -
- _Default:_ (last frame)
- _Remarks:_ Can be increased in a rerun.

### END <time (ps)>

The first time in picoseconds (`-e`).

- _Output:_ -
- _Default:_ (last frame)
- _Remarks:_ Can be decreased in a rerun.

### DT <time (ps)>

The time step to consider, i.e., only take frames that have (t - t\_begin)%dt == 0.

- _Output:_ -
- _Default:_ (all frames)
- _Remarks:_ Can be increased in a rerun.

### TRUNC_INTER <float (nm)>

The truncation threshold to turn interactions **on**.

- _Output:_ -
- _Default:_ TRUNC.
- _Remarks:_ This should be lower than TRUNC.

### TRUNC_INTER_HIGH <float (nm)>

The truncation threshold to turn interactions **off**.

- _Output:_ -
- _Default:_ `TRUNC_INTER`
- _Remarks:_ This should be higher than `TRUNC_INTER` but lower than `TRUNC`. The two cutoffs, `TRUNC_INTER` and `TRUNC_INTER_HIGH`, are meant to avoid interactions turning on and off only due to being close to one of the truncation thresholds (and thereby the _number of encounters_ being overestimated). The two thresholds can be set to be equal, though, if this "buffer zone" is undesirable.

### TRUNC_LIFETIME <float (0.00..1.00)>

Truncation of _lifetimes_ to consider interactions in the "interaction types" plot (if a `.pdb` file has been given).

- _Output:_ `aggregate/interaction_types.png`
- _Default:_ 0.50
- _Remarks:_ For example, `TRUNC_LIFETIME 0.5` would only plot interactions that are physically relevant and last at least 50% of the trajectory.

### DR_MODE <init|prev|both>

Which sort of differential contact maps should ConAn make?

- _Choices:_
  - `init`: the difference of each contact map with respect to the initial frame.
  - `prev`: the difference of each contact map with respect to the previous frame.
  - `both`: both :-)
- _Output:_ `frames/%05d\_dr.png`, `movies/ConAn/ConAn\_dr.mp4`
- _Default:_ (nothing) -> neither `init`, nor `prev`.
- _Remarks:_ `init` shows a general drift from the initial configuration. `prev` shows the change between consecutive frames, and is probably most useful when DT is quite large.

### TRUNC_DR <float (nm)>

Truncation of _differential contact maps_.

- _Output:_ `frames/%05d\_dr.png`, `movies/ConAn/ConAn\_dr.mp4`
- _Default:_ TRUNC
- _Remarks:_ If TRUNC\_DR is too low and/or DT is too high, the plots can be hard to see.

### RUN_MDMAT <yes|no>

Should ConAn run `gmx mdmat`?

- _Output:_ `dmf.xpm`
- _Default:_ yes
- _Remarks:_ If set to `no`, ConAn will look for a file `dmf.xpm` in the working directory and will ignore the keywords `TRAJ` and `STRUCT`.

### MATRICES <yes|no>

Should ConAn generate matrices? If set to `no`, ConAn will only generate aggregate data files.

- _Output:_ `matrices/%05d.dat`
- _Default:_ yes
- _Remarks:_ If set to `no`, ConAn will only generate aggregate data files, saving time and space. Of course, movies are also impossible in this case. In case you are interested in movies and plots, but want to save space, check out the keyword `CLEAN_MATRICES`.

### CLEAN_MATRICES <yes|no>

- _Output:_ `frames/%05d.png` (but nothing in `matrices/`)
- _Default:_ no
- _Remarks:_ If set to yes, ConAn will clean up the matrices (text files) right after using them for the frames (.png files). Execution will be slower but much less disk space will be used. Movies are still possible to make.

### NTERM <int>

N-terminus residue ID.

- _Output:_ -
- _Default:_ 1
- _Remarks:_ Note that residue IDs are always be contiguous but they will start at `NTERM`. You can still work with non-contiguous residue IDs, but only if you define domains (see `DOMAINS` keyword).

### PEARSON_TIME <yes|no>

The Pearson correlation with time will be plotted for each residue.

- _Output:_ `aggregate/pearson_time.png`
- _Default:_ no
- _Remarks:_ This can identify drifts in contacts.

### PEARSON_OBS <file name>

A file with a time dependence of some external observables. The first line of the file should be the number of observables N, the next N lines contain the titles of the observables, then the following lines contain N+1 columns time and the values of the observables). The times do not need to exactly match the ones from the frames but all the frames need to be present here.

- _Output:_ `aggregate/pearson_obs_%d.png`
- _Default:_ (not to do it)
- _Remarks:_ See also `IGNORE_OBS_TIME` below. An example file is in `example_input/area_gyrate_temperature_100ps.dat`.

### PEARSON_INTER <yes|no>

Should CONAN compute inter-residue cross-correlations?

- _Output:_ `aggregate/pearson_inter_residue.png`
- _Default:_ (not do to it)
- _Remark:_ If set to `yes`, also the clustering of residues will be performed based on this inter-residue cross-correlation (with d = 1-r^2).

### PCA_MAKE <int>

The number of principal components to analyze and project to.

- _Output:_ `pca/*`
- _Default:_ 0
- _Remarks:_ This analysis can take a while!

### PCA_READ <path>

Read in the principal components from another run and project this run onto it.

- _Output:_ `pca/*`
- _Default:_ (nothing)
- _Remarks:_ Note that the mean value of inter-residue distances will be different in different runs! A rigorous analysis would involve:
  1. Performing PCA on one group of frames (e.g., wild-type protein).
  2. Read in PC from 1., project frames from both groups (e.g., wild-type+mutant protein) onto these PCs.

### REREAD <path>

A path to existing matrices. This path must also contain a file called `index.dat`, specifying `NTERM`, `NRES`, and `TRUNC`. CONAN will read all matrices there. Matrices do not need to contain all entries, i.e., data can be sparse.

- _Output:_ -
- _Default:_ (not to do it)
- _Remarks:_ if you want to restart from a previous ConAn run, it is **faster** to go through the `dmf.xpm` file again (copy it/link it here and use `RUN_MDMAT no`)! Use this only if you have some special data file (COM distances, or some modified data?).

### ZOOM_LIST <file name>

A file name with a few interesting residue pairs to be "zoomed in" on. CONAN will plot the time evolution of the contact distance between each of these residues and the correlation with other contacts. Each line of the file should contain two numbers (residue IDs).

- _Output:_ `zoom/*`
- _Default:_ (not to do it)
- _Remarks:_ -

### MAKE_MOVIE <yes|no>

Should ConAn make a movie out of the frames?

- _Output:_ `frames/*`, `movies/*`
- _Default:_ yes
- _Remarks:_ In case something went wrong with the encoding, you can re-create them from `frames/`. Creating the frames can take some time.

### COORD_PDB <file name>

A pdb file for ConAn to read and understand the sequence of our protein.

- _Output:_ `aggregate/interaction\_types.png`, `aggregate/local\_interaction.{pdb|png}`
- _Default:_ (not to give it)
- _Remarks:_ This is needed for the "interaction type" plots and the "local interaction type" PDB, in which residues with shorter interactions will have low beta values.

### SHADOW_TOL <float>

Remove "shadowed contacts": consider only contacts (i,j) for which there is NO third residue k where

shadow\_tol\*(d\_ik + d\_jk) < d\_ij.

- _Output:_ (all output will be affected)
- _Default:_ (not to do it, i.e., `shadow_tol` = infinite)
- _Remarks:_ Any residues (i, j) for which there IS such a third residue k will be set to d_ij = `TRUNC`. Use with caution!

### IGNORE_OBS_TIME <yes|no>

Ignore the time variable in the observable file.
- _Output:_ -
- _Default:_ no
- _Remarks:_ This can be useful if the "trajectory" is actually a concatenated list of previous trajectories. In this case, however, **the number of observations and the number of frames in the contact map analysis must match exactly**! Also, the time column must still be there for consistency.

### K_RES_CLUSTERS <0|range of ints>

The number of clusters for residue-based cluster analysis of the protein. If not 0, the user should enter a range of values. Possible are formats as: `1 2 3`, `1-3`, or even `1-3 5 8-10` (which would mean 1 2 3 5 8 9 10). If 0, the number will be chosen interactively (when the dendrogram is displayed), using the same format as above.

- _Output:_ `cluster_res/*`
- _Default:_ (not to do it)
- _Remarks:_ See the description of the output for more details.

### K\_TRAJ\_CLUSTERS <0|range of ints>

The number of clusters for frame-based cluster analysis of the trajectory. ConAn will perform a k-medoid clustering of the trajectory, based on the RMSD between the contact maps.

- _Output:_ `cluster_trj/*`
- _Default:_ (not to do it)
- _Remarks:_ See `K_RES_CLUSTERS` for input format.

### DOMAINS <file name>

A file giving the domains of the protein. Not all residues must belong to a domain and technically one residue could belong to more than one domain (but the plot will be confusing). The file should contain 3 columns, the first two giving the first and last residue IDs and the third one the name of the given domain enclosed in quotes.

- _Output:_ `\*\*domains.png`
- _Default:_ (not to give it)
- _Remarks:_ Gnuplot syntax applies, for example, `{/Symbol a}` can be used for α. An example input file is in
`/example_input/domains.txt`

### STAGES <file name>

A file giving stages of the trajectory. ConAn will not perform a full analysis on each stage but rather plot the difference between the end and the beginning of each stage.

- _Output:_ `frames/stage_%05d_dr.png`
- _Default:_ (not to give it)
- _Remarks:_ These stages can be overlapping. The format of the stage file is the same as the one describing the domains. This can give a quick overview of what happened in the trajectory (for example, one can divide a longer trajectory into 10 equidistant "stages").

### ECONOMY <yes|no>

CONAN can be run in "economy mode", i.e., the data is not stored in memory in raw format but only aggregate quantities. This means, however, that only frames, matrices, movies, and aggregate plots will be printed. Use this if you have too many residues or too many frames.

### Specific keywords for asymmetric runs

Instead of a symmetric run (comparing intra-chain contacts), CONAN can be used to analyze the interface between two sets of residues, group X and group Y (to be shown in the X- and Y-axes of the matrices). The average number of inter-group interactions will be plotted and inserted in a PDB file of each group if given. CONAN will automatically switch to "asymmetric mode" if any of the below keywords are given. Note that all plots by CONAN are square formatted, so if the residue numbers between X and Y are very different, the contact maps can look distorted.

### START_X <int>

The starting residue of _group X_. 1-indexed, based on what is in the XPM file.

### NRES_X <int>

The number of residues in _group X_ (i.e., _group X_ will be `START_X` <= index `START_X` + `NRES_X` where index is 1-indexed.)

### NTERM_X <int>

The N-terminus index of _group X_ for plotting purposes and to read the sequence from the PDB (i.e., the residue IDs will be: `NTERM_X` <= resID < `NTERM_X + NRES_X`).

- _Default:_ equal to `START_X`

### DOMAINS_X <file name>

Domain file for plotting of _group X_. This is based on the `NTERM_X` index.

### COORD_PDB_X <file name>

PDB coordinate file for _group X_. The residue IDs in this file must correspond to `NTERM_X`.

### START_Y, NRES_Y, NTERM_Y, DOMAINS_Y, COORD_PDB_Y

As above. If only one of `COORD_PDB_X`, `COORD_PDB_Y` is given, the same PDB will be used for both groups X and Y.

### DIMER

"Dimer mode": if both groups X and Y are in fact different copies of the identical group (i.e., a homodimer), CONAN can take this into account for read-in and cluster analysis. This requires that the residue IDs be equal.

### Example for an asymmetric run

Suppose we are working with a PDB of 3000 residues but want to analyze the contacts formed between residues 1001..1050 and residues 2001..2050 while keeping this numbering but ignoring all other pairs. We can proceed as follows:

1. We create a special index file for these residues only.

  - In the `.xpm` file these appear as 1..50 and 51..100, respectively. The `.xpm` file will contain a symmetric matrix, but CONAN will ignore the "on-diagonal" parts.

2. We set:

```
START_X### 1
NRES_X    50
NTERM_X 1001
START_Y   51
NRES_Y    50
NTERM_Y 2001
```

Note that if we want to use a PDB file, conan.py will look for 1001..1050 and 2001..2050 in the residue ID column. 

## Specific keywords for comparative runs

The "comparative CONAN" script is a simple tool intended to give a quick comparison between two runs ("run A" and "run B") on similar systems. The two runs must have the same number of residues\* and it is recommended that they also have the same simulation time.

\* of course, if the residues are similar in number, the user him/herself can create an index file to create the same number of residues, by e.g., skipping over gaps in alignment.

The input supports the following keywords.

### RUN_A <path>

The location of "run A". **Compulsory.**

### RUN_B <path>

The location of "run B". **Compulsory.**

### GNUS_PATH <path>

The location of the gnuplot scripts. **Compulsory.**

### TITLE_A <string>

The title of "run A" between double quotes. Optional; default is "run A".

### TITLE_B <string>

See description of `TITLE_A`.

### R_CUT <float (nm)>

Consider only pairs of residues for which the average r_ij < r_cut for at least one of the two runs. Optional; default is considering all residue pairs. This is especially useful for changes in distances/RMSF, where large changes between large distances are possibly uninteresting. In particular, the RMSF of a pairwise distance can change drastically, depending on which side it is of the original truncation radius, so it makes sense to include some `R_CUT` < `TRUNC`.

### LIFE_CUT <float (nm)>

Consider only pairs of residues for which the interaction lifetime is larger than `life_cut` for at least one of the two runs. Optional, default is considering all residue pairs. Similarly to `R_CUT`, `LIFE_CUT` can restrict the analysis to "interesting" pairs. If both `R_CUT` and `LIFE_CUT` are given, then residue pairs must pass both tests to be considered.

Note: all "lifetime" values will be taken from the original analyses. The "comparative CONAN" script does not process any contact maps, just compares aggregate quantities it found in the folders!

## Output

CONAN will generate output in a series of folders.

### backup.{i}/

In case CONAN is run over an older CONAN run, all of the below folders will be moved into backup folders (`backup.1` is the oldest run, then `backup.2`, ...), including the `input/` folder, for users who prefer not to give meaningful names to input files. `dmf.xpm` is not copied into the backup folder to save space.

### movies/

- mp4: the animated map with the time-evolution of contacts.
- ConAn\_dr.init.mp4: the animated map with the time-evolution of the change in contacts compared to the initial frame [i.e., r\_ij(t) - r\_ij(t\_0)]
- ConAn\_dr.init.mp4: the animated map with the time-evolution of the change in contacts compared to the previous frame [i.e., r\_ij(t) – r\_ij(t - Δt)].
- (same files with .domains.mp4): animated maps with the domains shown.

_Hint:_ movies can be turned off with the keyword `MAKE_MOVIE no`.

### frames/

This folder contains the .png files with frames that form the movies mentioned above.

_Hint:_ the keyword `MAKE_MOVIE no` also turns off the generation of these frames.

### matrices/

This folder contains the raw data (i.e., distances) that form the frames before. The format is `<i> <j> <r_ij>` (or `i j dr_ij`). Note: there is an empty line between sections with different `i`'s, to respect **gnuplot** requirements.

_Hint:_ the keyword `MATRICES no` turns off the generation of these data files as well. The keyword `CLEAN_MATRICES yes` cleans up the raw data files but keeps the `.png` frames and will generate the movies.

### input/

The input file will be copied here for later reference.

### aggregate/

This folder has the most important plots and raw data files:

- `time_mdmat_rmsd.dat`: the time-evolution of the RMSD with respect to the initial frame and with respect to the previous frame. The format is: `t RMSD (t, t0) RMSD (t, t – Δt)`.
- RMSD (t, t0) compares each structure to the initial frame (i.e., it is the norm of the difference matrices `xxxxx_dr.init.dat`) and it is plotted in the file `rmsd_contact_map_initial.png`. RMSD (t, t – Δt)
compares each structure to the previous frame (i.e., it is the norm of the difference matrices `xxxxx_dr.prev.dat`) and it is plotted in the file `rmsd_contact_map_previous.png`. RMSD w.r.t. the initial frame can help diagnose (lack of) equilibration, while RMSD w.r.t. the previous frame (assuming Δt is rather large) can help identify key moments in the simulation.
- `mdmat_average_rmsf.dat` the average and standard deviation of each inter-residue distance. The format is: `i j <r_ij> σ_ij` The average contact map is plotted in `avg_mdmat.png` while the standard deviation in the contact map is plotted in `stddev_mdmat.png`. This latter plot is a pairwise alternative of the more commonly used RMSF.
- The file `timeline.dat` contains information on contact formation and loss. The format is:
`i j t_first t_last t_life t_max-t_0 N_meet`
where:
  - `t_first` shows the **first** time these two residues formed an interaction (if no such interaction has been made, the time of the last frame is used) This information is plotted in `timeline_first_encounter.png`. This can be important for classifying _folding_ trajectories.
  - `t_last`  is the **last** time these two residues formed an interaction (if no such interaction has been made, the time of the first frame is used). This information is plotted in `timeline\_last\_encounter.png`. This can be important for classifying unfolding trajectories.
  - `t_life` is the fraction of frames in which these two residues formed an interaction (between 0 and 1). This information is plotted in `interaction_lifetime.png`.
  - `t_max- t_0` is a constant column for plotting purposes only (the total time of the trajectory).
  - `N_meet` is the number of meetings of these two residues (under the assumption that Δt is small enough). This information is plotted in `num_encounter.png`. Furthermore, the average encounter time, `t_encounter = t_life *(t_max - t_0)/N_meet`, is plotted in: `avg_encounter.png`
  - A file named `interaction_types.png`, that shows a contact map in which the type of interaction (defined as hydrogen bonds, salt bridges or hydrophobic contacts) are in place. The interaction types are calculated for contacts having a lifetime longer (in percentage) than the one defined in the parameter file with the keyword `TRUNC_LIFETIME` (see below).
  - A series of files showing contact maps reporting the correlation coefficient between inter-residue distances and each of the given observable(s) provided by the user in a text file that CONAN can read (through the keyword `PEARSON_OBS` – see below). These files are called `pearson_time.png` and `pearson_obs_{i}.png`.
  - A plain-text file, `native_contacts.dat`, with the number of native and non-native contacts. Native contacts are simply any interactions that are below `TRUNC` in the first frame.

### cluster_res/{method}

This folder will contain the results of the cluster analysis of _residues_. Three different methods are available, `{method} = crosscorr`, `distance`, or `lifetime`. For each of these options, there will be a folder with:

- `linkage.txt`: the linkage created by Ward's clustering algorithm. (see `scipy.cluster.hierarchy.linkage`) for the description of the 4xN matrix.
- `dendrogram.png`: A dendrogram of residues (or the top 100 nodes, if there are more than 100 residues).

After the user has chosen a desired set of numbers of frames, there will be created a new folder for each number chosen. Under each of these folders, there will be three more files:

- `assignment.dat`: the assignment of each residue to a cluster.
- `frames_cluster_assignment.png`: the plotted version of the assignment data.
- `summary.txt`: summarizes the clustering result. Readable by the Chimera script.

### cluster_trj/

This folder will contain the results of the cluster analysis of _frames_ (trajectories).

- `linkage.txt`: as described in `cluster_res`.
- `rmsd_interframe.dat`: the raw data for interframe RMSD, with the format `t_i t_j RMSD(t_i, t_j)` (in ns and nm, respectively). Plotted in `rmsd_interframe.png`. If there are more than 1000 frames, CONAN will only print out approximately 1000 frames (N/( (int) N/1000)).

After the user has chosen a desired set of numbers of frames, there will be created a new folder for each number chosen. Under each of these folders `n` (showing the subdivision into `n` clusters), there will be another file and `n` folders:

- `assignment.dat`: the assignment of each frame to a cluster. Plotted in `frames_cluster_assignment.png` (color coded by similarity to the cluster center and showing the cluster centers).

In the subfolders `n/cluster0..n/cluster(n-1)`, one can find:

- The average contact map in `average_map.dat`, plotted in `average_map.png`.
- The contact map of the central frame plotted (as a link) in: `center.png` (ideally, this should be very similar to the average map).
- The contact map of all the frames in the cluster (linked) in `{index}.png`
- The interaction lifetime in `lifetime.dat`, plotted in `lifetime.png` (this is useful in case the clusters are very diffuse).

### pca/

The output of the principal component analysis, as:

- `pca_output.txt`: a file with the format `i ε_i Σ_i %_i`, where i is the PC index, ε\_i is its eigenvalue, Σ\_i is the running sum, and %\_i is the portion of the variance explained by all the PCs up to and including i, i.e., Σ\_i / Σ\_N (where N is the total number of PCs possible, not just the requested number). This information is plotted in `pca_variance.png`.
- `pc.{n}.dat`: the components of PC #n as `i j c_{ij}` where c\_{ij} is the (dimensionless) component of the PC. The PC is normalized (so that the squared sum of the coefficient is 1) and the sign is chosen so that the correlation of the projection with the time is always positive, i.e., residue pairs with positive components tend to move farther from each other according to the projection on this principal component in the main trajectory. This is plotted in `pca_{i}.png`.
- `pca_projections.dat`: the time-dependent projections on the principal components. The file contains n+1 columns, where the first column is time (in ps) and the next n columns are the projections on the requested principal components (in nm).

### zoom/

Analysis "zoomed in" on specific residue pairs.

- `dist_{i}_{j}.dat`: The time evolution of the distance between residues i and j and the status of the interaction between the two. The format is: `<t (ps)> <d_ij (nm)> <int_ij (0|1)>` where `int_ij` shows if the interaction is ON (1) or OFF (0). This is plotted in `dist_{i}_{j}.png`.
- `2d_correlate_{i}_{j}.dat`: The Pearson correlation between d\_ij and d\_kl, i.e., shows which residue pairs (k, l) change their distance together with i and j. This is plotted in `2d_correlate_{i}_{j}.png`.

### blocking/

"Blocking analysis" of the contact map to estimate its correlation time (see Flyvbjerg and Petersen, JCP 1989). This is done directly on the data by blocking each non-trivial time series (i.e., any residue pair that has a non-constant distance as a function of time) and normalizing the standard error as: ε\_tot=∑ε\_i^2N_res

In the limit of decorrelation of all time series, there will be a plateau in ε\_tot. ConAn estimates the correlation time as: T\_corr≈Δt⋅(ε\_platε\_0)^2

- `block_out.txt`: `i N_blocks^i ε^i ε(ε^i) T_corr`

where `i` is the number of blocking steps, `ε^i` is the standard error at the so-far blocked data, `ε(ε^i)` is the standard error on the standard error, and `T_corr` is the estimated correlation time.

The user needs to visualize and assess whether or not they see a plateau in the standard error. In case there is no such plateau, the correlation time computed by CONAN is a lower limit.

Hint: note that getting the autocorrelation time from PC projections is another way of defining this quantity (Ernst et al, JCP 2015).

## Common recipes

### Binary analysis

If you prefer a binary analysis, with contacts defined as any two residues within r\_0, just set `TRUNC <2*r0>` and `NLEVEL 2`. The factor 2 is necessary as the XPM file will be rounded. For example, if you consider inter-residue interactions to be anything that is higher than 0.5 nm, you can write:

```
TRUNC 1.0
NLEVEL 0.5
```

Any two residues within 0.5 nm will appear as d=0 nm and any two residues farther than this distance will appear as d=1 nm.

### Analyzing part of a trajectory

If you think you identified a part of a trajectory that you would like to concentrate on, you can re-run CONAN with changed `BEGIN` and `END` points and `RUN_MDMAT no`.

### Uniting several trajectories

`.xpm` files can be concatenated without issues, although time stamps could be duplicated. You can go around this with `PATCH_TIME yes`, which will make time stamps consecutive.

### Smaller than picosecond time steps

CONAN assumes the time steps between frames to be an integer multiple of 1 ps. If your data is more finely spaced than that, you can just use `PATCH TIME yes` and set `DT 1` which will unite all your frames with a fictional 1 ps time step.

### Comparing several parts of the same trajectory

If you want a detailed comparison between parts of a trajectory, one option is to run a whole CONAN first, then changing `BEGIN` and `END` in the input file as well as `RUN_MDMAT no`, keeping the outputs in separate subfolders, and then using `conan_comp.py` for a comparison.

### Using metrics other than smallest distance

One option for other metrics can be to just use a special group, for example, repeating the analysis for just the Ca atoms or an atom that you consider to represent the amino acid's/lipid's/etc position.

Otherwise, if you have another way of measuring distances and have plain-text outputs, you can transform them into the CONAN matrix format (check them out) and use the `REREAD` keyword.

### Coarse-grained input

CONAN can handle any input, as long as atoms in the same residue have the same residue ID (in fact, by changing the residue IDs in your reference file, you can trick CONAN into inter-domain contacts or any other criteria). If you need the `INPUT_PDB` file to be read in for the purposes of identifying physical interactions, be sure that the backbone atom is called either `BB` (as Martini does) or `CA` in the PDB file.

### PCA or cluster analysis on two populations

Just unite the two `.xpm` files (as mentioned above, they can be concatenated and they can still be read), use the `PATCH_TIME yes` command, run the PCA or cluster analysis, and later separate the output accordingly.
