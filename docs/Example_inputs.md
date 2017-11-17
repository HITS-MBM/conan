
# Examples for typical use cases
In the following, we will assume that you trajectory is `traj.xtc` and your structure file is `traj.tpr`, and your gnuplot scripts are in `~/conan/gnuplot_inputscripts`.
## Absolute minimum
The absolute minimum input file that CONAN can still run is:

    TRAJ traj.xtc
    COORD traj.tpr

This will use all of Gromacs standard values, for example, only `40` levels of distances, and a main truncation of 1.5 nm. No plots will be generated, only a series of plain-text matrices.
## Recommended first try
    TRAJ traj.xtc
    COORD traj.tpr
    NLEVEL 1001
    DT 1000
    ECONOMY yes
    TRUNC 1.0
    TRUNC_INTER 0.5
    GNUS_PATH ~/conan/gnuplot_inputscripts
CONAN will generate contact maps with a stride of 1 ns, a main cutoff of 1.0 nm, and defining interactions at 0.5 nm. When CONAN calls `gmx mdmat`, choose group 2 = Protein-H. We are running in "economy mode", i.e., the memory usage will be kept to a minimum.

The output will be a series of matrices (`matrices/*.dat`, plain text format), frames (`frames/*.png` ), and a movie (`movies/CONAN.mp4`). Furthermore, there will be a series of outputs in the `aggregate/` folder. See the manual for their description.

Consider adding the options `DOMAINS` (defining domains) and `DR_MODE` (making differential contact maps).
## "Which contacts form contacts together?"
This question could be answered in 6 different ways.

 1. If we mean "which contacts form together with subset X?" I.e., we have _a priori_ knowledge of some interesting residue pairs, we should give them as a "zoom in list" and CONAN will correlate contact formation/loss with respect to each pair. Keyword: `ZOOM`
 2. If we mean "which residues tend to gain/lose contacts together (cooperatively) and which residues tend to gain contacts when others lose them (competitively)?", check for inter-residue cross-correlation. Keyword: `PEARSON_INTER` .
 3. If we mean "which pairs of residues gained/lost contact *around time t*?", a good idea is to define stages of the residue. Keyword: `STAGES`.
 4. If we mean "which contacts have high covariance between them?", the answer could be principal component analysis. The first principal component can give us an idea of contacts forming/breaking together. Keyword: `PCA_MAKE`.
 5. Hierarchical clustering can help us identify subpopulations in which the contact maps are similar. Comparing clusters can then help us understand what happened. Keyword: `K_TRAJ_CLUSTERS`.
 6. Hierarchical clustering of *residues* can find subdomains of the protein (experimental). Keyword: `K_RES_CLUSTERS`.
 
Let's say we are not quite sure, so we will just try "everything." For points 1 and 3, we will need a separate file. For this illustration, let's suppose that pairs 1-10 and 5-25 are somehow special (perhaps they are evolutionarily conserved, or breaking them opens a pocket, or ...). For point 3, suppose we know that *something* happened around 5-6 ns in our trajectory (maybe we saw it in the trajectory, or maybe some interesting observable changed). The input file will be (assuming the `dmf.xpm` is still in the same folder):

    RUN_MDMAT no
    TRUNC_INTER 0.5
    GNUS_PATH ~/conan/gnuplot_inputscripts
    ZOOM_LIST zoom_list.txt
    PEARSON_INTER yes
    STAGES stage.txt
    PCA_MAKE 10
    K_TRAJ_CLUSTERS 1-5
    K_RES_CLUSTERS  1-5

Note that `NLEVEL` and `TRUNC` will be automatically detected. We chose 1-5 clusters for simplicity. The number of clusters can also be detected 
The "zoom list file" `zoom_list.txt` will be given as:

    #comments explaining why we chose them (this line is optional)
    #this line is also optional
    1 10
    5 25
   Check out the files `zoom/2d_correlate_1_10.png` and `zoom/2d_correlate_5_25.png` to see which pairs correlate with each of our "interesting" pairs.
 
   And the "stage list file" stage.txt will be:
   
    #optional line again... note that times are given in ps.
    5000 6000 "transition between A and B"
   Check out the file `frames/stage_0001_dr.png` for the change between 5 and 6 ns.
   
## Statistics!

CONAN can perform a few types of statistical analysis for you. In this example, we will correlate the evolution of the contact map with time and with observables.

    RUN_MDMAT no
    TRUNC_INTER 0.5
    GNUS_PATH ~/conan/gnuplot_inputscripts
    PEARSON_INTER yes
    PEARSON_TIME yes
    PEARSON_OBS obs.txt

Suppose we want to correlate inter-residue distances with the overall end-to-end distance. The observable file `obs.txt` will look like:

    1
    R_e
	   0 10.50
    1000 10.75
    2000  9.75
    3000  9.90
    <etc etc etc>
The inter-residue cross-correlation will be plotted in `aggregate/pearson_inter.png`, the correlation with time (showing drifts) will be in `aggregate/pearson_time.png`, and the correlation with end-to-end distance in `aggregate/pearson_obs_1.png` (in fact, pearson_data_obs_1.dat also shows the slope and the intercept).

## Uniting trajectories as an ensemble

There is no need to concatenate trajectories with `gmx trjcat`. You can unite the `dmf.xpm` files from subtrajectories (in bash, through for instance `for i in {1..10}; do cat run_${i}/dmf.xpm >> dmf.xpm; done`)  and analyze them. A typical example is:

    RUN_MDMAT no
    TRUNC_INTER 0.5
    GNUS_PATH ~/conan/gnuplot_inputscripts
    DT 1000
    PATCH_TIME yes
    OBS_IGNORE_TIME yes
    PEARSON_OBS obs.txt

Two keywords are interesting here:

 - `PATCH_TIME yes`: the time will increase in time monotonously in steps of `DT` (in this case, 1 ns). This eliminates the problem of duplicate times. Of course, the "time" variable will be somewhat abstract.
 - `IGNORE_OBS_TIME yes`: correlating with the observable will be done frame-by-frame. The observable file still needs a time column, but it will be ignored. This means that one can concatenate several observable files (except the header) and duplicate times will not be an issue.