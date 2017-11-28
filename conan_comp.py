#!/usr/bin/python3
# Import a bunch of stuff

import sys
import os
import math
import gc
import numpy as np
from scipy.stats import linregress
from scipy.stats import pearsonr
from scipy import var

# Define directory for source files

def write_matrix(outputf, time, nterm, nres, trunc, dict_frame, dr):
    if (not dr):
        outputf.write ("# TIME %d ps NTERM %d NRES %d\n"%(time, nterm, nres))
        outputf.write ("# i  j   r_ij\n")
        outputf.write ("# r_ij is the shortest distance between any atom of residues  i, j.\n")
        defaultr = trunc
    else:
        outputf.write ("# TIME %d ps NTERM %d NRES %d\n"%(time, nterm, nres))
        outputf.write ("# i  j   dr_ij\n")
        outputf.write ("# dr_ij is the change (compared to the previous frame)\n")
        outputf.write ("# in shortest distance between any atom of residues  i, j.\n")
        defaultr = 0.0
    
    for i in range(nterm, nterm+nres):
        for j in range(nterm, nterm+nres):
            pair = (i, j)
            r = dict_frame.get(pair, defaultr)
            outputf.write("%4d %4d %9.4f\n"%(pair+(r,)) )
        outputf.write("\n")
            
def combine_2d_files(dict_a, dict_b, pairs_filter, nterm, nres, name):
    file_triangle = open(name+".tri.dat", "w")
    file_difference = open(name+".diff.dat", "w")
    for n1 in range(nterm, nres+nterm):
        for n2 in range(nterm, nres+nterm):
            val_a = dict_a[(n1, n2)]
            val_b = dict_b[(n1, n2)]
            nval = len(dict_a[(n1, n2)])
            if (n1 > n2):
                val = val_a
            else:
                val = val_b
            diff = [val_b[i] - val_a[i] for i in range(nval)]
            file_triangle.write(("%5i %5i"+nval*"%15.6f"+"\n") %((n1, n2)+tuple(val)) )
            if ((n1, n2) in pairs_filter):
                file_difference.write(("%5i %5i"+nval*" %15.6f"+"\n") %((n1, n2)+tuple(diff)) )
            else:
                file_difference.write(("%5i %5i"+nval*" %15.6f"+"\n") %((n1, n2)+nval*(0.0, ) ) )
        file_triangle.write("\n")
        file_difference.write("\n")
    file_triangle.close()
    file_difference.close()

def get_2d(file):
    dict_file = dict()
    inputf = open(file)
    line = inputf.readline()
    while (line.startswith('#')):
        line = inputf.readline() # skip all comments...
    n1 = int(line.split()[0])
    n2 = int(line.split()[1])
    vals = tuple ([ float(x) for x in line.split()[2:] ])
    dict_file[(n1, n2)] = vals
    for line in inputf:
        if (len(line.split()) > 0):
            n1 = int(line.split()[0])
            n2 = int(line.split()[1])
            vals = tuple ([ float(x) for x in line.split()[2:] ])
            dict_file[(n1, n2)] = vals
    return dict_file

def get_1d(file):
    dict_file = dict()
    inputf = open(file)
    line = inputf.readline()
    while (line.startswith('#')):
        line = inputf.readline() # skip all comments...
    for line in inputf:
        if (len(line.split()) > 0):
            n = int(line.split()[0])
            vals = tuple ([ float(x) for x in line.split()[1:] ])
            dict_file[n] = vals
    return dict_file

# In[ ]:

# Function to check the input file for mdmat and run mdmat

if len(sys.argv)!=2:
    print('Usage: You need to provide a single text file in which specify all the options for creating a matrix!')
finput=open(sys.argv[1])

gnus_path = ''

r_cut = float('Inf')
life_cut = 0.0
domains = 0
title_a = "Run A"
title_b = "Run B"

for line in finput:
    if line.startswith('#') or len(line.split())<2:
        continue
    else:
        key = line.split()[0].upper()
        val = line.split()[1]
        if key == 'RUN_A':
            run_a = str(val)
        elif key == 'RUN_B':
            run_b = str(val)
        elif key == 'TITLE_A':
            title_a = line[line.find('"')+1:line.rfind('"')]
        elif key == 'TITLE_B':
            title_b = line[line.find('"')+1:line.rfind('"')]
        elif key == 'R_CUT':
            r_cut = float(val)
        elif key == 'LIFE_CUT':
            life_cut = float(val)
        elif key == 'GNUS_PATH':
            gnus_path = str(val)
        else:
            print("Cannot understand the keyword %s given in the input file." % key)

# Check if trajectories and/or coordinates are present in the input file..if not stop!

# check if we have the 2x three files in the two directories:

for i in [run_a, run_b]:
    for j in ['mdmat_average_rmsf.dat', 'timeline.dat']:
        file_name = i+'/aggregate/'+j
        if (not os.path.exists(file_name)):
            print("Required file does not exist!", file_name)
            sys.exit()

print("Found all 6 files!")

ok = True
dict_avg_rmsf_a = get_2d(run_a+'/aggregate/mdmat_average_rmsf.dat')
dict_avg_rmsf_b = get_2d(run_b+'/aggregate/mdmat_average_rmsf.dat')

mina = min([x[0] for x in dict_avg_rmsf_a.keys()])
maxa = max([x[0] for x in dict_avg_rmsf_a.keys()])
minb = min([x[0] for x in dict_avg_rmsf_b.keys()])
maxb = max([x[0] for x in dict_avg_rmsf_b.keys()])
if (mina == minb and maxa == maxb):
    nterm = mina
    nres = maxa - nterm + 1
    print("Found residues between %i and %i."%(mina,maxa))
else:
    print("Residue number mismatch. Good luck next time.")
    sys.exit()
pairs_filter = set([])
dict_avg_rmsf_a = get_2d(run_a+'/aggregate/mdmat_average_rmsf.dat')
dict_avg_rmsf_b = get_2d(run_b+'/aggregate/mdmat_average_rmsf.dat')
pairs_filter_life = set([])
dict_timeline_a = get_2d(run_a+'/aggregate/timeline.dat')
#for pair in dict_avg_rmsd_a:
#    if (dict_avg_rmsd_a[pair][1] <= r_cut and dict_timeline_a[pair][2] >= life_cut):
#        pairs_filter = pairs_filter | set([pair])
dict_timeline_b = get_2d(run_b+'/aggregate/timeline.dat')
#for pair in dict_avg_rmsd_a:
#    if (dict_avg_rmsd_b[pair][1] <= r_cut and dict_timeline_b[pair][2] >= life_cut):
#        pairs_filter = pairs_filter | set([pair])
pairs_filter = dict_avg_rmsf_a.keys()
maxt_a = max([x[0] for x in dict_timeline_a.values()])
maxt_b = max([x[0] for x in dict_timeline_b.values()])
maxr_a = max([x[0] for x in dict_avg_rmsf_a.values()])
maxr_b = max([x[0] for x in dict_avg_rmsf_b.values()])
maxs_a = max([x[1] for x in dict_avg_rmsf_a.values()])
maxs_b = max([x[1] for x in dict_avg_rmsf_b.values()])

maxt = max(maxt_a, maxt_b)
print("Maximum time:",maxt)
maxr = max(maxr_a, maxr_b)
print("Maximum distance:",maxr)
maxs = max(maxs_a, maxs_b)
print("Maximum RMSF:",maxs)

print("Files have been read in!")

combine_2d_files(dict_avg_rmsf_a, dict_avg_rmsf_b, pairs_filter, nterm, nres, 'average_rmsf')
combine_2d_files(dict_timeline_a, dict_timeline_b, pairs_filter, nterm, nres, 'timeline')

for i in [run_a, run_b]:
    file_name = i+'/domains.gnu'
    if (os.path.exists(file_name) and domains == 0):
        domains = 1
        os.system("cp " + file_name + " .")
        
print("Calling gnuplot...")
os.system('gnuplot -e "domains=%d;minres=%d;maxres=%d;maxt=%f;maxr=%f;maxs=%f;'%(domains,nterm,nterm-1+nres,maxt/1000,maxr,maxs)+"title_A='%s';title_B='%s'"%(title_a,title_b)+'" %s/compare.gnu'%(gnus_path))
print('gnuplot -e "domains=%d;minres=%d;maxres=%d;maxt=%f;maxr=%f;maxs=%f;'%(domains,nterm,nterm-1+nres,maxt/1000,maxr,maxs)+"title_A='%s';title_B='%s'"%(title_a,title_b)+'" %s/compare.gnu'%(gnus_path))
