#!/usr/bin/env python3
# Import a bunch of stuff
import time
import sys
import os
import math
import gc
import numpy as np
import matplotlib.pyplot as plt 
import itertools
from scipy.stats import pearsonr
from scipy.stats import sem
from scipy.stats import linregress
from scipy.spatial.distance import pdist
from scipy import var
from scipy.sparse.linalg import eigsh
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster

def weight_contacts(r, trunc):
# this function is used for the "weighted number of contacts" for the cross-correlation.
    return ((trunc - r) / trunc)
# 0.0 when r>rcut... 

def update_interaction(r, trunc_inter, trunc_inter_high, status, inter_mode):
# This function is used for the update of the status of an interaction
# Two modes are implemented at the moment:
# 0 = "history" - status can only be 0 or 1. when the distance is in the buffer, keep the last status.
# 1 = "linear" - status can be fractional, switching between 0 and 1 linearly when in the buffer.
    if (trunc_inter == trunc_inter_high):
        new_status = int(r < trunc_inter) # If the buffer doesn't exist, the old status is ignored.
    elif (inter_mode == 0):
        if (r < trunc_inter):
            new_status = 1
        elif (r > trunc_inter_high):
            new_status = 0
        else:
            new_status = status
    elif (inter_mode == 1):
        if (r < trunc_inter):
            new_status = 1
        elif (r > trunc_inter_high):
            new_status = 0
        else:
            new_status = (r - trunc_inter)/(trunc_inter_high - trunc_inter)
    return new_status
# new_status can be an integer (1 or 0) but also a real number between (0, 1).
# Add new modes if needed.
# At the moment, inter_mode = 0 is the only one accessible
# ("# of encounters" and "avg encounters" would be needlessly complicated).


def read_in_observables(inputf, time_vec):
# this function will read in the observables from the given file. It gets as input the (open) file
# as well as the available time frames. It returns the time vector (of frames that have both time information
# and observables) and several vectors of observables.
    time_vec_out = []
    frame_indices = []
    obs_vec = []
    titles  = []
    units = []
    ntitles = -1
    for line in inputf:
        if ( not line.startswith("#") and (len(line.split()) > 0) ):
            if (ntitles == -1):
                ntitles = int(line.split()[0])
            elif (len(titles)<ntitles):
                title = line
                titles.append(title)
            else:
                t = int(float(line.split()[0]))
                if (t in time_vec):
                    ind = time_vec.index(t)
                    time_vec_out.append(t)
                    frame_indices.append(ind)
                    if (not obs_vec):
                        obs_vec = [ [] for val in range(ntitles)]
                    for n in range (0, ntitles):
                        obs_vec[n].append(float(line.split()[n+1]))
    return frame_indices, time_vec_out, titles, obs_vec

def cond_index(i, j, n):
# the index in the condensed distance matrix of an element (i,j) from the upper triangle (i<j).
    return int(n*(n-1)/2 - (n-i)*(n-i-1)/2 + j - i - 1)

def remove_shadows(dict_frame, nterm, nres, tol):
# assuming a dictionary of a frame, we remove any contact between residue pair (i,j) where there is a
# "shortcut" through a third residue k; this is defined as r_ij > tol* (r_ik + r_jk)
    dict_copy = dict(dict_frame)
    for pair in dict_frame:
        i = pair[0]
        j = pair[1]
        d0 = dict_frame[pair]
        for k in range(nterm, nterm+nres):
            pair1 = tuple(sorted((i,k)))
            pair2 = tuple(sorted((j,k)))
            if (pair1 in dict_frame and pair2 in dict_frame):
                if ((dict_frame[pair1] + dict_frame[pair2])*tol < d0 ):
                    dict_copy.pop((i, j))
                    break
    return dict_copy

def read_in_observables_blind(inputf):
# this function ignores the time column from the observable file and just returns the vectors.
    obs_vec = []
    titles  = []
    units   = []
    ntitles = -1
    for line in inputf:
        if ( not line.startswith("#") and (len(line.split()) > 0) ):
            if (ntitles == -1):
                ntitles = int(line.split()[0])
            elif (len(titles)<ntitles):
                quotes = [i for i,j in enumerate(line) if j=='"']
                title = line[quotes[0]+1:quotes[1]]
                unit = line[quotes[2]+1:quotes[3]]
                titles.append(title)
                units.append(unit)
            else:
                if (not obs_vec):
                    obs_vec = [ [] for val in range(ntitles)]
                for n in range (0, ntitles):
                    obs_vec[n].append(float(line.split()[n+1]))
    return titles, units, obs_vec

def which(program):
# This function checks for the existence of an executable
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return ''


# This function will generate the necessary gnuplot preamble that will label the domains in the axes.
def change_bf_in_pdb (inputpdb, outputpdb, bf_dict):
    for line in inputpdb:
        lineoutput = line
        if (len(line.split()) > 5):
            if (line.split()[0]=='ATOM'):
                resid = int(line[22:26])
                bf = bf_dict.get(resid,0.0)
                lineoutput = line[:60] + "%6.2f"%bf + line[66:]
        outputpdb.write(lineoutput)

def interaction(str):
# This function identifies a possible interaction between two resiudes, assuming this happens
# through the side chains (but will ignore hydrophobic packing of lysine and other exotic effects).
    negative = ['E', 'D']
    positive = ['R', 'K']
    charged  = negative + positive
    hydrophobic = ['A', 'F', 'I', 'L', 'M', 'P', 'V', 'W', 'Y']
    donor    = ['C', 'H', 'N', 'Q', 'S', 'T', 'Y', 'R', 'K']
    acceptor = ['D', 'E', 'N', 'Q', 'S', 'T', 'Y']
    r1 = str[0]
    r2 = str[1]
    if (r1 in hydrophobic and r2 in hydrophobic):
        return 1 #hydrophobic
    elif (r1 in negative and r2 in positive) or (r1 in positive and r2 in negative):
        return 3 #salt bridge 
    elif ((r1 in donor and r2 in acceptor) or (r1 in acceptor and r2 in donor)):
        return 2 # hydrogen bond
    else:
        return 0 #nothing...
        
def read_zoom_list(zoomfile):
# This function reads in the "zoom list" with particularly interesting residue pairs.
    zoom_list = []
    for line in zoomfile:
        if ( not line.startswith("#") and (len(line.split()) > 0) ):
            i1          = int(line.split()[0])
            i2          = int(line.split()[1])
            pair = tuple(sorted( (i1, i2)))
            zoom_list.append(pair)
    return zoom_list

def zoom_on_pair(pair, rvec, time_vec, trunc, inter_low, inter_high, pairs_list, pairs_legend, vectors):
# This function executes the "zooming in" on one of the residue pairs on the list.
    inter = 0
    dfile = open("zoom/dist_%i_%i.dat"%pair, "w")
    pearson2dfile = open("zoom/2d_correlate_%i_%i.dat"%pair, "w")
    dfile.write("# t (ps) r (nm) interact?\n")
    for i in range(len(rvec)):
        r = rvec[i]
        if (r < inter_low):
            inter = 1
        elif (r > inter_high):
            inter = 0
        dfile.write("%9d %12.6f %3i\n"%(time_vec[i], rvec[i], inter))
    dfile.close()
    pearson2d_dict = dict()
    for i in range(len(pairs_list)):
        pair2 = pairs_list[i]
        rvec2 = vectors[i]
        rvalue, pvalue = pearsonr(rvec, rvec2)
        pearson2d_dict[pair2] = rvalue 
        pearson2d_dict[pair2[::-1]] = rvalue 
    for i in range(xterm, xterm+xres):
        for j in range(yterm, yterm+yres):
            pair2 = (i,j)
            x = pearson2d_dict.get(pair2, 0.0)
            pearson2dfile.write("%5i %5i %12.6f\n"%(i, j, x))
        pearson2dfile.write("\n")
    pearson2dfile.close()
    if (gnus_path): os.system("""gnuplot -e "maxz=1;domains=%i;inputfile='zoom/2d_correlate_%i_%i.dat' ; outputfile='zoom/2d_correlate_%i_%i.png'; title_time = 'Correlation with the contact distance between residues (%i, %i)'" %s/script_corr.gnu """%
                ((domains,) + 3*pair + (gnus_path,) ))
    dfile.close()
    if (gnus_path): os.system("""gnuplot -e "inter_high=%f;inter_low=%f;inputfile='zoom/dist_%i_%i.dat'; outputfile='zoom/dist_%i_%i.png'; title_plot='Time dependence of distance between residues (%i,%i)'"  %s/1d_zoom.gnu """%((inter_high,inter_low)+3*pair + (gnus_path,)))

def aggregate(rvec, trunc, inter_low, inter_high):
# This function "aggregates" the vector of inter-residue distance of one particular pair.
    dwells = []
    n = 0
    first = -1
    last = -1
    i = 0
    for r in rvec:
        if (n==0 and r<inter_low):
            n = 1
            if (len(dwells) == 0):
                first = i
        elif (n > 0):
            if (r > inter_high):
                dwells.append(n)
                n = 0
                last = i
            else:
                n += 1
                last = i
        i = i + 1
    if (n > 0):
        dwells.append(n)
        last = i-1
    if (len(dwells) == 0):
        last = 0
        first = i - 1
    n_tot = sum(dwells)
    avg = np.mean(rvec)
    std = np.std(rvec)
    return avg, std, first, last, n_tot, len(dwells)

def read_sequence(inputpdb):
#This is a function reading in all the residues from a PDB file (assuming the name is the fourth column).
    residue_dict = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','GLH':'E','PHE':'F','GLY':'G','HIS':'H','HIE':'H','HID':'H','HIP':'H','ILE':'I','LYS':'K','LYN':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'}
    list = []
    for line in inputpdb:
        if (len(line.split()) > 5):
            if (line.split()[0]=='ATOM' and line.split()[2] in ['CA', 'BB']):
                resid = int(line[22:26])
                resname = residue_dict.get(line.split()[3],'X')
                list.append((resid,resname))
    interact_pair_dict = dict()
    i = 0
    for r1 in list:
        for r2 in list:
            i = i+1
            s = r1[1]+r2[1]
            pair = (r1[0], r2[0])
            if (r1[0] != r2[0]):
                inter = interaction(s)
            else:
                inter = 0
            
            if (inter>0):
                interact_pair_dict[pair] = inter
    return interact_pair_dict

def read_sequence_asymm(inputpdb_x, inputpdb_y, xterm, xres, yterm, yres):
#This is an analogue to the function above but using two PDBs and identifying those interactions.
    residue_dict = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'}
    list_x = []
    for line in inputpdb_x:
        if (len(line.split()) > 5):
            if (line.split()[0]=='ATOM' and line.split()[2] in ['CA', 'BB']):
                resid = int(line[22:26])
                resname = residue_dict.get(line.split()[3],'X')
                list_x.append((resid,resname))
    list_y = []
    for line in inputpdb_y:
        if (len(line.split()) > 5):
            if (line.split()[0]=='ATOM' and line.split()[2] in ['CA', 'BB']):
                resid = int(line[22:26])
                resname = residue_dict.get(line.split()[3],'X')
                list_y.append((resid,resname))
    interact_pair_dict = dict()
    i = 0
    for r1 in list_x:
        for r2 in list_y:
            if (r1[0] in range(xterm, xterm+xres) and r2[0] in range(yterm, yterm+yres)):
                i = i+1
                s = r1[1]+r2[1]
                pair = (r1[0], r2[0])
                inter = interaction(s)
                if (inter>0):
                    interact_pair_dict[pair] = inter
    return interact_pair_dict

def read_stages(stagef):
#This is a function reading in the stage file...
    stages_list = []
    for line in stagef:
        if ( not line.startswith("#") and (len(line.split()) > 0) ):
            t1 = int(line.split()[0])
            t2 = int(line.split()[1])
            if (len(line.split()) > 2):
                stage_name = line[line.find('"')+1:line.rfind('"')]
                stages_list.append((t1, t2, stage_name))
            else:
                stages_list.append( (t1, t2, "Stage %d"%(len(stages_list)+1) ) )
    return stages_list

def create_tics(domainf, nterm, nres):
#This function creates the tics for the domains file for gnuplot.
    gnuf_1d = open("domains_1d.gnu","w")
    gnuf_2d = open("domains.gnu","w")
    rmin = nterm
    rmax = nterm+nres-1

    mtics = []
    tics  = []

    for line in domainf:
        if ( not line.startswith("#") and (len(line.split()) > 0) ):
            i1          = int(line.split()[0])
            i2          = int(line.split()[1])
            domain_name = line[line.find('"')+1:line.rfind('"')]
            tics.append((domain_name, 0.5 * (i1 + i2)))
            mtics       = mtics + [i1 - 0.5, i2 + 0.5]
    mtics_copy = set(mtics)
    mtics = set(mtics)

    for i in mtics_copy:
        if (i < rmin or i > rmax):
            mtics.remove(i)

    tics_str  = ", ".join(['"%s" %.1f'%i for i in tics])
    xtics_str = "set xtics (" + tics_str + ")\n"
    ytics_str = "set ytics (" + tics_str + ")\n"

    gnuf_2d.write(xtics_str)
    gnuf_2d.write(ytics_str)
    gnuf_1d.write(xtics_str)

    for i in mtics:
        gnuf_2d.write("set xtics add (%.1f 1) \n"%(i))
        gnuf_2d.write("set ytics add (%.1f 1) \n"%(i))
        gnuf_1d.write("set xtics add (%.1f 1) \n"%(i))
    gnuf_2d.close()
    gnuf_1d.close()

def create_tics_asymm(domainf_x, domainf_y, xterm, xres, yterm, yres):
    gnuf_1dx = open("domains_1d.x.gnu","w")
    gnuf_1dy = open("domains_1d.y.gnu","w")
    gnuf_2d = open("domains.gnu","w")
    rminx = xterm
    rmaxx = xterm+xres-1
    rminy = yterm
    rmaxy = yterm+yres-1

    mxtics = []
    xtics  = []
    mytics = []
    ytics  = []

    for line in domainf_x:
        if ( not line.startswith("#") and (len(line.split()) > 0) ):
            i1          = int(line.split()[0])
            i2          = int(line.split()[1])
            domain_name = line[line.find('"')+1:line.rfind('"')]
            xtics.append((domain_name, 0.5 * (i1 + i2)))
            mxtics       = mxtics + [i1 - 0.5, i2 + 0.5]
    for line in domainf_y:
        if ( not line.startswith("#") and (len(line.split()) > 0) ):
            i1          = int(line.split()[0])
            i2          = int(line.split()[1])
            domain_name = line[line.find('"')+1:line.rfind('"')]
            ytics.append((domain_name, 0.5 * (i1 + i2)))
            mytics      = mytics + [i1 - 0.5, i2 + 0.5]

    mxtics_copy = set(mxtics)
    mxtics = set(mxtics)

    for i in mxtics_copy:
        if (i < rminx or i > rmaxx):
            mxtics.remove(i)

    mytics_copy = set(mytics)
    mytics = set(mytics)

    for i in mytics_copy:
        if (i < rminy or i > rmaxy):
            mytics.remove(i)

    xtics_str  = ", ".join(['"%s" %.1f'%i for i in xtics])
    ytics_str  = ", ".join(['"%s" %.1f'%i for i in ytics])
    xtics_str = "set xtics (" + xtics_str + ")\n"
    ytics_str2d = "set ytics (" + ytics_str + ")\n"
    ytics_str1d = "set xtics (" + ytics_str + ")\n"

    gnuf_2d.write(xtics_str)
    gnuf_1dx.write(xtics_str)
    gnuf_2d.write(ytics_str2d)
    gnuf_1dy.write(ytics_str1d)

    for i in mxtics:
        gnuf_2d.write("set xtics add (%.1f 1) \n"%(i))
        gnuf_1dx.write("set xtics add (%.1f 1) \n"%(i))
    for i in mytics:
        gnuf_2d.write("set ytics add (%.1f 1) \n"%(i))
        gnuf_1dy.write("set xtics add (%.1f 1) \n"%(i))
    gnuf_2d.close()
    gnuf_1dx.close()
    gnuf_1dy.close()

def plot_pc (xterm, xres, yterm, yres, ii, vec, pairs, pairs_dict, gnus_path, domains):
    pc_file = open("pca/pc.%i.dat"%(ii+1), "w+")
    maxz = max(-min(vec), max(vec))
    pc_dict = dict()
    for pair in pairs:
        pc_dict[pair] = vec[pairs_legend[pair]]
        if (not asymm):
            pc_dict[pair[::-1]] = vec[pairs_legend[pair]]
    for i in range(xterm, xterm+xres):
        for j in range(yterm, yterm+yres):
            pair = i, j
            x = pc_dict.get(pair, 0.0)
            pc_file.write("%5i %5i %12.6f\n"%(i, j, x))
        pc_file.write("\n")
    pc_file.close()
    if (gnus_path): os.system("""gnuplot -e "domains=%i;i=%i;maxz=%f;input_file='pca/pc.%i.dat'" %s/script_pca.gnu """%(domains, ii+1, maxz, ii+1, gnus_path))

def project_pc(v,  time_vec, vectors, pairs_list, change_sign = True):
    file_projections = open("pca/pca_projections.dat", "w+")
    nframes = len(time_vec)
    n = np.shape(v)[1]
    projections = np.zeros((nframes, n))
    ordered_v = np.array(v)
    print("Projecting principal components...")
    for i in range(nframes):
        time = time_vec[i]
        dist_frame = [vec[i] for vec in vectors]
        projections[i] = np.array([ np.dot(v[:, j ], dist_frame) for j in range(n)])
    if (change_sign):
        for i in range(n):
            slope, intercept, r, p, std_err = linregress(projections[:, i], time_vec)
            if (r < 0):
                ordered_v[:, i] = -v[:, i]
                projections[:, i] *= -1
    for i in range(nframes):
        dist_frame = [vec[i] for vec in vectors]
        file_projections.write("%12i"%time_vec[i])
        file_projections.write((n*" %12.6f")%tuple(projections[i])+"\n")
    file_projections.close()
    return ordered_v

# This function reads the index from the first frame of the dmf.xpm file
def read_index(inputf):
    dict_legend = dict()
    line = inputf.readline()
    while (line.split()[0]!="static"):
        line = inputf.readline()
    line = inputf.readline()[1:-3]
    nres    = int(line.split()[0])
    nlevels = int(line.split()[2])
    nchar   = int(line.split()[3])
    value = 0.0
    chars   = []

    for i in range(nlevels):
        line = inputf.readline()
        char = line[1:].split()[0]
        chars.append(char)
        if (len(char) != nchar):
            print('warning! this character is of the wrong length:',char)
            print('carrying on but it is a good idea to check.')
    trunc = float(line.split()[-2][1:-1])
    dx = trunc/(nlevels-1)
    value = 0.0
    for i in range(nlevels-1):
       dict_legend[chars[i]]=value
       value = value + dx

    dict_legend[chars[nlevels-1]] = trunc #make sure the truncation threshold is the same as the highest value.
    return nres, nlevels, nchar, dict_legend, trunc

# Function to read the values and the title (time) of a frame

def read_next_time(inputf):
    line = inputf.readline()
    if not line:
       return float('-Inf')
    get_out = False
    while (not get_out):
        line = inputf.readline()
        if (not line):
            return float('-Inf')
            break
        get_out = len(line.split()) > 1
        if (get_out):
            get_out = get_out and (line.split()[1] == 'title:')
    title = line[line.find('"')+1:line.rfind('"')]
    time  = float(title[2:].split()[0])
    return time

def read_frame(inputf, nres, nlevels, nchar, dict_legend, nterm, trunc, asymm_data, dimer):
    line = inputf.readline()
    dict_frame = dict()
    
    if (asymm_data):
        start_x, start_y, xres, yres, xterm, yterm = asymm_data

    #skip all the stuff with /* comments...
    
    while (line.split()[0]!="static"):
        line = inputf.readline()
    
    #skip the legend - we know already...
  
    for i in range(nlevels+1):
        line = inputf.readline()

    #skip the axes...
    
    line = inputf.readline()
    
    while (not line.split()[0].startswith('"')):
        line = inputf.readline()
    n = nres-1
    
    if (not asymm_data):
        for i in range(nres):
            if (i < nres-1):
                linestr = line[1:-3]
            else:
                linestr = line[1:-2]
            z = [dict_legend[linestr[k:k+nchar]] for k in range(0, n * nchar, nchar)]
            for m in range(n):
                if (z[m] < trunc):
                    dict_frame[(m + nterm, n + nterm)] = z[m]
            line = inputf.readline()
            n=n-1
    else:
        for i in range(nres):
            if (i < nres-1):
                linestr = line[1:-3]
            else:
                linestr = line[1:-2]
            z = [dict_legend[linestr[k:k+nchar]] for k in range(0, nres * nchar, nchar)]
            for m in range(nres):
                if (z[m] < trunc and ( (m+1) in range(start_x, start_x+xres)) and ( (n+1) in range(start_y, start_y+yres)) ):
                    dict_frame[(m + 1 + xterm - start_x, n + 1 + yterm - start_y)] = z[m]
            line = inputf.readline()
            n=n-1
    
    if (dimer):
        contacts_above = 0.0
        contacts_below = 0.0
        for i, j in dict_frame:
            if (i < j):
                contacts_above += weight_contacts(dict_frame[(i, j)], trunc)
            elif (i > j):
                contacts_below += weight_contacts(dict_frame[(i, j)], trunc)
        if (contacts_above < contacts_below):
            dict_frame_transpose = {(j, i): dict_frame[(i, j)] for i, j in dict_frame}
            return dict_frame_transpose
        else:
            return dict_frame
    else:
        return dict_frame


def write_matrix(outputf, time, xterm, xres, yterm, yres, trunc, dict_frame, dr, asymm):
#This function prints out one matrix.
# dr is a boolean variable showing whether the "frame" is really a difference of frames.
    if (not dr):
        outputf.write ("# i  j   r_ij\n")
        outputf.write ("# r_ij is the shortest distance between any atom of residues  i, j.\n")
        defaultr = trunc
    else:
        outputf.write ("# i  j   dr_ij\n")
        outputf.write ("# dr_ij is the change (compared to the previous frame)\n")
        outputf.write ("# in shortest distance between any atom of residues  i, j.\n")
        defaultr = 0.0
    
    for i in range(xterm, xterm+xres):
        for j in range(yterm, yterm+yres):
            if (not asymm):
                if (i == j):
                    r = 0
                else:
                    if (i < j):
                        pair = (i, j)
                    else:
                        pair = (j, i)
                    r = dict_frame.get(pair, defaultr)
            else:
                r = dict_frame.get((i, j), defaultr)
            outputf.write("%4d %4d %9.4f\n"%(i, j, r) )
        outputf.write("\n")
            
def update_aggregate(r, time, trunc, trunc_inter, trunc_inter_high, agg, inter_mode = 0):
    first, last, n_int, n_enc, status, sumr, sumr2 = tuple(agg)
    sumr += r
    sumr2 += r*r
    new_status = update_interaction(r, trunc_inter, trunc_inter_high, status, inter_mode) # at the moment, "history" is the only one accessible.
    if (new_status):
        n_int += 1
        last = time
        if (not status): #switch from OFF to ON
            if (first == -1):
                first = time
            n_enc += 1
    return np.array([first, last, n_int, n_enc, new_status, sumr, sumr2])

def prepare_dimer_transpose(pairs_list, pairs_legend):
    npad = 0
    n_len = len(pairs_list)
    new_legend = pairs_legend
    map = []
    pad_map = []
    for i, pair in enumerate(pairs_list):
        pair_transpose = pair[::-1]
        if pair_transpose in pairs_list:
            j = pairs_legend[pair_transpose]
        else:
            j = n_len + npad
            pad_map.append(i)
            npad += 1
        map.append(j)
    map = map + pad_map
    return npad, map
    
def plot_frames(dict_frame, first_frame, prev_frame, pairs_list, dr_mode, xterm, xres, yterm, yres, trunc, matrices, clean_matrices, n, time, rmsd_perframe_file, asymm, dimer):
    if (matrices):
        outputf = open("matrices/%05d.dat"%n,"w")
        write_matrix(outputf, time, xterm, xres, yterm, yres, trunc, dict_frame, False, asymm)
        outputf.close()
        if (clean_matrices):
            if (gnus_path): os.system("""gnuplot -e "domains=%d;maxz=%f;label_str='Residue index';inputfile='matrices/%05d.dat';outputfile='frames/%05d.png';"""%(domains,trunc,n,n)+
                      """title_str='t = %7.3f (ns)';cblabel_str = 'Distance (nm)" %s/script_single.gnu"""%(time * 0.001, gnus_path))
            os.system("""rm matrices/%05d.dat"""%n)
    if (first_frame):
        diff_dict_0 = {pair: dict_frame.get(pair, trunc) - first_frame.get(pair, trunc) for pair in (set(dict_frame.keys()) | set(first_frame.keys()))}
        diff_dict_p = {pair: dict_frame.get(pair, trunc) -  prev_frame.get(pair, trunc) for pair in (set(dict_frame.keys()) | set( prev_frame.keys()))}
        rmsd_perframe_0 = rmsd_frame(dict_frame, first_frame, trunc, xres, yres, asymm, dimer)
        rmsd_perframe_p = rmsd_frame(dict_frame, prev_frame,  trunc, xres, yres, asymm, dimer)
        rmsd_perframe_file.write("%9.4f %9.4f %9.4f\n"%(0.001*time, rmsd_perframe_p, rmsd_perframe_0))
        if (matrices):
            if (dr_mode in [1, 3]):
                outputf_dr = open("matrices/%05d_dr.init.dat"%n,"w")
                write_matrix(outputf_dr, time, xterm, xres, yterm, yres, trunc, diff_dict_0, True, asymm)
                outputf_dr.close()
                if (clean_matrices):
                    if (gnus_path):
                        os.system('gnuplot -e "domains=0;maxz=%f;inputfile='%(trunc_dr)+"'matrices/%05d_dr.init.dat'"%n+";outputfile='frames/%05d_dr.init.png';title_time='Change up to t = %7.3f (ns)'"%(n,time)+'" %s/script_dr.gnu'%(gnus_path))
                        if (domains):
                            os.system('gnuplot -e "domains=1;maxz=%f;inputfile='%(trunc_dr)+"'matrices/%05d_dr.init.dat'"%n+";outputfile='frames/%05d_dr.init.domains.png';title_time='Change up to t = %7.3f (ns)'"%(n,time)+'" %s/script_dr.gnu'%(gnus_path))
                    os.system("rm matrices/%05d_dr.init.dat"%n)
            if (dr_mode in [2, 3]):
                outputf_dr = open("matrices/%05d_dr.prev.dat"%n,"w")
                write_matrix(outputf_dr, time, xterm, xres, yterm, yres, trunc, diff_dict_p, True, asymm)
                outputf_dr.close()
                if (clean_matrices):
                    if (gnus_path):
                        os.system('gnuplot -e "domains=0;maxz=%f;inputfile='%(trunc_dr)+"'matrices/%05d_dr.prev.dat'"%n+";outputfile='frames/%05d_dr.prev.png';title_time='Change at t = %7.3f (ns)'"%(n,time)+'" %s/script_dr.gnu'%(gnus_path))
                        if (domains):
                            os.system('gnuplot -e "domains=1;maxz=%f;inputfile='%(trunc_dr)+"'matrices/%05d_dr.prev.dat'"%n+";outputfile='frames/%05d_dr.prev.domains.png';title_time='Change at t = %7.3f (ns)'"%(n,time)+'" %s/script_dr.gnu'%(gnus_path))
                    os.system("rm matrices/%05d_dr.prev.dat"%n)
    print("done with frame corresponding to ",time,"ps",end="\r")
    
def read_process_frames(inputf, nres, nlevels, nchar, dict_legend, asymm_data, trunc, trunc_inter, trunc_inter_high, begin, end, dt, patch_time, gnus_path,
                        dr_mode, domains, pearson_inter, economy, dimer, reread, stages_list, interact_pair_dict, trunc_lifetime):
    outputf_timeline = open("aggregate/timeline.dat", "w")
    outputf_timeline.write ("#  i  j   t_first t_last t_total (# encounters)\n")
    outputf_timeline.write ("# first/last time residues i and j are within the interaction distance.\n")
    outputf_timeline.write ("# total: fraction of frames in which residues i and j are within interaction distance (maximum = 1.0).\n")
    outputf_avg_rmsf = open("aggregate/mdmat_average_rmsf.dat", "w")
    outputf_avg_rmsf.write ("#  i  j  r_ij^avg sigma_(rij)\n")
    outputf_avg_rmsf.write ("# average of the distance between residues i and j and the standard deviation (computed across frames.\n")
    native_file = open("aggregate/native_contacts.dat", "w")
    native_file.write("# time (ns)  #native #non-native\n")
    rmsd_perframe_file = open("aggregate/time_mdmat_rmsd.dat","w")
    rmsd_perframe_file.write("# time (ns)  RMSD_mdmat_previous RMSD_mdmat_total\n")
    if (stages_list):
        frames_stages = [ stage + (None, None) for stage in stages_list]
            
    asymm = bool(asymm_data)
    if (asymm):
        start_x, start_y, xres, yres, xterm, yterm = asymm_data
        interaction_map_x = [[] for x in range(xres)]
        interaction_map_y = [[] for x in range(yres)]
        total_interact_dict_x = dict()
        total_interact_dict_y = dict()
    else:
        start_x, start_y, xres, yres, xterm, yterm = (1, 1, nres, nres, nterm, nterm)
        interaction_map = [[] for x in range(nres)]
        local_interact_dict = dict()
    if (inputf):
        inputline = inputf.readline()
    prev_frame = None
    first_frame = None
    pairs_list = []
    time_vec = []
    native_list = []
    vectors = np.zeros((0, 0), dtype = np.float16)   # the main data structure. rows = pairs. columns = frames.
    emptiness = np.zeros((1, 0), dtype = np.float16) # the row that "newcomer pairs" will inherit.
    pairs_legend = dict()
    if (interact_pair_dict):
        interact_file = open("aggregate/interaction_types.dat", 'w')
    if (pearson_inter):
        residue_vectors = np.zeros([nres,0])
    aggregates = np.zeros([0, 7]) # aggregate values for all pairs
    empty_aggregate = np.array([-1, 0, 0, 0, 0, 0, 0]) #first, last, lifetime, N_encounter, status, sumr, sumr2
    n = 0
    if (reread):
        matrixfiles = []
        for i in (range(100000)):
            filename = "%s/%05i.dat"%(reread, i)
            if (os.path.exists(filename)):
                matrixfiles.append(filename)
        if (not begin):
            begin = 0
        if (not dt):
            dt = 1000
        if (matrixfiles == []):
            print("Reread folder is empty, will quit now.")
            sys.exit()
        else:
            print("Found %i matrices to reread from."%len(matrixfiles))
            time = begin
            end = float('Inf')
    else:        
        time = read_next_time(inputf)
        if (not dt):
            time_next = read_next_time(inputf)
            dt = time_next - time
            inputf.seek(0)
            inputline = inputf.readline()
            time = read_next_time(inputf)
            if (dt > float('-inf')):
                print("Detected dt = %i ps"%dt)
        time0= time
        if (not end):
            end = float('Inf')
        if (not begin):
            begin = time0
    if (asymm):
        total_interactionf_x = open("aggregate/total_interactions.x.dat","w")
        total_interactionf_x.write("%i\n"%xres)
        total_interactionf_y = open("aggregate/total_interactions.y.dat","w")
        total_interactionf_y.write("%i\n"%yres)
    else:
        local_interactionf = open("aggregate/local_interactions.dat","w")
        local_interactionf.write("%i\n"%nres)
    while (time > float('-Inf') and time <=end):
        if (time >= begin and ((time-begin) % dt == 0)):
            npairs = vectors.shape[0]
            frame_column = np.zeros([npairs, 1], dtype = np.float16)
            if (pearson_inter):
                residue_column = np.zeros([nres, 1])
            if (reread):
                inputf = open(matrixfiles[n])
                dict_frame = read_matrix(inputf, xterm, xres, yterm, yres, asymm, trunc)
                inputf.close()
            else:
                dict_frame = read_frame(inputf, nres, nlevels, nchar, dict_legend, nterm, trunc, asymm_data, dimer)
            if (stages_list):
                for i in range(len(frames_stages)):
                    stage = frames_stages[i]
                    if (not stage[3] and time >= stage[0]):
                        stage = stage[:3] + (dict_frame,) + (stage[4],)
                    if (not stage[4] and time >= stage[1]):
                        stage = stage[:4] + (dict_frame,)
                    frames_stages[i] = stage
            for i, pair in enumerate(pairs_list):
                if (not economy):
                    r = dict_frame.get(pair, trunc)
                    frame_column[i] = r
            for pair in (set(dict_frame.keys()) - set(pairs_list)):
                pairs_legend[pair] = len(pairs_list)
                i, j = pair
                if (asymm):
                    interaction_map_x[i-xterm].append(len(pairs_list))
                    interaction_map_y[j-yterm].append(len(pairs_list))
                else:
                    if (i < j-3): # exclude "trivial interactions"
                        interaction_map[i-nterm].append(len(pairs_list))
                        interaction_map[j-nterm].append(len(pairs_list))
                    pairs_legend[pair[::-1]] = len(pairs_list)
                aggregates = np.vstack([aggregates, empty_aggregate])
                pairs_list.append(pair)
                if (not economy):
                    vectors = np.vstack([vectors, emptiness])
                    frame_column = np.vstack( [   frame_column, np.array([np.float16(dict_frame[pair])])   ] )
            vectors    = np.hstack([vectors, frame_column])
            plot_frames(dict_frame, first_frame, prev_frame, pairs_list, dr_mode, xterm, xres, yterm, yres, trunc, matrices, clean_matrices, n, time, rmsd_perframe_file, asymm, dimer)
            if (dt == float('-inf')):
                print("Detected a single frame. Not very dynamic, but we will plot the distances and finish.")
                if (gnus_path): os.system("""gnuplot -e "domains=%d;maxz=%f;label_str='Residue index';inputfile='matrices/00000.dat';outputfile='frames/00000.png';"""%(domains,trunc)+
                          """title_str='Inter-residue distance';cblabel_str = 'Distance (nm)" %s/script_single.gnu"""%(gnus_path))
                sys.exit()
            if (not first_frame):
                first_frame = dict_frame
            prev_frame = dict_frame
            time_vec.append(time)
            for i, pair in enumerate(pairs_list):
                r = dict_frame.get(pair, trunc)
                agg = aggregates[i]
                aggregates[i] = update_aggregate(r, time, trunc, trunc_inter, trunc_inter_high, agg)
                if (len(time_vec) == 1 and aggregates[i][4] > 0 ):
                    native_list.append(i)
                if (pearson_inter):
                    weight = weight_contacts(r, trunc)
                    residue_column[pair[0] - nterm] += weight
                    residue_column[pair[1] - nterm] += weight
            total_contacts = sum(aggregates[:, 4])
            native         = sum(aggregates[native_list, 4])
            non_native     = total_contacts - native
            native_file.write("%9.4f %12.4f %12.4f\n"%(0.001*time, native, non_native))
            n += 1
            emptiness = trunc * np.ones([1, n], dtype = np.float16)
            if (pearson_inter):
                residue_vectors = np.hstack([residue_vectors, residue_column])
            empty_aggregate = np.array([-1, 0., 0., 0., 0., n*trunc, n*trunc**2])
        if (reread):
            if (n == len(matrixfiles)):
                time = float('-Inf')
            else:
                time = time + dt
        else:
            time_tmp = read_next_time(inputf)
            if (patch_time and time_tmp > float('-Inf')):
                time = time + dt
            else:
                time = time_tmp
    density_r2_metric = dict()
    if (pearson_inter):
        results_inter = []
        for i in range(nres):
            for j in range(i+1, nres):
                slope, intercept, r, pvalue, std_err = linregress(residue_vectors[i], residue_vectors[j])
                results_inter.append(r)
                density_r2_metric[(i, j)] = 1 - r**2
                density_r2_metric[(j, i)] = 1 - r**2
        pearson_2d_density_file = open("aggregate/pearson_map_density.dat", "w")
        for i in range(nres):
            for j in range(nres):
                if (i == j):
                    r = 1.0
                elif (i < j):
                    r = results_inter[cond_index(i, j, nres)]
                else:
                    r = results_inter[cond_index(j, i, nres)]
                pearson_2d_density_file.write("%5i %5i %15.6f\n"%(nterm + i, nterm + j, r))
            pearson_2d_density_file.write("\n")
        pearson_2d_density_file.close()
    print("")
    t0 = time_vec[0]
    tmax = time_vec[-1]
    tn = len(time_vec)
    avg_metric_dict0 = dict()
    life_metric_dict0 = dict()
    for i in range(xres):
        local_interactions = 0
        local_lifetime = 0.0
        for j in range(yres):
            if (i==j and not asymm):
                outputf_timeline.write("%4d %4d %12i %12i %12.6f %12.6f %5i\n"%(i+nterm, j+nterm, t0, tmax, 1.0, tmax-t0, 1))
                outputf_avg_rmsf.write("%4d %4d %9.4f %9.4f\n"%(i+nterm, j+nterm, 0.0, 0.0))
                continue
            pair = (i+xterm, j+yterm)
            if (pair in pairs_legend):
                agg = tuple(aggregates[pairs_legend[pair]])
            else:
                agg = empty_aggregate
            first, last, n_int, n_enc, status, sumr, sumr2 = agg
            if (first == -1):
                first = tmax
            avg = sumr/tn
            std = np.sqrt(sumr2/tn - avg**2)
            outputf_timeline.write("%4d %4d %12i %12i %12.6f %12.6f %5i\n"%(i+xterm, j+yterm, first, last, n_int/tn, tmax-t0, n_enc))
            if (interact_pair_dict):
                if (pair in interact_pair_dict and n_int/tn >= trunc_lifetime):
                    interact_file.write("%4d %4d %i\n"%(i+xterm, j+yterm, interact_pair_dict[pair]))
                else:
                    interact_file.write("%4d %4d %i\n"%(i+xterm, j+yterm, 0))
            if (n_int/tn > 0 and abs(i-j) > 3):
                local_interactions += 1
                local_lifetime += n_int/tn
            outputf_avg_rmsf.write("%4d %4d %9.4f %9.4f\n"%(i + xterm, j + yterm, avg, std))
            avg_metric_dict0[(i, j)] = avg
            life_metric_dict0[(i, j)] = 1 - n_int/tn
        outputf_timeline.write("\n")
        outputf_avg_rmsf.write("\n")
        if (interact_pair_dict):
            interact_file.write("\n")
        if (local_interactions > 0):
            avg_lifetime = 1.0*local_lifetime/local_interactions
        else:
            avg_lifetime = 0
    if (asymm):
        for i in range(xres):
            int_i = sum( x[2] for x in aggregates[interaction_map_x[i]])
            total_interactionf_x.write("%5i%12.6f\n"%(xterm + i, int_i/tn))
            total_interact_dict_x[xterm + i] = int_i/tn
        for j in range(yres):
            int_j = sum( x[2] for x in aggregates[interaction_map_y[j]])
            total_interactionf_y.write("%5i%12.6f\n"%(yterm + j, int_j/tn))
            total_interact_dict_y[yterm + j] = int_j/tn
        total_interactionf_x.close()
        total_interactionf_y.close()
    else:
        for i in range(nres):
            vec_i = np.array([x[2] for x in aggregates[interaction_map[i]]])
            int_i = sum(vec_i)/tn
            if (int_i > 0):
                int_i /= len([x for x in vec_i if x > 0])
            local_interactionf.write("%5i%12.6f\n"%(i, int_i))
            local_interact_dict[nterm + i] = int_i
        local_interactionf.close()
    outputf_timeline.close()
    outputf_avg_rmsf.close()
    native_file.close()

    if (stages_list):
        for n, stage in enumerate(frames_stages):
            if (stage[3]!=None and stage[4]!=None):
                print("Found stage:", stage[2])
            else:
                print("Did NOT find stage:", stage[2])
                continue
            stage_outf = open("matrices/stage_%04d_dr.dat"%(n+1), "w")
            title   = stage[2]
            frame_1 = stage[3]
            frame_2 = stage[4]
            for i in range(xterm, xterm+xres):
                for j in range(yterm, yterm+yres):
                    if (not asymm):
                        if (i < j):
                            pair = (i, j)
                        else:
                            pair = (j, i)
                    else:
                        pair = (i, j)
                    r0 = frame_1.get(pair, trunc)
                    r  = frame_2.get(pair, trunc)
                    dr = r-r0
                    stage_outf.write("%4d %4d %9.4f\n"%(i, j, dr))
                stage_outf.write("\n")
            stage_outf.close()
            if (gnus_path):
                os.system('gnuplot -e "domains=0;maxz=%f;inputfile='%(trunc_dr)+"'matrices/stage_%04d_dr.dat'"%(n+1)+";outputfile='frames/stage_%04d_dr.png';title_time='%s'"%(n+1,title)+'" %s/script_dr.gnu'%(gnus_path))
                if (domains):
                    os.system('gnuplot -e "domains=1;maxz=%f;inputfile='%(trunc_dr)+"'matrices/stage_%04d_dr.dat'"%(n+1)+";outputfile='frames/stage_%04d_dr.domains.png';title_time='%s'"%(n+1,title)+'" %s/script_dr.gnu'%(gnus_path))
    if (interact_pair_dict):
        interact_file.close()
        if (gnus_path): os.system('gnuplot -e "domains=%d" %s/interaction_types.gnu'%(domains,gnus_path))
    if (gnus_path): os.system('gnuplot -e "domains=%i" %s/aggregate.gnu'%(domains, gnus_path))
    rmsd_perframe_file.close()
    if (gnus_path): os.system('gnuplot -e "domains=%d" %s/1d_plots.gnu'%(domains,gnus_path))
    if (make_movie and matrices and not clean_matrices and gnus_path):
        print("Calling gnuplot for all frames... This may take a minute or two. (Aggregate plots are finished already - check them out.)")
        os.system('gnuplot -e "dr_mode=%d;domains=%d;maxz=%f;maxdr=%f;nframes=%d;begin=%d;dt=%d"'%(dr_mode,domains,trunc,trunc_dr,tn,begin,dt)
              +' %s/script_all.gnu'%(gnus_path))
    if (make_movie and gnus_path):
    # In this case we make a video using mencoder
        os.system('mencoder -vf scale=1200:1080 -ovc x264 -x264encopts bitrate=2000 -mf fps=24:type=png -oac copy -noskip -o movies/CONAN.mp4 "mf://frames/%05d.png"')
        if (dr_mode in [1, 3]):
            os.system('mencoder -vf scale=1200:1080 -ovc x264 -x264encopts bitrate=2000 -mf fps=24:type=png -oac copy -noskip -o movies/CONAN_dr.init.mp4 "mf://frames/%05d_dr.init.png"')
        if (dr_mode in [2, 3]):
            os.system('mencoder -vf scale=1200:1080 -ovc x264 -x264encopts bitrate=2000 -mf fps=24:type=png -oac copy -noskip -o movies/CONAN_dr.prev.mp4 "mf://frames/%05d_dr.prev.png"')
        if (domains):
            os.system('mencoder -vf scale=1200:1080 -ovc x264 -x264encopts bitrate=2000 -mf fps=24:type=png -oac copy -noskip -o movies/CONAN.domains.mp4 "mf://frames/%05d.domains.png"')
            if (dr_mode in [1, 3]):
                os.system('mencoder -vf scale=1200:1080 -ovc x264 -x264encopts bitrate=2000 -mf fps=24:type=png -oac copy -noskip -o movies/CONAN_dr.init.domains.mp4 "mf://frames/%05d_dr.init.domains.png"')
            if (dr_mode in [2, 3]):
                os.system('mencoder -vf scale=1200:1080 -ovc x264 -x264encopts bitrate=2000 -mf fps=24:type=png -oac copy -noskip -o movies/CONAN_dr.prev.domains.mp4 "mf://frames/%05d_dr.prev.domains.png"')
    print("number of frames=",len(time_vec))
    for i, vec in enumerate(vectors):
        if (len(time_vec)!=len(vec)):
            print("Length mismatch!", i, len(time_vec), len(vec), pairs_list[i])
            print(vec)
            sys.exit()
    if (not asymm):
        return begin, dt, time_vec, pairs_list, pairs_legend, vectors, avg_metric_dict0, life_metric_dict0, density_r2_metric, local_interact_dict
    else:
        return begin, dt, time_vec, pairs_list, pairs_legend, vectors, avg_metric_dict0, life_metric_dict0, density_r2_metric, total_interact_dict_x, total_interact_dict_y
    
    
def create_list_of_frames(directory, begin, end, dt):
    i = 1
    sfile = directory + "%05d.dat"%i
    while (not os.path.isfile(sfile)):
       i = i+1
# we found the first file in the subfolder...
    list_of_frames = []
    time = -1
    if (end == None):
        end = float('inf')
    if (begin == None):
        begin = float('inf')
    while (os.path.isfile(sfile) and (time <= end)):
        matrixf = open(sfile)
        firstline = matrixf.readline()
        time  = int(firstline.split()[2])
        nterm = int(firstline.split()[5])
        nres  = int(firstline.split()[7])
        if ( (time >= begin) and ((time-begin) % dt == 0) and (time <= end) ):
            list_of_frames.append((time, sfile))
        i = i+1
        matrixf.close()
        sfile = directory + "%05d.dat"%i
    return nterm, nres, list_of_frames

def read_matrix(inputf, xterm, xres, yterm, yres, asymm, trunc):
# This function reads in ONE matrix.
# inputf: an open (plain-text) file output by CONAN (or having the same structure).
# Empty lines will be automatically detected.
# Not all entries need be there.
# In case of an asymmetric run, the indices must correspond to the read in ones.
# For example, if xterm is 1..100 and yterm is 101..200, the function will expect the same numbers in the file.
# Any line about a pair that is NOT there will be just ignored.
    dict_frame = dict()
    for line in inputf:
        if (line.strip().startswith('#') or len(line.split()) == 0):
            continue
        n1 = int(line.split()[0])
        n2 = int(line.split()[1])
        r  = float(line.split()[2])
        if (n1 >= xterm and n1 < xterm + xres and n1 >= xterm and n1 < xterm + xres and r < trunc):
            if asymm:
                dict_frame[(n1,n2)] = r
            elif (n1 < n2):
                dict_frame[(n1,n2)] = r
    return dict_frame

def read_pca(pca_dir, xterm, xres, yterm, yres, time_vec, asymm, vectors_mean0, pairs_list):
    v = []
    i = 1
    sfile = pca_dir + "/pc.%i.dat"%i
    while (os.path.isfile(sfile)):
        pc_file = open(sfile)
        pc_dict = read_matrix(pc_file, xterm, xres, yterm, yres, asymm, 1.0)
        vec_pc = []
        sum = 0.0
        for pair in pairs_list:
            coeff = pc_dict[pair]
            vec_pc.append(coeff)
            sum += coeff**2
        print("Norm of PC # %i = %12.6f"%(i, math.sqrt(sum)))
        v.append(vec_pc)
        i = i+1
        sfile = pca_dir + "/pc.%i.dat"%i
    n = i - 1
    print("Reading in PC's.")
    print("I found %i principal components..."%n)
    v = np.array(v)
    v = np.transpose(v)
    project_pc(v, time_vec, vectors_mean0, pairs_list, change_sign = False)

def handle_backups():
    backup = 0
    dirs = ['frames', 'input', 'matrices', 'movies', 'aggregate', 'cluster_res', 'cluster_trj', 'pca', 'zoom', 'blocking']
    dirs_exist = []
    for dir in dirs:
        if (os.path.isdir(dir)):
            dirs_exist.append(dir)
    if (dirs_exist):
        backup = 1
        back_dir_name = "backup.1"
        i = 1
        while (os.path.isdir(back_dir_name)):
            i = i+1
            back_dir_name = "backup.%i"%i
        print("Moving directories containing old calculation to: %s ..."%back_dir_name)
        os.system('mkdir %s'%back_dir_name)
        for dir in dirs_exist:
            os.system('mv %s %s'%(dir, back_dir_name))
    return backup

def read_all_matrices(list_of_frames, nterm, nres, trunc):
    for frame in list_of_frames:
        time  = frame[0]
        sfile = frame[1]
        matrixf = open(sfile)
        dict_frame = read_matrix(matrixf, nterm, nres, trunc)
        matrixf.close()
        all_data.append(time, dict_frame)
    return all_data

def get_contact_density(dict_frame, nterm, nres, trunc):
    contact_density_dict = dict()
    for i in range(nterm, nterm+nres):
        contact_density_dict[i] = 0.0
    for pair in dict_frame:
        i = pair[0]
        j = pair[1]
        if (i<j):
            r = dict_frame[pair]
            contact_density_dict[i] += weight_contacts(r, trunc)
            contact_density_dict[j] += weight_contacts(r, trunc)
    return contact_density_dict

def density_dicts_analyze(contact_density_dicts, nterm, nres):
    time_density_file = open("aggregate/time_density.dat", "w")
    pearson_1d_time_density_file = open("aggregate/pearson_sequence_time_density.dat", "w")
    pearson_2d_density_file = open("aggregate/pearson_map_density.dat", "w")
    density_r2_metric = dict()
    for frame in contact_density_dicts:
        time = int(frame[0])
        contact_density_dict = frame[1]
        contact_density = 0.0
        for i in range(nterm, nterm+nres):
            contact_density += contact_density_dict[i]
        time_density_file.write("%12i %15.6f\n"%(time, contact_density))
    time_density_file.close()
    vec_t = [frame[0] for frame in contact_density_dicts]
    for i in range(nterm, nterm+nres):
        vec_i = [frame[1][i] for frame in contact_density_dicts]
        slope, intercept, r, pvalue, std_err = linregress(vec_i, vec_t)
        pearson_1d_time_density_file.write("%5i %15.6f\n"%(i, r))
        for j in range(nterm, nterm+nres):
            vec_j = [frame[1][j] for frame in contact_density_dicts]
            if (i==j):
                r = 1.0
                density_r2_metric[(i,j)] = 0.0
            else:
                slope, intercept, r, pvalue, std_err = linregress(vec_i, vec_j)
                density_r2_metric[(i,j)] = 1-r**2
                if (pvalue > 0.05):
                    r = 0
            pearson_2d_density_file.write("%5i %5i %15.6f\n"%(i, j, r))
        pearson_2d_density_file.write("\n")
    return density_r2_metric
            
def get_best_center(cluster, metric_cond, n):
    best_total_dist = float('Inf')
    best_center = 0
    for i in cluster:
        total_dist = 0
        for j in cluster:
            if (i < j):
                total_dist += metric_cond[cond_index(i, j, n)]
            elif (i > j):
                total_dist += metric_cond[cond_index(j, i, n)]
        if (total_dist < best_total_dist):
            best_total_dist = total_dist
            best_center = i
    return best_center

def assign_tiebreak(metric_cond, clusters, resid, n):
    best_cluster = -1
    best_distance = float('inf')
    for i in range(len(clusters)):
        avg_dist = 0.0
        cluster = clusters[i]
        for res in cluster:
            if (resid < res):
                avg_dist += metric_cond[cond_index(resid, res, n)]
            elif (res < resid):
                avg_dist += metric_cond[cond_index(res, resid, n)]
        avg_dist /= len(cluster)
        if (avg_dist < best_distance):
            best_distance = avg_dist
            best_cluster = i
    return i

def block_analysis(vectors, nres, dt):
    b = np.array(vectors)
    i = 0
    block_outf = open("blocking/block_out.txt", "w+")
    block_outf.write("#N_iter N_blocks stderr (nm)  stderr_err (nm) corr.time est. (ps)\n")
    while (len(b) >= 2):
        sems = sem(b)
        sem_avg = np.linalg.norm(sems)/nres
        sem_avg_err = sem_avg*(1.0 / np.sqrt(2*(len(b)-1)) )
        if (i == 0):
            sem_avg0 = sem_avg
        block_outf.write(" %4i %9i %12.3e %12.3e   %15.2f\n"%(i, len(b), sem_avg, sem_avg_err, dt*(sem_avg/sem_avg0)**2 ) )
        c = [ 0.5* (b[2*i] + b[2*i+1]) for i in range(int(len(b)/2)) ]
        b = np.array(c)
        i = i+1
    block_outf.close()

def assign_elements(n_elements, metric_cond, centers):
    clusters = [ [x] for x in centers ]
    nowhere = []
    cost_list = []
    for i in range(n_elements):
        if (i not in centers):
            distances = []
            for j in range(len(centers)):
                if (i < centers[j]):
                    distances.append(metric_cond[cond_index(i, centers[j], n_elements)])
                else:
                    distances.append(metric_cond[cond_index(centers[j], i, n_elements)])
            best_dist = min(distances)
            if (len(set(distances))>1):
                best_cluster = distances.index(best_dist)
                clusters[best_cluster].append(i)
            else:
                nowhere.append(i)
#these are residues with the same distance from every center, ie, "infinitely far". we do a lookup for them later, based on elements we *did* manage to assign.
            cost_list.append(best_dist)
        else:
            cost_list.append(0.0)
    
    for i in nowhere:
        best_cluster = assign_tiebreak(metric_cond, clusters, i, n_elements)
        clusters[best_cluster].append(i)
    
    return cost_list, clusters

def print_res_clusters(summary_file, nterm, nres, clusters, centers):
    for i in range(len(clusters)):
        cluster = clusters[i]
        center = centers[i]
        summary_file.write("Cluster %i with center at resid %i has # of residues: %i #chimera\n"%(i, center, len(cluster)))
        summary_file.write("VMD selection: \n")
        vmd_str="resid "
        chimera_str = "Cluster%i residues: ["%i
        cluster.sort()
        current_section = [cluster[0]]
        for el in cluster[1:]:
            if (el > current_section[-1]+1):
                if (len(current_section)==1):
                    vmd_str = vmd_str + "%i "%(nterm+current_section[-1])
                    chimera_str = chimera_str + "%i, "%(nterm+current_section[-1])
                else:
                    vmd_str = vmd_str + "%i to %i "%(nterm+current_section[0],nterm+current_section[-1])
                    chimera_str = chimera_str + "%i-%i, "%(nterm+current_section[0],nterm+current_section[-1])
                current_section=[el]
            else:
                current_section.append(el)
        if (len(current_section)==1):
            vmd_str = vmd_str + "%i\n"%(nterm+current_section[-1])
            chimera_str = chimera_str + "%i] #chimera\n"%(nterm+current_section[-1])
        else:
            vmd_str = vmd_str + "%i to %i\n"%(nterm+current_section[0],nterm+current_section[-1])
            chimera_str = chimera_str + "%i-%i] #chimera\n"%(nterm+current_section[0],nterm+current_section[-1])
        summary_file.write(vmd_str)
        summary_file.write(chimera_str)

# get a condensed distance matrix:
def get_condensed_distances(metric_dict, nterm, nres):
    condensed_dist = []
    for i in range(nres):
        for j in range(i+1, nres):
            pair = (i+nterm, j+nterm)
            condensed_dist.append(metric_dict[pair])
    return condensed_dist

def print_hierarchical_clusters(dir, vectors, trunc, trunc_inter, nframes, clusters, centers, time_vec, cost_list, xterm, xres, yterm, yres, asymm, type, matrices):
    cluster_assignment_file = open("%s/assignment.dat"%dir, "w")
    cluster_assignment = []
    for i in range(len(clusters)):
        for j in clusters[i]:
            cluster_assignment.append((j, i))
        if (type == 'trj'):
            os.system("mkdir %s/cluster%i"%(dir, i))
            avg_cluster_file = open("%s/cluster%i/average_map.dat"%(dir, i), "w")
            lifetime_file = open("%s/cluster%i/lifetime.dat"%(dir, i), "w")
            avg_dict = dict()
            n = 0
            if (matrices):
                for j in clusters[i]:
                    os.system("ln -s ../../../frames/%05i.png %s/cluster%i/%05i.png"%(j, dir, i, n))
                    if (j in centers):
                        os.system("ln -s ../../../frames/%05i.png %s/cluster%i/center.png"%(j, dir, i))
                    n = n + 1
            for res_i in range(xterm, xterm+xres):
                for res_j in range(yterm, yterm+yres):
                    if (res_i == res_j and not asymm):
                        mean_d = 0
                        interact = 100.0
                    elif ((res_i, res_j) in pairs_legend):
                        vec = vectors[pairs_legend[(res_i, res_j)], :][clusters[i]]
                        if (len(vec) > 0):
                            interact = sum([100.0 for x in vec if x < trunc_inter])/len(vec)
                        else:
                            interact = 0.0
                        mean_d = np.mean(vec)
                    else:
                        mean_d = trunc
                        interact = 0.0
                    avg_cluster_file.write("%4d %4d %10.6f\n"%(res_i, res_j, mean_d))
                    lifetime_file.write("%4d %4d %10.6f\n"%(res_i, res_j, interact))
                avg_cluster_file.write("\n")
                lifetime_file.write("\n")
            avg_cluster_file.close()
            lifetime_file.close()
            if (gnus_path):
                os.system("""gnuplot -e "domains=%d;maxz=%f;label_str='Residue index';inputfile='%s/cluster%i/average_map.dat';outputfile='%s/cluster%i/average_map.png';"""%(domains, trunc, dir, i, dir, i)+
                      """title_str='Average contact map for cluster %i';cblabel_str = 'Distance (nm)'" %s/script_single.gnu"""%(i, gnus_path))
                os.system("""gnuplot -e "domains=%d;maxz=0;label_str='Residue index';inputfile='%s/cluster%i/lifetime.dat';outputfile='%s/cluster%i/lifetime.png';"""%(domains, dir, i, dir, i)+
                      """title_str='Interaction lifetime for cluster %i';cblabel_str = 'Lifetime (%%)'" %s/script_single.gnu"""%(i, gnus_path))
        else:
            summary_file = open("%s/summary.txt"%(dir), "w")
            print_res_clusters(summary_file, nterm, nres, clusters, centers)
            summary_file.close()
    cluster_assignment.sort(key = lambda x: x[0])
    cluster_assignment_file.write("# element_ind cluster_ind center? dist\n")
    for ind, el in cluster_assignment:
        cres = 0
        if (ind in centers):
            cres = 1
        cluster_assignment_file.write("%5i %i %i %9.3f\n"%(ind, el, cres, cost_list[ind]) )
    cluster_assignment_file.close()

def parse_numbers(vals):
    ks = set([])
    numbers = vals.replace("-", " - ").split()
    a = 0
    i = 0
    while (i < len(numbers)):
        if (numbers[i] == '-'):
            if (i + 1 < len(numbers)):
                l = set(range(a+1, int(numbers[i+1])+1))
                i += 1
        else:
            l = set([ int(numbers[i]) ])
            a = int(numbers[i])
        i += 1
        ks |= l
    return list(ks)

def call_hierarchical(dir, vectors, trunc, trunc_inter, condensed_dist, n_elements, ksi, title, time_vec, type, xterm, xres, yterm, yres, asymm, matrices):
    Z = linkage(condensed_dist, 'ward')
    linkage_file = open("%s/linkage.txt"%dir, "w+")
    for row in Z:
        linkage_file.write("%6i %6i %12.6f %6i\n"%tuple(row))
    linkage_file.close()
    plt.figure(figsize=(25, 10))
    plt.xlabel('sample index')
    plt.ylabel('distance')
    font = {'size': 20}
    plt.rc('font', **font)
    if (n_elements > 100):
        dendrogram(Z, leaf_rotation=90., leaf_font_size=12., truncate_mode='lastp', p = 100, show_contracted = True)
    else:
        dendrogram(Z, leaf_rotation=90., leaf_font_size=12.)
    plt.savefig(dir+'/dendrogram.png')
    if (ksi == [0]):
        print("Please examine the dendrogram and give the desired number(s) of clusters. (Close the window to continue).")
        plt.show()
        kstr = input("Number of clusters? (multiple values are allowed, or a single interval as i-j).\n")
        ks = parse_numbers(kstr)
    else:
        ks = ksi
    for k in ks:
        cost_list = [0 for i in range(n_elements)]
        os.system("mkdir -p %s/%i"%(dir, k))
        centers = []
        assignments = fcluster(Z, k, criterion='maxclust')
        assignments = list(assignments)
        clusters = [ [] for x in range(k) ]
        for i in range(n_elements):
            clusters[assignments[i]-1].append(i)
        clusters.sort(key = lambda x: len(x), reverse=True)
        for j, cluster in enumerate(clusters):
            center = get_best_center(cluster, condensed_dist, n_elements)
            centers.append(center)
            for frame in cluster:
                d = 0
                if (frame < center):
                    d = condensed_dist[cond_index(frame, center, n_elements)]
                elif (center < frame):
                    d = condensed_dist[cond_index(center, frame, n_elements)]
                cost_list[frame] = d
        clusters_centers = [ (clusters[i], centers[i]) for i in range(k) ]
        clusters_centers.sort(key = lambda x: x[1])
        clusters = [ x[0] for x in clusters_centers]
        centers  = [ x[1] for x in clusters_centers]
        print_hierarchical_clusters("%s/%i/"%(dir, k), vectors, trunc, trunc_inter, nframes, clusters, centers, time_vec, cost_list, xterm, xres, yterm, yres, asymm, type, matrices)
        if (gnus_path): os.system("""gnuplot -e "inputdir='%s/%i'" %s/cluster_assignment.gnu"""%(dir, k, gnus_path))

def rmsd_frame(dict_a, dict_b, trunc, xres, yres, asymm, dimer):
    sum2 = 0.0
    for pair in (set(dict_a.keys()) | set(dict_b.keys())):
        ra = dict_a.get(pair, trunc)
        rb = dict_b.get(pair, trunc)
        sum2 += (ra - rb)**2
    if (dimer):
        dict_b_t = {(j, i): dict_b[(i, j)] for i, j in dict_b}
        sum2_t = 0.0
        for pair in (set(dict_a.keys()) | set(dict_b_t.keys())):
            ra = dict_a.get(pair, trunc)
            rb = dict_b_t.get(pair, trunc)
            sum2_t += (ra - rb)**2
        if (sum2_t < sum2):
            sum2 = sum2_t
    if (asymm):
        return math.sqrt(sum2/(xres*yres))
    else:
        return math.sqrt(2*sum2)/nres
# In[ ]:

# Function to check the input file for mdmat and run mdmat

def mdmat_run(program, traj, coord, nlevels, trunc, mean, begin, end, dt, indexf):
    command = program + " mdmat -frames -f %s -s %s"%(opts['traj'], opts['coord'])
    if (nlevels):
        command = command+" -nlevels %d"%nlevels
    if (trunc):
        command = command+" -t %f"%trunc
    if (mean):
        command = command+" -mean"
    if (begin):
        command = command+" -b %d"%begin
    if (end):
        command = command+" -e %d"%end
    if (dt):
        command = command+" -dt %d"%dt
    if (indexf):
        command = command+" -n %s"%indexf
    os.system(command)

def read_index_folder(reread):
    finput = open('%s/index.dat'%reread)
    asymm = 0
    nterm = 1
    trunc = 1.5
    nres = None
    for line in finput:
        if line.strip().startswith('#') or len(line.split())<2:
            continue
        else:
            key = line.split()[0].upper()
            val = line.split()[1]
            if key == 'ASYMM' and val=='yes':
                asymm = 1
            elif key == 'NTERM':
                nterm = int(val)
            elif key == 'NRES':
                nres = int(val)
            elif key == 'TRUNC':
                trunc = float(val)
    return asymm, nterm, nres, trunc
    
def read_options(finput, opts):
    opts['trunc_given'] = False
    for line in finput:
        if line.strip().startswith('#') or len(line.split())<2:
            continue
        else:
            key = line.split()[0].upper()
            val = line.split()[1]
            if key == 'TRAJ':
                opts['traj'] = str(val)
            elif key == 'COORD':
                opts['coord'] = str(val)
            elif key == 'NLEVEL':
                opts['nlevels'] = int(val)
            elif key == 'TRUNC':
                opts['trunc'] = float(val)
                opts['trunc_given'] = True
            elif key == 'PATCH_TIME':
                opts['patch_time'] = (val == 'yes')
            elif key == 'MEAN':
                opts['mean'] = (val=='yes')
            elif key == 'INDEX':
                opts['indexf'] = str(val)
            elif key == 'BEGIN':    # ps!!!
                opts['begin'] = int(val)
            elif key == 'END':      # ps!!!
                opts['end'] = int(val)
            elif key == 'DT':
                opts['dt']  = int(val)
            elif key == 'TRUNC_INTER':
                opts['trunc_inter'] = float(val)
            elif key == 'TRUNC_INTER_HIGH':
                opts['trunc_inter_high'] = float(val)
            elif key == 'TRUNC_LIFETIME':
                opts['trunc_lifetime'] = float(val)
            elif key == 'TRUNC_DR':
                opts['trunc_dr'] = float(val)
            elif key == 'DR_MODE':
                if (val== 'init'):
                    opts['dr_mode'] = 1
                elif (val== 'prev'):
                    opts['dr_mode'] = 2
                elif (val== 'both'):
                    opts['dr_mode'] = 3
            elif key == 'GNUS_PATH':
                opts['gnus_path'] = str(val)
            elif key == 'RUN_MDMAT':
                opts['run_mdmat'] = (val == 'yes')
            elif key == 'MATRICES':
                opts['matrices'] = (val == 'yes')
            elif key == 'CLEAN_MATRICES':
                opts['clean_matrices'] = (val == 'yes')
            elif key == 'NTERM':
                opts['nterm'] = int(val)
                opts['xterm'] = int(val)
                opts['yterm'] = int(val)
            elif key == 'START_X':
                opts['start_x'] = int(val)
                if ('xterm' not in opts):
                    opts['xterm'] = int(val)
                opts['asymm'] = 1
            elif key == 'START_Y':
                opts['start_y'] = int(val)
                if ('yterm' not in opts):
                    opts['yterm'] = int(val)
                opts['asymm'] = 1
            elif key == 'NRES_X':
                opts['xres'] = int(val)
                opts['asymm'] = True
            elif key == 'NRES_Y':
                opts['yres'] = int(val)
                opts['asymm'] = True
            elif key == 'NTERM_X':
                opts['xterm'] = int(val)
                opts['asymm'] = True
            elif key == 'NTERM_Y':
                opts['yterm'] = int(val)
                opts['asymm'] = 1
            elif key == 'PEARSON_TIME':
                opts['pearson_time'] = (val == 'yes')
            elif key == 'PEARSON_INTER':
                opts['pearson_inter'] = (val == 'yes')
            elif key == 'SHADOW_TOL':
                opts['shadow'] = True
                opts['shadow_tol'] = float(val)
            elif key == 'PCA_MAKE':
                opts['pca_make'] = int(val)
            elif key == 'PCA_READ':
                opts['pca_dir'] = str(val)
            elif key == 'PEARSON_OBS':
                opts['pearson_obs_str'] = str(val)
            elif key == 'REREAD':
                opts['reread'] = str(val)
            elif key == 'ZOOM_LIST':
                opts['zoom'] = True
                zoom_file_name = str(val)
                zoom_file = open(zoom_file_name)
                opts['zoom_list'] = read_zoom_list(zoom_file)
                zoom_file.close()
            elif key == 'MAKE_MOVIE':
                opts['make_movie'] = (val != 'no')
            elif key == 'BLOCKING':
                if (val != 'no'):
                    opts['blocking'] = True
            elif key == 'COORD_PDB':
                opts['coordpdb'] = True
                opts['pdbf'] = str(val)
            elif key == 'COORD_PDB_X':
                opts['coordpdb'] = True
                opts['pdbx'] = str(val)
            elif key == 'COORD_PDB_Y':
                opts['coordpdb'] = True
                opts['pdby'] = str(val)
            elif key == 'IGNORE_OBS_TIME' and (val != 'no'):
                opts['ignore_obs_time'] = True
            elif key == 'K_RES_CLUSTERS':
                vals = " ".join(line.split()[1:])
                ks = parse_numbers(vals)
                if (ks != [0]):
                    print("Will create the following number of clusters of residues:"," ".join(["%i"%k for k in ks]) )
                else:
                    print("Will choose the number of clusters of residues interactively.")
                opts['k_res_clusters'] = ks
                opts['ks_res'] = ks
            elif key == 'K_TRAJ_CLUSTERS':
                vals = " ".join(line.split()[1:])
                ks = parse_numbers(vals)
                if (ks != [0]):
                    print("Will create the following number of clusters of frames:"," ".join(["%i"%k for k in ks]) )
                else:
                    print("Will choose the number of clusters of frames interactively.")
                opts['k_traj_clusters'] = ks
                opts['ks_traj'] = ks
            elif key == 'ECONOMY':
                opts['economy'] = (val != 'no')
            elif key == 'DOMAINS':
                opts['domains'] = True
                opts['domainf'] = str(val)
            elif key == 'DOMAINS_X':
                opts['domains'] = True
                opts['domainf_xn'] = str(val)
            elif key == 'DOMAINS_Y':
                opts['domains'] = True
                opts['domainf_yn'] = str(val)
            elif key == 'STAGES':
                opts['stages'] = True
                opts['stagef'] = str(val)
            elif key == 'DIMER':
                opts['dimer'] = (val == 'yes')
            else:
                print("Cannot understand the keyword %s given in the input file." % key)
    if (opts['asymm']):
        if ('shadow' in opts):
            print("Removing shadows is impossible in an asymmetric run.")
            opts.pop('shadow')
        if ('k_clusters' in opts):
            print("Clustering residues is impossible in an asymmetric run.")
            opts.pop('k_clusters')
        if (opts.get('pearson_inter', False)):
            print("Inter-residue cross-correlation is impossible in an asymmetric run.")
            opts['pearson_inter'] = False
    if ('pca_make' in opts) and ('pca_dir' in opts):
        print("PCA can either be performed or PC's can be read in, but not both. We will attempt reading them in.")
        opts.pop('pca_make')
    if ('gnus_path' not in opts):
        print(" .gnu paths have not been given. NO PLOTS will be produced.")
        opts['gnus_path'] = None
    if (not opts['matrices'] and opts['make_movie']):
        print("Can't make a movie without matrices! We will turn off movie making.")
        opts['make_movie'] = False
    if (opts['dimer']):
        if (not opts['asymm'] or not (opts.get('xterm', 0) == opts.get('yterm', 1) and opts.get('xres', 0) == opts.get('yres', 1))):
            print("""'Dimer mode' requires 1. asymmetric mode and 2. the indices of the X- and Y- residues to be identical. We will turn off the option now.""")
            opts['dimer'] = False
    if (opts['economy']):
        if (opts['zoom']):
            print("Zooming onto residues is not available in economy mode.")
            opts.pop('zoom')
        if ('pca_make' in opts):
            print("PCA is not available in economy mode.")
            opts.pop('pca_make')
        if ('k_traj_clusters' in opts):
            print("Trajectory clustering is not available in economy mode.")
            opts.pop('k_traj_clusters')
        if ('pearson_time' in opts):
            print("Pearson correlation with time is not available in economy mode.")
            opts.pop('pearson_time')
    return opts

def adjust_truncs(trunc, trunc_dr, trunc_inter, trunc_inter_high):
    if (not trunc_dr):
        trunc_dr = trunc
        print("No distance cutoff given for differential plot, will use maximum cutoff:", trunc)
    if (not trunc_inter or trunc_inter == trunc):
        trunc_inter      = 0.999999 * trunc
        trunc_inter_high = 0.999999 * trunc
        print("No distance cutoff given for defining interactions, will use maximum cutoff:",trunc)
    elif (not trunc_inter_high or trunc_inter_high == trunc):
        print("No buffer value given for the interaction cutoff. We will use the same value for both.")
        trunc_inter_high = trunc_inter
    if (trunc_inter_high < trunc_inter):
        print("""It must hold that TRUNC_INTER_HIGH >= TRUNC_INTER. We will set them to be equal.""")
        trunc_inter_high = trunc_inter
    return trunc_dr, trunc_inter, trunc_inter_high

if __name__ == '__main__':

    program = which("gmx")
    if program == '':
        print("Gromacs executable was not found")
    else:
        print("Using gmx executable in %s" % program)
# Set default value if none given:
    opts = {'trunc' : None, 'mean' : True, 'run_mdmat' : True, 'nterm' : 1, 'matrices' : True, 'dr_mode' : 0, 'trunc_lifetime' : 0.5, 'asymm' : False, 'clean_matrices' : False, 'dimer': False,
        'make_movie': True, 'trunc_dr': None, 'trunc_inter': None, 'trunc_inter_high': None, 'domains' : 0, 'begin' : None, 'dt': None, 'end': None, 'patch_time': False, 'economy': False,
        'pearson_inter': False, 'reread': False, 'indexf': None, 'zoom': False}
    if len(sys.argv)!=2:
        print('Usage: You need to provide a single text file in which specify all the options for creating a matrix!')
    finput=open(sys.argv[1])
    
# read in the options...
    opts = read_options(finput, opts)
    globals().update(opts)
# just add all the read in options into the global namespace.
    asymm_data = None
    if (asymm):
        print("Asymmetric run found!")
        asymm_data = start_x, start_y, xres, yres, xterm, yterm
    if ('coordpdb' in opts and asymm):
        pdbf = opts.get('pdbf', '')
        pdbx = opts.get('pdbx', '')
        pdby = opts.get('pdby', '')
        if (pdbx + pdby == ''):
            print ("No separate PDBs for X- and Y-residues given, will use the same one.")
            pdbx = pdbf
            pdby = pdbf
        elif (pdbx==''):
            print ("No PDB for X-residues given, will use the same as for Y.")
            pdbx = pdby
        elif (pdby==''):
            print ("No PDB for Y-residues given, will use the same as for X.")
            pdby = pdbx
        inputpdb_x = open(pdbx)
        inputpdb_y = open(pdby)
        interact_pair_dict = read_sequence_asymm(inputpdb_x, inputpdb_y, xterm, xres, yterm, yres)
        inputpdb_x.close()
        inputpdb_y.close()
    elif ('coordpdb' in opts):
        inputpdb = open(pdbf)
        interact_pair_dict = read_sequence(inputpdb)
        inputpdb.close()
    else:
        interact_pair_dict = None
    trunc = opts['trunc']
    
    backup = handle_backups()
    trunc_input = trunc
    if (opts.get('run_mdmat', False) and not reread):
        mdmat_run(program, traj, coord, nlevels, trunc, mean, begin, end, dt, indexf)
    if (not reread):
        inputf = open('dmf.xpm')
        nres, nlevels, nchar, dict_legend, trunc = read_index(inputf)
        inputf.close()
        inputf = open('dmf.xpm')
    else:
        asymm, nterm, nres, trunc = read_index_folder(reread)
        print ("asymm, nterm, nres, trunc", asymm, nterm, nres, trunc)
        inputf = None
        nlevels = None
        nchar = None
        dict_legend = None
    os.system('mkdir -v frames matrices aggregate movies input')
    os.system('cp %s input'%sys.argv[1])
    if (not asymm):
        xres = nres
        yres = nres
        xterm = nterm
        yterm = nterm
    if (trunc != trunc_input and trunc_given):
        print("Warning! Truncation threshold given is different to the original one. Will use the given one!")
        print("Given:     %9.4f"%trunc_input)
        print("Detected:  %9.4f"%trunc)
        trunc = trunc_input
    trunc_dr, trunc_inter, trunc_inter_high = adjust_truncs(trunc, trunc_dr, trunc_inter, trunc_inter_high)
    print("Will read in the data now...")
    if (domains):
        if (not asymm):
            domain_file = open(domainf)
            create_tics(domain_file, nterm, nres)
        else:
            domainf_x = open(domainf_xn)
            domainf_y = open(domainf_yn)
            create_tics_asymm(domainf_x, domainf_y, xterm, xres, yterm, yres)

    stages_list = None
    if ('stages' in opts):
        stage_file = open(stagef)
        stages_list = read_stages(stage_file)

    if (not asymm):
        begin, dt_read, time_vec, pairs_list, pairs_legend, vectors, avg_metric_dict0, life_metric_dict0, density_r2_metric, local_interact_dict = read_process_frames(inputf, nres, nlevels, nchar, dict_legend, asymm_data, trunc, trunc_inter, trunc_inter_high,
                                                                                  begin, end, dt, patch_time, gnus_path, dr_mode, domains, pearson_inter, economy, dimer, reread, stages_list, interact_pair_dict, trunc_lifetime)
    else:
        begin, dt_read, time_vec, pairs_list, pairs_legend, vectors, avg_metric_dict0, life_metric_dict0, density_r2_metric, total_interact_dict_x, total_interact_dict_y = read_process_frames(inputf, nres, nlevels, nchar, dict_legend, asymm_data, trunc, trunc_inter, trunc_inter_high,
                                                                                  begin, end, dt, patch_time, gnus_path, dr_mode, domains, pearson_inter, economy, dimer, reread, stages_list, interact_pair_dict, trunc_lifetime)
    if ('dt' not in opts):
        dt = dt_read
    if (not reread):
        inputf.close()
    nframes  = len(time_vec)
    
    if (zoom):
        print("Zooming in on interesting residue pairs...")
        os.system("mkdir zoom")
        for pair in zoom_list:
            if (pair in pairs_list):
                rvec = vectors[pairs_legend[pair]]
            else:
                rvec = [trunc for i in time_vec]
                print("Warning, pair (%i, %i) in zoom list never formed an interaction!"%pair)
            zoom_on_pair(pair, rvec, time_vec, trunc, trunc_inter, trunc_inter_high, pairs_list, pairs_legend, vectors)
        
    if ('coordpdb' in opts):
        if (not asymm):
            inputpdb = open(pdbf)
            outputpdb = open("aggregate/local_interactions.pdb","w")
            change_bf_in_pdb(inputpdb, outputpdb, local_interact_dict)
            outputpdb.close()
        else:
            inputpdb_x = open(pdbx)
            outputpdb_x = open("aggregate/total_interactions.x.pdb","w")
            change_bf_in_pdb(inputpdb_x, outputpdb_x, total_interact_dict_x)
            inputpdb_x.close()
            outputpdb_x.close()
            inputpdb_y = open(pdby)
            outputpdb_y = open("aggregate/total_interactions.y.pdb","w")
            change_bf_in_pdb(inputpdb_y, outputpdb_y, total_interact_dict_y)
            inputpdb_y.close()
            outputpdb_y.close()
    if opts.get('pearson_inter', False) and gnus_path:
        os.system('gnuplot -e "domains=%i;maxz=%f;inputfile='%(domains,trunc_dr)+"'aggregate/pearson_map_density.dat'"+";outputfile='aggregate/pearson_inter_residue.png';title_time='Inter-residue Pearson correlation'"+'" %s/script_corr.gnu'%(gnus_path))
    
    if (gnus_path): os.system('gnuplot -e "domains=%d" %s/1d_plots.gnu'%(domains,gnus_path))
    
    
    pearsonerr = False
    if ('pearson_time' in opts):
        pearsonf = open("aggregate/pearson_data.dat","w")
        pearson_results = dict()
        for i in range(xterm, xterm + xres):
            for j in range(yterm, yterm + yres):
                if (i, j) in pairs_legend:
                    rvec = vectors[pairs_legend[i, j]]
                    pearsonerr = (len(rvec) != len(time_vec))
                    if (pearsonerr):
                        print("Pearson error! Wrong number of frames.", len(rvec), len(time_vec))
                        print(i, j)
                        print(pairs_legend[i, j])
                        #print(rvec)
                        print(pairs_legend)
                        break
                    slope, intercept, rvalue, pvalue, std_err = linregress(time_vec, rvec)
                else:
                    rvalue = 0.0
                    slope = 0.0
                pearson_results[i, j]   = rvalue
            print("Computing correlations for residue:", i, end="\r")
            if (pearsonerr):
                break
        if (not pearsonerr):
            for i in range(xterm, xterm + xres):
                for j in range(yterm, yterm + yres):
                    pearsonf.write("%4d %4d %10.6f\n"%(i,j,pearson_results[i, j]))
                pearsonf.write("\n")
            pearsonf.close()
            if (gnus_path): os.system('gnuplot -e "domains=%d;maxz=1.0;inputfile='%(domains) + "'aggregate/pearson_data.dat'"+";outputfile='aggregate/pearson_time.png';title_time='Pearson correlation with time for pairwise distances'"+'" %s/script_corr.gnu'%(gnus_path))
    print("")
    
    if ('pearson_obs_str' in opts and not pearsonerr):
        pears_obsf = open(pearson_obs_str)
        if (not opts.get('ignore_obs_time', False)):
            frame_indices, time_vec_out, titles, obs_vec = read_in_observables(pears_obsf, time_vec)
        else:
            titles, obs_vec = read_in_observables_blind(pears_obsf)
        pears_obsf.close()
        iobs = 0
        for obs in obs_vec:
            iobs = iobs + 1
            pearson_results = dict()
            obs_pearsonfname = "aggregate/pearson_data_obs_%d.dat"%iobs
            obs_pearsonf = open(obs_pearsonfname, "w")
            minslope = 0.0
            maxslope = 0.0
            for i in range(nterm, nterm+nres):
                for j in range(i, nterm+nres):
                    if (i, j) in pairs_legend:
                        vecij = vectors[pairs_legend[i, j]]
                        if (opts.get('ignore_obs_time', False)):
                            rvec = vecij
                        else:
                            rvec = vecij[frame_indices]
                        pearsonerr = (len(rvec) != len(obs))
                        if (pearsonerr):
                            print("Pearson error! Wrong number of frames."," len (rvec) = %i len(obs) = %i"%(len(rvec), len(obs)) )
                            break
                        slope, intercept, pears_r, pvalue, std_err = linregress(rvec, obs)
                    else:
                        pears_r = 0.0
                        slope = 0
                        std_err = 0
                    pearson_results[i, j]   = (pears_r, slope, std_err)
                if (pearsonerr):
                    break
            if (not pearsonerr):
                for i in range(nterm, nterm+nres):
                    for j in range(nterm, nterm+nres):
                        if (i < j):
                            pair = (i, j)
                        else:
                            pair = (j, i)
                        obs_pearsonf.write("%4d %4d %10.6f %10.6f %10.6f\n"%((i,j) + pearson_results[pair]))
                    obs_pearsonf.write("\n")
                obs_pearsonf.close()
                if (gnus_path): os.system('gnuplot -e "domains=%d;maxz=1.0;inputfile='%(domains)+"'%s'"%(obs_pearsonfname)
                      +";outputfile='aggregate/pearson_obs_%d.png';title_time='Pearson correlation with %s for pairwise distances'"%(iobs,titles[iobs-1])+'" %s/script_corr.gnu'%(gnus_path))
    
    if ('pca_make' in opts):
        os.system("mkdir pca")
        print("Starting principal component analysis...")
        pca_file = open("pca/pca_output.txt", "w")
        vectors_mean0 = np.zeros([0, nframes], dtype = np.float16)
        print("There are a total of %i non-constant pairs."%len(vectors))
        print("Removing means....")
        for vec in vectors:
            meanvec = np.full(len(vec), np.average(vec), dtype = np.float16)
            vectors_mean0 = np.vstack([vectors_mean0, np.array(vec) - meanvec])
        print("Building covariance matrix....")
        covmat = np.cov(vectors)
        print("Done.")
        sum_var = 0.0
        for i in range(len(pairs_list)):
            sum_var = sum_var + covmat[i, i]
        pca_file.write("# Total variance: %f\n"%sum_var)
        pca_file.write("# i (eigvalue i) (running sum to i) (explained variance)\n")
        print("Diagonalizing....")
        w, v = eigsh(covmat, k=pca_make)
        w = w[::-1]
        v = np.fliplr(v)
        partial_var = 0.0
        ordered_v = project_pc(v, time_vec, vectors_mean0, pairs_list)
        for i in range(pca_make):
            partial_var += w[i]
            pca_file.write("%3i %12.6f %12.6f %12.6f\n"%(i+1, w[i], partial_var, partial_var/sum_var))
            plot_pc(xterm, xres, yterm, yres, i, ordered_v[:, i], pairs_list, pairs_legend, gnus_path, domains)
        pca_file.close()
        if (gnus_path): os.system("""gnuplot %s/script_pca_eigs.gnu """%(gnus_path))
    
    if ('pca_dir' in opts):
        os.system("mkdir pca")
        vectors_mean0 = []
        print("There are a total of %i non-constant pairs."%len(vectors))
        for vec in vectors:
            meanvec = np.full(len(vec), np.average(vec), dtype = np.float16)
            vectors_mean0.append(np.array(vec) - meanvec)
        read_pca(pca_dir, xterm, xres, yterm, yres, time_vec, asymm, vectors_mean0, pairs_list)

    if ('blocking' in opts):
        os.system("mkdir blocking")
        print("Performing blocking analysis...")
        vecarray = np.transpose (np.array(vectors))
        block_analysis(vecarray, nres, dt)
    
    if ('k_res_clusters' in opts):
        os.system('mkdir -v cluster_res')
        print("Attempting k-medoids based on interaction lifetimes...")
        os.system('mkdir -v cluster_res/lifetime')
        life_cond = get_condensed_distances(life_metric_dict0, 0, xres)
        call_hierarchical("cluster_res/lifetime", vectors, trunc, trunc_inter, life_cond, xres, ks_res, "lifetime", None, 'res', xterm, xres, yterm, yres, asymm, matrices)
        print("Attempting k-medoids based on simple distances...")
        os.system('mkdir -v cluster_res/distance')
        dist_cond = get_condensed_distances(avg_metric_dict0, 0, xres)
        call_hierarchical("cluster_res/distance", vectors, trunc, trunc_inter, dist_cond, xres, ks_res, "distance", None, 'res', xterm, xres, yterm, yres, asymm, matrices)
        if (opts.get('pearson_inter', False)):
            print("Attempting k-medoids based on inter-residue cross-correlation...")
            os.system('mkdir -v cluster_res/crosscorr')
            crosscorr_cond = get_condensed_distances(density_r2_metric, 0, xres)
            call_hierarchical("cluster_res/crosscorr", vectors, trunc, trunc_inter, crosscorr_cond, xres, ks_res, "distance", None, 'res', xterm, xres, yterm, yres, asymm, matrices)
    
    if ('k_traj_clusters' in opts):
        print("Starting contact map-based cluster analysis")
        print("")
        os.system("mkdir -p cluster_trj")
        if (asymm):
            norm = np.sqrt(xres*yres)
        else:
            norm = np.sqrt(2)*nres
        if (dimer):
            print("Dimer mode is turned on. The cluster analysis could be unbearably slow.")
            npad, map = prepare_dimer_transpose(pairs_list, pairs_legend)
            padvec = trunc * np.ones(npad)
            rmsd_interframe_cond = np.zeros((nframes*(nframes-1)/2,), dtype = np.float16)
            k = 0
            for i in range(nframes):
                print("t=%10i ps (%3.0f%% progress)"%(time_vec[i], 100.0*(k/((nframes*(nframes-1)/2) )) ), end="\r")
                veci = np.hstack([vectors[:,i], padvec])
                for j in range(i+1, nframes):
                    vecj = vectors[:, j]
                    vecj = np.hstack([vecj, padvec])
                    vecj_t = vecj[map]
                    d = np.linalg.norm(veci/norm - vecj/norm)
                    d_t = np.linalg.norm(veci/norm - vecj_t/norm)
                    if (d_t < d):
                        d = d_t
                    rmsd_interframe_cond[k] = d
                    k += 1
        else:
           print("Building distance matrix...")
           rmsd_interframe_cond = pdist(vectors.transpose()/norm)
           print("Done.")
        
        rmsd_outputf = open("cluster_trj/rmsd_interframe.dat","w")
        print("")
        if (nframes > 1000):
            delta = int(nframes/1000)
            print("There are more than 1000 frames, we will limit the interframe plots to approx. 1000x1000.")
        else:
            delta = 1
        for i in range(0, nframes, delta):
            time_i = time_vec[i] * 0.001
            for j in range(0, nframes, delta):
                time_j = time_vec[j] * 0.001
                if (i == j):
                    rmsd_interframe = 0
                elif (i < j):
                    rmsd_interframe = rmsd_interframe_cond[cond_index(i, j, nframes)]
                else:
                    rmsd_interframe = rmsd_interframe_cond[cond_index(j, i, nframes)]
                rmsd_outputf.write("%8.3f %8.3f %12.6f\n"%(time_i, time_j, rmsd_interframe ))
            rmsd_outputf.write("\n")
        rmsd_outputf.close()
        if (gnus_path): os.system("""gnuplot -e "domains=0;maxz=0;label_str='Time (ns)';inputfile='cluster_trj/rmsd_interframe.dat';outputfile='cluster_trj/rmsd_interframe.png';""" +
                      """title_str='Interframe RMSD';cblabel_str = 'RMSD (nm)" %s/script_single.gnu"""%(gnus_path))
        print("hierarchical...")
        os.system('mkdir -p cluster_trj')
        call_hierarchical('cluster_trj', vectors, trunc, trunc_inter, rmsd_interframe_cond, nframes, ks_traj, 'Trajectory', time_vec, 'trj', xterm, xres, yterm, yres, asymm, matrices)
    