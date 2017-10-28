import os
import re
import sys
import colorsys
from codecs import encode, decode
import chimera
from chimera import runCommand as rc
from chimera import replyobj

# Function to generate HexCol that chimera can use

def get_N_HexCol(N=5):
	HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in range(N)]
	hex_out = []
	for rgb in HSV_tuples:
		rgb = map(lambda x: int(x*255),colorsys.hsv_to_rgb(*rgb))
		hex_out.append("".join(map(lambda x: chr(x).encode('hex'),rgb)))
	return hex_out

# Now open the file and get the number of cluster

fileinp = sys.argv[2]
with open(fileinp, 'r') as fp:
	for line in fp:
		if line.startswith('Total'):
			nclu = int(re.findall('\d+', line)[0])

# Set background as white
rc("set bgColor white")

colors = get_N_HexCol(nclu)
count=0

with open(fileinp, 'r') as fp:
	rc("open "+sys.argv[3]) # Open the pdb file
	rc("~disp #0")
	for line in fp:
		if re.match('Cluster \d+ has \d+ elements with center at resid \d+', line):
			center = int(re.findall(r'\d+$', line)[0]) # Force them into integer because that's what chimera wants!
			rc("represent sphere #0:%d" % (center))
			rc("disp #0:%d" % (center))
		if re.match(r'Cluster\d+ residues', line):
			ranges  = re.findall(r'\d+-\d+', line)
			singles = map(str.strip, re.findall(r'[ ]\d+[ ,\]]', line)) # Find single residues and remove spaces
			col=colors[count] # get the first colour out of the list that the function has created

			# Depict ranges of residues in a single cluster
			for i in ranges:
				first=int(re.findall('\d+', i)[0]) # Force them into integer because that's what chimera wants!
				last=int(re.findall('\d+', i)[1])  # Force them into integer because that's what chimera wants!
				rc("color #%s,r #0:%d-%d" % (col, first, last))
			count+=1
