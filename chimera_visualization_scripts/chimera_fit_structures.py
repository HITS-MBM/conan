import os
import chimera
from chimera import runCommand as rc
from chimera import replyobj

# Get filenames and iterate through them

conformers = [filename for filename in os.listdir('.') if filename.endswith('.gro')]
conformers_total = len(conformers)

# Get name of the reference structure
ref='cluster_center.gro'
rc("open "+ref) # This is our model #0

for filename in conformers:
	replyobj.status("Processing " + filename) # show what file is being worked on
	rc("open "+filename)
rc("mm #0 #1-%d" %(conformers_total)) # This fits all the structures found in the directory to the ref
