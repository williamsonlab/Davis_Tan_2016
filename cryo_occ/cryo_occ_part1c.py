###============================
### Chimera Occupancy Script Part 1c
### Version 1.0
### 12 January 2016
### By Yong Zi Tan & Joey Davis <joeydavis@gmail.com>
### 12 July 2016
###============================

from glob import glob
import os

######CHANGE ME TO THE APPROPRIATE THRESHOLD RANGES####
min=0.02
max=0.05
threshold_ref=0.02
TEST=False #true will just output commands, false will perform them
pixel_den = 1.31
#######################################################
num_steps=10
step_size = (max-min)/num_steps
all_files = glob('./03_cut/Cut_*')
command = 'mkdir 05_cutlist'
print command if TEST else os.system(command)

#iterate through thresholds
for i in range(int(num_steps)):
	threshold=min+i*step_size

	#make dir for each threshold (will hold vols for all maps, and a reference)
	new_dir = './04_vol_'+str(threshold) 
	command = 'mkdir '+new_dir
	print command if TEST else os.system('mkdir '+new_dir)

	#iterate through each masked density
	for f in all_files:
		new_file_name = new_dir+'/'+f.split('/')[-1]+'_vol'+str(threshold)+'.log' 
		
		#make log for each masked density (one log for each map file for each threshold)
		command = 'volume ' + f + ' ' + str(pixel_den) + ' thr='+str(threshold)+' > '+new_file_name
		print command if TEST else os.system(command)

	#make a listing of the log files for each theshold, separated by map/ref
	command = 'ls '+new_dir+'/*.log > '+'./05_cutlist/cut_volList_'+str(threshold)+'.txt' 
	print command if TEST else os.system(command)

#Determine the reference volume at desired threshold
ref_vol_dir = './04_volREF_'+str(threshold_ref)
command = 'mkdir ' + ref_vol_dir
print command if TEST else os.system(command)

for f in glob('./01_PDB_mrc/PDB_*.mrc'):
	ref_file_name = ref_vol_dir+'/'+f.split('/')[-1]+'_vol_'+str(threshold_ref)+'.log'
	command = 'volume ' + f + ' ' + str(pixel_den) + ' thr='+str(threshold_ref)+' > ' + ref_file_name
	print command if TEST else os.system(command)
command = 'ls '+ref_vol_dir+'/*.log > '+'./05_cutlist/cut_volListREF_'+str(threshold_ref)+'.txt' 
print command if TEST else os.system(command)

print 'completed part1C'

