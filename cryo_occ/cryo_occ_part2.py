###============================
### Chimera Occupancy Script Part 2
### Version 1.0
### 12 January 2016
### By Yong Zi Tan & Joey Davis <joeydavis@gmail.com>
### 12 July 2016
###============================
 
import os
from glob import glob

##############CHANGE ME BASED ON MODELS TO INSPECT#####################################
map_files = ['B_scaled_aligned.mrc', 'C_scaled_aligned.mrc', 'C1_scaled_aligned.mrc', 'C2_scaled_aligned.mrc', 'C3_scaled_aligned.mrc',
		'D_scaled_aligned.mrc', 'D1_scaled_aligned.mrc', 'D2_scaled_aligned.mrc', 'D3_scaled_aligned.mrc', 'D4_scaled_aligned.mrc',
		'E_scaled_aligned.mrc', 'E1_scaled_aligned.mrc', 'E2_scaled_aligned.mrc', 'E3_scaled_aligned.mrc', 'E4_scaled_aligned.mrc', 'E5_scaled_aligned.mrc']
ref_thresh = '0.02'
##############CHANGE ME BASED ON MODELS TO INSPECT#####################################

#define chains
chains = ["uL2","uL3","uL4","uL5","uL6","bL9","uL10","uL11","uL13","uL14","uL15","uL16",
"bL17","uL18","bL19","bL20","bL21","uL22","uL23","uL24","bL25","bL27","bL28","uL29",
"uL30","bL32","bL33","bL34","bL35","bL36","h1","h2","h3","h4","h5","h6","h7","h8","h9",
"h10","h11","h12","h13","h14","h16","h18","h19","h20","h21","h22","h23","h24","h25","h26",
"h27","h28","h29","h31","h32","h33","h34","h35","h36","h37","h38","h39","h40","h41","h42",
"h43","h44","h45","h46","h47","h48","h49","h50","h51","h52","h53","h54","h55","h56","h57",
"h58","h59","h60","h61","h62","h63","h64","h65","h66","h67","h68","h69","h70","h71","h72",
"h73","h74","h75","h76","h77","h78","h79","h80","h81","h82","h83","h84","h85","h86","h87",
"h88","h89","h90","h91","h92","h93","h94","h95","h96","h97","h98","h99","h100","h101",
"h102","h103","h104","h105","h106","h107","h108","h109","h110","h111","h112"]

ref_name = 'PDB_'+str(ref_thresh)#define ref name
cutlist_ref = glob('./05_cutlist/*REF_'+str(ref_thresh)+'.txt')[0]
print cutlist_ref

all_lists = glob('./05_cutlist/cut_volList_*.txt') #get all thresholds

for cutlist_thr in all_lists:# do the following for all thresholds
	threshold = cutlist_thr.split('_volList_')[1].split('.tx')[0] #get numeric threshold
	with open(cutlist_thr) as f: #get list of all volume files
	        volumes_files = [x.strip() for x in f.readlines()]
	with open(cutlist_ref) as f: #get list of all reference files
	        volumes_ref = [x.strip() for x in f.readlines()]

	output = open('calculated_volumes_thr_'+str(threshold)+'.txt',"w") #open output file

	#define vol dict, write model header
	vol_dict = {}
	output.write("#")
	for map_file in map_files:
        	vol_dict[map_file] = {}
        	output.write("\t%s" % map_file)
	vol_dict[ref_name] = {}
	output.write("\t%s" % ref_name)
	output.write("\n")

	#Pull out volumes from map files
	for vol_log in volumes_files:
        	element_name = vol_log.split("PDB_")[1].split('.mrc')[0]
		map_name = vol_log.split('_Bi_')[0].split('ut_')[1]
        	with open(vol_log) as f:
        	       lines = f.readlines()
        	data_line=lines[4].split()
        	vol_dict[map_name][element_name] = data_line[1]
	
	#Pull out volumes from reference
	for ref_vol_log in volumes_ref:
		element_name = ref_vol_log.split("PDB_")[1].split('.mrc')[0]
        	with open(ref_vol_log) as f:
        	       lines = f.readlines()
        	data_line=lines[4].split()
        	vol_dict[ref_name][element_name] = data_line[1]

	#Write out data for each chain
	for chain in chains:
      		output.write(chain)
      		for map_file in map_files:
			output.write("\t%s" % vol_dict[map_file][chain]) #write each experimental vol
		output.write("\t%s" % vol_dict[ref_name][chain]) #write the reference vol
		output.write("\n")
	output.close()
print 'completed part 2'
