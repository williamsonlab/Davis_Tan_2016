###============================
### Chimera Occupancy Script Part 1
### Version 1.0
### 12 January 2016
### By Yong Zi Tan & Joey Davis <joeydavis@gmail.com>
###============================

import os
from chimera import runCommand as rc # use 'rc' as shorthand for runCommand
from chimera import replyobj # for emitting status messages

###============================
### Parameters to Change
res = 6.50 # Resolution to filter the maps to
PDB1 = "./00_aligned_pdbs/proteinAnd5S.pdb" # Make sure PDB is aligned to ReferenceMap
PDB2 = "./00_aligned_pdbs/rnaSet1.pdb"
PDB3 = "./00_aligned_pdbs/rnaSet2.pdb"
ReferenceMap = "./00_source_maps/B_scaled_aligned.mrc" # One of the maps you would want to determine occupancy of

os.mkdir('./01_PDB_mrc')
###============================

chain_proteinAnd5S = ["u","v","w","y","x","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z","a","b","c","d","e","f"]
name_proteinAnd5S = ["h108","h109","h110","h111","h112","uL2","uL3","uL4","uL5","uL6","bL9","uL10","uL11","uL13","uL14","uL15","uL16","bL17","uL18","bL19","bL20","bL21","uL22","uL23","uL24","bL25","bL27","bL28","uL29","uL30","bL32","bL33","bL34","bL35","bL36"]
chain_rnaSet1 = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","P","R","S","T","U","V","W","X","Y","Z","a","b","c","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z","1","2","3","4","5","6","7","8","9"]
name_rnaSet1 = ["h1","h2","h3","h4","h5","h6","h7","h8","h9","h10","h11","h12","h13","h14","h16","h18","h19","h20","h21","h22","h23","h24","h25","h26","h27","h28","h29","h31","h32","h33","h34","h35","h36","h37","h38","h39","h40","h41","h42","h43","h44","h45","h46","h47","h48","h49","h50","h51","h52","h53","h54","h55","h56","h57","h58","h59","h60","h61"]
chain_rnaSet2 = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z","a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t"]
name_rnaSet2 = ["h62","h63","h64","h65","h66","h67","h68","h69","h70","h71","h72","h73","h74","h75","h76","h77","h78","h79","h80","h81","h82","h83","h84","h85","h86","h87","h88","h89","h90","h91","h92","h93","h94","h95","h96","h97","h98","h99","h100","h101","h102","h103","h104","h105","h106","h107"]

### Part 1: Using PDB model, make the 3D volumes of each chain, save them
### PDB1
rc("open " + ReferenceMap)
rc("open " + PDB1)
num1 = 1000
num2 = 2000
x = 0
for i in chain_proteinAnd5S:
	print "Processing %s" % i
	num1 += 1
	num2 += 1
	rc("sel :." + i) # select chain
	rc("molmap sel " + str(res) + " modelId " + str(num1)) # make chain into density
	rc("vop #" + str(num1) + " resample onGrid #0 modelID " + str(num2)) # make sure the mrc file has the same box size as the reference map
	rc("volume #" + str(num2) + " save ./01_PDB_mrc/PDB_" + name_proteinAnd5S[x] + ".mrc") # save map
	x += 1
rc("close all")

### PDB2
rc("open " + ReferenceMap)
rc("open " + PDB2)
num1 = 1000
num2 = 2000
x = 0
for i in chain_rnaSet1:
	print "Processing %s" % i
	num1 += 1
	num2 += 1
	rc("sel :." + i) # select chain
	rc("molmap sel " + str(res) + " modelId " + str(num1)) # make chain into density
	rc("vop #" + str(num1) + " resample onGrid #0 modelID " + str(num2)) # make sure the mrc file has the same box size as the reference map
	rc("volume #" + str(num2) + " save ./01_PDB_mrc/PDB_" + name_rnaSet1[x] + ".mrc") # save map
	x += 1
rc("close all")

### PDB3
rc("open " + ReferenceMap)
rc("open " + PDB3)
num1 = 1000
num2 = 2000
x = 0
for i in chain_rnaSet2:
	print "Processing %s" % i
	num1 += 1
	num2 += 1
	rc("sel :." + i) # select chain
	rc("molmap sel " + str(res) + " modelId " + str(num1)) # make chain into density
	rc("vop #" + str(num1) + " resample onGrid #0 modelID " + str(num2)) # make sure the mrc file has the same box size as the reference map
	rc("volume #" + str(num2) + " save ./01_PDB_mrc/PDB_" + name_rnaSet2[x] + ".mrc") # save map
	x += 1
rc("stop now")
