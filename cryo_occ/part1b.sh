#mkdir './02_mask'
#mkdir './03_cut'
#for file in ./01_PDB_mrc/PDB*; do 
#	suffix=${file##*/}; 
#	prefix="Bi_"; 
#	newfile="./02_mask/$prefix$suffix";
#	relion_mask_create --i "$file" --o "$newfile";
#done


####CHANGE ME TO THE CORRECT .MAP OR .MRC####
for file in ./00_source_maps/B*_scaled_aligned.mrc; do 
#############################################
	prefix=${file##*/};
	for BiPDB in ./02_mask/Bi_PDB*; do 
		suffix=${BiPDB##*/}; 
		newfile="./03_cut/Cut_"$prefix"_"$suffix; 
		proc3d "$file" "$newfile" maskfile="$BiPDB" &
	done; 
done
