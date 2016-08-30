
###cryo_occ consists of a series of scripts written by Joey Davis <joeydavis@gmail.com> and Yong-Zi Tan. They consist of the following:

1) cryo_occ_part0_rename_pdb.ipynb - a python notebook used to rename chains for the 50S subunit such that each element has a unique chain ID. The required input and resulting output files are also included

2) cryo_occ_part1a_chimera.py - a python script to be run in chimera. Used to generate 3D volumes for each element in the PDB file

3) cryo_occ_part1b.sh - a shell script used to generate masks

4) cryo_occ_part1c.py - a python script used to calculate volumes occupied as a function of the threshold chosen.

5) cryo_occ_part2.py - a python script used to generate tables of occupancy as a function of occupancy

6) cryo_occ_plot.ipynb - a python notebook used to plot occupancy vs. threshold curves for each element and each model analyzed.

7) cryo_occ_part0_diffmap.com - used to run diffmap to amplitude scale initial maps.

###Usage
1.	Amplitude scale maps to the lowest (worst) resolution map
    a.	  Get pixels/angstrom (open source maps in chimera – volume viewer – coordinates)
    b.	Edit diffmap.com for the relevant files and scale to the lowest res map you have
    c.	Repeat for each map (including the initial one)
2.	Align maps
    a.	Open maps
    b.	Fit in map (align to the same one you amplitude scaled to)
    c.	If there are negative values, use vop scale shift [amount[
    d.	vop #1 resample onGrid #0
    e.	Volume viewer to save as _scaled_aligned.mrc
3.	Align pdb file
    a.	Open 4ybb_grouped.py
    b.	Open scaled_aligned maps
    c.	Fit in map (align to the most complete structure)
    d.	Save pdb (check relative to and choose the base map from above) - $name_aligned_[base model]
4.	Make pdb reference files
    a.	Edit part1a_chimera.py to update the reference map, resolution, etc.
    b.	chimera --nogui --nostatus --script part1a_chimera.py
    c.	take uL2 and load in chimera – make sure it aligns and figure out a reasonable threshold for the upcoming comparisons
5.	Make mask files
    a.	Edit part1b.sh to update to .mrc or .map file
    b.	./part1b.sh
6.	Make cut files
    a.	Change thresholds to reasonable ranges based on scaled, aligned maps
    b.	Change the reference threshold [default is 0.05]
    c.	Change steps to 10
    d.	python part1c.py
7.	Save the tabulated occupancies
    a.	Change the map files to point to the scaled and aligned maps
    b.	Change the reference threshold 
    c.	./part2.py

