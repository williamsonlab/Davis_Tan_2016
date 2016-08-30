#!/usr/bin/perl

# nuccyl
# ======
# 
# Copyright (c) 2003-2005 Luca Jovine (luca.jovine@mssm.edu). All rights reserved. 
# 
# Permission to use, copy, modify, and distribute this software and its
# documentation for any purpose and without fee is hereby granted, provided
# that the above copyright notice appear in all copies and that both that
# copyright notice and this permission notice appear in supporting
# documentation, and that the name of Luca Jovine not be used in advertising
# or publicity pertaining to distribution of the software without specific,
# written prior permission. 
# 
# LUCA JOVINE DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING
# ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL
# LUCA JOVINE BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR
# ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER
# IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT
# OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE. 

$nuccyl_version = "1.5.2";
$nuccyl_revision = "1";
$nuccyl_date = "050304";

require "ctime.pl";

$switch = $ARGV[0];
($purine_match, $pyrimidine_match, $standard_nt_match, $modified_nt_match, $all_nt_match) = &load_nt_data;
if ($switch eq "-cyl_1") {
	&nuccyl_cyl_1;
} elsif ($switch eq "-cyl_2d") {
	&nuccyl_cyl_2d;
} elsif ($switch eq "-cyl_2r") {
	&nuccyl_cyl_2r;
} elsif ($switch eq "-cyl_3") {
	&nuccyl_cyl_3;
} elsif ($switch eq "-filled") {
	&nuccyl_filled;
} else {
	print "\n nuccyl $nuccyl_version - $nuccyl_date";
	print "\n =====================\n";
	print "\n Usage:\n";
	print "\n\tnuccyl.pl -cyl_1  <input_pdb_coordinate_file> <output_rnaview_coordinate_file> <output_nuccyl3_coordinate_file>\n";
	print "\n\tnuccyl.pl -cyl_2d <input_3dna_base_pair_file> <input_nuccyl3_coordinate_file>\n";
	print "\n\tnuccyl.pl -cyl_2r <input_rnaview_base_pair_file> <input_nuccyl3_coordinate_file>\n";
	print "\n\tnuccyl.pl -cyl_3  <input_nuccyl2_file> <input_pdb_coordinate_file> <fatom_number_start> <fatom> <fres> <fres_number_start> <object_option> [obj/no_obj]\n";
	print "\n\tnuccyl.pl -filled <input_pdb_coordinate_file> <output_pdb_coordinate_file> <output_pml_file> <object_option> [obj/one_obj] <pur_rgb_r> <pur_rgb_g> <pur_rgb_b> (<pyr_rgb_r> <pyr_rgb_g> <pyr_rgb_b>)\n\n";
}

sub bp_get_coords {
	$x = $y = $z = "";
	my ($chain,$type,$base,$atom) = @_ ;
	open (INPUT_PDB_FILE, "$input_pdb_file");
	$line = <INPUT_PDB_FILE>;
	while ($line ne "") {
		if (($line =~ /^ATOM\s{0,6}\d+\s+$atom\s+$type\s+$chain\s{0,3}$base\s+(\-{0,1}\d+\.\d+)\s+(\-{0,1}\d+\.\d+)\s+(\-{0,1}\d+\.\d+)/)||($line =~ /^HETATM\s+\d+\s+$atom\s+$type\s+$chain\s+$base\s+(\-{0,1}\d+\.\d+)\s+(\-{0,1}\d+\.\d+)\s+(\-{0,1}\d+\.\d+)/)){
			($x,$y,$z) = ($1,$2,$3);
			return ($x,$y,$z);
			last;
		}
		$line = <INPUT_PDB_FILE>;
	}
}

sub bp_get_coords_no_type {
	$x = $y = $z = "";
	my ($chain,$base,$atom) = @_ ;
	open (INPUT_PDB_FILE, "$input_pdb_file");
	$line = <INPUT_PDB_FILE>;
	while ($line ne "") {
		if (($line =~ /^ATOM\s{0,6}\d+\s+$atom\s+\w{1,4}\s+$chain\s{0,3}$base\s+(\-{0,1}\d+\.\d+)\s+(\-{0,1}\d+\.\d+)\s+(\-{0,1}\d+\.\d+)/)||($line =~ /^HETATM\s+\d+\s+$atom\s+\w{1,4}\s+$chain\s+$base\s+(\-{0,1}\d+\.\d+)\s+(\-{0,1}\d+\.\d+)\s+(\-{0,1}\d+\.\d+)/)){
			($x,$y,$z) = ($1,$2,$3);
			return ($x,$y,$z);
			last;
		}
		$line = <INPUT_PDB_FILE>;
	}
}

sub bp_intro {
	@_ = ($bp_number,$chain_1,$type_1,$base_1,$atom_1,$chain_2,$type_2,$base_2,$atom_2);
	print "\n\n Base pair \# $bp_number\n ------------------------------------------------------------------------------\n";
	print " nucleotide 1: $chain_1 $type_1 $base_1\n";
	print " nucleotide 2: $chain_2 $type_2 $base_2\n";
}

sub build_bp_list {
	local($in_file) = @_;
	$bp_list = "";
	open (IN_FILE, "$in_file");
	while (<IN_FILE>) {
		$bp_list .= $_;
	}
	close (IN_FILE);
	$bp_list;
}

sub build_nt_list {
	local($in_file) = @_;
	$nt_list="";
	open (IN_FILE, "$in_file");
	while (<IN_FILE>) {
		$line = $_;
		($pdb_atom_number,$pdb_atom_name,$pdb_residue_name,$pdb_chain,$pdb_residue_number,$pdb_x,$pdb_y,$pdb_z,$pdb_b_factor,$pdb_occ) = &pdb_scan($line);
 		if (($pdb_atom_number)&&($pdb_residue_name =~ /^(($all_nt_match))$/)) {
			unless ($nt_list =~ /$pdb_chain $pdb_residue_name $pdb_residue_number/) {
				$nt_id = "   " . $pdb_chain . " " . $pdb_residue_name . " " . $pdb_residue_number . "   ";
				$nt_list .= $nt_id;
			}
		}
	}
	@nt_list = split /   /, $nt_list;
	close (IN_FILE);
	@nt_list;
}

sub calculate_run_time {
	#
	local($time_start,$time_end)=@_;
	$time_difference=$time_end-$time_start;
	$seconds_temp=$time_difference;

	$minutes_temp=int($seconds_temp/60);
	$seconds_final=$seconds_temp-($minutes_temp*60);
	$hours_temp=int($minutes_temp/60);
	$minutes_final=$minutes_temp-($hours_temp*60);
	$days_temp=int($hours_temp/24);
	$hours_final=$hours_temp-($days_temp*24);
	$days_final=$days_temp;

	if($days_final==1){$ind_day="day"}else{$ind_day="days"};
	if($hours_final==1){$ind_hour="hour"}else{$ind_hour="hours"};
	if($minutes_final==1){$ind_min="minute"}else{$ind_min="minutes"};
	if($seconds_final==1){$ind_sec="second"}else{$ind_sec="seconds"};

	if(($days_final==0)&&($hours_final==0)&&($minutes_final==0)){
		print " Run time: $seconds_final $ind_sec\n";
	}elsif(($days_final==0)&&($hours_final==0)){
		print " Run time: $minutes_final $ind_min and $seconds_final $ind_sec\n";
	}elsif($days_final==0){
		print " Run time: $hours_final $ind_hour, $minutes_final $ind_min and $seconds_final $ind_sec\n";
	}else{
		print " Run time: $days_final $ind_day, $hours_final $ind_hour, $minutes_final $ind_min and $seconds_final $ind_sec\n";
	}
}

sub coordinate_middlepoint {
	#
	$x1 = $y1 = $z1 = $x2 = $y2 = $z2 = $x3_rounded = $y3_rounded = $z3_rounded = "";
	my ($x1,$y1,$z1,$x2,$y2,$z2) = @_ ;
	$x3 = ($x1+$x2)/2;
	$y3 = ($y1+$y2)/2;
	$z3 = ($z1+$z2)/2;
	$x3_rounded = sprintf("%.3f", $x3);
	$y3_rounded = sprintf("%.3f", $y3);
	$z3_rounded = sprintf("%.3f", $z3);
	return ($x3_rounded,$y3_rounded,$z3_rounded);
}

sub date_end {
	#
	$date_end=&ctime(time);
	$time_end=time();
	chop($date_end);
	$date_end=~s/BST //;
	$date_end=~s/  / /;
	print "\n $date_end: run finished\n";
	&calculate_run_time($time_start,$time_end);
	print "\n";
}

sub date_start {
	#
	$date_start = &ctime(time);
	$time_start = time();
	chop ($date_start);
	$date_start =~ s/BST //;
	$date_start =~ s/  / /;
	print "\n $date_start: run started\n";
}

sub load_modified_purine_data {
	#
	%modified_purine_replacements = (
		"  E" => "  A",
		"  R" => "  A",
		"  Y" => "  A",
		"1MA" => "  A",
		"2AR" => "  A",
		"6HA" => "  A",
		"12A" => "  A",
		"A23" => "  A",
		"AAM" => "  A",
		"ABM" => "  A",
		"ABR" => "  A",
		"ABS" => "  A",
		"ACP" => "  A",
		"ADN" => "  A",
		"ADP" => "  A",
		"AET" => "  A",
		"AMP" => "  A",
		"AP2" => "  A",
		"MA6" => "  A",
		"RMP" => "  A",
		"SMP" => "  A",
		"SRA" => "  A",
		"  P" => "  G",
		"  X" => "  G",
		" IG" => "  G",
		" LG" => "  G",
		" YG" => "  G",
		"1MG" => "  G",
		"2GP" => "  G",
		"2MG" => "  G",
		"3GP" => "  G",
		"5GP" => "  G",
		"6HG" => "  G",
		"7MG" => "  G",
		"8MG" => "  G",
		"8OG" => "  G",
		"BGM" => "  G",
		"DFG" => "  G",
		"DGT" => "  G",
		"FAG" => "  G",
		"G2S" => "  G",
		"G25" => "  G",
		"G7M" => "  G",
		"GDP" => "  G",
		"GMP" => "  G",
		"GMS" => "  G",
		"GN7" => "  G",
		"GSR" => "  G",
		"GSS" => "  G",
		"IGU" => "  G",
		"LGP" => "  G",
		"M2G" => "  G",
		"OMG" => "  G",
		"PGP" => "  G",
		"S6G" => "  G",
		"SGP" => "  G",
		"YYG" => "  G",
	);
	#$modified_purine_replacement_number = keys %modified_purine_replacements;
	@modified_purine_replacements = keys %modified_purine_replacements;
	#
	%modified_purines = %modified_purine_replacements;
	foreach $value (sort values %modified_purines) {
		$value =~ s/^\s+//;
	}
	#$modified_purine_number = keys %modified_purines;
	@modified_purines = values %modified_purines;
	#
	foreach $nt (@modified_purine_replacements) {
		$modified_purine_match = $modified_purine_match . "$nt|";
	}
	chop ($modified_purine_match);
	return $modified_purine_match;
}

sub load_modified_pyrimidine_data {
	#
	%modified_pyrimidine_replacements = (
		" CH" => "  C",
		" IC" => "  C",
		" LC" => "  C",
		"5CM" => "  C",
		"5IC" => "  C",
		"5MC" => "  C",
		"6HC" => "  C",
		"55C" => "  C",
		"C2P" => "  C",
		"C2S" => "  C",
		"C3P" => "  C",
		"C25" => "  C",
		"CAR" => "  C",
		"CBR" => "  C",
		"CCC" => "  C",
		"CDP" => "  C",
		"CMP" => "  C",
		"CMR" => "  C",
		"DCZ" => "  C",
		"DFC" => "  C",
		"DNR" => "  C",
		"EDC" => "  C",
		"GCK" => "  C",
		"I5C" => "  C",
		"IMC" => "  C",
		"MCY" => "  C",
		"OMC" => "  C",
		"5AT" => "  T",
		"6CT" => "  T",
		"6HT" => "  T",
		"ADT" => "  T",
		"ATD" => "  T",
		"DRT" => "  T",
		"MTR" => "  T",
		"SMT" => "  T",
		"T2S" => "  T",
		"T23" => "  T",
		"TAF" => "  T",
		"TDP" => "  T",
		"THX" => "  T",
		"TLB" => "  T",
		"TLC" => "  T",
		"TMP" => "  T",
		"TSP" => "  T",
		"2MU" => "  U",
		"3ME" => "  U",
		"4SU" => "  U",
		"5BU" => "  U",
		"5IU" => "  U",
		"5MD" => "  U",
		"5MU" => "  U",
		"70U" => "  U",
		"125" => "  U",
		"126" => "  U",
		"127" => "  U",
		"ADU" => "  U",
		"BRU" => "  U",
		"DHU" => "  U",
		"GMU" => "  U",
		"H2U" => "  U",
		"HEU" => "  U",
		"LHU" => "  U",
		"PSU" => "  U",
		"SSU" => "  U",
		"U2P" => "  U",
		"U3P" => "  U",
		"U5P" => "  U",
		"U8U" => "  U",
		"U25" => "  U",
		"UDP" => "  U",
		"UM3" => "  U",
		"UMP" => "  U",
		"UMS" => "  U",
	);
	#$modified_pyrimidine_replacement_number = keys %modified_pyrimidine_replacements;
	@modified_pyrimidine_replacements = keys %modified_pyrimidine_replacements;
	#
	%modified_pyrimidines = %modified_pyrimidine_replacements;
	foreach $value (sort values %modified_pyrimidines) {
		$value =~ s/^\s+//;
	}
	#$modified_pyrimidine_number = keys %modified_pyrimidines;
	@modified_pyrimidines = values %modified_pyrimidines;
	#
	foreach $nt (@modified_pyrimidine_replacements) {
		$modified_pyrimidine_match = $modified_pyrimidine_match . "$nt|";
	}
	chop ($modified_pyrimidine_match);
	return $modified_pyrimidine_match;
}

sub load_modified_unusual_nt_data {
	#
	%modified_unusual_nt_replacements = (
		"  N" => "XX1", # any 5'-monophosphate nucleotide [C5 H11 O7 P]
		" D3" => "XX2", # 1-(2-deoxy-beta-d-ribofuranosyl)-4-(3-benzamido) phenylimidazole [C21 H22 N3 O7 P1]
		"3DR" => "XX3",	# abasic dideoxyribose; 1p,2p-dideoxyribofuranose-5p-phosphate [C5 H11 O6 P1]
		"6MI" => "XX4", # 6-methyl-8-(2-deoxy-ribofuranosyl)isoxanthopteridine [C12 H16 N5 O8 P1]
		"D1P" => "XX5", # abasic deoxyribose; 2p-deoxy-ribofuranose-5p-phosphate [C5 H11 O7 P1]
		"DFT" => "XX6", # 1-(2-deoxyribofuranosyl)-2,4-difluoro-5-methyl- benzene-5pmonophosphate [C12 H15 O6 F2 P1]
		"DPY" => "XX7",	# 2-deoxyribofuranosyl-pyridine-2,6-dicarboxylic acid- 5p-monophosphate [C12 H14 N1 O10 P1]
		"DRP" => "XX8", # 2-deoxyribofuranosyl-pyridine-5p-monophosphate [C10 H14 N1 O6 P1]
		"FMP" => "XX9", # formycin-5p-monophosphate; formycin monophosphate [C10 H14 N5 O7 P1]
		"IRP" => "X10", # (1s)-1(9-deazahypoxanthin-9yl)1,4-dideoxy-1,4-imino-d- ribitol-5-phosphate [C11 H15 N4 O7 P1]
		"M1G" => "X11", # 3-(2-deoxy-beta-d-ribofuranosyl)-pyrido(5,6-a)-purine- 10-one-5p-monophosphate [C13 H14 N5 O7 P1]
		"MBZ" => "X12", # 1-(2-deoxyribofuranosyl)-4-methyl-benzoimidazole-5p- monophosphate [C13 H17 N2 O6 P1]
		"MEP" => "X13", # phosphoric acid mono-(4-methoxy-3-methyl-5-(5-methyl- 2,4-dioxo-3,4-dihydro-2h-pyrimidin-1-yl)-tetrahydro- furan-2-ylmethyl) ester [C12 H19 N2 O8 P1]
		"NP3" => "X14", # 1-(2-deoxy-ribofuranosyl)-1h-(3-nitro-pyrrol)-5p- phosphate [C9 H13 N2 O8 P1]
		"P5P" => "X15", # purine riboside-5p-monophosphate [C10 H13 N4 O7 P1]
		"PPU" => "X16", # puromycin-5p-monophosphate [C22 H30 N7 O8 P1]
		"PYP" => "X17", # 2p-deoxyribofuranosylpyrene-5p-monophosphate [C21 H19 O6 P1]
		"PYY" => "X18", # d-ribofuranosyl-benzene-5p-monophosphate [C11 H15 O7 P1]
		"T6A" => "X19", # n-(n-(9-b-d-ribofuranosylpurin-6-yl) carbamoyl)threonine-5p-monophosphate; n-(nebularin-6-ylcarbamoyl)-l-threonine-5p- monophosphate [C15 H20 N6 O11 P1]
		"TLN" => "X20", # ((1r,3r,4r,7s)-7-hydroxy-3-(thymin-1-yl)-2,5- dioxabicyclo(2.2.1)hept-1-yl)methyl dihydrogen phosphate dihydrogen phosphate [C11 H15 N2 O9 P1]
	);
	#$modified_unusual_nt_replacement_number = keys %modified_unusual_nt_replacements;
	@modified_unusual_nt_replacements = keys %modified_unusual_nt_replacements;
	#
	%modified_unusual_nts = %modified_unusual_nt_replacements;
	foreach $value (sort values %modified_unusual_nts) {
		$value =~ s/^\s+//;
	}
	#$modified_unusual_nt_number = keys %modified_unusual_nts;
	@modified_unusual_nts = values %modified_unusual_nts;
	#
	foreach $nt (@modified_unusual_nt_replacements) {
		$modified_unusual_nt_match = $modified_unusual_nt_match . "$nt|";
	}
	chop ($modified_unusual_nt_match);
	return $modified_unusual_nt_match;
}

sub load_nt_data {
	#
	$standard_purine_match = &load_standard_purine_data;
	$modified_purine_match = &load_modified_purine_data;
	$purine_match = "A|G|" . $standard_purine_match . "|" . $modified_purine_match;
	$purine_match =~ s/ //g;
	#
	$standard_pyrimidine_match = &load_standard_pyrimidine_data;
	$modified_pyrimidine_match = &load_modified_pyrimidine_data;
	$pyrimidine_match = "C|U|T|" . $standard_pyrimidine_match . "|" . $modified_pyrimidine_match;
	$pyrimidine_match =~ s/ //g;
	#
	%standard_nt_replacements = (%standard_purine_replacements, %standard_pyrimidine_replacements);
	%standard_nt = (%standard_purines, %standard_pyrimidines);
	#$standard_nt_number = keys %standard_nt;
	@standard_nt = keys %standard_nt;
	foreach $nt (@standard_nt) {
		$standard_nt_match = $standard_nt_match . "$nt|";
	}
	chop ($standard_nt_match);
	$standard_nt_match = "A|G|C|U|T|" . $standard_nt_match;
	$standard_nt_match =~ s/ //g;
	#
	%modified_nt_replacements = (%modified_purine_replacements, %modified_pyrimidine_replacements);
	%modified_nt = (%modified_purines, %modified_pyrimidines);
	#$modified_nt_number = keys %modified_nt;
	@modified_nt = keys %modified_nt;
	foreach $nt (@modified_nt) {
		$modified_nt_match = $modified_nt_match . "$nt|";
	}
	chop ($modified_nt_match);
	$modified_nt_match =~ s/ //g;
	#
	%all_nt_replacements = (%standard_nt_replacements, %modified_nt_replacements); 
	%all_nt = (%standard_nt, %modified_nt); 
	#$all_nt_number = keys %all_nt;
	@all_nt = keys %all_nt;
	foreach $nt (@all_nt) {
		$all_nt_match = $all_nt_match . "$nt|";
	}
	chop ($all_nt_match);
	$all_nt_match = "A|G|C|U|T|" . $all_nt_match;
	$all_nt_match =~ s/ //g;
	#
	$modified_unusual_nt_match = &load_modified_unusual_nt_data;
	#
	return ($purine_match, $pyrimidine_match, $standard_nt_match, $modified_nt_match, $all_nt_match);
}

sub load_standard_purine_data {
	#
	%standard_purine_replacements = (
		"ADE" => "  A", # adenosine-5p-monophosphate [C10 H14 N5 O7 P1]
		"GUA" => "  G", # guanosine-5p-monophosphate [C10 H14 N5 O8 P1]
	);
	#$standard_purine_replacement_number = keys %standard_purine_replacements;
	@standard_purine_replacements = keys %standard_purine_replacements;
	#
	%standard_purines = %standard_purine_replacements;
	foreach $value (sort values %standard_purines) {
		$value =~ s/^\s+//;
	}
	#$standard_purine_number = keys %standard_purines;
	@standard_purines = values %standard_purines;
	#
	foreach $nt (@standard_purine_replacements) {
		$standard_purine_match = $standard_purine_match . "$nt|";
	}
	chop ($standard_purine_match);
	return $standard_purine_match;
}

sub load_standard_pyrimidine_data {
	#
	%standard_pyrimidine_replacements = (
		"CYT" => "  C", # cytidine-5p-monophosphate [C9 H14 N3 O8 P1]
		"THY" => "  T", # thymidine-5p-monophosphate [ C10 H15 N2 O8 P1]
		"URA" => "  U", # uridine-5p-monophosphate [C9 H13 N2 O9 P1]
		"URI" => "  U", # uridine-5p-monophosphate [C9 H13 N2 O9 P1]
	);
	#$standard_pyrimidine_replacement_number = keys %standard_pyrimidine_replacements;
	@standard_pyrimidine_replacements = keys %standard_pyrimidine_replacements;
	#
	%standard_pyrimidines = %standard_pyrimidine_replacements;
	foreach $value (sort values %standard_pyrimidines) {
		$value =~ s/^\s+//;
	}
	#$standard_pyrimidine_number = keys %standard_pyrimidines;
	@standard_pyrimidines = values %standard_pyrimidines;
	#
	foreach $nt (@standard_pyrimidine_replacements) {
		$standard_pyrimidine_match = $standard_pyrimidine_match . "$nt|";
	}
	chop ($standard_pyrimidine_match);
	return $standard_pyrimidine_match;
}

sub nuccyl_cyl_1 {
	#
	# Example:
	#
	#	./nuccyl.pl -cyl_1 1EVV.pdb 1EVV_rnaview.pdb 1EVV_nuccyl3.pdb | tee nuccyl1.log
	#
	PRINT_INTRO:;
	print "\n nuccyl_cyl_1 $nuccyl_version - $nuccyl_date";
	print "\n ===========================\n";
	&date_start;
	
	INPUT_VARIABLE_SETUP:;
	$input_pdb_file = $ARGV[1];
	$output_rnaview_file = $ARGV[2];
	$output_nuccyl3_file = $ARGV[3];
	
	IO_SYNTAX_CHECK:;
	$syntax_1 = " Usage: nuccyl.pl -cyl_1 <input_pdb_coordinate_file> <output_rnaview_coordinate_file> <output_nuccyl3_coordinate_file>";
	unless ($input_pdb_file && $output_rnaview_file && $output_nuccyl3_file) {
		print "\n$syntax_1\n\n";
		exit;
	}
	
	IO_FILE_CHECK:;
	open (INPUT_PDB_FILE, "$input_pdb_file") || die "\n ERROR: Can't open input PDB file \"$input_pdb_file\"!\n\n";
	print "\n Input PDB coordinate file: $input_pdb_file\n";
	close (INPUT_PDB_FILE);
	if (-e "$output_rnaview_file") {
		unlink $output_rnaview_file;
	}
	if (-e "$output_nuccyl3_file") {
		unlink $output_nuccyl3_file;
	}

	GENERATE_FILE_FOR_RNAVIEW:;
	print "\n Generated PDB coordinate file for RNAView: ";
	open (INPUT_PDB_FILE, "$input_pdb_file");
	open (OUTPUT_RNAVIEW_FILE, ">$output_rnaview_file");
	$line = <INPUT_PDB_FILE>;
	RNAVIEW_LINE: while ($line ne "") {
		if ($line =~ /^(ATOM|HETATM)/) {
			foreach $key (sort keys %all_nt_replacements) {
				if ($line =~ s/^(ATOM|HETATM)(\s{0,6}\d{1,7}\s{0,}[\w|\*]+\s{0,})$key(\s{0,}\w\s{0,3}\w+\s{0,}\-{0,1}\d+\.\d+)/ATOM  $2$modified_nt_replacements{$key}$3/) {
					next RNAVIEW_LINE;
				}
			}
		}
		print OUTPUT_RNAVIEW_FILE $line;
		$line = <INPUT_PDB_FILE>;
	}
	close (OUTPUT_RNAVIEW_FILE);
	close (INPUT_PDB_FILE);
	print "$output_rnaview_file\n";
	
	GENERATE_FILE_FOR_NUCCYL3:;
	print "\n Generated PDB coordinate file for nuccyl3: ";
	open (INPUT_PDB_FILE, "$input_pdb_file");
	open (OUTPUT_NUCCYL3_FILE, ">$output_nuccyl3_file");
	$line = <INPUT_PDB_FILE>;
	NUCCYL3_LINE: while ($line ne "") {
		if ($line =~ /^HETATM/) {
			foreach $key (sort keys %modified_nt) {
				if ($line =~ s/^HETATM(\s{0,4}\d{1,7}\s{0,}[\w|\*]+\s{0,}$key\s{0,}\w\s{0,3}\w+\s{0,}\-{0,1}\d+\.\d+)/ATOM  $1/) {
					next NUCCYL3_LINE;
				}
			}
		}
		print OUTPUT_NUCCYL3_FILE $line;
		$line = <INPUT_PDB_FILE>;
	}
	close (OUTPUT_NUCCYL3_FILE);
	close (INPUT_PDB_FILE);
	print "$output_nuccyl3_file\n";

	&date_end;
	
	EXIT_1:;
	exit;
}

sub nuccyl_cyl_2d {
	#
	# Example:
	#
	#	./nuccyl.pl -cyl_2d dna.inp dna_nuccyl3.pdb | tee nuccyl2.log
	#
	PRINT_INTRO:;
	print "\n nuccyl_cyl_2d $nuccyl_version - $nuccyl_date";
	print "\n ============================\n";
	&date_start;
	
	INPUT_VARIABLE_SETUP:;
	$input_3dna_file = $ARGV[1];
	$input_nuccyl3_file_coord = $ARGV[2];
	$output_nuccyl3_file = "nuccyl2.out";
		
	IO_SYNTAX_CHECK:;
	$syntax_2d = " Usage: nuccyl.pl -cyl_2d <input_3dna_base_pair_file> <input_nuccyl3_coordinate_file>";
	unless (($input_3dna_file)&&($input_nuccyl3_file_coord)) {
		print "\n$syntax_2d\n\n";
		exit;
	}
	
	IO_FILE_CHECK:;
	open (INPUT_3DNA_FILE, "$input_3dna_file") || die "\n ERROR: Can't open input file \"$input_3dna_file\"!\n\n";
	print "\n Input 3DNA base pair file: $input_3dna_file\n";
	open (INPUT_NUCCYL3_FILE_COORD, "$input_nuccyl3_file_coord") || die "\n ERROR: Can't open input file \"$input_nuccyl3_file_coord\"!\n\n";
	print "\n Input nuccyl3 coordinate file: $input_nuccyl3_file_coord\n";
	close (INPUT_3DNA_FILE);
	if (-e "$output_nuccyl3_file") {
		unlink $output_nuccyl3_file;
	}

	EXTRACT_BP_INFO_3DNA:;
	print "\n Generated file for nuccyl3: ";
	open (INPUT_3DNA_FILE, "$input_3dna_file");
	open (OUTPUT_NUCCYL3_FILE, ">$output_nuccyl3_file");
	$line = <INPUT_3DNA_FILE>;
	while ($line ne "") {
		if ($line =~ /\s+(.)\:\.+(\d+).*\]($all_nt_match)\-{0,5}($all_nt_match).*\.(\d+)\_\:(.)/) {
			($chain_1,$base_1,$type_1,$type_2,$base_2,$chain_2) = ($1,$2,$3,$4,$5,$6);
			++$bp_number;
			if ($type_1 =~ /^(($purine_match))/) {
				$atom_1 = "N1";
			} elsif ($type_1 =~ /^(($pyrimidine_match))/) {
				$atom_1 = "N3";
			}
			if ($type_2 =~ /^(($purine_match))/) {
				$atom_2 = "N1";
			} elsif ($type_2 =~ /^(($pyrimidine_match))/) {
				$atom_2 = "N3";
			}
			# this checks whether any current standard residue name has been (temporarily) replacing a modified nucleotide name and, if so, it replaces it with the latter:
			open (IN_FILE, "$input_nuccyl3_file_coord");
			while (<IN_FILE>) {
				$line = $_;
				if ($line =~ /^ATOM\s{0,6}\d{1,7}\s{0,}[\w|\*]+\s{0,}(\w{1,4})\s{0,}$chain_1\s{0,3}$base_1\s{0,}\-{0,1}\d+\.\d+/) {
					$type_1_full=$1;
					last;
				}
			}
			while (<IN_FILE>) {
			$line = $_;
				if ($line =~ /^ATOM\s{0,6}\d{1,7}\s{0,}[\w|\*]+\s{0,}(\w{1,4})\s{0,}$chain_2\s{0,3}$base_2\s{0,}\-{0,1}\d+\.\d+/) {
					$type_2_full=$1;
					last;
				}
			}
			close (IN_FILE);
			#
			select(OUTPUT_NUCCYL3_FILE);
			$~ = 'BP_FORMAT';
	        write;
			$line = <INPUT_3DNA_FILE>;
		} else {
			$line = <INPUT_3DNA_FILE>;
		}
	}
	close (INPUT_3DNA_FILE);
	close (OUTPUT_NUCCYL3_FILE);
	select(STDOUT);
	print "$output_nuccyl3_file\n";

	FIND_UNPAIRED_NUCLEOTIDES_3DNA:;
	open (OUTPUT_NUCCYL3_FILE, ">>$output_nuccyl3_file");
	$bp_list = &build_bp_list($output_nuccyl3_file);
	@nt_list = &build_nt_list($input_nuccyl3_file_coord);
	$unpaired_nt_number = 0;
	foreach $nt_list (@nt_list) {
		if ($nt_list) {
			($pdb_chain, $pdb_residue_name, $pdb_residue_number) = split / /, $nt_list;
			unless ($bp_list =~ /\s$pdb_chain\s+$pdb_residue_name\s+$pdb_residue_number\s/) {
				$unpaired_nt_number++;
				# this checks whether current standard residue name has been (temporarily) replacing a modified nucleotide name and, if so, it replaces it with the latter:
				open (IN_FILE, "$input_nuccyl3_file_coord");
				while (<IN_FILE>) {
					$line = $_;
 					if ($line =~ /^ATOM\s{0,6}\d{1,7}\s{0,}[\w|\*]+\s{0,}(\w{1,4})\s{0,}$pdb_chain\s{0,3}$pdb_residue_number\s{0,}\-{0,1}\d+\.\d+/) {
						$pdb_residue_name_full=$1;
						last;
					}
				}
				close (IN_FILE);
				#
				if ($pdb_residue_name =~ /^(($purine_match))$/) {
					select(OUTPUT_NUCCYL3_FILE);
					$~ = 'NT_FORMAT_UNPAIRED_PURINES';
	        		write;
				} elsif ($pdb_residue_name =~ /^(($pyrimidine_match))$/) {
					select(OUTPUT_NUCCYL3_FILE);
					$~ = 'NT_FORMAT_UNPAIRED_PYRIMIDINES';
	        		write;
				}	
			}
		}
	}
	select (STDOUT);
	close (OUTPUT_NUCCYL3_FILE);

	&date_end;
	
	EXIT_2D:;
	exit;
}

sub nuccyl_cyl_2r {
	#
	# Example:
	#
	#	./nuccyl.pl -cyl_2r 1EVV_rnaview.pdb_sort.out 1EVV_nuccyl3.pdb | tee nuccyl2.log
	#
	PRINT_INTRO:;
	print "\n nuccyl_cyl_2r $nuccyl_version - $nuccyl_date";
	print "\n ============================\n";
	&date_start;
	
	INPUT_VARIABLE_SETUP:;
	$input_rnaview_file = $ARGV[1];
	$input_nuccyl3_file_coord = $ARGV[2];
	$output_nuccyl3_file = "nuccyl2.out";
	$tmp_file = "nuccyl2r.tmp";
	
	IO_SYNTAX_CHECK:;
	$syntax_2r = " Usage: nuccyl.pl -cyl_2r <input_rnaview_base_pair_file> <input_nuccyl3_coordinate_file>";
	unless (($input_rnaview_file)&&($input_nuccyl3_file_coord)) {
		print "\n$syntax_2r\n\n";
		exit;
	}
	
	IO_FILE_CHECK:;
	open (INPUT_RNAVIEW_FILE, "$input_rnaview_file") || die "\n ERROR: Can't open input file \"$input_rnaview_file\"!\n\n";
	print "\n Input RNAView base pair file: $input_rnaview_file\n";
	open (INPUT_NUCCYL3_FILE_COORD, "$input_nuccyl3_file_coord") || die "\n ERROR: Can't open input file \"$input_nuccyl3_file_coord\"!\n\n";
	print "\n Input nuccyl3 coordinate file: $input_nuccyl3_file_coord\n";
	close (INPUT_RNAVIEW_FILE);
	if (-e "$output_nuccyl3_file") {
		unlink $output_nuccyl3_file;
	}
	if (-e "$tmp_file") {
		unlink $tmp_file;
	}
	
	SORT_RNAVIEW_FILE:;
	system("sort $input_rnaview_file > $tmp_file");
	
	EXTRACT_BP_INFO_RNAVIEW:;
	print "\n Generated file for nuccyl3: ";
	open (TMP_FILE, "$tmp_file");
	open (OUTPUT_NUCCYL3_FILE, ">$output_nuccyl3_file");
	$line = <TMP_FILE>;
	while ($line ne "") {
		if ($line =~ /^\s{0,}\d+_\d+\,\ (.)\:\s{0,}(\d+)\ (\w)\-(\w)\s{0,}(\d+)\ (.)\:/) {
			($chain_1,$base_1,$type_1,$type_2,$base_2,$chain_2) = ($1,$2,$3,$4,$5,$6);
			++$bp_number;
			if ($type_1 =~ /^(($purine_match))/) {
				$atom_1 = "N1";
			} elsif ($type_1 =~ /^(($pyrimidine_match))/) {
				$atom_1 = "N3";
			}
			if ($type_2 =~ /^(($purine_match))/) {
				$atom_2 = "N1";
			} elsif ($type_2 =~ /^(($pyrimidine_match))/) {
				$atom_2 = "N3";
			}
			# this checks whether any current standard residue name has been (temporarily) replacing a modified nucleotide name and, if so, it replaces it with the latter:
			open (IN_FILE, "$input_nuccyl3_file_coord");
			while (<IN_FILE>) {
				$line = $_;
				if ($line =~ /^ATOM\s{0,6}\d{1,7}\s{0,}[\w|\*]+\s{0,}(\w{1,4})\s{0,}$chain_1\s{0,3}$base_1\s{0,}\-{0,1}\d+\.\d+/) {
					$type_1_full=$1;
					last;
				}
			}
			while (<IN_FILE>) {
			$line = $_;
				if ($line =~ /^ATOM\s{0,6}\d{1,7}\s{0,}[\w|\*]+\s{0,}(\w{1,4})\s{0,}$chain_2\s{0,3}$base_2\s{0,}\-{0,1}\d+\.\d+/) {
					$type_2_full=$1;
					last;
				}
			}
			close (IN_FILE);
			#
			select(OUTPUT_NUCCYL3_FILE);
			$~ = 'BP_FORMAT';
	        write;
			$line = <TMP_FILE>;
		} else {
			$line = <TMP_FILE>;
		}
	}
	close (TMP_FILE);
	close (OUTPUT_NUCCYL3_FILE);
	unlink $tmp_file;
	select(STDOUT);
	print "$output_nuccyl3_file\n";

	FIND_UNPAIRED_NUCLEOTIDES_RNAVIEW:;
	open (OUTPUT_NUCCYL3_FILE, ">>$output_nuccyl3_file");
	$bp_list = &build_bp_list($output_nuccyl3_file);
	@nt_list = &build_nt_list($input_nuccyl3_file_coord);
	$unpaired_nt_number = 0;
	foreach $nt_list (@nt_list) {
		if ($nt_list) {
			($pdb_chain, $pdb_residue_name, $pdb_residue_number) = split / /, $nt_list;
			unless ($bp_list =~ /\s$pdb_chain\s+$pdb_residue_name\s+$pdb_residue_number\s/) {
				$unpaired_nt_number++;
				# this checks whether current standard residue name has been (temporarily) replacing a modified nucleotide name and, if so, it replaces it with the latter:
				open (IN_FILE, "$input_nuccyl3_file_coord");
				while (<IN_FILE>) {
					$line = $_;
 					if ($line =~ /^ATOM\s{0,6}\d{1,7}\s{0,}[\w|\*]+\s{0,}(\w{1,4})\s{0,}$pdb_chain\s{0,3}$pdb_residue_number\s{0,}\-{0,1}\d+\.\d+/) {
						$pdb_residue_name_full=$1;
						last;
					}
				}
				close (IN_FILE);
				#
				if ($pdb_residue_name =~ /^(($purine_match))$/) {
					select(OUTPUT_NUCCYL3_FILE);
					$~ = 'NT_FORMAT_UNPAIRED_PURINES';
	        		write;
				} elsif ($pdb_residue_name =~ /^(($pyrimidine_match))$/) {
					select(OUTPUT_NUCCYL3_FILE);
					$~ = 'NT_FORMAT_UNPAIRED_PYRIMIDINES';
	        		write;
				}	
			}
		}
	}
	select (STDOUT);
	close (OUTPUT_NUCCYL3_FILE);
	
	&date_end;
	
	EXIT_2R:;
	exit;
}

sub nuccyl_cyl_3 {
	#
	# Example:
	#
	#	./nuccyl.pl -cyl_3 nuccyl2.out_edit 1EVV_nuccyl3.pdb 5000 B X 500 obj | tee nuccyl3.log
	#
	PRINT_INTRO:;
	print "\n nuccyl_cyl_3 $nuccyl_version - $nuccyl_date";
	print "\n ===========================\n";
	&date_start;
	
	INPUT_VARIABLE_SETUP:;
	$output_nuccyl2_file = $ARGV[1];
	$input_pdb_file = $ARGV[2];
	$fatom_number_start = $ARGV[3];
	$fatom = $ARGV[4];
	$fres = $ARGV[5];
	$fres_number_start = $ARGV[6];
	$obj_option = $ARGV[7];
	$focc = "1.00";
	$fb = "60.00";
	$output_err_file = "nuccyl3.problems";
	$output_pdb_file = "nuccyl3.pdb";
	$output_pml_file = "nuccyl3.pml";
	
	IO_SYNTAX_CHECK:;
	$syntax_3 = " Usage: nuccyl.pl -cyl_3 <input_nuccyl2_file> <input_pdb_coordinate_file> <fatom_number_start> <fatom> <fres> <fres_number_start> <object_option> [obj/no_obj]";
	unless ($output_nuccyl2_file && $input_pdb_file && $fatom_number_start && $fatom && $fres && $fres_number_start && $obj_option) {
		die "\n$syntax_3\n\n";
	}
	unless (($fatom_number_start =~ /^\d{1,5}$/) && ($fatom_number_start >= 1)){
		die "\n$syntax_3\n\n ERROR: fatom_number_start must be an integer between 1 and 99999!\n\n";
	}
	unless ($fatom =~ /^[A-Z]{1}$/){
		die "\n$syntax_3\n\n ERROR: fatom must be a single uppercase letter!\n\n";
	}
	unless ($fres =~ /^[A-Z]{1,4}$/){
		die "\n$syntax_3\n\n ERROR: fres must consist of 1-4 uppercase letters!\n\n";
	}
	unless (($fres_number_start =~ /^\d{1,4}$/) && ($fres_number_start >= 1)){
		die "\n$syntax_3\n\n ERROR: fres_number_start must be an integer between 1 and 9999!\n\n";
	}
	unless (($obj_option eq "obj") || ($obj_option eq "no_obj")){
		die "\n$syntax_3\n\n ERROR: object_option must be set to either \"obj\" or \"no_obj\"!\n\n";
	}

	IO_FILE_CHECK:;
	open (OUTPUT_NUCCYL2_FILE, "$output_nuccyl2_file") || die "\n ERROR: Can't open input file \"$output_nuccyl2_file\"!\n\n";
	close (OUTPUT_NUCCYL2_FILE);
	print "\n Input nuccyl2 file: $output_nuccyl2_file\n";
	open (INPUT_PDB_FILE, "$input_pdb_file") || die "\n ERROR: Can't open input file \"$input_pdb_file\"!\n\n";
	close (INPUT_PDB_FILE);
	print "\n Input PDB file: $input_pdb_file\n";
	if (-e "$output_err_file") {
		unlink $output_err_file;
	}
	if (-e "$output_pdb_file") {
		unlink $output_pdb_file;
	}
	if (-e "$output_pml_file") {
		unlink $output_pml_file;
	}
	
	CALCULATE_CYLINDER_COORDINATES:;
	open (OUTPUT_NUCCYL2_FILE, "$output_nuccyl2_file");
	open (OUTPUT_ERR_FILE, ">>$output_err_file");
	open (OUTPUT_PDB_FILE, ">>$output_pdb_file");
	open (OUTPUT_PML_FILE, ">>$output_pml_file");
	$line = <OUTPUT_NUCCYL2_FILE>;
	while ($line ne "") {
		#
		# 1. process base pairs
		# ---------------------
		#
		# the "^\s+" at the beginning of the regular expression makes sure we only read active base pairs:
		if ($line =~ /^pair_(\d+)\s+\:\s+(\w)\s+(\w+)\s+(\d+)\s+(\w+)\s+\-\s+(\w)\s+(\w+)\s+(\d+)\s+(\w+)/) {
			($bp_number,$chain_1,$type_1,$base_1,$atom_1,$chain_2,$type_2,$base_2,$atom_2) = ($1,$2,$3,$4,$5,$6,$7,$8,$9);
			&bp_intro($bp_number,$chain_1,$type_1,$base_1,$atom_1,$chain_2,$type_2,$base_2,$atom_2);
			#
			# reset all coordinates:
			$x_a = $y_a = $z_a = $x_ap = $y_ap = $z_ap = "";
			$x_b = $y_b = $z_b = $x_bp = $y_bp = $z_bp = "";
			$x_c = $y_c = $z_c = $x_cp = $y_cp = $z_cp = "";
			$x_base = $y_base = $z_base = "";
			$x_backbone_1 = $y_backbone_1 = $z_backbone_1 = "";
			$x_backbone_2 = $y_backbone_2 = $z_backbone_2 = "";
			#
			# get coordinates for the base atom of nt 1:
			($x_a,$y_a,$z_a) = &bp_get_coords($chain_1,$type_1,$base_1,$atom_1);
			if ($x_a && $y_a && $z_a) {
				print " \|\n coordinate A  ($atom_1)   =>\t$x_a\t$y_a\t$z_a\n";
			} else {
				print " \|\n coordinate A  ($atom_1)   =>\t*** WARNING - ATOM NOT FOUND! SKIPPING BASE PAIR $bp_number ($chain_1 $type_1 $base_1 - $chain_2 $type_2 $base_2)... ***\n";
				print OUTPUT_ERR_FILE " coordinate A  ($atom_1)   =>\t*** WARNING - ATOM NOT FOUND! SKIPPING BASE PAIR $bp_number ($chain_1 $type_1 $base_1 - $chain_2 $type_2 $base_2)... ***\n\n";
				$line = <OUTPUT_NUCCYL2_FILE>;
				next;
			}
			#
			# get coordinates for the base atom of nt 2:
			($x_ap,$y_ap,$z_ap) = &bp_get_coords($chain_2,$type_2,$base_2,$atom_2);
			if ($x_ap && $y_ap && $z_ap) {
				print " coordinate A\' ($atom_2)   =>\t$x_ap\t$y_ap\t$z_ap\n";
			} else {
				print " coordinate A\' ($atom_2)   =>\t*** WARNING - ATOM NOT FOUND! SKIPPING BASE PAIR $bp_number ($chain_1 $type_1 $base_1 - $chain_2 $type_2 $base_2)... ***\n";
				print OUTPUT_ERR_FILE " coordinate A\' ($atom_2)   =>\t*** WARNING - ATOM NOT FOUND! SKIPPING BASE PAIR $bp_number ($chain_1 $type_1 $base_1 - $chain_2 $type_2 $base_2)... ***\n\n";
				$line = <OUTPUT_NUCCYL2_FILE>;
				next;
			}			
			#
			# calculate middle point for bases (= coordinates for cylinder base end):
			($x_base,$y_base,$z_base) = &coordinate_middlepoint($x_a,$y_a,$z_a,$x_ap,$y_ap,$z_ap);
			print " => base coordinates for cylinder   =>\t$x_base\t$y_base\t$z_base\n"; 
			#
			# ///
			#
			# get coordinates for 5' P atom of nt 1:
			$p_5prime = "P";
			($x_b,$y_b,$z_b) = &bp_get_coords($chain_1,$type_1,$base_1,$p_5prime);
			# if 5' P exists, use it:
			if ($x_b && $y_b && $z_b) {
				print " \|\n coordinate B  ($p_5prime)    =>\t$x_b\t$y_b\t$z_b\n";
			# otherwise, use the C5* of nt 1 as a replacement:
			} else {
				print " \|\n coordinate B  ($p_5prime)    =>\t*** WARNING - $chain_1$base_1 $p_5prime ATOM NOT FOUND! TRYING TO USE $chain_1$base_1 C5\*... ***\n";
				print OUTPUT_ERR_FILE " coordinate B  ($p_5prime)    =>\t*** WARNING - $chain_1$base_1 $p_5prime ATOM NOT FOUND! TRYING TO USE $chain_1$base_1 C5\*... ***\n";
				$p_5prime = 'C5\*';
				($x_b,$y_b,$z_b) = &bp_get_coords($chain_1,$type_1,$base_1,$p_5prime);
				if ($x_b && $y_b && $z_b) {
					print " coordinate B  (C5\*)  =>\t$x_b\t$y_b\t$z_b\n";
					print OUTPUT_ERR_FILE " PROBLEM BYPASSED: coordinate B  (C5\*)  =>\t$x_b\t$y_b\t$z_b\n\n";
				} else {
					print " coordinate B  (C5\*)  =>\t*** WARNING - $chain_1$base_1 C5\* ATOM NOT FOUND! SKIPPING BASE PAIR $bp_number ($chain_1 $type_1 $base_1 - $chain_2 $type_2 $base_2)... ***\n";
					print OUTPUT_ERR_FILE " coordinate B  (C5\*)  =>\t*** WARNING - $chain_1$base_1 C5\* ATOM NOT FOUND! SKIPPING BASE PAIR $bp_number ($chain_1 $type_1 $base_1 - $chain_2 $type_2 $base_2)... ***\n\n";
					$line = <OUTPUT_NUCCYL2_FILE>;
					next;
				}
			}
			#
			# get coordinates for the 3' P atom of nt 1 (= P of nt following nt 1):
			$base_next = ($base_1 + 1);
			$p_3prime = "P";
			($x_c,$y_c,$z_c) = &bp_get_coords_no_type($chain_1,$base_next,$p_3prime);
			# if P of following nt exists, use it:
			if ($x_c && $y_c && $z_c) {
				print " coordinate C  ($p_3prime+1)  =>\t$x_c\t$y_c\t$z_c\n";
			# otherwise, use the P atom of the current nt itself as a replacement (this is the very 3' end of a PyMOL nucleic acid cartoon):
			} else {
				print " coordinate C  ($p_3prime+1)  =>\t*** WARNING - $chain_1$base_next $p_3prime ATOM NOT FOUND! TRYING TO USE $chain_1$base_1 $p_3prime... ***\n";
				print OUTPUT_ERR_FILE " coordinate C  ($p_3prime+1)  =>\t*** WARNING - $chain_1$base_next $p_3prime ATOM NOT FOUND! TRYING TO USE $chain_1$base_1 $p_3prime... ***\n";
				$p_3prime = 'P'; # alternatively, one could use O3* ($p_3prime = 'O3\*';), but then the base cylinder would not be connected to the backbone (PyMOL cartoon ends at the P atom)
				($x_c,$y_c,$z_c) = &bp_get_coords($chain_1,$type_1,$base_1,$p_3prime);
				if ($x_c && $y_c && $z_c) {
					print " coordinate C  ($p_3prime)    =>\t$x_c\t$y_c\t$z_c\n";
					print OUTPUT_ERR_FILE " PROBLEM BYPASSED: coordinate C  ($p_3prime)    =>\t$x_c\t$y_c\t$z_c\n\n";
				} else {
					print " coordinate C  ($p_3prime)    =>\t*** WARNING - $chain_1$base_1 $p_3prime ATOM NOT FOUND! SKIPPING BASE PAIR $bp_number ($chain_1 $type_1 $base_1 - $chain_2 $type_2 $base_2)... ***\n";
					print OUTPUT_ERR_FILE " coordinate C  ($p_3prime)    =>\t*** WARNING - $chain_1$base_1 $p_3prime ATOM NOT FOUND! SKIPPING BASE PAIR $bp_number ($chain_1 $type_1 $base_1 - $chain_2 $type_2 $base_2)... ***\n\n";
					$line = <OUTPUT_NUCCYL2_FILE>;
					next;
				}
			}			
			#
			# calculate middle point for backbone of nt 1 (= cylinder backbone end coordinates (1)):
			($x_backbone_1,$y_backbone_1,$z_backbone_1) = &coordinate_middlepoint($x_b,$y_b,$z_b,$x_c,$y_c,$z_c);
			print " => backbone coordinates for cylinder (1)  =>\t$x_backbone_1\t$y_backbone_1\t$z_backbone_1\n"; 
			#
			# ///
			#
			# get coordinates for 5' P atom of nt 2:
			$p_5prime = "P";
			($x_bp,$y_bp,$z_bp) = &bp_get_coords($chain_2,$type_2,$base_2,$p_5prime);
			# if 5' P exists, use it:
			if ($x_bp && $y_bp && $z_bp) {
				print " \|\n coordinate B\' ($p_5prime)    =>\t$x_bp\t$y_bp\t$z_bp\n";
			# otherwise, use the C5* of nt 2 as a replacement:
			} else {
				print " \|\n coordinate B\' ($p_5prime)    =>\t*** WARNING - $chain_2$base_2 $p_5prime ATOM NOT FOUND! TRYING TO USE $chain_2$base_2 C5\*... ***\n";
				print OUTPUT_ERR_FILE " coordinate B\' ($p_5prime)    =>\t*** WARNING - $chain_2$base_2 $p_5prime ATOM NOT FOUND! TRYING TO USE $chain_2$base_2 C5\*... ***\n";
				$p_5prime = 'C5\*';
				($x_bp,$y_bp,$z_bp) = &bp_get_coords($chain_2,$type_2,$base_2,$p_5prime);
				if ($x_bp && $y_bp && $z_bp) {
					print " coordinate B\' (C5\*)  =>\t$x_bp\t$y_bp\t$z_bp\n";
					print OUTPUT_ERR_FILE " PROBLEM BYPASSED: coordinate B\' (C5\*)  =>\t$x_bp\t$y_bp\t$z_bp\n\n";
				} else {
					print " coordinate B\' (C5\*)  =>\t*** WARNING - $chain_2$base_2 C5\* ATOM NOT FOUND! SKIPPING BASE PAIR $bp_number ($chain_1 $type_1 $base_1 - $chain_2 $type_2 $base_2)... ***\n";
					print OUTPUT_ERR_FILE " coordinate B\' (C5\*)  =>\t*** WARNING - $chain_2$base_2 C5\* ATOM NOT FOUND! SKIPPING BASE PAIR $bp_number ($chain_1 $type_1 $base_1 - $chain_2 $type_2 $base_2)... ***\n\n";
					$line = <OUTPUT_NUCCYL2_FILE>;
					next;
				}
			}
			#
			# get coordinates for the 3' P atom of nt 2 (= P of nt following nt 2):
			$base_next = ($base_2 + 1);
			$p_3prime = "P";
			($x_cp,$y_cp,$z_cp) = &bp_get_coords_no_type($chain_2,$base_next,$p_3prime);
			# if P of following nt exists, use it:
			if ($x_cp && $y_cp && $z_cp) {
				print " coordinate C' ($p_3prime+1)  =>\t$x_cp\t$y_cp\t$z_cp\n";
			# otherwise, use the P atom of the current nt itself as a replacement (this is the very 3' end of a PyMOL nucleic acid cartoon):
			} else {
				print " coordinate C' ($p_3prime+1)  =>\t*** WARNING - $chain_2$base_next $p_3prime ATOM NOT FOUND! TRYING TO USE $chain_2$base_2 $p_3prime... ***\n";
				print OUTPUT_ERR_FILE " coordinate C' ($p_3prime+1)  =>\t*** WARNING - $chain_2$base_next $p_3prime ATOM NOT FOUND! TRYING TO USE $chain_2$base_2 $p_3prime... ***\n";
				$p_3prime = 'P'; # alternatively, one could use O3* ($p_3prime = 'O3\*';), but then the base cylinder would not be connected to the backbone (PyMOL cartoon ends at the P atom)
				($x_cp,$y_cp,$z_cp) = &bp_get_coords($chain_2,$type_2,$base_2,$p_3prime);
				if ($x_cp && $y_cp && $z_cp) {
					print " coordinate C' ($p_3prime)    =>\t$x_cp\t$y_cp\t$z_cp\n";
					print OUTPUT_ERR_FILE " PROBLEM BYPASSED: coordinate C' ($p_3prime)    =>\t$x_cp\t$y_cp\t$z_cp\n\n";
				} else {
					print " coordinate C' ($p_3prime)    =>\t*** WARNING - $chain_2$base_2 $p_3prime ATOM NOT FOUND! SKIPPING BASE PAIR $bp_number ($chain_1 $type_1 $base_1 - $chain_2 $type_2 $base_2)... ***\n";
					print OUTPUT_ERR_FILE " coordinate C' ($p_3prime)    =>\t*** WARNING - $chain_2$base_2 $p_3prime ATOM NOT FOUND! SKIPPING BASE PAIR $bp_number ($chain_1 $type_1 $base_1 - $chain_2 $type_2 $base_2)... ***\n\n";
					$line = <OUTPUT_NUCCYL2_FILE>;
					next;
				}
			}			
			#
			# calculate middle point for backbone of nt 2 (= cylinder backbone end coordinates (2)):
			($x_backbone_2,$y_backbone_2,$z_backbone_2) = &coordinate_middlepoint($x_bp,$y_bp,$z_bp,$x_cp,$y_cp,$z_cp);
			print " => backbone coordinates for cylinder (2)  =>\t$x_backbone_2\t$y_backbone_2\t$z_backbone_2\n \|\n"; 
			#
			# ///
			#
			# at this point, for each active base pair, we need to generate 3 pseudo atoms with these coordinates:
			#
			#	$x_bp			$y_bp			$z_bp
			#	$x_backbone_1	$y_backbone_1	$z_backbone_1
			#	$x_backbone_2	$y_backbone_2	$z_backbone_2
			#
			select(STDOUT);
			print " ";
			$~ = 'FPDB_FORMAT_BASE';
			write;
			select(OUTPUT_PDB_FILE);
			$~ = 'FPDB_FORMAT_BASE';
			write;
			++$fatom_number_start;
			++$fres_number_start;
			#
			# ///
			#
			select(STDOUT);
			print " ";
			$~ = 'FPDB_FORMAT_BACKBONE_1';
			write;
			select(OUTPUT_PDB_FILE);
			$~ = 'FPDB_FORMAT_BACKBONE_1';
			write;
			++$fatom_number_start;
			++$fres_number_start;
			#
			# ///
			#
			select(STDOUT);
			print " ";
			$~ = 'FPDB_FORMAT_BACKBONE_2';
			write;
			select(OUTPUT_PDB_FILE);
			$~ = 'FPDB_FORMAT_BACKBONE_2';
			write;
			++$fatom_number_start;
			++$fres_number_start;
			#
			# ///
			#
			# finally, we need to output commands so that PyMOL will connect:
			#
			#	bp - backbone_1
			#   bp - backbone_2
			#
			select(STDOUT);
			print " \|\n";
			$bond_1 = ($fres_number_start - 3);
			$bond_2 = ($fres_number_start - 2);
			$bond_3 = ($fres_number_start - 1);
			#
			$object_1 = $chain_1 . "_" . $type_1 . $base_1;
			print " create $object_1, ((resi $bond_1 or resi $bond_2) and resn $fres)\n";
			if ($obj_option eq "obj") {
				print OUTPUT_PML_FILE " create $object_1, ((resi $bond_1 or resi $bond_2) and resn $fres)\n";
			} elsif ($obj_option eq "no_obj") {
				print OUTPUT_PML_FILE "\# create $object_1, ((resi $bond_1 or resi $bond_2) and resn $fres)\n";
			}
			$object_2 = $chain_2 . "_" . $type_2 . $base_2;
			print " create $object_2, ((resi $bond_1 or resi $bond_3) and resn $fres)\n";
			if ($obj_option eq "obj") {
				print OUTPUT_PML_FILE " create $object_2, ((resi $bond_1 or resi $bond_3) and resn $fres)\n";
			} elsif ($obj_option eq "no_obj") {
				print OUTPUT_PML_FILE "\# create $object_2, ((resi $bond_1 or resi $bond_3) and resn $fres)\n";
			}
			#
			print " bond $object_1///$bond_1/$fatom, $object_1///$bond_2/$fatom\n";
			if ($obj_option eq "obj") {
				print OUTPUT_PML_FILE " bond $object_1///$bond_1/$fatom, $object_1///$bond_2/$fatom\n";
			} elsif ($obj_option eq "no_obj") {
				print OUTPUT_PML_FILE " bond ////$bond_1/$fatom, ////$bond_2/$fatom\n";
			}
			print " bond $object_2///$bond_1/$fatom, $object_2///$bond_3/$fatom\n";
			if ($obj_option eq "obj") {
				print OUTPUT_PML_FILE " bond $object_2///$bond_1/$fatom, $object_2///$bond_3/$fatom\n";
			} elsif ($obj_option eq "no_obj") {
				print OUTPUT_PML_FILE " bond ////$bond_1/$fatom, ////$bond_3/$fatom\n";
			}
			#
			$line = <OUTPUT_NUCCYL2_FILE>;
		#
		# 2. process unpaired nucleotides
		# -------------------------------
		#
	} elsif ($line =~ /^unpaired_(\d+)\s+\:\s+(\w)\s+(\w+)\s+(\d+)\s+(\w+)/) { 
			($unpaired_number,$chain_1,$type_1,$base_1,$atom_1) = ($1,$2,$3,$4,$5);
			&unpaired_intro($unpaired_number,$chain_1,$type_1,$base_1,$atom_1);
			#
			# reset all coordinates:
			$x_a = $y_a = $z_a = "";
			$x_b = $y_b = $z_b = "";
			$x_c = $y_c = $z_c = "";
			$x_backbone_1 = $y_backbone_1 = $z_backbone_1 = "";
			#
			# get coordinates for the base atom:
			($x_a,$y_a,$z_a) = &bp_get_coords($chain_1,$type_1,$base_1,$atom_1);
			if ($x_a && $y_a && $z_a) {
				print " \|\n coordinate A  ($atom_1)   =>\t$x_a\t$y_a\t$z_a\n";
			} else {
				print " \|\n coordinate A  ($atom_1)   =>\t*** WARNING - ATOM NOT FOUND! SKIPPING NUCLEOTIDE $chain_1 $type_1 $base_1... ***\n";
				print OUTPUT_ERR_FILE " coordinate A  ($atom_1)   =>\t*** WARNING - ATOM NOT FOUND! SKIPPING NUCLEOTIDE $chain_1 $type_1 $base_1... ***\n\n";
				$line = <OUTPUT_NUCCYL2_FILE>;
				next;
			}
			#
			# in this case, cylinder base end coordinates equal those of the specified base atom itself:
			$x_base = $x_a;
			$y_base = $y_a;
			$z_base = $z_a;
			print " => base coordinates for cylinder   =>\t$x_base\t$y_base\t$z_base\n"; 
			#
			# ///
			#
			# get coordinates for the 5' P atom:
			$p_5prime = "P";
			($x_b,$y_b,$z_b) = &bp_get_coords($chain_1,$type_1,$base_1,$p_5prime);
			# if 5' P exists, use it:
			if ($x_b && $y_b && $z_b) {
				print " \|\n coordinate B  ($p_5prime)    =>\t$x_b\t$y_b\t$z_b\n";
			# otherwise, use the C5* as a replacement:
			} else {
				print " \|\n coordinate B  ($p_5prime)    =>\t*** WARNING - $chain_1$base_1 $p_5prime ATOM NOT FOUND! TRYING TO USE $chain_1$base_1 C5\*... ***\n";
				print OUTPUT_ERR_FILE " coordinate B  ($p_5prime)    =>\t*** WARNING - $chain_1$base_1 $p_5prime ATOM NOT FOUND! TRYING TO USE $chain_1$base_1 C5\*... ***\n";
				$p_5prime = 'C5\*';
				($x_b,$y_b,$z_b) = &bp_get_coords($chain_1,$type_1,$base_1,$p_5prime);
				if ($x_b && $y_b && $z_b) {
					print " coordinate B  (C5\*)  =>\t$x_b\t$y_b\t$z_b\n";
					print OUTPUT_ERR_FILE " PROBLEM BYPASSED: coordinate B  (C5\*)  =>\t$x_b\t$y_b\t$z_b\n\n";
				} else {
					print " coordinate B  (C5\*)  =>\t*** WARNING - $chain_1$base_1 C5\* ATOM NOT FOUND! SKIPPING BASE PAIR $bp_number ($chain_1 $type_1 $base_1 - $chain_2 $type_2 $base_2)... ***\n";
					print OUTPUT_ERR_FILE " coordinate B  (C5\*)  =>\t*** WARNING - $chain_1$base_1 C5\* ATOM NOT FOUND! SKIPPING BASE PAIR $bp_number ($chain_1 $type_1 $base_1 - $chain_2 $type_2 $base_2)... ***\n\n";
					$line = <OUTPUT_NUCCYL2_FILE>;
					next;
				}
			}
			#
			# get coordinates for the 3' P atom (= P of nt following current nt):
			$base_next = ($base_1 + 1);
			$p_3prime = "P";
			($x_c,$y_c,$z_c) = &bp_get_coords_no_type($chain_1,$base_next,$p_3prime);
			# if P of following nt exists, use it:
			if ($x_c && $y_c && $z_c) {
				print " coordinate C  ($p_3prime+1)  =>\t$x_c\t$y_c\t$z_c\n";
			# otherwise, use the P atom of the current nt itself as a replacement (this is the very 3' end of a PyMOL nucleic acid cartoon):
			} else {
				print " coordinate C  ($p_3prime+1)  =>\t*** WARNING - $chain_1$base_next $p_3prime ATOM NOT FOUND! TRYING TO USE $chain_1$base_1 $p_3prime... ***\n";
				print OUTPUT_ERR_FILE " coordinate C  ($p_3prime+1)  =>\t*** WARNING - $chain_1$base_next $p_3prime ATOM NOT FOUND! TRYING TO USE $chain_1$base_1 $p_3prime... ***\n";
				$p_3prime = 'P'; # alternatively, one could use O3* ($p_3prime = 'O3\*';), but then the base cylinder would not be connected to the backbone (PyMOL cartoon ends at the P atom)
				($x_c,$y_c,$z_c) = &bp_get_coords($chain_1,$type_1,$base_1,$p_3prime);
				if ($x_c && $y_c && $z_c) {
					print " coordinate C  ($p_3prime)    =>\t$x_c\t$y_c\t$z_c\n";
					print OUTPUT_ERR_FILE " PROBLEM BYPASSED: coordinate C  ($p_3prime)    =>\t$x_c\t$y_c\t$z_c\n\n";
				} else {
					print " coordinate C  ($p_3prime)    =>\t*** WARNING - $chain_1$base_1 $p_3prime ATOM NOT FOUND! SKIPPING BASE PAIR $bp_number ($chain_1 $type_1 $base_1 - $chain_2 $type_2 $base_2)... ***\n";
					print OUTPUT_ERR_FILE " coordinate C  ($p_3prime)    =>\t*** WARNING - $chain_1$base_1 $p_3prime ATOM NOT FOUND! SKIPPING BASE PAIR $bp_number ($chain_1 $type_1 $base_1 - $chain_2 $type_2 $base_2)... ***\n\n";
					$line = <OUTPUT_NUCCYL2_FILE>;
					next;
				}
			}			
			#
			# calculate middle point for nt backbone (= cylinder backbone end coordinates):
			($x_backbone_1,$y_backbone_1,$z_backbone_1) = &coordinate_middlepoint($x_b,$y_b,$z_b,$x_c,$y_c,$z_c);
			print " => backbone coordinates for cylinder (1)  =>\t$x_backbone_1\t$y_backbone_1\t$z_backbone_1\n \|\n"; 
			#
			# ///
			#
			# at this point, for each active unpaired nucleotide, we need to generate 2 pseudo atoms with these coordinates:
			#
			#	$x_bp			$y_bp			$z_bp
			#	$x_backbone_1	$y_backbone_1	$z_backbone_1
			#
			select(STDOUT);
			print " ";
			$~ = 'FPDB_FORMAT_BASE';
			write;
			select(OUTPUT_PDB_FILE);
			$~ = 'FPDB_FORMAT_BASE';
			write;
			++$fatom_number_start;
			++$fres_number_start;
			#
			# ///
			#
			select(STDOUT);
			print " ";
			$~ = 'FPDB_FORMAT_BACKBONE_1';
			write;
			select(OUTPUT_PDB_FILE);
			$~ = 'FPDB_FORMAT_BACKBONE_1';
			write;
			++$fatom_number_start;
			++$fres_number_start;
			#
			# ///
			#
			# finally, we need to output commands so that PyMOL will connect:
			#
			#	bp - backbone_1
			#
			select(STDOUT);
			print " \|\n";
			$bond_1 = ($fres_number_start - 2);
			$bond_2 = ($fres_number_start - 1);
			#
			$object_1 = $chain_1 . "_" . $type_1 . $base_1;
			print " create $object_1, ((resi $bond_1 or resi $bond_2) and resn $fres)\n";
			if ($obj_option eq "obj") {
				print OUTPUT_PML_FILE " create $object_1, ((resi $bond_1 or resi $bond_2) and resn $fres)\n";
			} elsif ($obj_option eq "no_obj") {
				print OUTPUT_PML_FILE "\# create $object_1, ((resi $bond_1 or resi $bond_2) and resn $fres)\n";
			}
			#
			print " bond $object_1///$bond_1/$fatom, $object_1///$bond_2/$fatom\n";
			if ($obj_option eq "obj") {
				print OUTPUT_PML_FILE " bond $object_1///$bond_1/$fatom, $object_1///$bond_2/$fatom\n";
			} elsif ($obj_option eq "no_obj") {
				print OUTPUT_PML_FILE " bond ////$bond_1/$fatom, ////$bond_2/$fatom\n";
			}
			#
			$line = <OUTPUT_NUCCYL2_FILE>;
		} else {
			$line = <OUTPUT_NUCCYL2_FILE>;
		}
	}
	close (OUTPUT_NUCCYL2_FILE);
	close (OUTPUT_ERR_FILE);
	close (OUTPUT_PDB_FILE);
	close (OUTPUT_PML_FILE);
	print "\n\n Generated PDB file for PyMOL: nuccyl3.pdb\n";
	print "\n Generated command (.pml) file for PyMOL: nuccyl3.pml\n";

	&date_end;
	
	EXIT_3:;
	exit;
}

sub nuccyl_filled {
	#
	# Example:
	#
	#	./nuccyl.pl -filled 1EVV.pdb 1EVV_pymol.pdb 1EVV_filled.pml no_obj 0 0 255 255 0 0 | tee 1EVV_filled.log
	#
	PRINT_INTRO:;
	print "\n nuccyl_filled $nuccyl_version - $nuccyl_date";
	print "\n ============================\n";
	&date_start;

	INPUT_VARIABLE_SETUP:;
	$input_pdb_file = $ARGV[1];
	$output_pdb_file = $ARGV[2];
	$output_pml_file = $ARGV[3];
	$obj_option = $ARGV[4];
	$pur_rgb_r = $ARGV[5];
	$pur_rgb_g = $ARGV[6];
	$pur_rgb_b = $ARGV[7];
	$pyr_rgb_r = $ARGV[8];
	$pyr_rgb_g = $ARGV[9];
	$pyr_rgb_b = $ARGV[10];
	
	IO_SYNTAX_CHECK:;
	$syntax_filled = " Usage: nuccyl.pl -filled <input_pdb_coordinate_file> <output_pdb_coordinate_file> <output_pml_file> <object_option> [obj/no_obj] <pur_rgb_r> <pur_rgb_g> <pur_rgb_b> (<pyr_rgb_r> <pyr_rgb_g> <pyr_rgb_b>)";
	unless (($input_pdb_file)&&($output_pdb_file)&&($output_pml_file)&&($obj_option)&&($pur_rgb_r =~ /\w+/)&&($pur_rgb_g =~ /\w+/)&&($pur_rgb_b =~ /\w+/)) {
		print "\n$syntax_filled\n\n";
		exit;
	}
	unless (($obj_option eq "obj") || ($obj_option eq "no_obj")){
		die "\n$syntax_filled\n\n ERROR: object_option must be set to either \"obj\" or \"no_obj\"!\n\n";
	}
	unless (($pyr_rgb_r)||($pyr_rgb_r == 0)) {
		$pyr_rgb_r = $pur_rgb_r;
	}
	unless (($pyr_rgb_g)||($pyr_rgb_g == 0)) {
		$pyr_rgb_g = $pur_rgb_g;
	}
	unless (($pyr_rgb_b)||($pyr_rgb_b == 0)) {
		$pyr_rgb_b = $pur_rgb_b;
	}
	unless (($pur_rgb_r =~ /^\d{1,3}$/)&&($pur_rgb_r >= 0)&&($pur_rgb_r <= 255)) {
		die "\n$syntax_filled\n\n ERROR: pur_rgb_r must be an integer between 0 and 255!\n\n";
	}
	unless (($pur_rgb_g =~ /^\d{1,3}$/)&&($pur_rgb_g >= 0)&&($pur_rgb_g <= 255)) {
		die "\n$syntax_filled\n\n ERROR: pur_rgb_g must be an integer between 0 and 255!\n\n";
	}
	unless (($pur_rgb_b =~ /^\d{1,3}$/)&&($pur_rgb_b >= 0)&&($pur_rgb_b <= 255)) {
		die "\n$syntax_filled\n\n ERROR: pur_rgb_b must be an integer between 0 and 255!\n\n";
	}

	unless (($pyr_rgb_r =~ /^\d{1,3}$/)&&($pyr_rgb_r >= 0)&&($pyr_rgb_r <= 255)) {
		die "\n$syntax_filled\n\n ERROR: pyr_rgb_r must be an integer between 0 and 255!\n\n";
	}
	unless (($pyr_rgb_g =~ /^\d{1,3}$/)&&($pyr_rgb_g >= 0)&&($pyr_rgb_g <= 255)) {
		die "\n$syntax_filled\n\n ERROR: pyr_rgb_g must be an integer between 0 and 255!\n\n";
	}
	unless (($pyr_rgb_b =~ /^\d{1,3}$/)&&($pyr_rgb_b >= 0)&&($pyr_rgb_b <= 255)) {
		die "\n$syntax_filled\n\n ERROR: pyr_rgb_b must be an integer between 0 and 255!\n\n";
	}
		
	IO_FILE_CHECK:;
	open (INPUT_PDB_FILE, "$input_pdb_file") || die "\n ERROR: Can't open input PDB file \"$input_pdb_file\"!\n\n";
	print "\n Input PDB coordinate file: $input_pdb_file\n";
	close (INPUT_PDB_FILE);
	if (-e "$output_pdb_file") {
		unlink $output_pdb_file;
	}
	if (-e "$output_pml_file") {
		unlink $output_pml_file;
	}

	GENERATE_OUTPUT_PDB_FILE:;
	print "\n Generated PDB coordinate file for PyMOL: ";
	open (INPUT_PDB_FILE, "$input_pdb_file");
	open (OUTPUT_PDB_FILE, ">$output_pdb_file");
	$line = <INPUT_PDB_FILE>;
	PYMOL_LINE: while ($line ne "") {
		if ($line =~ /^HETATM/) {
			foreach $key (sort keys %modified_nt) {
				if ($line =~ s/^HETATM(\s{0,4}\d{1,7}\s{0,}[\w|\*]+\s{0,}$key\s{0,}\w\s{0,3}\w+\s{0,}\-{0,1}\d+\.\d+)/ATOM  $1/) {
					next PYMOL_LINE;
				}
		
			}
		}
		print OUTPUT_PDB_FILE $line;
		$line = <INPUT_PDB_FILE>;
	}
	close (OUTPUT_PDB_FILE);
	close (INPUT_PDB_FILE);
	print "$output_pdb_file\n";
	
	CONVERT_RGB_VALUES:;
	$pur_rgb_r = ($pur_rgb_r/255);
	$pur_rgb_r = sprintf("%.2f", $pur_rgb_r);
	$pur_rgb_g = ($pur_rgb_g/255);
	$pur_rgb_g = sprintf("%.2f", $pur_rgb_g);
	$pur_rgb_b = ($pur_rgb_b/255);
	$pur_rgb_b = sprintf("%.2f", $pur_rgb_b);
	$pyr_rgb_r = ($pyr_rgb_r/255);
	$pyr_rgb_r = sprintf("%.2f", $pyr_rgb_r);
	$pyr_rgb_g = ($pyr_rgb_g/255);
	$pyr_rgb_g = sprintf("%.2f", $pyr_rgb_g);
	$pyr_rgb_b = ($pyr_rgb_b/255);
	$pyr_rgb_b = sprintf("%.2f", $pyr_rgb_b);
	
	CREATE_FILLED_BASE_CGO_OBJECT_INFORMATION:;
	print "\n Generated command (.pml) file for PyMOL: ";
	open (OUTPUT_PML_FILE, ">$output_pml_file");
	print OUTPUT_PML_FILE "from pymol.cgo import \*\nfrom pymol import cmd\n\n";
	@nt_list = &build_nt_list($input_pdb_file);
	$object_load = "";
	$object_load_1 = "";
	$object_load_2 = "";
	$object_name_1 = "";
	$object_name_2 = "";
	$object_name_total_pur = "";
	$object_name_total_pyr = "";
	$obj_printout_1 = "";
	$obj_printout_2 = "";
	foreach $nt_list (@nt_list) {
		$vertex_list_1 = "";
		$vertex_list_2 = "";
		if ($nt_list) {
			($pdb_chain, $pdb_residue_name, $pdb_residue_number) = split / /, $nt_list;
			#
			# 1. process purines:
			if ($pdb_residue_name =~ /^(($purine_match))$/) {
				#print "\n$pdb_chain, $pdb_residue_name, $pdb_residue_number -> PURINE\n"; # DEBUG
				#
				# first purine ring:
				open (IN_FILE, "$input_pdb_file");
				while (<IN_FILE>) {
					$line = $_;
					if ($line =~ /^(ATOM|HETATM)\s{0,6}(\d{1,7})\s{0,}(N1|C2|N3|C4|C5|C6)\s{0,}$pdb_residue_name\s{0,}$pdb_chain\s{0,3}$pdb_residue_number\s{0,}(\-{0,1}\d+\.\d+)\s{0,}(\-{0,1}\d+\.\d+)\s{0,}(\-{0,1}\d+\.\d+)\s{0,}(\d{1,3}\.\d{2})\s{0,}(\d{1,3}\.\d{2})/) {
						($pdb_atom_name,$pdb_x,$pdb_y,$pdb_z) = ($3,$4,$5,$6);
						$vertex_list_1 = $vertex_list_1 . " VERTEX, $pdb_x, $pdb_y, $pdb_z,";
						#chop $line; # DEBUG
						#print "\t$pdb_atom_name\t$pdb_x\t$pdb_y\t$pdb_z\t( $line )\n"; # DEBUG
					}
				}
				close (IN_FILE);
				$obj_name = $pdb_chain . "_" . $pdb_residue_name . "_" . $pdb_residue_number;
				$obj_name_1 = $pdb_chain . "_" . $pdb_residue_name . "_" . $pdb_residue_number . "_1";
				$obj_printout_1 = "$obj_name_1 = [ BEGIN, TRIANGLE_FAN, COLOR, $pur_rgb_r, $pur_rgb_g, $pur_rgb_b,$vertex_list_1 END]\n";
				print OUTPUT_PML_FILE "$obj_printout_1";
				#
				# second purine ring:
				open (IN_FILE, "$input_pdb_file");
				while (<IN_FILE>) {
					$line = $_;
					if ($line =~ /^(ATOM|HETATM)\s{0,6}(\d{1,7})\s{0,}(C4|C5|N7|C8|N9)\s{0,}$pdb_residue_name\s{0,}$pdb_chain\s{0,3}$pdb_residue_number\s{0,}(\-{0,1}\d+\.\d+)\s{0,}(\-{0,1}\d+\.\d+)\s{0,}(\-{0,1}\d+\.\d+)\s{0,}(\d{1,3}\.\d{2})\s{0,}(\d{1,3}\.\d{2})/) {
						($pdb_atom_name,$pdb_x,$pdb_y,$pdb_z) = ($3,$4,$5,$6);
						$vertex_list_2 = $vertex_list_2 . " VERTEX, $pdb_x, $pdb_y, $pdb_z,";
						#chop $line; # DEBUG
						#print "\t$pdb_atom_name\t$pdb_x\t$pdb_y\t$pdb_z\t( $line )\n"; # DEBUG
					}
				}
				close (IN_FILE);
				$obj_name_2 = $pdb_chain . "_" . $pdb_residue_name . "_" . $pdb_residue_number . "_2";
				$obj_printout_2 = "$obj_name_2 = [ BEGIN, TRIANGLE_FAN, COLOR, $pur_rgb_r, $pur_rgb_g, $pur_rgb_b,$vertex_list_2 END]\n";
				#
				if ($obj_option eq "obj") {
					print OUTPUT_PML_FILE "$obj_printout_2";
					$obj_load = "cmd.load_cgo($obj_name,\'$obj_name\')\n";
					print OUTPUT_PML_FILE "$obj_name = $obj_name_1 + $obj_name_2\n$obj_load\n";
				} elsif ($obj_option eq "no_obj") {
					print OUTPUT_PML_FILE "$obj_printout_2";
					print OUTPUT_PML_FILE "$obj_name = $obj_name_1 + $obj_name_2\n\n";
					$object_name_total_pur = $object_name_total_pur . " + " . $obj_name;
				}
			#
			# 2. process pyrimidines:
			} elsif ($pdb_residue_name =~ /^(($pyrimidine_match))$/) {
				#print "\n$pdb_chain, $pdb_residue_name, $pdb_residue_number -> PYRIMIDINE\n"; # DEBUG
				open (IN_FILE, "$input_pdb_file");
				while (<IN_FILE>) {
					$line = $_;
					if ($line =~ /^(ATOM|HETATM)\s{0,6}(\d{1,7})\s{0,}(N1|C2|N3|C4|C5|C6)\s{0,}$pdb_residue_name\s{0,}$pdb_chain\s{0,3}$pdb_residue_number\s{0,}(\-{0,1}\d+\.\d+)\s{0,}(\-{0,1}\d+\.\d+)\s{0,}(\-{0,1}\d+\.\d+)\s{0,}(\d{1,3}\.\d{2})\s{0,}(\d{1,3}\.\d{2})/) {
						($pdb_atom_name,$pdb_x,$pdb_y,$pdb_z) = ($3,$4,$5,$6);
						$vertex_list_1 = $vertex_list_1 . " VERTEX, $pdb_x, $pdb_y, $pdb_z,";
						#chop $line; # DEBUG
						#print "\t$pdb_atom_name\t$pdb_x\t$pdb_y\t$pdb_z\t( $line )\n"; # DEBUG
					}
				}
				close (IN_FILE);
				$obj_name_1 = $pdb_chain . "_" . $pdb_residue_name . "_" . $pdb_residue_number;
				$obj_printout_1 = "$obj_name_1 = [ BEGIN, TRIANGLE_FAN, COLOR, $pyr_rgb_r, $pyr_rgb_g, $pyr_rgb_b,$vertex_list_1 END]\n";
				$obj_load_1 = "cmd.load_cgo($obj_name_1,\'$obj_name_1\')\n";
				if ($obj_option eq "obj") {
					print OUTPUT_PML_FILE "$obj_printout_1$obj_load_1\n";
				} elsif ($obj_option eq "no_obj") {
					print OUTPUT_PML_FILE "$obj_printout_1\n";
					$object_name_total_pyr = $object_name_total_pyr . " + " . $obj_name_1;
				}

			}	
		}
	}
	if ($obj_option eq "no_obj") {
		$object_name_total_pur =~ s/^ \+ /filled_pur = /;
		$obj_load_total_pur = "cmd.load_cgo(filled_pur,\'filled_pur')\n";
		print OUTPUT_PML_FILE "$object_name_total_pur\n$obj_load_total_pur\n";
		$object_name_total_pyr =~ s/^ \+ /filled_pyr = /;
		$obj_load_total_pyr = "cmd.load_cgo(filled_pyr,\'filled_pyr')\n";
		print OUTPUT_PML_FILE "$object_name_total_pyr\n$obj_load_total_pyr\n";
	}
	close (OUTPUT_PML_FILE);
	
	print "$output_pml_file\n";
	
	&date_end;

	EXIT_FILLED:;
	exit;
}

sub pdb_scan {
	local($line)=@_;
	if ($line =~ /^(ATOM|HETATM)\s{0,6}(\d{1,7})\s{0,}([\w|\*]+)\s{0,}(\w{1,4})\s{0,}(\w)\s{0,3}(\w+)\s{0,}(\-{0,1}\d+\.\d+)\s{0,}(\-{0,1}\d+\.\d+)\s{0,}(\-{0,1}\d+\.\d+)\s{0,}(\d{1,3}\.\d{2})\s{0,}(\d{1,3}\.\d{2})/) {
		($pdb_atom_number,$pdb_atom_name,$pdb_residue_name,$pdb_chain,$pdb_residue_number,$pdb_x,$pdb_y,$pdb_z,$pdb_b_factor,$pdb_occ) = ($2,$3,$4,$5,$6,$7,$8,$9,$10,$11);
	}
}

sub unpaired_intro {
	@_ = ($unpaired_number,$chain_1,$type_1,$base_1,$atom_1);
	print "\n\n Unpaired nucleotide \# $unpaired_number\n ------------------------------------------------------------------------------\n";
	print " nucleotide: $chain_1 $type_1 $base_1\n";
}

format BP_FORMAT =
pair_@<<<<    :  @>>>>  @>>>>  @>>>>  @>>>>      -  @>>>>  @>>>>  @>>>>  @>>>>
$bp_number, $chain_1, $type_1_full, $base_1, $atom_1, $chain_2, $type_2_full, $base_2, $atom_2
.

format FPDB_FORMAT_BASE =
HETATM@>>>>  @@>>>>>  @>>>    @>>>>>>>@>>>>>>>@>>>>>>>  @>>> @>>>>           @
$fatom_number_start, $fatom, $fres, $fres_number_start, $x_base, $y_base, $z_base, $focc, $fb, $fatom
.

format FPDB_FORMAT_BACKBONE_1 =
HETATM@>>>>  @@>>>>>  @>>>    @>>>>>>>@>>>>>>>@>>>>>>>  @>>> @>>>>           @
$fatom_number_start, $fatom, $fres, $fres_number_start, $x_backbone_1, $y_backbone_1, $z_backbone_1, $focc, $fb, $fatom
.

format FPDB_FORMAT_BACKBONE_2 =
HETATM@>>>>  @@>>>>>  @>>>    @>>>>>>>@>>>>>>>@>>>>>>>  @>>> @>>>>           @
$fatom_number_start, $fatom, $fres, $fres_number_start, $x_backbone_2, $y_backbone_2, $z_backbone_2, $focc, $fb, $fatom
.

format NT_FORMAT_UNPAIRED_PURINES =
unpaired_@<<<<:  @>>>>  @>>>>  @>>>>     N1
$unpaired_nt_number, $pdb_chain, $pdb_residue_name_full, $pdb_residue_number
.

format NT_FORMAT_UNPAIRED_PYRIMIDINES =
unpaired_@<<<<:  @>>>>  @>>>>  @>>>>     N3
$unpaired_nt_number, $pdb_chain, $pdb_residue_name_full, $pdb_residue_number
.
