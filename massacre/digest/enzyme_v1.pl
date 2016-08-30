#!/usr/bin/perl

#
# enzyme_vX.pl
#
# Usage:
# ./enzyme_vX.pl --enzyme=enzymename --miss=X < fasta_list.txt
#
# X = number of possible consecutive missed enzymatic cleavage sites
# fasta_list.txt is a list of fasta formatted sequences
#
# Filenames should be formatted:
# something.protein.fasta
# Where "protein" is used in the output to specify the source protein of the peptides
#
# S1.S1.fasta
# S2.S2.fasta
# etc
#
# 30S.S1.fasta
# 30S.S2.fasta
# etc
#
# Both give the same result
#
# --enzyme
#   trypsin -> trypsin digest for proteins
#   t1 -> RNAse T1 for RNA

use Getopt::Long;
require "$ENV{'MASSACRE_PATH'}/msc/mikesubs.pl";

$max_miss_cut = 4;
$enzyme = "trypsin";

&GetOptions("miss=i" => \$max_miss_cut, "enzyme=s" => \$enzyme);

if($max_miss_cut == -1)
{
	die "Error: Must specify maximum number of missed trypsin cleavages: miss=X\n";
}

$enzyme = lc($enzyme);
if($enzyme ne "trypsin" && $enzyme ne "t1")
{
	die "Error: --enzyme=(trypsin|t1)\n";
}

print "#ENZYME $enzyme\n";
print "#MISS $max_miss_cut\n";

print "protein seq startres endres missed\n";

#
# Takes a list of .fasta file as input
#
while(defined($file=<STDIN>))
{
	chomp($file);

	# Read the .fasta file into an array and put the sequence into $seq
	@filecontents = readfile($file);
	$header = shift(@filecontents);
	$seq = join '', @filecontents;

	# This removes any white space
	@temp_seq_array = split ' ', $seq;
	$seq = join '', @temp_seq_array;

	# Often sequences terminate with an asterix (*) - this removes that
	if(substr($seq, -1, 1) eq "*")
	{
		substr($seq, -1, 1) = "";
	}
	
	#
	# Set up for filenames of format
	# $basename.$protein.fasta
	# $protein is important as it is output later to identify the source of the peptide
	#
	@file_array = split /\./, $file;
	$basename = $file_array[0];	
	$protein = $file_array[1];
			
	@all_array = ();
	#
	# These enzymatic cleavage subroutines are responsible for populating @all_array
	# This is simply an array that contains the indices of the cleavage residues
	# Consider trypsin cleaving after R and K unless it's followed by P
	#
	# 0123456789012345
	# AAARAAAKAAARPAAA
	#    *   *   
	#
	# @all_array = (3,7) The R at index 11 is ignored due to P at index 12
	#
	if($enzyme eq "trypsin")
	{
		do_trypsin();
	}
	elsif($enzyme eq "t1")
	{
		do_t1();
	}
	
	#
	# If the last residue was not already identified as a cleavage site, then we add in a pseudo-cleavage site.
	#
	$length = length($seq);
	$last_index = $length - 1;
	if($all_array[$#all_array] != $last_index)
	{
		$elem = push(@all_array, $last_index);
	}
	
	#
	# Now add a Pseudo-Cleavage Site to the very beginning, before the actual sequence
	#
	$elem = unshift(@all_array, "-1");
	#print "elem: $elem\n";
	#print "@all_array\n";

	
	#
	# Clear the Hash where Peptides are Stored
	#
	%peptide_hash = ();
	
	#
	# $num_trypsin_site is the actual number of internal trypsin sites, so not including the start/end pseudo-sites
	# Even if the last residue is R/K, this still counts as a pseudo-site since it doesn't do anything
	#
	# If we ask for more missed cleavages than sites, scale back the number of missed cleavages
	#
	$num_trypsin_site = $#all_array - 1;
	if($num_trypsin_site < $max_miss_cut)
	{
		$peptide_miss_cut = $num_trypsin_site;
	}
	else
	{
		$peptide_miss_cut = $max_miss_cut;
	}
	
	#
	# Get all the peptides
	#
	
	# Cycle through all possible number of consecutive missed cleavages
	for($ctr = 0; $ctr <= $peptide_miss_cut; $ctr++)
	{
		# Generate peptides containing $ctr consecutive missed cleavages
		for($ctr2 = 0; $ctr2 < ($#all_array - $ctr); $ctr2++)
		{
			#
			# Get the actual peptide sequence
			#
			$peptide = substr($seq, $all_array[$ctr2] + 1,  $all_array[$ctr2 + 1 + $ctr] - $all_array[$ctr2]);
			
			#
			# Join up "peptide sequence", "starting residue", "ending residue", "number missed cleavages"
			#
			$peptide = join ' ', $peptide, $all_array[$ctr2] + 2, $all_array[$ctr2 + 1 + $ctr] + 1, $ctr;

			#
			# Put the peptide/start/end in the hash
			#
			$peptide_hash{$peptide} += 1;
		}
	}
	
	#digest(@all_array);
	
	#
	# Now we should output all the peptides for this particular protein
	# This is where $protein is important to identify the source of the peptide
	#
	foreach $key (keys %peptide_hash)
	{
		print "$protein $key\n";
	}
}

sub do_trypsin
{
	#
	# Trypsin cuts after Arg(R) or Lys(K) except when followed by Pro(P)
	#
	# Find the cut locations, put them in @all_array, cut locations indexed at 0 and cut after the referenced amino acid
	#
	@r_array = ();
	@k_array = ();
	@all_array = ();
	
	$r_num = 0;
	$pos = -1;
	while( ($pos = index($seq, "R", $pos) ) > -1)
	{
		if(substr($seq, $pos+1, 1) ne "P")
		{
			$r_array[$r_num] = $pos;
			$r_num++;
		}
		$pos++;
	}
	
	$k_num = 0;
	$pos = -1;
	while( ($pos = index($seq, "K", $pos) ) > -1)
	{
		if(substr($seq, $pos+1, 1) ne "P")
		{
			$k_array[$k_num] = $pos;
			$k_num++;
		}
		$pos++;
	}
	
	#print "@r_array - @k_array\n";
	
	#
	# Tack @k_array onto the end of @r_array
	#
	$elem = push(@r_array, @k_array);
	#print "elem: $elem\n";
	#print "@r_array\n";
	@all_array = sort numerically @r_array;
	#print "@all_array\n";
}

sub do_t1
{
	# Cut after G
	@g_array = ();
	@all_array = ();
	
	$g_num = 0;
	$pos = -1;
	while( ($pos = index($seq, "G", $pos) ) > -1)
	{
		$g_array[$g_num] = $pos;
		$g_num++;
		$pos++;
	}
	
	@all_array = @g_array;
}

sub numerically
{
	$a <=> $b;
}
