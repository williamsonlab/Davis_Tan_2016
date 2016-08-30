#!/usr/bin/perl

use Getopt::Long;

###################
# masses_vX.pl #
###################
#
# Take a consolidated peptide list as input
# Input Header:
# seq mod missed num_copy digest_loc protein startres endres digest_notes
#
#
#
# Set up the elemental masses, as well as a few special groups
# a) N terminus gets an additional hydrogen (NH2 - uncharged)
# b) C terminus gets an additional oxygen and hydrogen (COOH - uncharged)
#
# Nterm: NH2 -> +1(H)
# Cterm: COOH -> +1(O) +1(H)
#
# Elemental Masses From:
# (1) http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=html&isotype=some
# (2) http://www.sisweb.com/referenc/source/exactmaa.htm
#
# NEW: Masses from CRC Handbook 2008, 88th Edition

##########
# MASSES #
##########
$c_mass = 12.0000000;
$h_mass = 1.00782503207;
$n_mass = 14.0030740048;
$o_mass = 15.99491461956;
$s_mass = 31.97207100; # Differs slightly vs (2) -> 31.972072
$p_mass = 30.97376163;

$h2_mass = 2.0141017778;
$c13_mass = 13.0033548378;
$n15_mass = 15.0001088982;

$h_delta = $h2_mass - $h_mass;
$c_delta = $c13_mass - $c_mass;
$n_delta = $n15_mass - $n_mass;

$delta{"H"} = $h_delta;
$delta{"C"} = $c_delta;
$delta{"N"} = $n_delta;

$c_avg = 12.0107358967645;
$h_avg = 1.00794075389926;
$n_avg = 14.006703211812;
$o_avg = 15.9994049271305;
$s_avg = 32.064787234926;
$p_avg = 30.97376163;
$n993_avg = 14.9931296539462;

#$na_mass = 22.98976967; # Note: Only one Na isotope exists, so this is a monoisotopic mass
#$k_mass = 38.9637069; # Note: Multiple K isotopes exist, monoisotopic mass used, average mass = 39.0982

###########
# Options #
###########
# These defaults apply to proteins and are replaced if $rna == 1
$modmassfile = "$ENV{'MASSACRE_PATH'}/digest/mod_mass.txt";
$aamassfile = "$ENV{'MASSACRE_PATH'}/digest/aa_mass.txt";
$aaelemfile = "$ENV{'MASSACRE_PATH'}/digest/aa_elem.txt";
# These are user-specified to override the automatic protein/RNA files
$aamassfile_user = "";
$aaelemfile_user = "";

$mz_low = 300; # Lowest ion mass to output
$mz_high = 2000; # Highest ion mass to output
$min_charge = 1; # Minimum ion charge to output
$max_charge = 6; # Maximum ion charge to output

# Define the type of spike (n15mass) to be calculated
# 0 = 100% 15N Spike, the standard thing (E coli)
# 1 = Fancy, for yeast AA spikes, etc
# if($spike_type == 1) -> need to read in a second variable, $spike_def
$spike_type = 0;
$spike_def = "";

$output_avg = 0;

# 0 for protein
# 1 for RNA (-'ve ion mode)
# 2 for RNA (+'ve ion mode)
$rna = 0; 
$cyclic = 0; # 1 for cyclic phosphates in RNA

$pseudo = 0; # Set this to 1 via massacre for pseudo digests to generate id of -1

&GetOptions(
	"minz=i" => \$min_charge,
	"z=i" => \$max_charge,
	"mzlow=i" => \$mz_low,
	"mzhigh=i" => \$mz_high,
	"rna=i" => \$rna,
	"cyclic=i" => \$cyclic,
	"output_avg=i" => \$output_avg,
	"aamassfile=s" => \$aamassfile_user,
	"aaelemfile=s" => \$aaelemfile_user,
	"modmassfile=s" => \$modmassfile,
	"spike_type=i" => \$spike_type,
	"spike_def=s" => \$spike_def,
	"pseudo=i" => \$pseudo
);

if($spike_type && $spike_def eq "")
{
	die "Error, if --spike_type==1 then --spike_def must be defined\n";
}

#
# Default amino acid library files are defined
# If $rna == 1, these are replaced by nucleotide files
#
if($rna >= 1)
{
	$aamassfile = "$ENV{'MASSACRE_PATH'}/digest/nuc_mass.txt";
	$aaelemfile = "$ENV{'MASSACRE_PATH'}/digest/nuc_elem.txt";
}

#
# The user can choose to override these defaults libraries
# They must specify both a mass library (for general mass calculation) and an element library for fancy spike calculations
#
if($aamassfile_user ne "" && $aaelemfile_user ne "")
{
	$aamassfile = $aamassfile_user;
	$aaelemfile = $aaelemfile_user;
}
elsif($aamassfile_user ne "")
{
	die "Error, must specify both --aaelemfile=something and --aamassfile=something\n";
}
elsif($aaelemfile_user ne "")
{
	die "Error, must specify both --aaelemfile=something and --aamassfile=something\n";
}

#
# This array just contains every amino acid/nucleotide defined
#
@aalib = ();

#
# Read in amino acid masses file
# We read in the headerline by default (aa n14mass n15mass)
# Since this is a hash it won't mess things up, and will still function properly if the header is omitted
#
open AAMASS, "$aamassfile" or die "Can't open $aamassfile\n";
while(defined($input=<AAMASS>))
{
	chomp($input);
	@array = split ' ', $input;
	$aa_mass{$array[0]} = $array[1];
	$aa_mass_n15{$array[0]} = $array[2];
	$aa_avg{$array[0]} = $array[3];
	$aa_n993{$array[0]} = $array[4];
	push(@aalib, $array[0]);
}
close AAMASS;

#
# Read in modifications mass file as per above
#
open MODMASS, "$modmassfile" or die "Can't open $modmassfile\n";
while(defined($input=<MODMASS>))
{
	chomp($input);
	@array = split ' ', $input;
	$mod_mass{$array[0]} = $array[1];
	$mod_mass_n15{$array[0]} = $array[2];
	$mod_avg{$array[0]} = $array[3];
	$mod_n993{$array[0]} = $array[4];
}
close MODMASS;

#
# Read in elemental composition of amino acids used for fancy spike
#
open AA_ELEM, "$aaelemfile" or die "Can't open $aaelemfile\n";
$header = <AA_ELEM>;
while(defined($input=<AA_ELEM>))
{
	chomp($input);
	@array = split ' ', $input;
	$aa = shift(@array);
	
	$elem{$aa}{"C"} = shift(@array);
	$elem{$aa}{"H"} = shift(@array);
	$elem{$aa}{"N"} = shift(@array);
	$elem{$aa}{"O"} = shift(@array);
	$elem{$aa}{"S"} = shift(@array);
	$elem{$aa}{"P"} = shift(@array);
}
close AA_ELEM;

#
# Add an additional wildcard amino acid/nucleotide to the library for fancy spike purposes
#
push(@aalib, "X");

#
# %aahash can be used to query the presence of any given letter in the library
#
for($ctr = 0; $ctr <= $#aalib; $ctr++)
{
	$aahash{$aalib[$ctr]} = 1;
}

#
# Parse Spike Definition
#
@spike_array1 = split /\,/, $spike_def;
for($ctr1 = 0; $ctr1 <= $#spike_array1; $ctr1++)
{
	@spike_array2 = split ' ', $spike_array1[$ctr1];
	
	# First element is the amino acid
	$aa = shift(@spike_array2);
	
	# If RNA, translate
	# NO LONGER NECESSARY, RNA has a dedicated nucleotide library
	#if($rna)
	#{
	#	$aa = $aa2nuc{$aa};
	#}
	
	if(!defined($aahash{$aa}))
	{
		die "Error, residue |$aa| not recognized\n";
	}
	
	# Next array elements are pairs of Element_Isotope/Number
	while($#spike_array2 >= 0)
	{
		$element = shift(@spike_array2);
		$number = shift(@spike_array2);
		
		if($aa ne "X")
		{
			#$spike{$aa} += $number*$delta{$element};
			
			if($number eq "*")
			{
				$spike{$aa} += $elem{$aa}{$element}*$delta{$element};
			}
			else
			{
				$spike{$aa} += $number*$delta{$element};
			}
		}
		else
		{
			# Just do <, ignores X as last element in @aalib
			for($ctr2 = 0; $ctr2 < $#aalib; $ctr2++)
			{
				# Star means all elements of that type in that amino acid
				if($number eq "*")
				{
					$spike{"X"}{$aalib[$ctr2]} += $elem{$aalib[$ctr2]}{$element}*$delta{$element};
				}
				# You can also specify N elements for every amino acid
				# eg X N 3 means 3 heavy nitrogen in every amino acid
				# Hard to imagine a situation where this would occur
				else
				{
					$spike{"X"}{$aalib[$ctr2]} += $number*$delta{$element};
				}
			}
		}
	}
}

#
# Print a header including the various parameters used in the digest
#
print "#MZ_LOW $mz_low\n";
print "#MZ_HIGH $mz_high\n";

print "#RNA $rna\n";
print "#CYCLIC $cyclic\n";

print "#MINCHARGE $min_charge\n";
print "#MAXCHARGE $max_charge\n";

print "#OUTPUTAVG $output_avg\n";

print "#SPIKE_TYPE $spike_type\n";
if($spike_type)
{
	print "#SPIKE_DEF $spike_def\n";
}
if($aamassfile_user ne "")
{
	print "#AAMASSFILE $aamassfile\n";
}
if($aaelemfile_user ne "")
{
	print "#AAELEMFILE $aaelemfile\n";
}
if($modmassfile ne "$ENV{'MASSACRE_PATH'}/digest/mod_mass.txt")
{
	print "#MODMASSFILE $modmassfile\n";	
}



# Protein N and C Termini
$nh2_mass = $h_mass;
$cooh_mass = $o_mass + $h_mass;

# RNA 5' and 3' Termini
$p5_mass = $o_mass + $h_mass;
$p3_mass = 0.0;

# Protein N and C Termini AVG
$nh2_avg = $h_avg;
$cooh_avg = $o_avg + $h_avg;

# RNA 5' and 3' Termini AVG
$p5_avg = $o_avg + $h_avg;
$p3_avg = 0.0;

# Regular Phosphate Terminus:
#   3' PO4(-2) or PO4H(-1)
#   2' OH
# Cyclic phosphate terminus
#   PO4(-1)
#
# So for cyclic:
#   1) Lose OH
#   2) Charge is +1 compared to regular version
#   3) |Max Charge| is one less
if($cyclic == 1)
{
	$p3_mass = -1 * $h_mass + -1 * $o_mass;
	$p3_avg = -1 * $h_avg + -1 * $o_avg;
	#print "p3 $p3_mass\n";
}

#
# Carbamidomethylation is modification of cysteine by iodoacetamide
# Reaction described here:
# http://prospector.ucsf.edu/prospector/4.0.7/html/instruct/allman.htm
# -SH -> S(CH2)(CO)(NH2)
# Additional Mass = C2H3NO ~ 24 + 3 + 14 + 16 ~ 57
# http://en.wikipedia.org/wiki/Iodoacetamide
#
$carbamido_mass = 2*$c_mass + 3*$h_mass + $n_mass + $o_mass;
$carbamido_avg = 2*$c_avg + 3*$h_avg + $n_avg + $o_avg;
$carbamido_n993 = 2*$c_avg + 3*$h_avg + $n993_avg + $o_avg;

#print "1H 2H 3H 4H 5H 6H 7H 8H 9H 10H 11H 12H 13H 14H 15H 1Na 2Na 1K 2K\n";
#print "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 1 2 1 2\n";

$line_ctr = 0;

@peptide_array = (); # array containing all the unique peptide lines to be output
@all_masses = (); # array containing all masses - 2D, with first index being charge

#
# Start Reading in Peptide File from STDIN
# This is the consolidated peptide list, with only one entry for each unique peptide
# Multiple copies of a given peptide are listed in the third array element ([2])
#
LINE: while(defined($input=<STDIN>))
{
	chomp($input);

	# Skip comments
	if(substr($input, 0, 1) eq "#")
	{
		print "$input\n";
		next LINE;
	}

	# seq mod missed num_copy digest_loc protein startres endres digest_notes
	@array = split ' ', $input;
	
	#
	# Process header line from input
	#
	if($array[0] eq "seq")
	{
		for($ctr = 0; $ctr <= $#array; $ctr++)
		{
			$name2col{$array[$ctr]} = $ctr;
		}
		
		next;
	}
	
	$seq = $array[$name2col{"seq"}];
	$mod = $array[$name2col{"mod"}];
	$missed = $array[$name2col{"missed"}];
	$num_copy = $array[$name2col{"num_copy"}];
	$digest_loc = $array[$name2col{"digest_loc"}];
	$protein = $array[$name2col{"protein"}];
	$startres = $array[$name2col{"startres"}];
	$endres = $array[$name2col{"endres"}];
	$digest_notes = $array[$name2col{"digest_notes"}];
	#$locationstring = $array[$name2col{"loc"}];
	#$permutation = $array[$name2col{"permutation"}];
	


	#
	# Start to calculate masses - there are four separate things calculated
	# 1) Monoisotopic unlabeled mass
	#		$pep_mass_base
	# 2) Monoisotopic spike mass (by default this is 100% 15N)
	#		$pep_mass_n15_base
	# 3) Average mass for unlabeled/natural abundance
	#		$pep_avg_base
	# 4) Average mass for 99.3% nitrogen labeling
	#		$pep_n993_base
	#
	if($rna == 0)
	{
		$pep_mass_base = $nh2_mass + $cooh_mass;
		# Not affected by 15N or residue-labeling spikes, so is the same as $pep_mass_base
		$pep_mass_n15_base = $nh2_mass + $cooh_mass;
		$pep_avg_base = $nh2_avg + $cooh_avg;
		$pep_n993_base = $nh2_avg + $cooh_avg;
	}
	elsif($rna >= 1)
	{
		$pep_mass_base = $p5_mass + $p3_mass;
		$pep_mass_n15_base = $p5_mass + $p3_mass;
		$pep_avg_base = $p5_avg + $p3_avg;
		$pep_n993_base = $p5_avg + $p3_avg;
	}

	$seqmod = join '', $seq, $mod;
	$seqmod =~ s/\-//g;
	
	# Add in the amino acids to the base peptide masses
	for($ctr = 0; $ctr < length($seqmod); $ctr++)
	{
		$aa = substr($seqmod, $ctr, 1);
		
		# Standard 15N Spike, just read in 15N masses from data structures
		if($spike_type == 0)
		{
			if(defined($aa_mass{$aa}) && defined($mod_mass{$aa}))
			{
				die "Error, |$aa| defined as both residue and modification\n";
			}
			
			# Search for the mass in the amino acids data structures
			if(defined($aa_mass{$aa}))
			{
				$pep_mass_base += $aa_mass{$aa};
				$pep_mass_n15_base += $aa_mass_n15{$aa};
				$pep_avg_base += $aa_avg{$aa};
				$pep_n993_base += $aa_n993{$aa};
			}
			# Search for mass in modifications
			elsif(defined($mod_mass{$aa}))
			{
				$pep_mass_base += $mod_mass{$aa};
				$pep_mass_n15_base += $mod_mass_n15{$aa};
				$pep_avg_base += $mod_avg{$aa};
				$pep_n993_base += $mod_n993{$aa};
			}
			else
			{
				die "Error, could not find mass for |$aa|\n";
			}
		}
		# Complex Spike, Labeled Amino Acids or the like
		# Leave pep_avg and pep_n993 alone
		# Only modify pep_mass_n15 behavior
		elsif($spike_type == 1)
		{
			if(defined($aa_mass{$aa}) && defined($mod_mass{$aa}))
			{
				die "Error, |$aa| defined as both residue and modification\n";
			}
			
			# Search for the mass in the amino acids data structures
			if(defined($aa_mass{$aa}))
			{
				$pep_mass_base += $aa_mass{$aa};
				$pep_mass_n15_base += $aa_mass{$aa}; # Start by adding 14N mass
				$pep_avg_base += $aa_avg{$aa};
				$pep_n993_base += $aa_n993{$aa};
			}
			# Search for mass in modifications
			elsif(defined($mod_mass{$aa}))
			{
				$pep_mass_base += $mod_mass{$aa};
				$pep_mass_n15_base += $mod_mass{$aa}; # Start by adding 14N mass
				$pep_avg_base += $mod_avg{$aa};
				$pep_n993_base += $mod_n993{$aa};
			}
			else
			{
				die "Error, could not find mass for |$aa|\n";
			}
			
			if(defined($spike{$aa}))
			{
				$pep_mass_n15_base += $spike{$aa}; # Add based on defined rules
			}
			elsif(defined($spike{"X"}))
			{
				$pep_mass_n15_base += $spike{"X"}{$aa}; # Add in default labeling
			}
		}
		else
		{
			die "Error, --spike_type=$spike_type unrecognized\n";
		}
	}

	
	# OLD
	#@pep_mass = ( $pep_mass_base, $pep_mass_base, $pep_mass_base, $pep_mass_base );
	#@pep_mass_n15 = ( $pep_mass_n15_base, $pep_mass_n15_base, $pep_mass_n15_base, $pep_mass_n15_base );
	#@pep_avg = ( $pep_avg_base, $pep_avg_base, $pep_avg_base, $pep_avg_base ); 
	#@pep_n993 = ( $pep_n993_base, $pep_n993_base, $pep_n993_base, $pep_n993_base );
	# OLD

	$pep_mass = $pep_mass_base;
	$pep_mass_n15 = $pep_mass_n15_base;
	$pep_avg = $pep_avg_base;
	$pep_n993 = $pep_n993_base;

	
	#
	# Add on hydrogen masses and calculate actual ion charges
	# Protein masses are defined as neutral so leave them alone
	#
	
	#
	# If an RNA oligo, first bring the oligo to neutral charge
	# Each nucleotide residue has a -1 charge, except the 3' which has -2
	#
	# Isn't this kind of ridiculous? Why not just use a neutral-charged library?
	#
	if($rna >= 1)
	{	
		# Cyclic phosphate terminus has -1 charge, not -2
		if($cyclic == 1)
		{
			$start_charge = abs(-1*length($seq));
		}
		else
		{
			$start_charge = abs(-1*length($seq) - 1);
		}

		# Add on hydrogens to bring the oligo to neutral
		$pep_mass += $start_charge*$h_mass;
		$pep_mass_n15 += $start_charge*$h_mass;
		$pep_avg += $start_charge*$h_avg;
		$pep_n993 += $start_charge*$h_avg;
	}
	
	CHARGE: for($ctr = $min_charge; $ctr <= $max_charge; $ctr++)
	{
		$outseq = $seq;
	
		# Treat proteins and RNA (+'ve ion mode) the same way
		if($rna == 0 || $rna == 2)
		{
			# If the min charge is say +3, we have to add 3 H the first time
			if($ctr == $min_charge)
			{
				$pep_mass += $ctr*$h_mass;
				$pep_mass_n15 += $ctr*$h_mass;
				$pep_avg += $ctr*$h_avg;
				$pep_n993 += $ctr*$h_avg;
			}
			else
			{
				$pep_mass += $h_mass;
				$pep_mass_n15 += $h_mass;
				$pep_avg += $h_avg;
				$pep_n993 += $h_avg;
			}
			
			$mz = $pep_mass / $ctr;
			$mz_n15 = $pep_mass_n15 / $ctr;
			$mz_avg = $pep_avg / $ctr;
			$mz_n993 = $pep_n993 / $ctr;
		
			if($ctr == 1)
			{
				$iontype = "M+H";
			}
			else
			{
				$iontype = join '', "M+", $ctr, "H";
			}
		}
		elsif($rna == 1)
		{
			# Only proceed if $ctr (ie current -'ve charge) <= $start_charge
			# We cannot remove hydrogens
			if($ctr > $start_charge)
			{
				next CHARGE;
			}
	
			# If the min charge is say -3, we have to remove 3 H the first time
			if($ctr == $min_charge)
			{
				$pep_mass -= $ctr*$h_mass;
				$pep_mass_n15 -= $ctr*$h_mass;
				$pep_avg -= $ctr*$h_avg;
				$pep_n993 -= $ctr*$h_avg;
			}
			else
			{
				$pep_mass -= $h_mass;
				$pep_mass_n15 -= $h_mass;
				$pep_avg -= $h_avg;
				$pep_n993 -= $h_avg;
			}
			
			$mz = $pep_mass / $ctr;
			$mz_n15 = $pep_mass_n15 / $ctr;			
			$mz_avg = $pep_avg / $ctr;
			$mz_n993 = $pep_n993 / $ctr;

			if($ctr == 1)
			{
				$iontype = "M-H";
			}
			else
			{
				$iontype = join '', "M-", $ctr, "H";
			}
			
			$num_h = $start_charge - $ctr;
			for($ctr_h = 0; $ctr_h < $num_h; $ctr_h++)
			{
				$outseq = join '', $outseq, "H";
			}
		}


		$peptide_line = join ' ', $mz, $mz_n15, $protein, $startres, $endres, $iontype, $ctr, $missed, $outseq, $mod, $digest_loc, $digest_notes;

		if($output_avg == 1)
		{
			$peptide_line = join ' ', $peptide_line, $mz_avg, $mz_n993;
		}

		push(@peptide_array, $peptide_line);
		push(@{ $all_masses[$ctr] }, $mz);
		push(@{ $all_masses[$ctr] }, $mz_n15);
	}
}


#
#$peptide_line[$ctr2] = join ' ', $mz[$ctr2], $mz_n15[$ctr2], @loc_array, $iontype, $ctr, $missed, $outseq, $mod[$ctr2];			
#$peptide_line[$ctr2] = join ' ', $peptide_line[$ctr2], $mz_avg[$ctr2], $mz_n993[$ctr2];
#
# @peptide_array is full of $peptide_line entries
# mz mzN15 loc iontype ctr missed seq mod (mzavg mzavgN15)

# @sort_array is sorted by mz(N14)
@sort_array = sort bymz @peptide_array;

#@sort_array_n15 = sort bymzn15 @peptide_array;

# @all_masses[charge][] contains both n14 and n15 masses
# Two dimensional array with the first being charge
for($ctr = $min_charge; $ctr <= $max_charge; $ctr++)
{
	@{ $sort_masses[$ctr] } = sort bymz2 @{ $all_masses[$ctr] };

	for($ctr2 = 0; $ctr2 <= $#{ $sort_masses[$ctr] }; $ctr2++)
	{
		#print "$ctr $sort_masses[$ctr][$ctr2]\n";
		
		#
		# We make an index of the sort_masses array
		# All peptides are unique, so all masses *might* be unique
		# But maybe not with N15 and various modifications
		# Also, sequence permutations
		# Even so, the index still works as we take +/- 1 for proximity below, which in the case of redundancy gives a proximity of 0
		#
		$sort_hash[$ctr]{ $sort_masses[$ctr][$ctr2] } = $ctr2;
	}
}



#
# Filter for m/z threshold
# Calculate proximity values based on ion types for N14 and N15 masses
# Print
#
print "digestionid digestpepid n14mass n15mass protein startres endres iontype charge missed seq mod digest_loc digest_notes ";
if($output_avg == 1)
{
	print "n14avgmass n993avgmass ";
}
print "n14prox n15prox\n";

$num_ion_out = 0;
$num_pep_out = 0;
for($ctr = 0; $ctr <= $#sort_array; $ctr++)
{
	@array = split ' ', $sort_array[$ctr];

	#
	# Do m/z filtering here, N14 mass >= 300 and N15 mass <= 6000
	# Probably not the most efficient way, but things seem fast enough not to worry
	#
	$iontype = $array[5];
	$n14mass = $array[0];
	$n15mass = $array[1];
	$pep = $array[8];
	if($iontype eq "M+H")
	{
		$charge = 1;
	}
	else
	{
		# Split on the + sign and remove the H from the end
		# This will still work if charge is >= 10
		@array2 = split /\+/, $iontype;
		$charge = $array2[1];
		substr($charge, -1, 1) = "";
	}
	
	if($n14mass >= $mz_low && $n15mass <= $mz_high)
	{
		# Seek lower m/z with the same iontype for proximity
		$n14_lowval = abs($sort_masses[$charge][ $sort_hash[$charge]{$n14mass} - 1 ] - $n14mass);
		$n15_lowval = abs($sort_masses[$charge][ $sort_hash[$charge]{$n15mass} - 1 ] - $n15mass);
		
		# Seek higher m/z with the same iontype for proximity
		$n14_highval = abs($sort_masses[$charge][ $sort_hash[$charge]{$n14mass} + 1 ] - $n14mass);
		$n15_highval = abs($sort_masses[$charge][ $sort_hash[$charge]{$n15mass} + 1 ] - $n15mass);
		
		if($n14_lowval <= $n14_highval)
		{
			$n14_prox = $n14_lowval;
		}
		else
		{
			$n14_prox = $n14_highval;
		}
		
		if($n15_lowval <= $n15_highval)
		{
			$n15_prox = $n15_lowval;
		}
		else
		{
			$n15_prox = $n15_highval;
		}
		
		$num_ion_out++;
		
		if(defined($pep_hash{$pep}))
		{
		}
		else
		{
			$num_pep_out++;
			$pep_hash{$pep} = $num_pep_out;
		}
		
		if($pseudo)
		{
			printf("-1 -1 %s %.6f %.6f\n", $sort_array[$ctr], $n14_prox, $n15_prox);
		}
		else
		{
			printf("%d %d %s %.6f %.6f\n", $num_ion_out, $pep_hash{$pep}, $sort_array[$ctr], $n14_prox, $n15_prox);
		}
	}
}

sub bymz
{
	@array_a = split ' ', $a;
	@array_b = split ' ', $b;
	
	$array_a[0] <=> $array_b[0];
}

sub bymzn15
{
	@array_a = split ' ', $a;
	@array_b = split ' ', $b;
	
	$array_a[1] <=> $array_b[1];
}

sub bymz2
{
	$a <=> $b;
}
