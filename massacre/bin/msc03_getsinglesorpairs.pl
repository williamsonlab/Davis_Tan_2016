#!/usr/bin/perl

use Getopt::Long;
use IO::Handle;
require "$ENV{'MASSACRE_PATH'}/msc/mikesubs.pl";

# Note! This is the same mode as used in msc01 - it helps determine the approach used here
# mode = 0 -> Both N14 and N15 matches
# mode = 1 -> N14 Only
# mode = 2 -> N15 Only
$matchmode = 0;

# 0 -> No peptide redundancy
# 1 -> Peptide redundancy in same protein okay
# 2 -> Output all
$pairmode = 0;

# Peaks must be within this range (in minutes) to count as a pair
$rt_thresh = 0.2;

&GetOptions("input=s" => \$inputfile, "output=s" => \$outputfile, "matchmode=i" => \$matchmode, "pairmode=i" => \$pairmode, "rt_thresh=f" => \$rt_thresh);

@input_array = readfile($inputfile);
@header_array = split ' ', shift(@input_array);
$header_line = join ' ', @header_array;
for($ctr = 0; $ctr <= $#header_array; $ctr++)
{
	$name2col{ $header_array[$ctr] } = $ctr;
	$col2name{ $ctr } = $header_array[$ctr];
}

$num_feat = 0;
$num_mult = 0;
open OUTPUT, ">$outputfile" or die "Can't open $outputfile\n";
print OUTPUT "$header_line\n";

#if($pairmode != 0 && $pairmode != 2)
if($pairmode != 2 && $pairmode != 1 && $pairmode != 0) # 6/10/2009 -> Redundancy Stuff
{
	print "Error - pairmode '$pairmode' currently unsupported\n";
	die;
}

# N14 or N15 Only, Eliminate All Redundancies
if(($matchmode == 1 || $matchmode == 2) && $pairmode == 0)
{
	if($matchmode == 1)
	{
		$matchtype = "N14";
	}
	else
	{
		$matchtype = "N15";
	}
	print "Searching for non-redundant $matchtype only matches\n";
	STDOUT->autoflush(1);

	for($ctr = 0; $ctr <= $#input_array; $ctr++)
	{
		@array = split ' ', $input_array[$ctr];
		# Ensure only matches of the appropriate type (N14 or N15) are output
		# This should be done in msc01/step1, this is just a safeguard
		if($matchtype eq $array[ $name2col{"matchtype"} ])
		{
			$ailid = $array[$name2col{"ailid"}];
			push(@{$ailhash{$ailid}}, $input_array[$ctr]);
		}
	}
	
	$numfeat = 0;
	# For each key, only output ones which have arrays of length 1 ($#=0)
	foreach $ailid (keys %ailhash)
	{
		if($#{$ailhash{$ailid}} == 0)
		{
			print OUTPUT "$ailhash{$ailid}[0]\n";
			$numfeat++;
		}
	}
	
	print "Keeping $numfeat $matchtype only matches\n";
}
# N14 or N15 Only, Eliminate Redundancies EXCEPT Within a Single Protein
if(($matchmode == 1 || $matchmode == 2) && $pairmode == 1)
{
	if($matchmode == 1)
	{
		$matchtype = "N14";
	}
	else
	{
		$matchtype = "N15";
	}
	print "Searching for protein-redundant $matchtype only matches\n";
	STDOUT->autoflush(1);

	for($ctr = 0; $ctr <= $#input_array; $ctr++)
	{
		@array = split ' ', $input_array[$ctr];
		# Ensure only matches of the appropriate type (N14 or N15) are output
		# This should be done in msc01/step1, this is just a safeguard
		if($matchtype eq $array[ $name2col{"matchtype"} ])
		{
			$ailid = $array[$name2col{"ailid"}];
			push(@{$ailhash{$ailid}}, $input_array[$ctr]);
		}
	}
	
	$numfeat = 0;
	# For each key, only output ones which have arrays of length 1 ($#=0)
	foreach $ailid (keys %ailhash)
	{
		if($#{$ailhash{$ailid}} == 0)
		{
			print OUTPUT "$ailhash{$ailid}[0]\n";
			$numfeat++;
		}
		else
		{

			@array = split ' ', $ailhash{$ailid}[0];
			$protein_0 = $array[$name2col{"protein"}];
		
			$rflag = 0;
			for($ctr = 1; $ctr <= $#{$ailhash{$ailid}}; $ctr++)
			{
				@array = split ' ', $ailhash{$ailid}[$ctr];
				$protein = $array[$name2col{"protein"}];
				if($protein ne $protein_0)
				{
					$rflag = 1;
				}
			}
			
			# All Same Protein
			if($rflag == 0)
			{
				for($ctr = 0; $ctr <= $#{$ailhash{$ailid}}; $ctr++)
				{
					print OUTPUT "$ailhash{$ailid}[$ctr]\n";
					$numfeat++;
				}
			}
		}
	}
	
	print "Keeping $numfeat $matchtype only matches\n";
}
# N14 or N15 only, Output everything
if(($matchmode == 1 || $matchmode == 2) && $pairmode == 2)
{
	if($matchmode == 1)
	{
		$matchtype = "N14";
	}
	else
	{
		$matchtype = "N15";
	}
	print "Searching for all $matchtype only matches\n";
	STDOUT->autoflush(1);

	for($ctr = 0; $ctr <= $#input_array; $ctr++)
	{
		@array = split ' ', $input_array[$ctr];
		# Ensure only matches of the appropriate type (N14 or N15) are output
		# This should be done in msc01/step1, this is just a safeguard
		if($matchtype eq $array[ $name2col{"matchtype"} ])
		{
			print OUTPUT "$input_array[$ctr]\n";
		}
	}
	
	print "Keeping $ctr $matchtype only matches\n";
}

# N14/N15 Pairs - Eliminate All Redundancies
if($matchmode == 0 && $pairmode == 0)
{
	print "Searching for non-redundant N14/N15 feature pairs\n";
	print "Using RT Pair Threshold of $rt_thresh minutes\n";
	STDOUT->autoflush(1);

	# First, get all pairs
	# Go through list of matches and build up hashes
	for($ctr = 0; $ctr <= $#input_array; $ctr++)
	{
		@array = split ' ', $input_array[$ctr];
		$matchtype = $array[ $name2col{"matchtype"} ];
		$rt = $array[ $name2col{"rt"} ];
		$seq = $array[ $name2col{"seq"} ];
		$mod = $array[ $name2col{"mod"} ];
		$charge = $array[ $name2col{"charge"} ];
		
		$seqmod = join '', $seq, $mod;
		
		# The seqmod (sequence+modifications) and charge uniquely identify the ion
		$key = join '_', $seqmod, $charge;
	
		if($matchtype eq "N14")
		{
			push @{ $n14_line_hash{$key} }, $input_array[$ctr];
			push @{ $n14_rt_hash{$key} }, $rt;
		}
		elsif($matchtype eq "N15")
		{
			push @{ $n15_line_hash{$key} }, $input_array[$ctr];
			push @{ $n15_rt_hash{$key} }, $rt;
		}
		else
		{
			print "Error: matchtype $matchtype not recognized\n";
			die;
		}
	}	
	
	$initial_feat = 0;
	foreach $key (keys %n14_line_hash)
	{
		# For a given key, compare all the N14 matches to N15 matches by RT
		for($ctr1 = 0; $ctr1 <= $#{ $n14_line_hash{$key} }; $ctr1++)
		{
			for($ctr2 = 0; $ctr2 <= $#{ $n15_line_hash{$key} }; $ctr2++)
			{
				$rtdiff = abs($n14_rt_hash{$key}[$ctr1] - $n15_rt_hash{$key}[$ctr2]);
				if($rtdiff <= $rt_thresh)
				{
					# Everything matches, so output the pair
					$initial_pairs[$initial_feat][0] = $key;
					$initial_pairs[$initial_feat][1] = $ctr1;
					$initial_pairs[$initial_feat][2] = $ctr2;
					#print OUTPUT "$n14_line_hash{$key}[$ctr1]\n";
					#print OUTPUT "$n15_line_hash{$key}[$ctr2]\n";	
					$initial_feat++;
				}
			}
		}
	}

	# Set all redun flags to 0 to start
	for($ctr1 = 0; $ctr1 < $initial_feat; $ctr1++)
	{
		$redun[$ctr1] = 0;
	}
	
	# Now, do all vs. all comparison of pairs, and eliminate ones with shared peaks
	for($ctr1 = 0; $ctr1 < $initial_feat; $ctr1++)
	{
		@a1_n14 = split ' ', $n14_line_hash{$initial_pairs[$ctr1][0]}[$initial_pairs[$ctr1][1]];
		@a1_n15 = split ' ', $n15_line_hash{$initial_pairs[$ctr1][0]}[$initial_pairs[$ctr1][2]];

		$ailid_1_n14 = $a1_n14[$name2col{"ailid"}];
		$ailid_1_n15 = $a1_n15[$name2col{"ailid"}];

		for($ctr2 = $ctr1 + 1; $ctr2 < $initial_feat; $ctr2++)
		{
			@a2_n14 = split ' ', $n14_line_hash{$initial_pairs[$ctr2][0]}[$initial_pairs[$ctr2][1]];
			@a2_n15 = split ' ', $n15_line_hash{$initial_pairs[$ctr2][0]}[$initial_pairs[$ctr2][2]];
			
			$ailid_2_n14 = $a2_n14[$name2col{"ailid"}];
			$ailid_2_n15 = $a2_n15[$name2col{"ailid"}];
						
			if($ailid_1_n14 == $ailid_2_n14 ||
				$ailid_1_n15 == $ailid_2_n15 ||
				$ailid_1_n14 == $ailid_2_n15 ||
				$ailid_1_n15 == $ailid_2_n14)
			{
				$redun[$ctr1] = 1;
				$redun[$ctr2] = 1;
			}
		}
	}
	
	$num_feat = 0;
	# Now output the ones that are actually non-redundant
	for($ctr1 = 0; $ctr1 < $initial_feat; $ctr1++)
	{
		if($redun[$ctr1] == 1)
		{
			next;
		}
		
		print OUTPUT "$n14_line_hash{$initial_pairs[$ctr1][0]}[$initial_pairs[$ctr1][1]]\n";
		print OUTPUT "$n15_line_hash{$initial_pairs[$ctr1][0]}[$initial_pairs[$ctr1][2]]\n";
		$num_feat++;		
	}
	
	print "Found $num_feat N14/N15 feature pairs\n";
	STDOUT->autoflush(1);	
}
# N14/N15 Pairs - Eliminate Redundancies EXCEPT Within a Single Protein
# This is tricky.  For N14/N15 Only, you have a single peak to deal with at any time
# In this case you can have networks of peaks connected by pairs
#
#    X   X     X   X  <- Peaks
# 1  0---0
# 2  0---0
# 3      0-----0
# 4  0---------0
# 5            0---0
# In this example, eliminating all redundancies would eliminate all of these as well.
# Even if these are the same protein, we can't really keep them here, as they match to different sets of peaks and this will interfere
#
#    X   X
# 1  0---0
# 2  0---0
# This case might be reasonable to keep
#
# Just do a pairwise comparison, and flag as redundant except when it's the same protein and both peaks match
# For the upper case, 1 and 2 are still flagged redundant based on interacting with 3 etc
# For the bottom case, they are not
if($matchmode == 0 && $pairmode == 1)
{
	print "Searching for protein-redundant N14/N15 feature pairs\n";
	print "Using RT Pair Threshold of $rt_thresh minutes\n";
	STDOUT->autoflush(1);
	
	# First, get all pairs
	# Go through list of matches and build up hashes
	for($ctr = 0; $ctr <= $#input_array; $ctr++)
	{
		@array = split ' ', $input_array[$ctr];
		$matchtype = $array[ $name2col{"matchtype"} ];
		$rt = $array[ $name2col{"rt"} ];
		$seq = $array[ $name2col{"seq"} ];
		$mod = $array[ $name2col{"mod"} ];
		$charge = $array[ $name2col{"charge"} ];
		
		$seqmod = join '', $seq, $mod;
		
		# The seqmod (sequence+modifications) and charge uniquely identify the ion
		$key = join '_', $seqmod, $charge;
	
		if($matchtype eq "N14")
		{
			push @{ $n14_line_hash{$key} }, $input_array[$ctr];
			push @{ $n14_rt_hash{$key} }, $rt;
		}
		elsif($matchtype eq "N15")
		{
			push @{ $n15_line_hash{$key} }, $input_array[$ctr];
			push @{ $n15_rt_hash{$key} }, $rt;
		}
		else
		{
			print "Error: matchtype $matchtype not recognized\n";
			die;
		}
	}	
	
	$initial_feat = 0;
	foreach $key (keys %n14_line_hash)
	{
		# For a given key, compare all the N14 matches to N15 matches by RT
		for($ctr1 = 0; $ctr1 <= $#{ $n14_line_hash{$key} }; $ctr1++)
		{
			for($ctr2 = 0; $ctr2 <= $#{ $n15_line_hash{$key} }; $ctr2++)
			{
				$rtdiff = abs($n14_rt_hash{$key}[$ctr1] - $n15_rt_hash{$key}[$ctr2]);
				if($rtdiff <= $rt_thresh)
				{
					# Everything matches, so output the pair
					$initial_pairs[$initial_feat][0] = $key;
					$initial_pairs[$initial_feat][1] = $ctr1;
					$initial_pairs[$initial_feat][2] = $ctr2;
					#print OUTPUT "$n14_line_hash{$key}[$ctr1]\n";
					#print OUTPUT "$n15_line_hash{$key}[$ctr2]\n";	
					$initial_feat++;
				}
			}
		}
	}

	# Set all redun flags to 0 to start
	for($ctr1 = 0; $ctr1 < $initial_feat; $ctr1++)
	{
		$redun[$ctr1] = 0;
	}
	
	# Now, do all vs. all comparison of pairs, and eliminate ones with shared peaks
	for($ctr1 = 0; $ctr1 < $initial_feat; $ctr1++)
	{
		@a1_n14 = split ' ', $n14_line_hash{$initial_pairs[$ctr1][0]}[$initial_pairs[$ctr1][1]];
		@a1_n15 = split ' ', $n15_line_hash{$initial_pairs[$ctr1][0]}[$initial_pairs[$ctr1][2]];

		$ailid_1_n14 = $a1_n14[$name2col{"ailid"}];
		$ailid_1_n15 = $a1_n15[$name2col{"ailid"}];
		$protein1 = $a1_n14[$name2col{"protein"}];

		for($ctr2 = $ctr1 + 1; $ctr2 < $initial_feat; $ctr2++)
		{
			@a2_n14 = split ' ', $n14_line_hash{$initial_pairs[$ctr2][0]}[$initial_pairs[$ctr2][1]];
			@a2_n15 = split ' ', $n15_line_hash{$initial_pairs[$ctr2][0]}[$initial_pairs[$ctr2][2]];
			
			$ailid_2_n14 = $a2_n14[$name2col{"ailid"}];
			$ailid_2_n15 = $a2_n15[$name2col{"ailid"}];
			$protein2 = $a2_n14[$name2col{"protein"}];
			
			if($protein1 eq $protein2 && $ailid_1_n14 == $ailid_2_n14 && $ailid_1_n15 == $ailid_2_n15)
			{
				# Same protein, and both peaks match, keep
			}
			elsif($ailid_1_n14 == $ailid_2_n14 ||
				$ailid_1_n15 == $ailid_2_n15 ||
				$ailid_1_n14 == $ailid_2_n15 ||
				$ailid_1_n15 == $ailid_2_n14)
			{
				$redun[$ctr1] = 1;
				$redun[$ctr2] = 1;
			}
		}
	}
	
	$num_feat = 0;
	# Now output the ones that are actually non-redundant
	for($ctr1 = 0; $ctr1 < $initial_feat; $ctr1++)
	{
		if($redun[$ctr1] == 1)
		{
			next;
		}
		
		print OUTPUT "$n14_line_hash{$initial_pairs[$ctr1][0]}[$initial_pairs[$ctr1][1]]\n";
		print OUTPUT "$n15_line_hash{$initial_pairs[$ctr1][0]}[$initial_pairs[$ctr1][2]]\n";
		$num_feat++;		
	}
	
	print "Found $num_feat N14/N15 feature pairs\n";
	STDOUT->autoflush(1);	
}
# N14/N15 Pairs - Output all pairs, regardless of redundancies
if($matchmode == 0 && $pairmode == 2)
{
	print "Searching for all N14/N15 feature pairs\n";
	print "Using RT Pair Threshold of $rt_thresh minutes\n";
	STDOUT->autoflush(1);

	# Go through list of matches and build up hashes
	for($ctr = 0; $ctr <= $#input_array; $ctr++)
	{
		@array = split ' ', $input_array[$ctr];
		$matchtype = $array[ $name2col{"matchtype"} ];
		$rt = $array[ $name2col{"rt"} ];
		$seq = $array[ $name2col{"seq"} ];
		$mod = $array[ $name2col{"mod"} ];
		$charge = $array[ $name2col{"charge"} ];
		
		$seqmod = join '', $seq, $mod;
		
		# The seqmod (sequence+modifications) and charge uniquely identify the ion
		$key = join '_', $seqmod, $charge;
	
		if($matchtype eq "N14")
		{
			push @{ $n14_line_hash{$key} }, $input_array[$ctr];
			push @{ $n14_rt_hash{$key} }, $rt;
		}
		elsif($matchtype eq "N15")
		{
			push @{ $n15_line_hash{$key} }, $input_array[$ctr];
			push @{ $n15_rt_hash{$key} }, $rt;
		}
		else
		{
			print "Error: matchtype $matchtype not recognized\n";
			die;
		}
	}	
	
	foreach $key (keys %n14_line_hash)
	{
		# For a given key, compare all the N14 matches to N15 matches by RT
		for($ctr1 = 0; $ctr1 <= $#{ $n14_line_hash{$key} }; $ctr1++)
		{
			for($ctr2 = 0; $ctr2 <= $#{ $n15_line_hash{$key} }; $ctr2++)
			{
				$rtdiff = abs($n14_rt_hash{$key}[$ctr1] - $n15_rt_hash{$key}[$ctr2]);
				if($rtdiff <= $rt_thresh)
				{
					# Everything matches, so output the pair
					print OUTPUT "$n14_line_hash{$key}[$ctr1]\n";
					print OUTPUT "$n15_line_hash{$key}[$ctr2]\n";	
					$num_feat++;
				}
			}
		}
	}

	print "Found $num_feat N14/N15 feature pairs\n";
	STDOUT->autoflush(1);
}

close OUTPUT;
