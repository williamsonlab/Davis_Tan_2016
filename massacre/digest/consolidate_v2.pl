#!/usr/bin/perl

#
# consolidate_vX.pl
#
# Consolidates the output of modify_vX.pl based on alphabetized sequence
# The alphabetization is important primarily for RNA.  It means that AUGC and GCAU are consolidated as the same peptide since they have an identical mass
# For proteins this condition is much less common but still possible.
# When no redundancy exists, or when redundancy is of identical sequence, then the sequence and modification strings are left intact.
# When the redundancy occurs between permutations of the same set of letters, then the sequence and mod strings are alphabetized and the permutation column is set to 1
#

while(defined($input=<STDIN>))
{
	chomp($input);
	
	if(substr($input, 0, 1) eq "#")
	{
		print "$input\n";
		next;
	}
	
	@array = split ' ', $input;

	if(lc($array[0]) eq "protein")
	{
		for($ctr = 0; $ctr <= $#array; $ctr++)
		{
			$name2col{$array[$ctr]} = $ctr;
		}
		next;
	}
	
	$protein = $array[$name2col{"protein"}];
	$seq = $array[$name2col{"seq"}];
	$startres = $array[$name2col{"startres"}];
	$endres = $array[$name2col{"endres"}];
	$missed = $array[$name2col{"missed"}];
	$mod = $array[$name2col{"mod"}];
	
	$seqmod = join '', $seq, $mod;
	
	$seqmod =~ s/\-//g;
	$seqmod = string_sort($seqmod);
	#@seqmod_array = split '', $seqmod;
	#@seqmod_array = sort @seqmod_array;
	#$seqmod = join '', @seqmod_array;
	
	#print "$seq $mod $seqmod\n";
	
	push(@{$data{$seqmod}}, $input);
}

#print "seq mod missed num_copy loc permutation\n";
print "seq mod missed num_copy digest_loc protein startres endres digest_notes\n";

$unique_id = 0;

KEY: foreach $key (keys %data)
{	
	@notes = ();
	
	push(@notes, $unique_id);
	$unique_id++;
	
	$num_copy = $#{$data{$key}} + 1;
	
	if($num_copy == 1)
	{
		@array = $data{$key}[0];
	}
	
	$location = "";
	
	@protein_array = ();
	@startres_array = ();
	@endres_array = ();
	@seq_array = ();
	@mod_array = ();
	@missed_array = ();
	
	#
	# Basically just extract the relevant information and push it into some arrays
	# These arrays are easily queried with array_compare()
	#
	for($ctr = 0; $ctr <= $#{$data{$key}}; $ctr++)
	{
		#print "CTR: $ctr\n";
		@array = split ' ', $data{$key}[$ctr];

		$protein = $array[$name2col{"protein"}];
		$startres = $array[$name2col{"startres"}];
		$endres = $array[$name2col{"endres"}];
		$seq = $array[$name2col{"seq"}]; # Common
		$missed = $array[$name2col{"missed"}]; # Common
		$mod = $array[$name2col{"mod"}]; # Common
		
		push(@protein_array, $array[$name2col{"protein"}]);
		push(@startres_array, $array[$name2col{"startres"}]);
		push(@endres_array, $array[$name2col{"endres"}]);
		push(@seq_array, $array[$name2col{"seq"}]);
		push(@mod_array, $array[$name2col{"mod"}]);
		push(@missed_array, $array[$name2col{"missed"}]);
		
		$location_sub = join ',', $protein, $startres, $endres;
		
		if($ctr == 0)
		{
			$location = join '', $location, $location_sub;			
		}
		else
		{
			$location = join ',', $location, $location_sub;			
		}
	}
	
	$same_protein = 0;
	$same_startres = 0;
	$same_endres = 0;
	$same_seq = 0;
	$same_mod = 0;
	$same_missed = 0;
	
	#
	# Since many of our MS measurements are on the per-protein level, having redundant peptides from the same protein isn't a big deal, and so we keep this protein.
	#
	if(array_compare(@protein_array))
	{
		$same_protein = 1;
	}
	else
	{		
		#
		# If the proteins are not all the same we need to consider the S4A/S4B case that's so common in yeast (different isoforms). Check for these and generate an S4X entry if necessary.
		#		
		$ab = 0;

		if($num_copy == 2)
		{
			$proteina = $protein_array[0];
			$proteinb = $protein_array[1];
			$ida = substr($proteina, -1, 1);
			$idb = substr($proteinb, -1, 1);
			substr($proteina, -1, 1) = "";
			substr($proteinb, -1, 1) = "";
			
			if($proteina eq $proteinb)
			{
				if($ida eq "A" && $idb eq "B")
				{
					$ab = 1;
					push(@output_notes, "AB");
				}
				elsif($ida eq "B" && $idb eq "A")
				{
					$ab = 1;
					push(@output_notes, "AB");
				}
			}
		}
		
		if($ab == 1)
		{
			$same_protein = 2;
			$protein = join '', $proteina, "X";
		}
		else
		{
			# Wasn't an A/B issue, so use XXX to denote different proteins
			$protein = "XXX";
		}
	}
	
	#
	# If EITHER the startres or endres are different, set both to 0
	# It's hard (impossible?) to imagine a situation where it's possible that one is the same and one is different since we are querying a single peptide of defined length here.
	#
	# One interesting possibility is different proteins (XXX) but the same residue numbers. This is allowed and in longer peptides likely indicative of coincidence or isoforms that do not adopt an A/B nomenclature. In shorter peptides (AR or something) it's probably just a coincidence.
	#
	if(array_compare(@startres_array))
	{
		$same_startres = 1;
	}
	else
	{
		$startres = 0;
		$endres = 0;
	}
	
	if(array_compare(@endres_array))
	{
		$same_endres = 1;
	}
	else
	{
		$startres = 0;
		$endres = 0;
	}
	
	#
	# If the mod strings are not all equal they must be a permutation so alphabetize
	#
	if(array_compare(@mod_array))
	{
		$same_mod = 1;
	}
	else
	{
		push(@notes, "MP");
		$mod = string_sort($mod);
	}
	
	#
	# If the seq strings are not all equal they must be a permutation so alphabetize
	#
	if(array_compare(@seq_array))
	{
		$same_seq = 1;
	}
	else
	{
		push(@notes, "SP");
		$seq = string_sort($seq);
		$mod = string_sort($mod);
	}
	
	#
	# It's possible that the same peptide has a different number of missed cleavages - a terminal Lys in one protein will not indicate a missed cleavage but a central Lys in another protein would
	#
	if(array_compare(@missed_array))
	{
		$same_missed = 1;
	}
	else
	{
		$missed = -1;
	}
	
	$notes_output = join '_', @notes;
	if($notes_output eq "")
	{
		$notes_output = "-";
	}
	
	print "$seq $mod $missed $num_copy $location $protein $startres $endres $notes_output\n";
}

sub string_sort
{
	my @vars = @_;
	my $string = shift(@vars);
	my @array = split '', $string;
	@array = sort @array;
	return join '', @array;
}

sub array_compare
{
	my @array = @_;
	
	my $base_value = shift(@array);
	
	my $sameflag = 1;
	
	for($ctr = 0; $ctr <= $#array; $ctr++)
	{
		if($array[$ctr] ne $base_value)
		{
			$sameflag = 0;
		}
	}
	
	return $sameflag;
}