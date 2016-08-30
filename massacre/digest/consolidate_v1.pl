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

print "seq mod missed num_copy loc permutation\n";

foreach $key (keys %data)
{
	$permutation = 0;
	
	$num_copy = $#{$data{$key}} + 1;
	
	$location = "";
	
	for($ctr = 0; $ctr <= $#{$data{$key}}; $ctr++)
	{
		@array = split ' ', $data{$key}[$ctr];
		
		$protein = $array[$name2col{"protein"}];
		$startres = $array[$name2col{"startres"}];
		$endres = $array[$name2col{"endres"}];
		$seq = $array[$name2col{"seq"}]; # Common
		$missed = $array[$name2col{"missed"}]; # Common
		$mod = $array[$name2col{"mod"}]; # Common

		$location_sub = join ',', $protein, $startres, $endres;
		
		if($ctr == 0)
		{
			$location = join '', $location, $location_sub;			
		}
		else
		{
			$location = join ',', $location, $location_sub;			

			# Compare the sequence and mod strings
			# Recall consolidation is based on residue content, not actual sequence
			# Meaning GANR == AAGR
			# So now what we do is detect that permutations are involved, and alphabetize
			if($seq ne $prev_seq || $mod ne $prev_mod)
			{
				$seq = string_sort($seq);
				$mod = string_sort($mod);
				
				$permutation = 1;
			}
		}
		
		$prev_seq = $seq;
		$prev_mod = $mod;
	}
	
	# Since modifications now use a delimiter, a sequence of length 1 has a mod string of ""
	# modify_v3.pl (and above) have been fixed so that this should no longer be a problem
	if($mod eq "")
	{
		$mod = "-";
	}
	
	print "$seq $mod $missed $num_copy $location $permutation\n";
}

sub string_sort
{
	my @vars = @_;
	my $string = shift(@vars);
	my @array = split '', $string;
	@array = sort @array;
	return join '', @array;
}