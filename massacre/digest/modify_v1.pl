#!/usr/bin/perl

#
# modify_vX.pl
#
# ./modify_vX.pl --modmethod=X --modfile=STRING --carbamido=Y --rna=Z
#
# modmethod 0 - only unmodified peptides are output, ignore modfile
# modmethod 1 - only modified peptides are output
# modmethod 2 - unmodified and modified versions are both output
#
# modfile is a file containing a list of modifications
#
# Protein AA Mod
# S5      1  a
# S12     88 b
# S12     88 d
#
# As you can see, multiple modifications to a single residue are supported
#
# carbamido 0 - ignore carbamido modified cysteines
# carbamido 1 - include only carbamido modified cysteines
# carbamido 2 - include carbamido unmodified and modified cysteines
#
# rna 0 - assume protein
# rna 1 - assume RNA, so ignore the carbamido 
#
# OUTPUT
# protein seq startres endres missed mod\n";


use Getopt::Long;
require "$ENV{'MASSACRE_PATH'}/msc/mikesubs.pl";

$modmethod = 0;
$modfile = "";
$rna = 0;
$carbamido = 1;

&GetOptions("modmethod=i" => \$modmethod, "modfile=s" => \$modfile, "carbamido=i" => \$carbamido, "rna=i" => \$rna);

if($modmethod != 0 && $modfile eq "")
{
	die "Error, must specify --modfile=SOMETIME\n";
}

print "#CARBAMIDO $carbamido\n";
print "#MODMETHOD $modmethod\n";
print "#MODFILE $modfile\n";
print "#RNA $rna\n";
print "protein seq startres endres missed mod\n";

#
# Read in the modifications specified in $modfile
#
if($modmethod > 0)
{
	@modfilearray = readfile($modfile);
	@array = split ' ', $modfilearray[0];
	if(lc($array[0]) eq "protein")
	{
		$header = shift(@modfilearray);
	}
	for($ctr = 0; $ctr <= $#modfilearray; $ctr++)
	{
		@array = split ' ', $modfilearray[$ctr];
		$protein = $array[0];
		$aa = $array[1];
		$modification = $array[2];
		
		#print "$protein -> $aa -> $modification\n";
		
		push(@{$mod_data{$protein}[$aa]},$modification); # Hash -> Array -> Array
	}
	
}

#
# Read in the output from enzyme_vX.pl
# protein seq startres endres missed
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
		$header = $input;
		next;
	}
	
	$protein = shift(@array);
	$seq = shift(@array);
	$startres = shift(@array);
	$endres = shift(@array);
	$missed = shift(@array);
	
	@blank_array = ();
	@carb_array = ();
	@mod_array = ();
	@both_array = ();
	$nullstring = "";
	
	@seq_array = split '', $seq;
	
	for($ctr = $startres; $ctr <= $endres; $ctr++)
	{
		$residue = shift(@seq_array);
		
		$bothstring = "";
		$modstring = "";
		$modstring = join '', @{$mod_data{$protein}[$ctr]};
		
		push(@blank_array, $nullstring);
		push(@mod_array, $modstring);

		if($residue eq "C" && $rna == 0)
		{
			push(@carb_array, "c");
			$bothstring = join '', $modstring, "c";
			push(@both_array, $bothstring);
		}
		else
		{
			push(@carb_array, $nullstring);
			push(@both_array, $modstring);
		}
	}
	
	$blank = join '-', @blank_array;
	$carb = join '-', @carb_array;
	$mod = join '-', @mod_array;
	$both = join '-', @both_array;
	
	@output_array = ();
	
	# mod crb
	#  0   0  -> $blank
	#  0   1  -> $carb
	#  0   2  -> $blank $carb
	#  1   0  -> $mod
	#  1   1  -> $both
	#  1   2  -> $mod $both
	#  2   0  -> $blank $mod
	#  2   1  -> $carb $both
	#  2   2  -> $blank $carb $mod $both
	
	if($modmethod == 0 && $carbamido == 0)
	{
		push(@output_array, $blank);
	}
	elsif($modmethod == 0 && $carbamido == 1)
	{
		push(@output_array, $carb);
	}
	elsif($modmethod == 0 && $carbamido == 2)
	{
		push(@output_array, $blank);
		push(@output_array, $carb);
	}
	elsif($modmethod == 1 && $carbamido == 0)
	{
		push(@output_array, $mod);
	}
	elsif($modmethod == 1 && $carbamido == 1)
	{
		push(@output_array, $both);
	}
	elsif($modmethod == 1 && $carbamido == 2)
	{
		push(@output_array, $mod);
		push(@output_array, $both);
	}
	elsif($modmethod == 2 && $carbamido == 0)
	{
		push(@output_array, $blank);
		push(@output_array, $mod);
	}
	elsif($modmethod == 2 && $carbamido == 1)
	{
		push(@output_array, $carb);
		push(@output_array, $both);
	}
	elsif($modmethod == 2 && $carbamido == 2)
	{
		push(@output_array, $blank);
		push(@output_array, $carb);
		push(@output_array, $mod);
		push(@output_array, $both);
	}
	
	for($ctr = 0; $ctr <= $#output_array; $ctr++)
	{
		# Always output the first one
		if($ctr == 0)
		{
			print "$input $output_array[$ctr]\n";
		}
		# Only output the subsequent ones if they are new
		else
		{
			$different = 1;
			for($ctr2 = 0; $ctr2 < $ctr; $ctr2++)
			{
				if($output_array[$ctr] eq $output_array[$ctr2])
				{
					$different = 0;
				}
			}
			
			if($different == 1)
			{
				print "$input $output_array[$ctr]\n";
			}
		}
	}
	
	#print "$input|$mod|$blank|\n";
}
