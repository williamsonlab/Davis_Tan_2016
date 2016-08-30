#!/usr/bin/perl

use Getopt::Long;
require "$ENV{'MASSACRE_PATH'}/msc/mikesubs.pl";
require "$ENV{'MASSACRE_PATH'}/msc/mikevars.pl";

$frac_lab_type = 0;
$blanks = 0;
$setf = 0;

&GetOptions("data=s" => \$datatype, "flab=i" => \$frac_lab_type, "blanks=i" => \$blanks, "setf=i" => \$setf);

if($datatype eq "30S")
{
	%pro2loc = %pro2loc30S;
}
elsif($datatype eq "50S")
{
	%pro2loc = %pro2loc50S;
}
elsif($datatype eq "70S")
{
	%pro2loc = %pro2loc70S;
}
elsif($datatype eq "30Schlamy")
{
	%pro2loc = %pro2loc30Schlamy;
}
elsif($datatype eq "40Shuman")
{
	%pro2loc = %pro2loc40Shuman;
}
elsif($datatype eq "60Shuman")
{
	%pro2loc = %pro2loc60Shuman;
}
elsif($datatype eq "80Shuman")
{
	%pro2loc = %pro2loc80Shuman;
}
elsif($datatype eq "40Syeast")
{
	%pro2loc = %pro2loc40Syeast;
}
elsif($datatype eq "60Syeast")
{
	%pro2loc = %pro2loc60Syeast;
}
elsif($datatype eq "80Syeast")
{
	%pro2loc = %pro2loc80Syeast;
}
elsif($datatype eq "RNA")
{
	%pro2loc = %pro2locRNA;
}
elsif($datatype eq "proteomic")
{
	# Don't use %pro2loc, one will be generated dynamically later
}
else
{
	die "Datatype: |$datatype| (--data=something) unrecognized\n";
}

@input_array = <STDIN>;
chomp(@input_array);
@header_array = split /\,/, shift(@input_array);
for($ctr = 0; $ctr <= $#header_array; $ctr++)
{
	$name2col{ lc($header_array[$ctr]) } = $ctr;
	$col2name{ $ctr } = $header_array[$ctr];
}

if(defined($name2col{"frac_lab"}))
{
	$frac_lab_flag = 1;
}
else
{
	$frac_lab_flag = 0;
}

# First go through pro2loc and establish the highest loc
$highloc = 0;
foreach $protein (keys %pro2loc)
{
	$loc = $pro2loc{$protein};
	if($loc > $highloc)
	{
		$highloc = $loc;
	}
}

# Go through contents of _iso.csv file line by line
for($ctr = 0; $ctr <= $#input_array; $ctr++)
{
	@array = split /\,/, $input_array[$ctr];
	$protein = $array[ $name2col{'protein'} ];
	if($frac_lab_flag == 1)
	{
		$frac_lab = $array[ $name2col{'frac_lab'} ];
	}
	else
	{
		if($frac_lab_type == 0)
		{
			$frac_lab = $array[ $name2col{'amp_l'} ] / ( $array[ $name2col{'amp_u'} ] + $array[ $name2col{'amp_l'} ] );
		}
		elsif($frac_lab_type == 1)
		{
			$frac_lab = ($array[ $name2col{'amp_u'} ] + $array[ $name2col{'amp_l'} ]) / ($array[ $name2col{'amp_u'} ] + $array[ $name2col{'amp_l'} ] + $array[ $name2col{'amp_f'} ]);
		}
	}

	# This overrides previous frac lab definitions
	if($setf == 1)
	{
		# Frac Lab: L/(U+L)
		$frac_lab = $array[ $name2col{'amp_l'} ] / ( $array[ $name2col{'amp_u'} ] + $array[ $name2col{'amp_l'} ] );
	}
	elsif($setf == 2)
	{
		# Protein Level: (U+L)/(U+L+F)
		$frac_lab = ($array[ $name2col{'amp_u'} ] + $array[ $name2col{'amp_l'} ]) / ($array[ $name2col{'amp_u'} ] + $array[ $name2col{'amp_l'} ] + $array[ $name2col{'amp_f'} ]);
	}
	elsif($setf == 3)
	{
		# Protein Level: (L)/(U+S)
		$frac_lab = $array[ $name2col{'amp_l'} ] / ( $array[ $name2col{'amp_u'} ] + $array[ $name2col{'amp_s'} ] );
	}

	if(defined($pro2loc{$protein}))
	{
		$loc = $pro2loc{$protein};
	}
	else # unrecognized protein, so tack it onto the hash and proceed
	{
		$highloc += 1;
		$pro2loc{$protein} = $highloc;
		$loc = $pro2loc{$protein};
	}
	
	# Keep track of which proteins are in the .csv data
	# We use $loc to keep track of this, rather than protein name, for below reasons
	if(defined($present_proteins{$loc}))
	{
		$present_proteins{$loc} += 1;
	}
	else
	{
		$present_proteins{$loc} = 1;
	}
	
	# We use $loc and not $protein since if multiple protein names map to the same location
	# ...we will want to calculate the avg/SD together
	push( @{ $frac_hash{$loc} }, $frac_lab );

	# We now build $loc2pro directly from pro2loc
	# In the case of multiple proteins matching to one location, this will be overwritten
	# So the end value is simply the last protein name for that location which happens to occur
	#$loc2pro{$loc} = $protein;
}

foreach $protein (keys %pro2loc)
{
	$loc2pro{ $pro2loc{$protein} } = $protein;
}

foreach $loc (keys %frac_hash)
{
	$sum = 0;
	$vals = "";
	for($ctr = 0; $ctr <= $#{ $frac_hash{$loc} }; $ctr++)
	{
		$sum += $frac_hash{$loc}[$ctr];
	}
	$avg = $sum / ($#{ $frac_hash{$loc} } + 1);
	
	$var = 0;
	for($ctr = 0; $ctr <= $#{ $frac_hash{$loc} }; $ctr++)
	{
		$var += ($frac_hash{$loc}[$ctr] - $avg)**2;
		
		# Convert to nicer looking format
		$frac_hash{$loc}[$ctr] = sprintf("%6.4f", $frac_hash{$loc}[$ctr]);
	}
	
	if($#{ $frac_hash{$loc} } == 0)
	{
		$sd = 0;
	}
	else
	{
		$sd = sqrt($var/$#{ $frac_hash{$loc} });
	}

	$vals = join ' ', @{ $frac_hash{$loc} };

	$output_line = sprintf("%10s %6.4f %6.4f %5d %5d", $loc2pro{$loc}, $avg, $sd, $loc, $#{ $frac_hash{$loc} }+1);
	$output_line = join ' ', $output_line, $vals;

	push(@output, $output_line);
}

# Output proteins which have no data associated with them
if($blanks == 1)
{
	foreach $loc (keys %loc2pro)
	{
		if(!defined($frac_hash{$loc}))
		{
			$output_line = sprintf("%10s %6.4f %6.4f %5d %5d", $loc2pro{$loc}, 0, 0, $loc, 0);
			push(@output, $output_line);
		}
	}
}

@output = sort byloc @output;

print "Protein  flab   +/- loc nval vals\n";
for($ctr = 0; $ctr <= $#output; $ctr++)
{
	print "$output[$ctr]\n";
}

sub byloc
{
	@aa = split ' ', $a;
	@ab = split ' ', $b;
	
	$aa[3] <=> $ab[3];
}
