#!/usr/bin/perl

use Getopt::Long;

# Put data into a histogram format
# This is an attempt at something more general than "bins.pl"

$numbins = 25;      # The number of bins
$binsize = 1;       # The size of the bins
$graphwidth = 70;   # The width of the plotted histogram
$col = 1;			# By default, use the first column of data
$absflag = 1;       # By default, take the absolute value - negatives aren't much fun in this script
$header = 0;		# By default, assume input file has no header
$print_toobig = 0;	# By default, don't print $numtoobig in initial data stream
$csv = 0;			# By default, whitespace-delimited data
$colname = "";		# By default, specify column number, set this to override
$offset = 0;		# By default, do not offset the data
$xval = 0;			# By default, do not include minimum x values in the output

&GetOptions("numbins=i" => \$numbins, "binsize=f" => \$binsize, "graphwidth:i" => \$graphwidth, "col:i" => \$col, "header:i" => \$header, "toobig:i" => \$print_toobig, "csv:i" => \$csv, "colname:s" => \$colname, "offset:f" => \$offset, "xval:i" => \$xval);

$most = 0;          # The size of the largest bin (Don't Change!)

# Initialize the array that holds the counts
for ($ctr = 0; $ctr < $numbins; $ctr++)
{
  $counts[$ctr] = 0;
  $normcounts[$ctr] = 0;
}

if($header == 1)
{
	$headerline = <STDIN>;
	if($csv == 0)
	{
		@harray = split ' ', $headerline;
	}
	elsif($csv == 1)
	{
		@harray = split /\,/, $headerline;
	}
	else
	{
		die "Error, --csv=0 or --csv=1\n";
	}
	
	for($ctr = 0; $ctr <= $#harray; $ctr++)
	{
		$name2col{$harray[$ctr]} = $ctr;
	}
}

if($colname ne "")
{
	if($header == 0)
	{
		die "Error, --header=1 if --colname defined\n";
	}
	
	$col = $name2col{$colname};
	$col += 1; # Convert from index to real column number
}

$numdata = 0;
$numtoobig = 0;
# Put the data in the bins
while(defined($input=<STDIN>))
{
	chomp($input);

	if($csv == 0)
	{
		@array = split ' ', $input;
	}
	elsif($csv == 1)
	{
		@array = split /\,/, $input;
	}
	else
	{
		die "Error, --csv=0 or --csv=1\n";
	}

	$input = $array[$col - 1];

	$input += $offset;

	if($absflag == 1)
	{
		$input = abs($input);
	}

	#####
	# This is data or application specific stuff here
	#####

	#####
	# End data-specific stuff
	#####

	$whichbin = $input / $binsize;

	$counts[$whichbin]++;

	$data[$numdata] = $input;

	$numdata++;
  
	#
	# Some things will be bigger than the specified range - place them in a special bin
	#
	if($whichbin >= $numbins)
	{
		$numtoobig++;
	}
}

# Print out the results
for ($ctr = 0; $ctr < $numbins; $ctr++)
{
	$column = $ctr * $binsize - $offset;

	printf "$counts[$ctr]";

	if($xval)
	{
		print " $column";
	}
	
	print "\n";
}

if($print_toobig == 1)
{
	print "$numtoobig\n";
}

for ($ctr = 0; $ctr < $numbins; $ctr++)
{
	if($counts[$ctr] > $most)
	{
		$most = $counts[$ctr];

		$modebin = $ctr + 1;
	}
}

print "#----------\n";
# Print out the Parameters Used:
print "# Numbins: $numbins\n";
print "# Binsize: $binsize\n";
print "# GrWidth: $graphwidth\n";


for ($ctr = 0; $ctr < $numbins; $ctr++)
{
	$normcounts[$ctr] = ($counts[$ctr] / $most) * $graphwidth;
	print "# $normcounts[$ctr]\n";
}

print "#----------\n";

for ($ctr = 0; $ctr < $numbins; $ctr++)
{
	$column = $ctr * $binsize - $offset;
	$column_end = ($ctr+1) * $binsize;

	# This is problem specific code
	#$column -= 90;

	printf "#%6.1f -< %-6.1f |", $column, $column_end;

	for($ctrb = 0; $ctrb < $normcounts[$ctr]; $ctrb++)
	{
		print "*";
	}

	print "\n";
}

print "#----------\n";
print "#Number Exceeding Data Range: $numtoobig\n";
print "#----------\n";

#####
# Now do some statistics
#####
$sum = 0;

for($ctr = 0; $ctr < $numdata; $ctr++)
{
	$sum += $data[$ctr];
}

$average = $sum / $numdata;

$variance_sum = 0;

for($ctr = 0; $ctr < $numdata; $ctr++)
{
	$variance = ($data[$ctr] - $average)*($data[$ctr] - $average);

	$variance_sum += $variance;
}

$sd = sqrt($variance_sum/($numdata-1));  

printf "# Average: %.3f\n", $average;
printf "# Stand D: %.3f\n", $sd;
printf "# Numdata: %d\n", $numdata;
printf "# ModeBin: %d (Indexed from 1)\n", $modebin;
