#!/usr/bin/perl

# Using a relative OD and a frac_lab, determine the pool size based on fit to curves
# Determine both a single pool size and a protein/intermediate pool size

use Getopt::Long;

$od = 1.5;
$frac = 0.37;
$pools = 1; # Num pool

&GetOptions("od=f" => \$od, "frac=f" => \$frac, "pools=i" => \$pools);

$dir = "$ENV{'MASSACRE_PATH'}/bin/pool/bal_od_fifo";

$hostname = `hostname`;
chomp($hostname);

if($hostname ne "urocyclon")
{
	$dir = "$ENV{'MASSACRE_PATH'}/bin/pool/bal_od_fifo";
}

# Open up the default file just to do a one-time OD matching
open FILE, "$dir/bal_000_000.txt" or die "Can't open bal_000_000.txt\n";
$ctr = 0;
$last_diff = 1000;
while(defined($input=<FILE>))
{
	chomp($input);
	@array = split ' ', $input;
	$diff = abs($od - $array[0]);
	if($diff > $last_diff)
	{
		$loc = $ctr - 1;
		last;
	}	
	$ctr++;
	$last_diff = $diff;
}
close FILE;

#print "loc: $loc\n";

# Single pool
for($ctr1 = 0; $ctr1 <= 50; $ctr1++)
{
	$pool1 = sprintf("%03d", $ctr1);

	$file = join '', $dir, "/", "bal_", $pool1, "_000.txt";
	open FILE, "$file" or die "Can't open $file\n";
	
	for($lctr = 0; $lctr <= $loc; $lctr++)
	{
		$input = <FILE>;
	}
	close FILE;
	chomp($input);
	@array = split ' ', $input;
	
	$frac_file = $array[2];

	$diff = abs($frac_file - $frac);

	if($ctr1 == 0)
	{
		$best_pool = $ctr1;
		$best_diff = $diff;
	}
	else
	{
		if($diff < $best_diff)
		{
			$best_pool = $ctr1;
			$best_diff = $diff;
		}
	}
}

print "Single Pool: $best_pool Diff: $best_diff\n";

if($pools <= 1)
{
	exit;
}

# Two pools
for($ctr1 = 0; $ctr1 <= 50; $ctr1++) # Total sum
{
	
	for($ctrp = 0; $ctrp <= $ctr1; $ctrp++)
	{
		$ctri = $ctr1 - $ctrp;

		$pool1 = sprintf("%03d", $ctrp);

		$pool2 = sprintf("%03d", $ctri);
		
		$file = join '', $dir, "/", "bal_", $pool1, "_", $pool2, ".txt";
		#print "$file\n";
		open FILE, "$file" or die "Can't open $file\n";
		
		@farray = <FILE>;
		$input = $farray[$loc];
		close FILE;
		
#		for($lctr = 0; $lctr <= $loc; $lctr++)
#		{
#			$input = <FILE>;
#		}
#		close FILE;
		chomp($input);
		@array = split ' ', $input;
		
		$frac_file = $array[2];
		
		$diff = abs($frac_file - $frac);
		
		if($ctr1 == 0 && $ctr2 == 0)
		{
			$best_pool1 = $ctrp;
			$best_pool2 = $ctri;
			$best_diff = $diff;
		}
		else
		{
			if($diff < $best_diff)
			{
				$best_pool1 = $ctrp;
				$best_pool2 = $ctri;
				$best_diff = $diff;
			}
		}
	}
}

print "Double Pool: $best_pool1 $best_pool2 Diff: $best_diff\n";
