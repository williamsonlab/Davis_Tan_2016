#!/usr/bin/perl

$header = <STDIN>;
@array = split /\,/, $header;
$xmax = $array[0];
$ymax = $array[1];
$mzmin = $array[2];
$mzmax = $array[3];

$dmz = $mzmax - $mzmin;

while(defined($input=<STDIN>))
{
	@array = split /\,/, $input;
	if($#array != 2)
	{
		next;
	}
	
	$x = $array[0];
	$rt = $array[1];
	$int = $array[2];
	
	$mz = $mzmin + ($x/$xmax)*$dmz;
	$int = 100*($int);
	
	$data[10*$mz][$rt] += $int;
}

print "X ";
for($ctr2 = 10*$mzmin; $ctr2 <= 10*$mzmax; $ctr2++)
{
	$mzval = $ctr2 / 10;
	print "$mzval ";
}
print "\n";

for($ctr1 = 0; $ctr1 < $ymax; $ctr1++)
{
	print "$ctr1 ";
	
	for($ctr2 = 10*$mzmin; $ctr2 <= 10*$mzmax; $ctr2++)
	{
		if(defined($data[$ctr2][$ctr1]))
		{
			print "$data[$ctr2][$ctr1] ";
		}
		else
		{
			print "0 ";
		}
	}
	print "\n";
}