#!/usr/bin/perl

use Getopt::Long;
require "$ENV{'MASSACRE_PATH'}/msc/mikesubs.pl";

&GetOptions("isoin=s" => \$isoinfile, "csv=s" => \$csvfile);

@isoin = readfile($isoinfile);
@csv = readfile($csvfile);

$isoin_header = shift(@isoin);
$csv_header = shift(@csv);

@isoin_header_array = split ' ', $isoin_header;
@csv_header_array = split /\,/, $csv_header;

for($ctr = 0; $ctr <= $#isoin_header_array; $ctr++)
{
	$isoin_name2col{ $isoin_header_array[$ctr] } = $ctr;
}
for($ctr = 0; $ctr <= $#csv_header_array; $ctr++)
{
	$csv_name2col{ $csv_header_array[$ctr] } = $ctr;
}

for($ctr = 0; $ctr <= $#csv; $ctr++)
{
	@array = split /\,/, $csv[$ctr];
	$isoid = $array[ $csv_name2col{"isoid"} ];
	$goodids{$isoid} = 1;
}

print "$isoin_header\n";
for($ctr = 0; $ctr <= $#isoin; $ctr++)
{
	@array = split ' ', $isoin[$ctr];
	$isoid = $array[ $isoin_name2col{"isoid"} ];
	if($goodids{$isoid} == 1)
	{
		print "$isoin[$ctr]\n";
	}
}
