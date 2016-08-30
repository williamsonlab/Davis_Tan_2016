#!/usr/bin/perl

use Getopt::Long;
use IO::Handle;
require "$ENV{'MASSACRE_PATH'}/msc/mikesubs.pl";

&GetOptions("input=s" => \$inputfile, "output=s" => \$outputfile, "rt_min=f" => \$rt_min, "rt_max=f" => \$rt_max);

print "Applying Retention Time Filter\n";
print "Removing Matches: < $rt_min or > $rt_max (minutes)\n";
STDOUT->autoflush(1);

@input_array = readfile($inputfile);
@header_array = split ' ', shift(@input_array);
$header_line = join ' ', @header_array;
for($ctr = 0; $ctr <= $#header_array; $ctr++)
{
	$name2col{ $header_array[$ctr] } = $ctr;
	$col2name{ $ctr } = $header_array[$ctr];
}

open OUTPUT, ">$outputfile" or die "Can't open $outputfile\n";

print OUTPUT "$header_line\n";

$num_removed = 0;
$num_kept = 0;
for($ctr = 0; $ctr <= $#input_array; $ctr++)
{
	@array = split ' ', $input_array[$ctr];

	$rt = $array[ $name2col{ "rt" } ];
	if($rt < $rt_min || $rt > $rt_max)
	{
		$num_removed++;
		next;
	}
	
	$num_kept++;
	print OUTPUT "$input_array[$ctr]\n";
}
close OUTPUT;

print "Removed $num_removed matches\n";
print "Kept $num_kept matches\n";
STDOUT->autoflush(1);
