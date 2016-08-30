#!/usr/bin/perl

use Getopt::Long;

@fracs = ("0.75", "0.50", "0.25");

&GetOptions("file=s" => \$file, "data=s" => \$data);

if($file eq "")
{
	die "Error, must define file... --file=somefile.csv\n";
}

open FILE, "$file" or die "Can't open $file\n";
@file = <FILE>;
close FILE;

chomp(@file);
$header = shift(@file);
@array = split /\,/, $header;
for($ctr = 0; $ctr <= $#array; $ctr++)
{
	$name2col{$array[$ctr]} = $ctr;
}


`$ENV{'MASSACRE_PATH'}/bin/massiso2plotdata.pl --data=$data < $file > $file.1.00.stats`;
@linearray = `$ENV{'MASSACRE_PATH'}/bin/histo.pl --numbins=10 --binsize=0.1 --col=3 --header=1 < $file.1.00.stats | egrep "Average|Numdata"`;
@linearray0 = split ' ', $linearray[0];
@linearray1 = split ' ', $linearray[1];
$average_deviation = $linearray0[2];
$numdata = $linearray1[2];

@linearray = `$ENV{'MASSACRE_PATH'}/bin/histo.pl --numbins=5 --binsize=1 --col=5 --header=1 --toobig=1 < $file.1.00.stats | head -6`;

print "$average_deviation original 1.00 $numdata";
for($ctr = 1; $ctr <= 5; $ctr++)
{
	printf(" %2d", $linearray[$ctr]);
}
print "\n";

@outdata = ();

foreach $name (keys %name2col)
{
	if(substr($name, 0, 6) ne "avg_wr")
	{
		next;
	}
	
	$col = $name2col{$name};

	@sort = sort bycol @file;

	for($ctr = 0; $ctr <= $#fracs; $ctr++)
	{
		$outfile = join '.', $file, $name, $fracs[$ctr];
		open OUT, ">$outfile" or die "Can't open $outfile\n";
		$num = sprintf("%.0f", ($fracs[$ctr] * ($#sort + 1)));
		
		print OUT "$header\n";
		for($ctr2 = 0; $ctr2 < $num; $ctr2++)
		{
			print OUT "$sort[$ctr2]\n";
		}
		
		close OUT;
		
		`$ENV{'MASSACRE_PATH'}/bin/massiso2plotdata.pl --data=$data < $outfile > $outfile.stats`;

		@linearray = `$ENV{'MASSACRE_PATH'}/bin/histo.pl --numbins=10 --binsize=0.1 --col=3 --header=1 < $outfile.stats | egrep "Average|Numdata"`;
		@linearray0 = split ' ', $linearray[0];
		@linearray1 = split ' ', $linearray[1];
		$average_deviation = $linearray0[2];
		$numdata = $linearray1[2];
		
		@linearray = `$ENV{'MASSACRE_PATH'}/bin/histo.pl --numbins=5 --binsize=1 --col=5 --header=1 --toobig=1 < $outfile.stats | head -6`;
		
		$outline = join ' ', $average_deviation, $name, $fracs[$ctr], $numdata;
		for($ctr2 = 1; $ctr2 <= 5; $ctr2++)
		{
			$outline = join ' ', $outline, sprintf("%2d", $linearray[$ctr2]);
		}
		push(@outdata, $outline);
		#print "$average_deviation $name $fracs[$ctr] $numdata\n";
	}
}

@outdata = sort by_frac_sd @outdata;

for($ctr = 0; $ctr <= $#outdata; $ctr++)
{
	print "$outdata[$ctr]\n";
}

sub bycol
{
	@aa = split /\,/, $a;
	@bb = split /\,/, $b;
	
	$aa[$col] <=> $bb[$col];
}

sub by_frac_sd
{
	@aa = split ' ', $a;
	@bb = split ' ', $b;
	
	$aa[2] <=> $bb[2]
		or
	$aa[0] <=> $bb[0];
}
