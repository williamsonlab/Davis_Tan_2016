#!/usr/bin/perl -w

#use strict;
use Graphics::AquaTerm ':all';
use Getopt::Long;
use IO::Handle;
require "$ENV{'MASSACRE_PATH'}/msc/mikesubs.pl";
require "$ENV{'MASSACRE_PATH'}/msc/mikevars.pl";

#############
# Variables #
#############
my $pi = 3.14159;

my $xint = 20;
my $yrange = 500;

my $xoff = 125;
my $yoff = 125;

my @black = (0, 0, 0);
my @white = (1, 1, 1);
my @red = (1, 0, 0);
my @green = (0, 1, 0);
my @blue = (0, 0, 1);
my @g25 = (0.25, 0.25, 0.25);
my @g50 = (0.5, 0.5, 0.5);
my @g75 = (0.75, 0.75, 0.75);

my $title = "";
my $xtitle = "";
my $ytitle = "";

my $numfontsize = 12;
my $ticfontsize = 18;
my $labelfontsize = 18;
my $titlefontsize = 24;

&GetOptions("data=s" => \$datatype, "title=s" => \$title, "xtitle=s" => \$xtitle, "ytitle=s" => \$ytitle, "ticfontsize=i" => \$ticfontsize, "labelfontsize=i" => \$labelfontsize, "titlefontsize=i" => \$titlefontsize, "numfontsize=i" => \$numfontsize);

if($datatype eq "30S")
{
	%pro2loc = %pro2loc30S;
	%loc2lab = %loc2lab30S;
	$xint = 30;
}
elsif($datatype eq "50S")
{
	%pro2loc = %pro2loc50S;
	%loc2lab = %loc2lab50S;
}
elsif($datatype eq "70S")
{
	%pro2loc = %pro2loc70S;
	%loc2lab = %loc2lab70S;
}
elsif($datatype eq "30Schlamy")
{
	%pro2loc = %pro2loc30Schlamy;
	%loc2lab = %loc2lab30Schlamy;
}
else
{
	die "Datatype: |$datatype| (--data=something) unrecognized\n";
}

# This one is dynamically generated
foreach $protein (keys %pro2loc)
{
	#print "$protein\n";
	$loc = $pro2loc{$protein};
	push( @{ $loc2pro{$loc} }, $protein );
}

@locs = ();
foreach $loc (keys %loc2pro)
{
	push(@locs, $loc);
}

@locs = sort bynum @locs;

# Here we add an additional $xint so that the last protein is not on the right axis
$plot_width = $locs[$#locs]*$xint + $xint + 2*$xoff;
$plot_height = $yrange + 2*$yoff;

#####################
# Read in Data File #
#####################
@file = <STDIN>;
chomp(@file);
$header = shift(@file);
for($ctr = 0; $ctr <= $#file; $ctr++)
{
	@array = split ' ', $file[$ctr];
	$protein = shift(@array);
	$avg = shift(@array);
	$sd = shift(@array);
	$loc = shift(@array);
	$nval = shift(@array);
	for($ctr2 = 0; $ctr2 <= $#array; $ctr2++)
	{
		push( @{ $vals{$protein} }, $array[$ctr2] );
	}
	$stats{$protein}[0] = $avg;
	$stats{$protein}[1] = $sd;
	$stats{$protein}[2] = $nval;
}

aqtInit();
aqtOpenPlot(1);
aqtSetPlotSize($plot_width, $plot_height);
aqtSetBackgroundColor(@white);
aqtSetPlotTitle($title);

DrawAxes();
DoPoints();
DoAvgSD();
DrawTitles();

aqtRenderPlot();

sub DrawAxes
{
	@xarray = ( $xoff, $plot_width - $xoff, $plot_width - $xoff, $xoff, $xoff );
	@yarray = ( $yoff, $yoff, $plot_height - $yoff, $plot_height - $yoff, $yoff );

	aqtSetColor(@black);
	aqtSetLinewidth(2);
	aqtAddPolyline(\@xarray, \@yarray);
	
	@yarray = ( $yoff - 20, $yoff - 20, $yoff, $yoff, $yoff - 20 );
	aqtAddPolyline(\@xarray, \@yarray);
	
	aqtSetFontsize($ticfontsize);
	aqtSetFontname("Helvetica");
	
	# Tics and Labels
	for($ctr = 0; $ctr <= $#locs; $ctr++)
	{
		if($locs[$ctr] == 0)
		{
			next;
		}

		$xloc = $xoff + $locs[$ctr]*$xint;
		@xarray = ($xloc, $xloc);
		@yarray = ($yoff - 20, $yoff - 30);
		aqtAddPolyline(\@xarray, \@yarray);


		$label = $loc2lab{$locs[$ctr]};
		aqtAddLabel($label, $xloc, $yoff - 35, 90, 2);
	}
	
	for($ctr = 0; $ctr <= 10; $ctr++)
	{
		$frac = $ctr / 10;
		$label = sprintf("%3.1f", $frac);
		$yloc = $yoff + $frac*$yrange;
		@xarray = ($xoff - 10, $xoff);
		@yarray = ($yloc, $yloc);
		aqtAddPolyline(\@xarray, \@yarray);
		aqtAddLabel($label, $xoff - 15, $yloc, 0, 2);
	}
}

sub DoPoints
{
	foreach $protein (keys %vals)
	{
		$loc = $pro2loc{$protein};
		if($loc == 0)
		{
			next;
		}
		$xloc = $xoff + $loc*$xint;
		for($ctr = 0; $ctr <= $#{ $vals{$protein} }; $ctr++)
		{
			$frac_lab = $vals{$protein}[$ctr];
			$yloc = $yoff + $frac_lab*$yrange;
			DrawPoint();
		}
	}
}

sub DoAvgSD
{
	foreach $protein (keys %stats)
	{
		$loc = $pro2loc{$protein};
		if($loc == 0)
		{
			next;
		}
		$xloc = $xoff + $loc*$xint;
		$avg = $stats{$protein}[0];
		$sd = $stats{$protein}[1];
		$nval = $stats{$protein}[2];
		DrawAvgSD();
	}
}

sub DrawPoint
{
	my $np = 90;
	my $radius = (1/3)*$xint;
	my $lctr;
	my $x;
	my $y;
	my @xarray = ();
	my @yarray = ();
	
	for($lctr = 0; $lctr < $np; $lctr++)
	{
		$theta = $pi * (($lctr/($np/90)) / 180);
		$x = $radius * cos($theta);
		$y = $radius * sin($theta);
		
		$xarray[$lctr] = $xloc + $x;
		$yarray[$lctr] = $yloc + $y;
		
		$xarray[$lctr + 1*$np] = $xloc - $y;
		$yarray[$lctr + 1*$np] = $yloc + $x;
		$xarray[$lctr + 2*$np] = $xloc - $x;
		$yarray[$lctr + 2*$np] = $yloc - $y;
		$xarray[$lctr + 3*$np] = $xloc + $y;
		$yarray[$lctr + 3*$np] = $yloc - $x;		
	}

	push(@xarray, $xloc + $radius);
	push(@yarray, $yloc);

	aqtSetColor( @g75 );
	aqtSetLinewidth(2);
	aqtAddPolyline(\@xarray, \@yarray);
}

sub DrawAvgSD
{
	my $radius = (1/3)*$xint;
	my $yloc = $yoff + $avg*$yrange;
	
	aqtSetColor(@g25);
	aqtSetLinewidth(3);
	
	if($nval >= 2)
	{
		@xarray = ( $xloc - $radius, $xloc + $radius );
		@yarray = ( $yloc, $yloc );
		aqtAddPolyline(\@xarray, \@yarray);
		@yarray = ( $yloc + $sd*$yrange, $yloc + $sd*$yrange );
		aqtAddPolyline(\@xarray, \@yarray);
		@yarray = ( $yloc - $sd*$yrange, $yloc - $sd*$yrange );
		aqtAddPolyline(\@xarray, \@yarray);	
		@xarray = ( $xloc, $xloc);
		@yarray = ( $yloc - $sd*$yrange, $yloc + $sd*$yrange );
		aqtAddPolyline(\@xarray, \@yarray);
	}
	
	aqtSetColor(@g50);
	aqtSetFontsize($numfontsize);
	aqtAddLabel($nval, $xloc, $yoff - 10, 0, 1);
}

sub DrawTitles
{
	aqtSetColor(@black);
	aqtSetFontsize($labelfontsize);
	aqtSetFontname("Helvetica");
	if($xtitle ne "")
	{
		$xloc = $plot_width / 2;
		$yloc = $xoff - 105;
		aqtAddLabel($xtitle, $xloc, $yloc, 0, 1);
	}
	if($ytitle ne "")
	{
		$xloc = $yoff - 75;
		$yloc = $plot_height / 2;
		aqtAddLabel($ytitle, $xloc, $yloc, 90, 1);
	}
	aqtSetFontsize($titlefontsize);
	if($title ne "")
	{
		$xloc = $plot_width / 2;
		$yloc = $plot_height - 50;
		aqtAddLabel($title, $xloc, $yloc, 0, 1);
	}
}


sub bynum
{
	$a <=> $b;
}
