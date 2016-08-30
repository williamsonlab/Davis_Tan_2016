#!/usr/bin/perl -w

use Getopt::Long;

# Auto-Detect what type of input file is given
# a) Standard .stats/.data file
#		-look for "Protein" as the first field of the first line
#		-continue to parse as a single panel entry
#		-allow multiple datasets in a single file?
# b) Anything else is interpreted as configuration file, as per below

# Sections of Input File
# 1) Global Options
# 2) Define Panels
# 3) Define Datasets

# >variables
# data 30S
# title This is some awesome data
# xtitle Various Proteins
# ytitle Fraction Whatever
# >panels ymin ymax height
# 0.0 1.0 500
# 0.0 1.0 500
# 0.0 1.0 500
# >datasets file panel showpoints color(0 = black, 1+ = colors) ycut1 cutcolor1 ... ycutn cutcolorn
# somefile1.txt 1 1 1
# somefile2.txt 1 1 2
# somefile3.txt 2 1 3
# somefile4.txt 3 1 4
# >lines panel yvalue solid/dash color(0 = g75, 1+ = colors)
# 1 0.2 0 1
# 2 0.3 1 2
# 3 0.4 1 3
# >end

use Getopt::Long;
use Graphics::AquaTerm ':all';
require "$ENV{'MASSACRE_PATH'}/msc/mikesubs.pl";
require "$ENV{'MASSACRE_PATH'}/msc/mikevars.pl";

# Statics
my $pi = 3.14159265359;
my $loff = 125;
my $roff = 50;
my $toff = 75;
my $boff = 125;
my $xint = 30;
my $num_data = 0; # This is tracked automatically by the script (files, not points)

# Colors
my @black = (0, 0, 0);
my @white = (1, 1, 1);
my @g25 = (0.25, 0.25, 0.25);
my @g50 = (0.5, 0.5, 0.5);
my @g75 = (0.75, 0.75, 0.75);

my @red = (1, 0, 0);
my @red2 = (1, 0.75, 0.75);

my @green = (0, 0.75, 0);
my @green2 = (0.75, 1, 0.75);

my @blue = (0, 0, 1);
my @blue2 = (0.75, 0.75, 1);

my @tangerine = (1, 0.5, 0);
my @tangerine2 = (1, 0.9, 0.75);

my @strawberry = (1, 0, 0.5);
my @strawberry2 = (1, 0.75, 0.9);

my @plum = (0.5, 0.0, 0.5);
my @plum2 = (0.75, 0.5, 0.75);

my @asparagus = (0.5, 0.5, 0);
my @asparagus2 = (0.75, 0.75, 0.5);

my @colors1 = ([@black], [@red], [@green], [@blue], [@tangerine], [@strawberry], [@plum], [@asparagus]);
my @colors2 = ([@g75], [@red2], [@green2], [@blue2], [@tangerine2], [@strawberry2], [@plum2], [@asparagus2]);

# Global Options
my $numpfontsize = 12;
my $axislabelfontsize = 18;
my $xaxisticfontsize = 18;
my $yaxisticfontsize = 18;
my $titlefontsize = 24;
my $fontname = "Helvetica";
my $xtitle = "Protein";
my $ytitle = "Fraction Labeled";
my $title = "My Awesome Data";
my $datatype = "";

# Per Panel Options
my $ymin = 0;
my @ymins = ();
my $ymax = 1;
my @ymaxs = ();
my $height = 500;
my @heights = ();

# Per Datafile Options
my $datafile = "stdin";
my @datafiles = ();
my $panel = 0;
my @panels = ();
my $showp = 1;
my @showps = ();
my $datacolor = 0;
my @datacolors = ();
my $ycut = 0.8;
my @ycuts = ();
my $ycutcolor = 0;
my @ycutcolors = ();

# Per Line Options
my @lines = ();
my @linepanels = ();
my @linestyles = ();
my $linecolor = 0;
my @linecolors = ();

# Read in command line options
&GetOptions(
	"ymin=f" => \$ymin,
	"ymax=f" => \$ymax,
	"height=f" => \$height,
	"data=s" => \$datatype,
	"title=s" => \$title
	);

my @input_array = <STDIN>;
chomp(@input_array);

@array = split ' ', $input_array[0];

# Read in a single dataset in the old style
if($array[0] eq "Protein")
{
	# Set plot to a single panel with a single dataset
	$num_panel = 1;
	$num_data = 1;
	
	# Setup arrays from single file input
	push(@ymins, $ymin);
	push(@ymaxs, $ymax);
	push(@heights, $height);
	push(@datafiles, $datafile);
	push(@panels, $panel);
	push(@showps, $showp);
	push(@linecolors, $linecolor);
	push(@datacolors, $datacolor);
	push(@{ $ycuts[0] }, $ycut);
	push(@{ $ycutcolors[0] }, $ycutcolor);
}
# parse a new style input file
else
{
	INPUT_LINE: for($ictr = 0; $ictr <= $#input_array; $ictr++)
	{
		if(substr($input_array[$ictr], 0, 1) eq ">")
		{
			#print "main $input_array[$ictr]\n";
			$entry_type = substr($input_array[$ictr], 0, 4);
			#print "Processing $entry_type\n";
			
			for($ictr2 = $ictr+1; $ictr2 <= $#input_array; $ictr2++)
			{
				if(substr($input_array[$ictr2], 0, 1) eq ">")
				{
					print "$input_array[$ictr2]\n";
					$ictr = $ictr2 - 1;
					last;
				}

				@array = split ' ', $input_array[$ictr2];
				#print "sub ------- $input_array[$ictr2]\n";
				
				if($entry_type eq ">var")
				{
					$variable = shift(@array);
					
					if($variable eq "data")
					{
						$datatype = shift(@array);
					}
					elsif($variable eq "xtitle")
					{
						$xtitle = join ' ', @array;
					}
					elsif($variable eq "ytitle")
					{
						$ytitle = join ' ', @array;
					}
					elsif($variable eq "title")
					{
						$title = join ' ', @array;
					}
					elsif($variable eq "xint")
					{
						$xint = shift(@array);
					}
					elsif($variable eq "xticfontsize")
					{
						$xaxisticfontsize = shift(@array);
					}
					elsif($variable eq "yticfontsize")
					{
						$yaxisticfontsize = shift(@array);
					}
					elsif($variable eq "labelfontsize")
					{
						$axislabelfontsize = shift(@array);
					}
					elsif($variable eq "titlefontsize")
					{
						$titlefontsize = shift(@array);
					}
					else
					{
						print "WARNING: No variable $variable found\n";
					}
				}
				if($entry_type eq ">pan")
				{
					$num_panel++;
					push(@ymins, shift(@array));
					push(@ymaxs, shift(@array));
					push(@heights, shift(@array));
				}
				if($entry_type eq ">dat")
				{
					push(@datafiles, shift(@array));
					# Users index panels starting from 1, so adjust to 0 indexing
					$panel_num = shift(@array);
					$panel_num -= 1;
					push(@panels, $panel_num);
					push(@showps, shift(@array));
					push(@datacolors, shift(@array));
					$nycut = ($#array + 1) / 2;
					for($ycutctr = 0; $ycutctr < $nycut; $ycutctr++)
					{
						#print "Hello! $#datafiles\n";
						#print "$#{ $ycuts[$#datafiles]}\n";
						$ycut = shift(@array);
						push(@{ $ycuts[$#datafiles] }, $ycut);
						$ycutcolor = shift(@array);
						push(@{ $ycutcolors[$#datafiles] }, $ycutcolor);
						#print "$#{ $ycuts[$#datafiles]}\n";
					}
				}
				if($entry_type eq ">lin")
				{
					$panel_num = shift(@array);
					$panel_num -= 1;
					push(@linepanels, $panel_num);
					push(@lines, shift(@array));
					push(@linestyles, shift(@array));
					push(@linecolors, shift(@array));
				}
				if($entry_type eq ">end")
				{
					last INPUT_LINE;
				}
			}
		}
	}
}

if($datatype eq "30S")
{
	%pro2loc = %pro2loc30S;
	#$xint = 30;
}
elsif($datatype eq "30Splot")
{
	%pro2loc = %pro2loc30Splot;
	#$xint = 30;
}
elsif($datatype eq "50S")
{
	%pro2loc = %pro2loc50S;
}
elsif($datatype eq "50Splot")
{
	%pro2loc = %pro2loc50Splot;
}
elsif($datatype eq "70S")
{
	%pro2loc = %pro2loc70S;
}
elsif($datatype eq "30Schlamy")
{
	%pro2loc = %pro2loc30Schlamy;
}
else
{
	die "Datatype: |$datatype| (--data=something) unrecognized\n";
}

# Dynamically generate %loc2pro which is used for x-axis labels
# Set up exclusions for proteins with multiple possible labels (ie S20 vs S20L26)
foreach $protein (keys %pro2loc)
{
	if(($datatype eq "30S" || $datatype eq "30Splot") && $protein eq "S20L26")
	{
		next;
	}

	$loc = $pro2loc{$protein};
	$loc2pro{$loc} = $protein;	
}

@locs = ();
foreach $loc (keys %loc2pro)
{
	push(@locs, $loc);
}

@locs = sort bynum @locs;

# Setup the plot itself
$plot_width = ($locs[$#locs] + 1)*$xint + $loff + $roff;
$plot_height = $boff + $toff;
for($ctr = 0; $ctr <= $#heights; $ctr++)
{
	$plot_height += $heights[$ctr];
}

#print "$plot_width, $plot_height\n";

aqtInit();
aqtOpenPlot(1);
aqtSetPlotSize($plot_width, $plot_height);
aqtSetBackgroundColor(@white);
aqtSetFontname($fontname);

DrawLines();
DrawFrame();
DrawTitles();

for($dctr = 0; $dctr <= $#datafiles; $dctr++)
{
	@data = ();
	if($datafiles[$dctr] eq "stdin")
	{
		@data = @input_array;
	}
	else
	{
		open DFILE, "$datafiles[$dctr]" or die "Can't open $datafiles[$dctr]\n";
		@data = <DFILE>;
		close DFILE;
	}
	
	chomp(@data);
	
	for($dctr2 = 0; $dctr2 <= $#data; $dctr2++)
	{
		@array = split ' ', $data[$dctr2];
		if($array[0] eq "Protein")
		{
			if($dctr2 > 0)
			{
				$num_data++;
			}
			next;
		}
		
		$protein = shift(@array);
		if(defined($pro2loc{$protein}))
		{
			$loc = $pro2loc{$protein};
		}
		else
		{
			$loc = 0;
		}
		#print "$protein -> $loc\n";
		$xloc = $loff + $loc*$xint;
		$avg = shift(@array);
		$sd = shift(@array);
		$junkloc = shift(@array);
		$nump = shift(@array);
		
		if($showps[$dctr])
		{
			for($pointctr = 0; $pointctr < $nump; $pointctr++)
			{
				$val = shift(@array);
				SetupPoint();
			}
		}
		
		DrawAvgSD();
	}
	
	$num_data++;
}

aqtRenderPlot();

sub DrawLines()
{	
	aqtSetLinewidth(2);
	
	$xa = $loff;
	$xb = $plot_width - $roff;
	@xarray = ($xa, $xb);
	
	for($lctr = 0; $lctr <= $#lines; $lctr++)
	{
		if($lines[$lctr] < $ymins[$linepanels[$lctr]] || $lines[$lctr] > $ymaxs[$linepanels[$lctr]])
		{
			next;
		}
	
		aqtSetColor(@{$colors2[$linecolors[$lctr]]});	
	
		$ya = $boff;
		for($pctr2 = $linepanels[$lctr] + 1; $pctr2 < $num_panel; $pctr2++)
		{
			$ya += $heights[$pctr2];
		}
		
		$ya += $heights[$linepanels[$lctr]] * ($lines[$lctr] - $ymins[$linepanels[$lctr]]) / ($ymaxs[$linepanels[$lctr]] - $ymins[$linepanels[$lctr]]);
		@yarray = ($ya, $ya);

		if($linestyles[$lctr]) # Dashed
		{
			@lineDash = (2, 2);
			aqtSetLinestylePattern(\@lineDash, 0);
		}
		else # Solid
		{
			@lineDash = (0, 0);
			aqtSetLinestylePattern(\@lineDash, 0);
		}

		aqtAddPolyline(\@xarray, \@yarray);
	}
}

sub DrawFrame()
{
	@lineDash = (0, 0);
	aqtSetLinestylePattern(\@lineDash, 0);
	aqtSetColor(@black);
	aqtSetLinewidth(2);
	
	$xa = $loff;
	$xb = $plot_width - $roff;
	@xarray = ($xa, $xb, $xb, $xa, $xa);

	# Box for each panel
	for($pctr = 0; $pctr < $num_panel; $pctr++)
	{
		$ya = $boff;
		for($pctr2 = ($pctr+1); $pctr2 < $num_panel; $pctr2++)
		{
			$ya += $heights[$pctr2];
		}
		$yb = $ya + $heights[$pctr];
		@yarray = ($ya, $ya, $yb, $yb, $ya);

		aqtAddPolyline(\@xarray, \@yarray);
	}

	# Tick marks and labels
	for($pctr = 0; $pctr < $num_panel; $pctr++)
	{	
		# x-axis
		aqtSetFontsize($xaxisticfontsize);
		$ya = $boff;
		for($pctr2 = ($pctr+1); $pctr2 < $num_panel; $pctr2++)
		{
			$ya += $heights[$pctr2];
		}
		$yb = $ya - 10;
		@yarray = ($ya, $yb);
		for($pctr2 = 0; $pctr2 <= $#locs; $pctr2++)
		{
			if($locs[$pctr2] == 0)
			{
				next;
			}
			$xa = $loff + $locs[$pctr2]*$xint;
			@xarray = ($xa, $xa);
			aqtAddPolyline(\@xarray, \@yarray);
			
			if($pctr < ($num_panel-1))
			{
				next;
			}
			$label = $loc2pro{$locs[$pctr2]};
			aqtAddLabel($label, $xa, $ya - 15, 90, 2);
		}
		
		# y-axis
		aqtSetFontsize($yaxisticfontsize);
		$xa = $loff;
		$xb = $loff - 10;
		@xarray = ($xa, $xb);
		for($pctr2 = $ymins[$pctr]; $pctr2 <= $ymaxs[$pctr]; $pctr2+=0.1)
		{
			$yb = $ya + $heights[$pctr]*($pctr2-$ymins[$pctr])/($ymaxs[$pctr]-$ymins[$pctr]);
			@yarray = ($yb, $yb);
			aqtAddPolyline(\@xarray, \@yarray);
			
			if(sprintf("%.1f", $pctr2) eq sprintf("%.1f", $ymaxs[$pctr]) && $pctr > 0)
			{
				next;
			}
			$label = sprintf("%3.1f", $pctr2);
			aqtAddLabel($label, $loff-15, $yb, 0, 2);
		}
	}
	
}

sub DrawTitles()
{
	aqtSetColor(@black);
	aqtSetFontsize($axislabelfontsize);
	
	if($xtitle ne "")
	{
		$xloc = $loff + ($plot_width - $loff - $roff)/2;
		$yloc = $boff - 75;
		aqtAddLabel($xtitle, $xloc, $yloc, 0, 1);
	}
	if($ytitle ne "")
	{
		$xloc = $loff - 75;
		$yloc = $boff + ($plot_height - $boff - $toff)/2;
		aqtAddLabel($ytitle, $xloc, $yloc, 90, 1);
	}
	
	if($title ne "")
	{
		aqtSetFontsize($titlefontsize);
		$xloc = $loff + ($plot_width - $loff - $roff)/2;
		$yloc = $plot_height - 50;
		aqtAddLabel($title, $xloc, $yloc, 0, 1);
	}
}

sub DrawAvgSD()
{
	aqtSetLinewidth(3);

	if($loc == 0)
	{
		return;
	}

	if($avg > $ymaxs[$panels[$dctr]] || $avg < $ymins[$panels[$dctr]])
	{
		return;
	}

	aqtSetColor(@{$colors1[$datacolors[$dctr]]});

	# Now override color based on cutoff
	for($ycutctr = 0; $ycutctr <= $#{ $ycuts[$dctr] }; $ycutctr++)
	{
		if($avg < $ycuts[$dctr][$ycutctr])
		{
			aqtSetColor(@{$colors1[   $ycutcolors[$dctr][$ycutctr]  ]});
		}
	}

	my $radius = (1/3)*$xint;

	$ya = $boff; # Baseline of the panel in question
	for($pctr2 = ($panels[$dctr]+1); $pctr2 < $num_panel; $pctr2++)
	{
		$ya += $heights[$pctr2];
	}
	
	$yb = $ya + $heights[$panels[$dctr]]*($avg-$ymins[$panels[$dctr]])/($ymaxs[$panels[$dctr]]-$ymins[$panels[$dctr]]);
	$yc = $ya + $heights[$panels[$dctr]]*(($avg+$sd)-$ymins[$panels[$dctr]])/($ymaxs[$panels[$dctr]]-$ymins[$panels[$dctr]]);
	$yd = $ya + $heights[$panels[$dctr]]*(($avg-$sd)-$ymins[$panels[$dctr]])/($ymaxs[$panels[$dctr]]-$ymins[$panels[$dctr]]);

	$xa = $xloc - $radius;
	$xb = $xloc + $radius;
	
	if($nump == 0)
	{
		return;
	}
	elsif($nump == 1)
	{
		$xpoint = $xloc;
		$ypoint = $yb;
		DrawPoint();
		return;
	}
	
	@xarray = ($xa, $xb);
	@yarray = ($yb, $yb);
	aqtAddPolyline(\@xarray, \@yarray);
	@yarray = ($yc, $yc);
	aqtAddPolyline(\@xarray, \@yarray);
	@yarray = ($yd, $yd);
	aqtAddPolyline(\@xarray, \@yarray);
	
	@xarray = ($xloc, $xloc);
	@yarray = ($yc, $yd);
	
	aqtAddPolyline(\@xarray, \@yarray);
}

sub SetupPoint
{
	if($loc == 0)
	{
		return;
	}

	if($val > $ymaxs[$panels[$dctr]] || $val < $ymins[$panels[$dctr]])
	{
		return;
	}

	$xpoint = $xloc;
	$ypoint = $boff; # Baseline of the panel in question
	for($pctr2 = ($panels[$dctr]+1); $pctr2 < $num_panel; $pctr2++)
	{
		$ypoint += $heights[$pctr2];
	}
	
	$ypoint = $ypoint + $heights[$panels[$dctr]]*($val-$ymins[$panels[$dctr]])/($ymaxs[$panels[$dctr]]-$ymins[$panels[$dctr]]);
	
	aqtSetColor(@{$colors2[$datacolors[$dctr]]});
	DrawPoint();
}

sub DrawPoint()
{
	my $np = 90;
	my $radius = (1/3)*$xint;
	my $npctr;
	my $x;
	my $y;
	my @xparray = ();
	my @yparray = ();
	my $theta;
	
	for($npctr = 0; $npctr < $np; $npctr++)
	{
		$theta = $pi * (($npctr/($np/90)) / 180);
		$x = $radius * cos($theta);
		$y = $radius * sin($theta);
		
		$xparray[$npctr] = $xpoint + $x;
		$yparray[$npctr] = $ypoint + $y;
		
		$xparray[$npctr + 1*$np] = $xpoint - $y;
		$yparray[$npctr + 1*$np] = $ypoint + $x;
		$xparray[$npctr + 2*$np] = $xpoint - $x;
		$yparray[$npctr + 2*$np] = $ypoint - $y;
		$xparray[$npctr + 3*$np] = $xpoint + $y;
		$yparray[$npctr + 3*$np] = $ypoint - $x;	
	}
	
	push(@xparray, $xpoint + $radius);
	push(@yparray, $ypoint);
	
	aqtAddPolyline(\@xparray, \@yparray);
}

sub bynum
{
	$a <=> $b;
}
