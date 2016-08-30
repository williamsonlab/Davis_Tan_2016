#!/usr/bin/perl

use Getopt::Long;
use IO::Handle;
use GD;
require "$ENV{'MASSACRE_PATH'}/msc/mikesubs.pl";

#
# 1) Capture the command from massacre.pl and pass it along to msc04_makepeaks_2blu.pl
# 2) Initiate the above script on bluefish
#    -copy the massive .csv file directly to the temp dir
# 3) Everything completes as normal on bluefish except:
#    a) plot files are not made
#    b) highres contour plots are not made
# 4) 3a and 3b are done locally from returned text files
#
#

$t_start = time();
$pid = $$;
$pwd = `pwd`;
chomp($pwd);
$hostname = `hostname`;
chomp($hostname);

# Note! This is the same mode as used in msc01 - it helps determine the approach used here (ie read in one or two lines at a time)
# mode = 0 -> Both N14 and N15 matches
# mode = 1 -> N14 Only
# mode = 2 -> N15 Only
$matchmode = 0;

# Extract this many units on either side of the N14/N15 monoisotopic peaks
$mz_expand = 4;
# Extract an additional this many units beyond the N15 monoisotopic peak
$mz_iso_expand = 0;
# Include this much in minutes above/below the average retention time
$rt_expand = 0.1;
$rt_expand_lo = 0.0;
$rt_expand_hi = 0.0;
# Threshold for summing points
$mz_bin_tol = 0.005;
# Do we filter based on...
$filter_abundance = 1;
$filter_rt = 0;
$filter_mz = 0;
$abundance_mincut = 2.5;
$rt_mincut = 5;
$rt_maxcut = 50;
$mz_mincut = 300;
$mz_maxcut = 1300;

$blueid = "";
$machine = "garibaldi00";

$keepcon = 0; # Set this to 1 to keep the raw contour plot data
#$spectra3D = 0; # Set this to 1 to extract 3D (non-summed) minispectra
$extractstyle = 0; # This replaces and extends $spectra3D. 0 = mz_bin_tol, 1 binned, 2 = 3D
$extractbin = 0.05; # Size of extraction bins in m/z domain

$ailtype = 0; # This is used to help calculate mz_bin_tol in 2blu

$gnuplot = "gnuplot"; # urocyclon + xtopher both work fine now

&GetOptions("input=s" => \$inputfile, "output=s" => \$outputfile, "csv=s" => \$csvfile, "dir=s" => \$dir, "matchmode=i" => \$matchmode, "filter_abundance=i" => \$filter_abundance, "filter_mz=i" => \$filter_mz, "filter_rt=i" => \$filter_rt, "abundance_min=f" => \$abundance_mincut, "rt_min=f" => \$rt_mincut, "rt_max=f" => \$rt_maxcut, "mz_min=f" => \$mz_mincut, "mz_max=f" => \$mz_maxcut, "id=s" => \$sample_id, "plot=s" => \$plotfile, "blueid=s" => \$blueid, "sortstyle=i" => \$sortstyle, "slice_per_pt=f" => \$slice_per_pt, "mz_per_pt=f" => \$mz_per_pt, "rt_expand=f" => \$rt_expand, "rt_expand_lo=f" => \$rt_expand_lo, "rt_expand_hi=f" => \$rt_expand_hi, "special_range=i" => \$special_range, "special_range_residues=s" => \$special_range_residues, "special_range_values=s" => \$special_range_values, "machine=s" => \$machine, "keepcon=i" => \$keepcon, "ailtype=i" => \$ailtype, "extractstyle=i" => \$extractstyle, "extractbin=f" => \$extractbin);


$qsubfile = join '.', "doblu_m4", $pid;
print "Generating $machine qsub file: $qsubfile\n";
STDOUT->autoflush(1);

open TEMPIN, "$inputfile" or die "Can't open $inputfile\n";
@tempinarray = <TEMPIN>;
if($matchmode == 1 || $matchmode == 2)
{
	$num_spectra = $#tempinarray + 1;
}
else
{
	$num_spectra = ($#tempinarray + 1) / 2;
}
# 1/31/2011 - Removing this $num_spectra - 300 business
# Presumably this was just an offset for walltime calculation
# This should provide a boost to walltime that scales better for low num_spectra
# Having some issues again with orbitrap data
#$num_spectra = $num_spectra - 300;
if($num_spectra < 0)
{
	$num_spectra = 0;
}

#$walltime = 60 + 0.125 * $num_spectra; # Add 0.125 minutes for every spectrum
#$walltime = 60 + 0.25 * $num_spectra; # Updated 3/10/2008, Joan being terminated
#$walltime = 60 + 0.35 * $num_spectra; # Updated 4/17/2008, Gabi being terminated
#$walltime = 60 + 0.40 * $num_spectra; # Updated 5/19/2010, mzdata being terminated
#$walltime = 60 + 0.60 * $num_spectra; # Updated 2/24/2012, Anna being terminated
#$walltime = 60 + 1 * $num_spectra; # Updated 2/26/2012, Anna being terminated
#$walltime = 60 + 2 * $num_spectra; # Updated 2/28/2012, Elizabeth/Tony's TTOF samples being terminated
$walltime = 60 + 10 * $num_spectra; # Updated 3/31/2012, Anna STILL being terminated; garibaldi also running slow
# Should improve speed of the SHIFT loop?

$hours = sprintf("%d", $walltime/60);
$mins  = sprintf("%d", $walltime - 60*$hours);

# VP; make sure that maximum walltime request is under 1000 hours.
if ($hours >= 1000) { ($hours, $mins) = (999, 0); }

`mkdir $dir`;

# QSUB Header Stuff
open QSUB, ">$qsubfile" or die "Error: Can't open $qsubfile\n";
print QSUB "#PBS -l walltime=$hours:$mins:00\n";
print QSUB "#PBS -l cput=$hours:$mins:00\n";
print QSUB "#PBS -l mem=10gb\n";

# Change to the current working dir
print QSUB "cd \$PBS_O_WORKDIR\n\n";

$rofile = join '', $pwd, "/rofile.", $pid;

# Actually run makepeaks - this should run normally(?)
print QSUB "\$MASSACRE_PATH/msc/msc04_makepeaks_2blu.pl --input=$inputfile --output=$outputfile --dir=$dir --csv=$csvfile --matchmode=$matchmode --filter_abundance=$filter_abundance --filter_mz=$filter_mz --filter_rt=$filter_rt --abundance_min=$abundance_mincut --mz_min=$mz_mincut --mz_max=$mz_maxcut --rt_min=$rt_mincut --rt_max=$rt_maxcut --id=$sample_id --plot=$plotfile --pid=$pid --sortstyle=$sortstyle --slice_per_pt=$slice_per_pt --mz_per_pt=$mz_per_pt --rt_expand=$rt_expand --rt_expand_lo=$rt_expand_lo --rt_expand_hi=$rt_expand_hi --special_range=$special_range --special_range_residues=$special_range_residues --special_range_values=$special_range_values --rofile=$rofile --blueid=$blueid --hostname=$hostname --ailtype=$ailtype --extractstyle=$extractstyle --extractbin=$extractbin\n";

# Do some cleanup
print QSUB "/bin/rm -f `basename $csvfile .mzML.gz`.mzML \n";

close QSUB;

$blueline = `qsub $qsubfile`;
print "$blueline";
STDOUT->autoflush(1);

@array = split /\./, $blueline;
$jobid = $array[0];

print "Waiting for job to finish...\n";
STDOUT->autoflush(1);

`touch $rofile`;
open ROFILE, "$rofile" or die "Can't open $rofile\n";

sleep (60);

WAIT: while(1) {
    while(defined($roline=<ROFILE>)) {
        print "$roline";
        STDOUT->autoflush(1);
    }

    my $qstat_out = `qstat $jobid 2>&1`;
    if ($qstat_out =~ /Unknown Job Id/) {
        print "Remove job $jobid finished.\n" and last WAIT;
    } else { sleep 180; }
}

$blueline = `cat doblu_m4.$pid.e$jobid`;
print "Start $machine eout\n";
print "$blueline";
print "End $machine eout\n";
STDOUT->autoflush(1);
#$blueline = `/bin/rm -f doblu_m4.$pid.e$jobid`;
#print "$blueline";
#$blueline = `/bin/rm -f doblu_m4.$pid.o$jobid`;
#print "$blueline";
#$blueline = `/bin/rm -f doblu_m4.$pid`;
#print "$blueline";
#STDOUT->autoflush(1);


open BLUE, "doblu_m4.$pid.e$jobid" or die "Can't open doblu_m4.$pid.e$jobid\n";
while(defined($line=<BLUE>))
{
	print "$line";
}
close BLUE;
STDOUT->autoflush(1);

open BLUO, "doblu_m4.$pid.o$jobid" or die "Can't open doblu_m4.$pid.o$jobid\n";
while(defined($line=<BLUO>))
{
	print "$line";
}
close BLUO;
STDOUT->autoflush(1);

# Turns out there is little use for these plots
# I will leave this if statement in in case they are desired in the future
if(0 == 1)
{
	#
	# Now generate the plots - we have .txt files and .con files
	#
	print "Generating plots\n";
	STDOUT->autoflush(1);

	open TXTLIST, "txtlist.$pid" or die "Can't open txtlist.$pid\n";
	while(defined($txtline=<TXTLIST>))
	{
		chomp($txtline);

		$pngfile = $txtline;
		substr($pngfile, -3, 3) = "png";

		#print "Plotting $pngfile\n";
		#print "Plotting $dir$spectra[$spec_ctr]{'label'}.png\n";
		$gnuplotfile = join '', "gnuplot.tmp.in.", $pid;
		open GNU, ">$gnuplotfile" or die "Can't open gnuplot.tmp.in.$pid\n";
	#	print GNU "set term png medium size 800,600\n"; # Does not work with /sw/bin/gnuplot, use below
	#	print GNU "set xtics out\n";
		print GNU "set term png picsize 800 600\n";
		#print GNU "set output \"$dir$spectra[$spec_ctr]{'label'}.png\"\n";
		print GNU "set output \"$pngfile\"\n";
	#	print GNU "plot \"$dir$spectra[$spec_ctr]{'label'}.txt\" title \"$dir$spectra[$spec_ctr]{'label'}.txt\" with lines\n";
		print GNU "plot \"$txtline\" title \"$txtline\" with lines\n";
		close GNU;

		`$gnuplot < $gnuplotfile`;
		`rm $gnuplotfile`;
	}
	close TXTLIST;
}

print "Generating local contour plots\n";
STDOUT->autoflush(1);

open CONLIST, "conlist.$pid" or die "Can't open conlist.$pid\n";
while(defined($conline=<CONLIST>))
{
	chomp($conline);

	$confile = $conline;
	substr($confile, -3, 3) = "contour.png";

	open CON, "$conline" or die "Can't open $conline\n";
	@conarray = <CON>;
	close CON;
	chomp(@conarray);
	@bounds = split /\,/, shift(@conarray);
	@rect = split /\,/, pop(@conarray);

	$contour_plot_mini = new GD::Image($bounds[0], $bounds[1]);
	for($ctr = 0; $ctr <= $#conarray; $ctr++)
	{
		@xyz = split /\,/, $conarray[$ctr];
		$color = $contour_plot_mini->colorResolve(255-$xyz[2]*255, 255-$xyz[2]*255, 255-$xyz[2]*255);
		$contour_plot_mini->setPixel($xyz[0], $xyz[1], $color);
	}
	$red = $contour_plot_mini->colorResolve(255, 0, 0);
	$contour_plot_mini->rectangle($rect[0], $rect[1], $rect[2], $rect[3], $red);

	if(defined($bounds[2])) # Only proceed if this is defined, ie new data
	{
		$mz_min_plot = $bounds[2];
		$mz_max_plot = $bounds[3];

		for($ctr = sprintf("%d", $mz_min_plot); $ctr <= $mz_max_plot; $ctr++)
		{
			$x = $bounds[0] * ($ctr - $mz_min_plot)/($mz_max_plot - $mz_min_plot);

			# The integer mass is within the red bounds, so plot a tic
			if($x >= $rect[0] && $x <= $rect[2])
			{
				$contour_plot_mini->line($x, $rect[1], $x, $rect[1]+5, $red);
				$contour_plot_mini->string(gdMediumBoldFont, $x, $rect[1]+5, $ctr, $red);
			}
		}
	}

	open CONPNG, ">$confile" or die "Can't open $confile\n";
	binmode CONPNG;
	print CONPNG $contour_plot_mini->png;
	close CONPNG;

	# Now remove the .con file, as it is no longer needed
	# Do this below because it was messing things up for double peaks
	#`rm $conline`;
}
close CONLIST;

if($keepcon == 0)
{
	open CONLIST, "conlist.$pid" or die "Can't open conlist.$pid\n";
	while(defined($conline=<CONLIST>))
	{
		chomp($conline);

		`/bin/rm -f $conline`;
	}
	close CONLIST;
}

#
# Clean Up Local Files
#
`/bin/rm -f doblu_m4.*`;
`/bin/rm -f txtlist.$pid`;
`/bin/rm -f conlist.$pid`;
`/bin/rm -f $rofile`;
#`/bin/rm -f newconlist.$pid`;

$t_end = time();

$t_elapsed = $t_end - $t_start;
$t_elapsed = $t_elapsed / 60;
printf("Elapsed Time: %.1f minutes\n", $t_elapsed);
