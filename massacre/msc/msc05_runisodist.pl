#!/usr/bin/perl

use Getopt::Long;
use IO::Handle;
require "$ENV{'MASSACRE_PATH'}/msc/mikesubs.pl";

# 1) Parse the rather large input file into the standard isodist input format
# 2) Copy all needed files to bluefish
# 3) Issue the qsub command on bluefish to start the job
#    a) Dynamic walltime calculation based on the number of fits to be done
# 4) Monitor the job for completion
# 5) Retrieve the data
# 6) Do gnuplot stuff locally? - advantage of being able to control gnuplot installation easily

$pid = $$; # use this as a unique ID - as long as one person isn't running on two machines, should be fine
$pwd = `pwd`;
chomp($pwd);
$hostname = `hostname`;
chomp($hostname);

$gnuplot = "gnuplot"; # urocyclon + xtopher both work fine now

$fit_model = 3;
$niter = 4;
$sig_global = 100.0;
$baseline = 1.0;
$ul_amp = 1.2;
$l_amp = 0.5;
$offset = 0.01;
$gw = 0.0003;
$frac_lab = 0.82;
$nround = 5;
$nfitpar = 6;
$fitschedule = "1 0 0 1 1 0 1 0 0 0 0 1 1 0 1 1 1 1 1 1 1 0 0 1 1 1 1 1 1 1";

$machine = "garibaldi00";

&GetOptions("input=s" => \$inputfile, "output=s" => \$outputfile, "fit_model=i" => \$fit_model, "baseline=f" => \$baseline, "ul_amp=f" => \$ul_amp, "l_amp=f" => \$l_amp, "offset=f" => \$offset, "gw=f" => \$gw, "frac_lab=f" => \$frac_lab, "nround=i" => \$nround, "nfitpar=i" => \$nfitpar, "fitschedule=s" => \$fitschedule, "dir=s" => \$dir, "fitsdir=s" => \$fitsdir, "blueid=s" => \$blueid);

# Need to make sure this directory exists, in case it has not been created yet
`mkdir $fitsdir`;

print "Reading $inputfile\n";
STDOUT->autoflush(1);
@input_array = readfile($inputfile);
@header_array = split ' ', shift(@input_array);
$header_line = join ' ', @header_array;
$header_line_csv = join ',', @header_array;
for($ctr = 0; $ctr <= $#header_array; $ctr++)
{
	$name2col{ $header_array[$ctr] } = $ctr;
	$col2name{ $ctr } = $header_array[$ctr];
}

$tempbatchfile = join '.', $inputfile, $pid;
print "Generating isodist batch input file: $tempbatchfile\n";
STDOUT->autoflush(1);
open TEMPBATCH, ">$tempbatchfile" or die "Error: Can't open $tempbatchfile\n";

for($ctr = 0; $ctr <= $#input_array; $ctr++)
{
	@array = split ' ', $input_array[$ctr];
	$file = $array[ $name2col{"file"} ];
	# strip off the original $dir from the front, and replace with $fitsdir
	@file_array = split /\//, $file;
	$file_array[0] = $fitsdir;
	$file = join "/", @file_array;

	$charge = $array[ $name2col{"charge"} ];
	$seq = $array[ $name2col{"seqmod"} ]; # Start using this once isodist can handle all modifications
	#$seq = $array[ $name2col{"seq"} ];

	print TEMPBATCH "$seq $charge $file\n";
}

close TEMPBATCH;

# Template for the batch config file that we will generate here
#####
#fitit (specify routine to run)
#          15N_pulse_short.batch (batchfile)
#          3 (fit_model 1=Natural Abundance, 2=Fractional Label, 3=Unlabeled + Labeled)
#          4 (niter -> number of iterations?)
#          100.0 (sig_global)
#          1.0 (baseline)
#          1.2 (unlabeled_amp)
#          0.5 (labeled_amp)
#          0.01 (offset)
#          0.0003 (gaussian_width)
#          0.82 (fraction_labeled)
#5 (nround -> number of fitting rounds, below)
#6 (nfitpar -> number of parameters to fit, below)
#B OFF GW UL_AMP LAB_AMP FRAC (fit schedule - next lines indicate whether to fit each variable on that pass
#1 0 0 1 1 0
#1 0 0 0 0 1
#1 0 1 1 1 1
#1 1 1 0 0 1
#1 1 1 1 1 1
#####


$tempinfile = join '.', "isodist_batch.in", $pid;
print "Generating isodist config file: $tempinfile\n";
STDOUT->autoflush(1);
open TEMPIN, ">$tempinfile" or die "Error: Can't open $tempinfile\n";

print TEMPIN "fitit\n";
print TEMPIN "          $tempbatchfile\n";
print TEMPIN "          $fit_model\n";
print TEMPIN "          $niter\n";
#print TEMPIN "          $sig_global\n";
printf(TEMPIN "          %.6f\n", $sig_global);
#print TEMPIN "          $baseline\n";
printf(TEMPIN "          %.6f\n", $baseline);
#print TEMPIN "          $ul_amp\n";
printf(TEMPIN "          %.6f\n", $ul_amp);
#print TEMPIN "          $l_amp\n";
printf(TEMPIN "          %.6f\n", $l_amp);
#print TEMPIN "          $offset\n";
printf(TEMPIN "          %.6f\n", $offset);
#print TEMPIN "          $gw\n";
printf(TEMPIN "          %.6f\n", $gw);
#print TEMPIN "          $frac_lab\n";
printf(TEMPIN "          %.6f\n", $frac_lab);
print TEMPIN "$nround\n";
print TEMPIN "$nfitpar\n";
print TEMPIN "B OFF GW UL_AMP LAB_AMP FRAC\n";

@array = split ' ', $fitschedule;
for($ctr1 = 0; $ctr1 < $nround; $ctr1++)
{
	for($ctr2 = 0; $ctr2 < $nfitpar; $ctr2++)
	{
		print TEMPIN "$array[$ctr1*$nfitpar + $ctr2] ";
	}
	print TEMPIN "\n";
}

close TEMPIN;

$npeptide = $#input_array + 1;
$min_per_pep = 0.5;
$time = $npeptide * $min_per_pep;
$hours = sprintf("%d", $time/60);
$mins = sprintf("%d", $time - 60*$hours);
if($hours == 0 && $mins == 0)
{
	$mins = 5;
}

$qsubfile = join '.', "doblu_m5", $pid;
print "Generating $machine qsub file: $qsubfile\n";
STDOUT->autoflush(1);

open QSUB, ">$qsubfile" or die "Error: Can't open $qsubfile\n";
print QSUB "#PBS -l walltime=$hours:$mins:00\n";
print QSUB "#PBS -l cput=$hours:$mins:00\n\n";

print QSUB "cd \$PBS_O_WORKDIR\n\n";

print QSUB "mkdir $fitsdir\n";

close QSUB;

print "Copying files to $machine\n";
STDOUT->autoflush(1);

$blueline = `mkdir $pid/`;
print "$blueline";
STDOUT->autoflush(1);

print "Running job on $machine\n";
STDOUT->autoflush(1);

$blueline = `qsub $pid/$qsubfile`;
print "$blueline";

@array = split /\./, $blueline;
$jobid = $array[0];

print "Waiting for job to finish...\n";
STDOUT->autoflush(1);

$sleeptime = 0;

sleep (60);

WAIT: while(1) {
    my $qstat_out = `qstat $jobid 2>&1`;

    if ($qstat_out =~ /Unknown Job Id/) {
        print "Remote job $jobid finished.\n" and last WAIT;
    } else {
        $sleeptime_min = $sleeptime / 60;
        printf("... %d minutes on $machine\n", $sleeptime_min);
        STDOUT->autoflush(1);
    }

    sleep (300); # Wait 5 minutes
    $sleeptime += 300;
}

#
# Once we determine job has finished, copy files back to our local directory
# ( this is all now done in the QSUB script )
#print "Copying files back to local machine\n";
#STDOUT->autoflush(1);

#$blueline = `scp -r \"$blueid\@bluefish.scripps.edu:$pid/$dir/*.fit\" $dir/`;
#print "$blueline";
#$blueline = `scp -r \"$blueid\@bluefish.scripps.edu:$pid/$dir/*.dat\" $dir/`;
#print "$blueline";
#$blueline = `scp $blueid\@bluefish.scripps.edu:$pid/$inputfile.csv.$pid .`;
#print "$blueline";
#STDOUT->autoflush(1);

#
# Now we need to generate some plots using gnuplot
#
#print "Generating plots\n";
#STDOUT->autoflush(1);
#for($ctr = 0; $ctr <= $#input_array; $ctr++)
#{
#	@array = split ' ', $input_array[$ctr];
#	$file = $array[ $name2col{"file"} ];
#
#	$fitfile = $file;
#	substr($fitfile, -3, 3) = "fit";
#	$datfile = $file;
#	substr($datfile, -3, 3) = "dat";
#	$pngfile = $file;
#	substr($pngfile, -3, 3) = "fit.png";
#
#	$gnuplotfile = join '.', "gnuplot", $pid;
#	open GNUPLOT, ">$gnuplotfile" or die "Error: Can't open $gnuplotfile\n";
#	print GNUPLOT "set term png\n";
#	print GNUPLOT "set output \"$pngfile\"\n";
#	print GNUPLOT "plot \"$datfile\" with points pointtype 7 pointsize 0.7, ";
#	print GNUPLOT "\"$fitfile\" with lines\n";
#	close GNUPLOT;
#
#	$gnuline = `$gnuplot < $gnuplotfile`;
#	print "$gnuline";
#	STDOUT->autoflush(1);
#}

#
# Consolidate the generated .csv file with our input file into the final output file
#
print "Consolidating output files and generating plots\n";
STDOUT->autoflush(1);

$isooutfile = join '.', $inputfile, "csv", $pid;
@input_array2 = readfile($isooutfile);
@header_array2 = split /\,/, shift(@input_array2);
for($ctr = 0; $ctr <= $#header_array2; $ctr++)
{
	# Need to add in this extra step to trim spurious whitespace from the isodist output fields
	@column_name_array = split ' ', $header_array2[$ctr];
	$column_name = $column_name_array[0];

	#$name2col2{ $header_array2[$ctr] } = $ctr;
	#$col2name2{ $ctr } = $header_array2[$ctr];
	$name2col2{ $column_name } = $ctr;
	$col2name2{ $ctr } = $column_name;
}

#@iso_keep = ( "mw", "tim", "chisq", "b", "b_err", "off", "off_err", "gw", "gw_err", "ul_amp", "ul_amp_err", "l_amp", "l_amp_err", "frac_n", "frac_n_err", "frac_lab", "frac_lab_err", "symb", "mz" );
@iso_keep = ( "file", "protein", "pep", "mf", "mw", "z_charge", "tim", "chisq", "b", "b_err", "ul_amp", "ul_amp_err", "l_amp", "l_amp_err", "off", "off_err", "gw", "gw_err", "frac_n", "frac_n_err", "frac_lab", "frac_lab_err", "symb", "mz" );
@iso_keep_print = ( "isofile", "isoprotein", "isopep", "isomf", "mw", "isoz_charge", "tim", "chisq", "b", "b_err", "ul_amp", "ul_amp_err", "l_amp", "l_amp_err", "off", "off_err", "gw", "gw_err", "frac_n", "frac_n_err", "frac_lab", "frac_lab_err", "symb", "mz" );
# file protein pep mf z_charge -> for now keep entire contents of old isodist output
$iso_keep_line = join ' ', @iso_keep;
$iso_keep_line_csv = join ',', @iso_keep;
$iso_keep_line_csv_print = join ',', @iso_keep_print;

open OUT, ">$outputfile" or die "Error: Can't open $outputfile\n";
#print OUT "$header_line $iso_keep_line\n";
#print OUT "$header_line_csv,$iso_keep_line_csv,comment\n";
print OUT "$iso_keep_line_csv_print,$header_line_csv,comment\n";
ISOSEARCH: for($ctr = 0; $ctr <= $#input_array; $ctr++)
{
	#
	# For each line from the .isoin file, we need to search and see if there actually is an output from isodist
	# isodist will skip input lines for various reasons (ie bogus peptides)
	#
	#print OUT "$input_array[$ctr] ";
	@array = split ' ', $input_array[$ctr]; # This is the .isoin
	$outline_csv = join ',', @array;
	#print OUT "$outline_csv,";
	$file = $array[ $name2col{"file"} ];
	substr($file, -4, 4) = ""; # Trim off the .txt extension
	# replace the $dir with the new $fitsdir
	@file_array = split /\//, $file;
	$file_array[0] = $fitsdir;
	$file = join "/", @file_array;

	for($ctr2 = 0; $ctr2 <= $#input_array2; $ctr2++)
	{
		@array2 = split /\,/, $input_array2[$ctr2];

		$isofile = $array2[ $name2col2{"file"} ];

		#if($array2[ $name2col2{"file"} ] eq $file)
		if($isofile eq $file)
		{
			#print OUT "$outline_csv,";
			for($ctr3 = 0; $ctr3 <= $#iso_keep; $ctr3++)
			{
				#print OUT "$array2[ $name2col2{ $iso_keep[$ctr3] } ] ";
				print OUT "$array2[ $name2col2{ $iso_keep[$ctr3] } ],";
			}
			print OUT "$outline_csv"; # For now tack this all onto end
			print OUT "\n";

			#
			# Now do the gnuplot stuff here
			#
			$fitfile = join '', $isofile, ".fit";
			$datfile = join '', $isofile, ".dat";
			$pngfile = join '', $isofile, ".fit.png";

			$gnuplotfile = join '.', "gnuplot", $pid;
			open GNUPLOT, ">$gnuplotfile" or die "Error: Can't open $gnuplotfile\n";

				#print GNUPLOT "set term png picsize 800 600\n";
				# This notation does not work with /sw/bin/gnuplot, use below
				print GNUPLOT "set term png medium size 800,600\n";
				print GNUPLOT "set xtics out\n";

			print GNUPLOT "set output \"$pngfile\"\n";

#				print GNUPLOT "plot \"$datfile\" with points pointtype 4 pointsize 1.5, ";
				print GNUPLOT "plot \"$datfile\" with points pointtype 2 pointsize 1.0, ";

			print GNUPLOT "\"$fitfile\" with lines\n";
			close GNUPLOT;

			$gnuline = `$gnuplot < $gnuplotfile`;
			print "$gnuline";
			STDOUT->autoflush(1);

			# Remove .fit and .dat to save space - not needed for now
			`/bin/rm -f $datfile`;
			`/bin/rm -f $fitfile`;

			next ISOSEARCH;
		}
	}
}
close OUT;

#
# Clean up some of the files both locally and on bluefish
#
print "Cleaning up files\n";
STDOUT->autoflush(1);

$blueline = `/bin/rm -f $pid/`;
print "$blueline";
$blueline = `cat doblu_m5.$pid.e$jobid`;
print "Start $machine eout\n";
print "$blueline";
print "End $machine eout\n";
STDOUT->autoflush(1);
$blueline = `/bin/rm -rf doblu_m5.$pid.e$jobid`;
$blueline = `/bin/rm -rf doblu_m5.$pid.o$jobid`;
print "$blueline";
STDOUT->autoflush(1);

$localline = `/bin/rm -f *.$pid`;
print "$localline";
STDOUT->autoflush(1);
