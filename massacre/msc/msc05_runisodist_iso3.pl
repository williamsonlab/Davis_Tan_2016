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

$baseline = 1.0;
$offset = 0.01;
$gw = 0.0003;
$niter = 5;
$keepfit = 0;
$baseline_fit_type = 0;
$ailtype = 0;

$isodist_binary = "isodist3";

&GetOptions("input=s" => \$inputfile, "output=s" => \$outputfile, "baseline=f" => \$baseline, "offset=f" => \$offset, "gw=f" => \$gw, "dir=s" => \$dir, "fitsdir=s" => \$fitsdir, "blueid=s" => \$blueid, "model=s" => \$input_model, "atoms=s" => \$input_atoms, "isodist_binary=s" => \$isodist_binary, "machine=s" => \$machine, "keepfit=i" => \$keepfit, "baseline_fit_type=i" => \$baseline_fit_type, "ailtype=i" => \$ailtype);

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

$num_in = $#input_array + 1; # Number of input lines from .isoin

####################
# isodist_v3 notes #
####################
# In addition to the batch file, which is unchanged, there are now 3 files:
# 1) isodist.in - similar to old version
# 2) atom definitions file
# 3) residue definitions file
#
# We can use the same naming scheme (isodist_batch.in.$pid) for #1
# For #2, we should really be able to use one global atom definition file
# For #3, we will have several pregenerated files available from a dropdown menu
# There will also be the option to specify a new file so users can add new models without altering the code of massacre in any way


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

# Template for the new v3 config file
#####
#fitit       = program options:  fitit tryit
#15N_pulse_short.batch  = batchfile:  file containing peptides, chgs, peaks
#/Users/jrwill/prog/isodist_8_10/atom_definitions.txt  = atomfile
#/Users/jrwill/prog/isodist_8_10/res_definitions.txt = resfile
#5         = niter   # of iteractions for each round of least squares
#100.0 = sigma   std deviation of noise (currently read, but not used)
#1.0   fixed       = B  :  fixed value for baseline, held constant  (Default = 1.0)
#0.01  variable    = OFF:  initial guess for accuracy offset  (Default = 0.01)
#0.003 variable    = GW :  initial guess for gaussian width (Use 0.0003 to 0.003)
#####

#$tempinfile = join '.', "isodist3.in", $pid;
$tempinfile = "isodist3.in";
#$tempatomsfile = join '.', "isodist3.atoms", $pid;
$tempatomsfile = "isodist3.atoms";
#$tempmodelfile = join '.', "isodist3.model", $pid;
$tempmodelfile = "isodist3.model";
print "Generating isodist config file: $tempinfile\n";
STDOUT->autoflush(1);
open TEMPIN, ">$tempinfile" or die "Error: Can't open $tempinfile\n";

if("isodist" eq "v2")
{
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
}

print TEMPIN "fitit\n";
print TEMPIN "$tempbatchfile\n";
print TEMPIN "isodist3.atoms\n";
print TEMPIN "isodist3.model\n";
print TEMPIN "$niter = niter\n";
print TEMPIN "100.0 = sigma\n";
if($baseline_fit_type == 0)
{
	printf(TEMPIN "%.6f auto = B\n", $baseline);
}
elsif($baseline_fit_type == 1)
{
	printf(TEMPIN "%.6f fixed = B\n", $baseline);
}
printf(TEMPIN "%.6f variable = OFF\n", $offset);
printf(TEMPIN "%.6f variable = GW\n", $gw);

close TEMPIN;

# Copy atoms and model files to the local directory with proper name
`cp $input_atoms $tempatomsfile`;
`cp $input_model $tempmodelfile`;

$npeptide = $#input_array + 1;
if($ailtype == 5)
{
	$min_per_pep = 1.0; # Orbitrap/Sequest data is taking a long time 2/24/2011
}
else
{
	$min_per_pep = 0.5;
}
$time = 15 + ($npeptide * $min_per_pep);
$hours = sprintf("%d", $time/60);
$mins = sprintf("%d", $time - 60*$hours);
if($hours == 0 && $mins == 0)
{
	$mins = 5;
}

$qsubfile = join '.', "doblu_m5", $pid;
print "Generating $machine qsub file: $qsubfile\n";
print "Requesting walltime of $hours hours and $mins minutes\n";
STDOUT->autoflush(1);

open QSUB, ">$qsubfile" or die "Error: Can't open $qsubfile\n";
print QSUB "#PBS -l walltime=$hours:$mins:00\n";
print QSUB "#PBS -l cput=$hours:$mins:00\n\n";

print QSUB "cd \$PBS_O_WORKDIR\n\n";

#print QSUB "mkdir $fitsdir\n";
print QSUB "$isodist_binary isodist3.in\n\n";
print QSUB "/bin/mv -f $tempbatchfile.csv $inputfile.csv.$pid\n";
close QSUB;

#print "Copying files to $machine\n";
#STDOUT->autoflush(1);

#$blueline = `mkdir $pid/`;
#print "$blueline";
#STDOUT->autoflush(1);

print "Running job on $machine\n";
STDOUT->autoflush(1);

$blueline = `qsub $qsubfile`;
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
	$header_array2[$ctr] = $column_name_array[0];

	#$name2col2{ $header_array2[$ctr] } = $ctr;
	#$col2name2{ $ctr } = $header_array2[$ctr];
	$name2col2{ $column_name } = $ctr;
	$col2name2{ $ctr } = $column_name;
}

$num_fit = $#input_array2 + 1; # Number of lines fit in the .batch.csv output from isodist

#@iso_keep = ( "mw", "tim", "chisq", "b", "b_err", "off", "off_err", "gw", "gw_err", "ul_amp", "ul_amp_err", "l_amp", "l_amp_err", "frac_n", "frac_n_err", "frac_lab", "frac_lab_err", "symb", "mz" );
@iso_keep = ( "file", "protein", "pep", "mf", "mw", "z_charge", "tim", "chisq", "b", "b_err", "ul_amp", "ul_amp_err", "l_amp", "l_amp_err", "off", "off_err", "gw", "gw_err", "frac_n", "frac_n_err", "frac_lab", "frac_lab_err", "symb", "mz" );
@iso_keep_print = ( "isofile", "isoprotein", "isopep", "isomf", "mw", "isoz_charge", "tim", "chisq", "b", "b_err", "ul_amp", "ul_amp_err", "l_amp", "l_amp_err", "off", "off_err", "gw", "gw_err", "frac_n", "frac_n_err", "frac_lab", "frac_lab_err", "symb", "mz" );
# file protein pep mf z_charge -> for now keep entire contents of old isodist output

# Now since isodist output columns are model-dependent, we just keep them all
@iso_keep = @header_array2;
# Change a few for naming conflicts
@iso_keep_print = @iso_keep;
$iso_keep_print[0] = "isofile";
$iso_keep_print[1] = "isoprotein";
$iso_keep_print[2] = "isopep";
$iso_keep_print[4] = "isoz_charge";

$iso_keep_line = join ' ', @iso_keep;
$iso_keep_line_csv = join ',', @iso_keep;
$iso_keep_line_csv_print = join ',', @iso_keep_print;

open OUT, ">$outputfile" or die "Error: Can't open $outputfile\n";
#print OUT "$header_line $iso_keep_line\n";
#print OUT "$header_line_csv,$iso_keep_line_csv,comment\n";
print OUT "$iso_keep_line_csv_print,$header_line_csv";

if($isodist_binary eq "isodist3_rsd_ks")
{
	print OUT ",ks_dstat,ks_pval";
}

print OUT ",comment\n";

$num_out = 0; # Number of fits matched to the .isoin

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
			$num_out++;

			#print OUT "$outline_csv,";
			for($ctr3 = 0; $ctr3 <= $#iso_keep; $ctr3++)
			{
				#print OUT "$array2[ $name2col2{ $iso_keep[$ctr3] } ] ";
				print OUT "$array2[ $name2col2{ $iso_keep[$ctr3] } ],";
			}
			print OUT "$outline_csv"; # For now tack this all onto end

			#
			# Now do the gnuplot stuff here
			# Need to set these files up before calc_ks
			#
			$fitfile = join '', $isofile, ".fit";
			$datfile = join '', $isofile, ".dat";
			$pngfile = join '', $isofile, ".fit.png";

			if($isodist_binary eq "isodist3_rsd_ks")
			{
				# Calculate Kolmogorov-Smirnov Test
				$dstat = calc_ks($datfile);
				print OUT ",$dstat";
			}

			print OUT "\n";

			$gnuplotfile = join '.', "gnuplot", $pid;
			open GNUPLOT, ">$gnuplotfile" or die "Error: Can't open $gnuplotfile\n";

			#print GNUPLOT "set term png picsize 800 600\n";
			# This notation does not work with /sw/bin/gnuplot, use below
			print GNUPLOT "set term png medium size 800,600\n";
			print GNUPLOT "set xtics out\n";

			print GNUPLOT "set output \"$pngfile\"\n";
			#print GNUPLOT "set y2tics border\n";

			print GNUPLOT "plot \"$datfile\" using 1:2 with points pointtype 2 pointsize 1.0, ";

			print GNUPLOT "\"$fitfile\" with lines, ";
			#print GNUPLOT "\"$datfile\" using 1:3 axes x1y2 with lines\n";
			print GNUPLOT "\"$datfile\" using 1:3 with lines\n";
			close GNUPLOT;

			$gnuline = `$gnuplot < $gnuplotfile`;
			print "$gnuline";
			STDOUT->autoflush(1);

			# Remove .fit and .dat to save space - not needed for now
			if($keepfit == 0)
			{
				`/bin/rm -f $fitfile`;
				`/bin/rm -f $datfile`;
			}

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

$blueline = `/bin/rm -rf $pid/`;
print "$blueline";
$blueline = `cat doblu_m5.$pid.e$jobid`;
print "Start $machine eout\n";
print "$blueline";
print "End $machine eout\n";
STDOUT->autoflush(1);
$blueline = `/bin/rm -r doblu_m5.$pid.e$jobid`;
$blueline = `/bin/rm -r doblu_m5.$pid.o$jobid`;
print "$blueline";
STDOUT->autoflush(1);

$localline = `/bin/rm -f *.$pid`;
print "$localline";
STDOUT->autoflush(1);

if($num_in > 0)
{
	$percent_out = 100 * ($num_out/$num_in);
}
else
{
	$percent_out = 0;
}
$percent_out = sprintf("%.1f", $percent_out);
print "$num_in entries in the .isoin\n";
print "$num_fit entries were fit\n";
print "$num_out fit entries matched the .isoin ($percent_out%)\n";

sub calc_ks
{
	my @params = @_;
	my $ksdat = shift(@params);
	open KSDAT, "$ksdat" or die "Can't open $ksdat\n";
	my $line;
	my @dat = ();
	my @fit = ();

	while(defined($line=<KSDAT>))
	{
		chomp($line);
		my @array = split /\,/, $line;
		my $mz = shift(@array);
		my $datval = shift(@array);
		push(@dat, $datval);
		my $rsd = shift(@array);
		my $fitval = shift(@array);
		push(@fit, $fitval);
	}
	close KSDAT;

	open KSTEMPDAT, ">$pid.kstempdat" or die "Can't open $pid.kstempdat\n";
	for(my $ctr = 0; $ctr <= $#dat; $ctr++)
	{
		print KSTEMPDAT "$dat[$ctr]\n";
	}
	close KSTEMPDAT;

	open KSTEMPFIT, ">$pid.kstempfit" or die "Can't open $pid.kstempfit\n";
	for(my $ctr = 0; $ctr <= $#fit; $ctr++)
	{
		print KSTEMPFIT "$fit[$ctr]\n";
	}
	close KSTEMPFIT;

	open KSTEMPR, ">$pid.kstempR" or die "Can't open $pid.kstempR\n";
	print KSTEMPR "x_fit <- scan(\"$pid.kstempfit\")\n";
	print KSTEMPR "x_dat <- scan(\"$pid.kstempdat\")\n";
	print KSTEMPR "ks.test(x_fit,x_dat)\n";

	`R CMD BATCH $pid.kstempR`;

	open KSTEMPROUT, "$pid.kstempR.Rout" or die "Can't open $pid.kstempR.Rout\n";
	my @rarray = <KSTEMPROUT>;
	close KSTEMPROUT;

	my $dstat;
	my $pval;

	for(my $ctr = 0; $ctr <= $#rarray; $ctr++)
	{
		chomp($rarray[$ctr]);
		my @array = split ' ', $rarray[$ctr];

		if($array[0] eq "D")
		{
			$dstat = $array[2];
			$pval = $array[$#array];
		}
	}

	my $returnval = join '', $dstat, $pval;

	`/bin/rm -f $pid.kstempR`;
	`/bin/rm -f $pid.kstempR.Rout`;
	`/bin/rm -f $pid.kstempfit`;
	`/bin/rm -f $pid.kstempdat`;

	return $returnval;
}

sub old_calc_ks
{
	my @params = @_;
	my $ksdat = shift(@params);
	open KSDAT, "$ksdat" or die "Can't open $ksdat\n";
	my $line;
	my @dat = ();
	my $datsum = 0;
	my @fit = ();
	my $fitsum = 0;

	while(defined($line=<KSDAT>))
	{
		chomp($line);
		my @array = split /\,/, $line;
		my $mz = shift(@array);
		my $datval = shift(@array);
		push(@dat, $datval);
		$datsum += $datval;
		my $rsd = shift(@array);
		my $fitval = shift(@array);
		push(@fit, $fitval);
		$fitsum += $fitval;
	}
	close KSDAT;

	my $dstat = 0;
	my $nval = 0;
	my $datrunsum = 0;
	my $fitrunsum = 0;

	$nval = $#dat + 1;

	for(my $ctr = 0; $ctr <= $#dat; $ctr++)
	{
		$datrunsum += $dat[$ctr];
		$fitrunsum += $fit[$ctr];

		my $datfrac = $datrunsum / $datsum;
		my $fitfrac = $fitrunsum / $fitsum;

		my $tempstat = abs($datfrac - $fitfrac);
		if($tempstat > $dstat)
		{
			$dstat = $tempstat;
		}
	}

	my $returnval = join ',', $dstat, $nval;

	return $returnval;
}
