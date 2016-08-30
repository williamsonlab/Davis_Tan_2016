#!/usr/bin/perl

# MS MS scan RT values are not always at the center of the elution profile
# Fix this

use Getopt::Long;
use IO::Handle;
use Text::ParseWords;

$pid = $$;
$pwd = `pwd`;
chomp($pwd);
$hostname = `hostname`;
chomp($hostname);
$gnuplot = "gnuplot";
$qsubfile = join '.', "doblu_m0", $pid;
$rofile = join '', $pwd, "/rofile.", $pid;

$has_rt = -1; # Flag to determine if RT is embedded, or scan number is used. If scan number is used, it still needs to be output in place of RT

require "$ENV{'MASSACRE_PATH'}/msc/mikesubs.pl";
require "$ENV{'MASSACRE_PATH'}/msc/mikevars.pl";

$mascotfile = "";
$outputfile = "";
$blueid = "";
$csvfile = "";
$dir = "";
$machine = "";
$msms_proteinkeep = 2;
$ailtype = "";

&GetOptions("ail=s" => \$mascotfile, "output=s" => \$outputfile, "blueid=s" => \$blueid, "csvfile=s" => \$csvfile, "dir=s" => \$dir, "machine=s" => \$machine, "msms_proteinkeep=i" => \$msms_proteinkeep, "ailtype=i" => \$ailtype);

$tempfile = join '.', $mascotfile, "rttemp", $pid;

if($mascotfile eq "" || $outputfile eq "" || $blueid eq "" || $csvfile eq "" || $dir eq "" || $machine eq "" || $ailtype eq "")
{
	print "ERROR: Must define --ail=something --output=something --blueid=something --csvfile=something --dir=something --machine=something --ailtype=something\n";
	STDOUT->autoflush(1);
	die;
}

if($msms_proteinkeep == 0)
{
	print "Keeping only ribosomal proteins\n";
	STDOUT->autoflush(1);
}
elsif($msms_proteinkeep == 1)
{
	print "Keeping only previously defined proteins\n";
	STDOUT->autoflush(1);
}
elsif($msms_proteinkeep == 2)
{
	print "Keeping all proteins\n";
	STDOUT->autoflush(1);
}

`mkdir $dir`;

if($ailtype == 4)
{
	parse_mascot();
}
elsif($ailtype == 5)
{
	parse_sequest();
}
else
{
	print "ailtype=|$ailtype| unknown\n";
	STDOUT->autoflush(1);
	die;
}

do_qsub_file();

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

WAIT: while (1) {
    while (defined ($roline = <ROFILE>)) {
        print "$roline" and STDOUT->autoflush(1);
    }

    my $qstat_out = `qstat $jobid 2>&1`;

    if ($qstat_out =~ /Unknown Job Id/) {
        print "Remote job $jobid finished.\n" and last WAIT;
    } else { sleep 180; }
}

# Print remote errors
$blueline = `cat doblu_m0.$pid.e$jobid`;
print "Start $machine eout\n";
print "$blueline";
print "End $machine eout\n";
STDOUT->autoflush(1);

# Remove some remote files
$blueline = `/bin/rm -f doblu_m0.$pid.e$jobid`;
print "$blueline";
$blueline = `/bin/rm -f doblu_m0.$pid.o$jobid`;
print "$blueline";
$blueline = `/bin/rm -f doblu_m0.$pid`;
print "$blueline";
STDOUT->autoflush(1);

# Remove some local files
`rm $rofile`;
`rm $qsubfile`;

open OUTPUT, ">$outputfile" or die "Can't open $outputfile\n";
if($ailtype == 4)
{
	print OUTPUT "$mascotheader,rt,corr_rt,rtcorrfile,rt_max_int\n";
}
elsif($ailtype == 5)
{
	print OUTPUT "$sequestheader\trt\tcorr_rt\trtcorrfile\trt_max_int\n";
}

print "Generating RT correction plots\n";
STDOUT->autoflush(1);
###########################
# Do Plotting and Fitting #
###########################
# The order of the mascot file is not preserved by msc00_adjustmsmsrt_blu.pl
# Basically, things are output in RT order
# Unfortunately, this destroys the integrity of the file because the information in the mascot output is sparse, and subsequent lines depend upon previous lines for things like the protein accession number and description
# So we need to resort $tempfile.2 based on the initial order
# Luckily, initial line numbers are embedded as the first element of the $key
open TMP2, "$tempfile.2" or die "Can't open $tempfile.2\n";
@tmp2_array = <TMP2>;
close TMP2;
@tmp2_array = sort by_key_number @tmp2_array;

#open TMP2, "$tempfile.2" or die "Can't open $tempfile.2\n";
#while(defined($input=<TMP2>))
while(defined($input=shift(@tmp2_array)))
{
	chomp($input);
	@array = split ' ', $input;
	$key = shift(@array);
	$rt = shift(@array);
	$max_int = shift(@array);

	open GNU, ">m0.gnu.$pid" or die "Can't open m0.gnu.$pid\n";
	print GNU "f(x) = A*exp(-B*(x-C)**2)\n";
	print GNU "A = $max_int\n";
	print GNU "B = 50\n"; # Determined from a test - has to do with width
	print GNU "C = $rt\n";
	print GNU "set fit logfile \"m0.gnu.log.$pid\"\n";
	print GNU "fit f(x) \"$dir/$key.rt\" via A,B,C\n";
	close GNU;

	`$gnuplot < m0.gnu.$pid >& m0.gnu.tmp`;

	open GNUTMP, "m0.gnu.log.$pid" or die "Can't open m0.gnu.log.$pid\n";
	@gnuarray = <GNUTMP>;
	close GNUTMP;
	`rm m0.gnu.log.$pid`;

	for($ctr = 0; $ctr <= $#gnuarray; $ctr++)
	{
		@array = split ' ', $gnuarray[$ctr];
		if($array[0] eq "A" && $array[1] eq "=")
		{
			$a = $array[2];
		}
		if($array[0] eq "B" && $array[1] eq "=")
		{
			$b = $array[2];
		}
		if($array[0] eq "C" && $array[1] eq "=")
		{
			$c = $array[2];
		}
	}

	$new_rt = $c;

	open GNU, ">m0.gnu.$pid" or die "Can't open m0.gnu.$pid\n";
	print GNU "f(x) = A*exp(-B*(x-C)**2)\n";
	print GNU "A = $a\n";
	print GNU "B = $b\n";
	print GNU "C = $c\n";
	print GNU "set arrow from $rt, graph 0 to $rt, graph 1 nohead ls 1\n";
	print GNU "set arrow from $new_rt, graph 0 to $new_rt, graph 1 nohead ls 2\n";
	print GNU "set label \"$rt\" at graph 0.975, 0.06 right tc rgbcolor \"#FF0000\"\n";
	print GNU "set label \"$new_rt\" at graph 0.975, 0.03 right tc rgbcolor \"#00FF00\"\n";
	print GNU "set term png\n";
	print GNU "set output \"$dir/$key.png\"\n";
	print GNU "plot \"$dir/$key.rt\" with lines, f(x)\n";
	close GNU;

	`$gnuplot < m0.gnu.$pid >& m0.gnu.tmp`;

	`rm m0.gnu.$pid`;
	`rm m0.gnu.tmp`;

	if($ailtype == 4)
	{
		print OUTPUT "$mascothash{$key},$rt,$new_rt,$dir/$key.png,$max_int\n";
	}
	elsif($ailtype == 5)
	{
		print OUTPUT "$sequesthash{$key}\t$rt\t$new_rt\t$dir/$key.png,$max_int\n";
	}
}

close OUTPUT;

# Remove some of the temp files
`rm $tempfile`;
`rm $tempfile.2`;


################
################
################
################
################


################
# do_squb_file #
################
sub do_qsub_file
{
	$num_spectra = $count;
	$walltime = 60 + 0.35 * $num_spectra;
	$hours = sprintf("%d", $walltime/60);
	$mins = sprintf("%d", $walltime-60*$hours);

	open QSUB, ">$qsubfile" or die "Error: Can't open $qsubfile\n";
	print QSUB "#PBS -l walltime=$hours:$mins:00\n";
	print QSUB "#PBS -l cput=$hours:$mins:00\n";

	print QSUB "cd \$PBS_O_WORKDIR\n\n";

	print QSUB "msc00_adjustmsmsrt_blu.pl --input=$tempfile --csv=$csvfile --dir=$dir --blueid=$blueid --hostname=$hostname --rofile=$rofile --ailtype=$ailtype  --has_rt=$has_rt\n";

	print QSUB "/bin/rm -f `basename $csvfile .mzML.gz`.mzML \n";
}

##################
# parse_mascot() #
##################
# Generate @mascotdata
#	$mascotdata[n][0] = key
#	$mascotdata[n][1] = charge
#	$mascotdata[n][2] = m/z
#	$mascotdata[n][3] = RT
#	$mascotdata[n][4] = original file line
#
# WTF @mascotdata does not seem to be used
# %mascothash stores all the values from the file for output
# Also need to generate $tempfile which is sent to msc00_adjustmsmsrt_blu.pl
#
sub parse_mascot
{
	my @filearray = readfile($mascotfile);

	open TEMP, ">$tempfile" or die "Can't open $tempfile\n";

	while(1)
	{
		my $line = shift(@filearray);
		my @array = quotewords(",", 0, $line);

		if($#filearray == 0 || $#filearray == -1)
		{
			print "ERROR: Does not appear to be a mascot export file\n";
			STDOUT->autoflush(1);
			die;
		}

		%mascotname2col = ();
		for(my $ctr = 0; $ctr <= $#array; $ctr++)
		{
			$mascotname2col{$array[$ctr]} = $ctr;
		}

		# Detect start of actual data by header line
		# This is important as there are tons of junk lines preceding the data
		if(defined($mascotname2col{"pep_exp_mz"}) && defined($mascotname2col{"pep_scan_title"}))
		{
			$mascotheader = $line;
			last;
		}
	}

	$count = 0;

	my $mpro;

	# Organize all entries into a hash (%mascotraw) sorted by ion (protein_peptide_charge)
	$were_noacc2pro_peptides = 0;
	for(my $ctr = 0; $ctr <= $#filearray; $ctr++)
	{
		my $mseq;
		my $mcrg;
		my $mrt;
		my $mmz;

		my @array = quotewords(",", 0, $filearray[$ctr]);

		if(substr($array[0], 0, 15) eq "Peptide matches")
		{
			last;
		}

		# If the prot_acc field is blank, it means a continuation of the previous protein
		# So only define a new protein if the field contains a value
		if($array[$mascotname2col{"prot_acc"}] ne "")
		{
			$mpro = $array[$mascotname2col{"prot_acc"}];

			# Replace underscores by dashes since underscores are used as delimiters elsewhere (that's probably not the best choice)
			$mpro =~ s/_/-/g;

			$mdesc = $array[$mascotname2col{"prot_desc"}];
		}

		$mseq = $array[$mascotname2col{"pep_seq"}];
		$mcrg = $array[$mascotname2col{"pep_exp_z"}];
		$mmz = $array[$mascotname2col{"pep_exp_mz"}];

		# pep_scan_title
        #
        # MFG Files generated by Qualitative Analysis result in mascot file that have an entry similar to this in the pep_scan_title
		# 1_30S_1to10dil_060409_cap_MSMS.d, MS/MS of 439.7410431 2+ at 33.7646166666667 mins
		# So RT is the penultimate array element of @array2, and is easily extracted in minutes
		my @array2 = split ' ', $array[$mascotname2col{"pep_scan_title"}];

		if($array2[$#array2] eq "mins")
		{
			$mrt = $array2[$#array2-1];
			$has_rt = 1;
		}
        #
        # QTOF .d files converted to .mgf using pwiz (or /jrwill_lab/bin/d2mgf.pl) generate files that look like this
        #BEGIN IONS
        #TITLE=scanId=190470
        #RTINSECONDS=190.463
        #PEPMASS=402.17558288574219 3177.2652859687805
        #CHARGE=2
        #82.42181236 11.19999981
        #84.04106262 50.01922989
        #84.07439832 62.96982574
        #95.08340327 74.30882263
        #
        # So in this case the pep_scan_title in the mascot export file is scanId=190470, where this actually corresponds to the retention time in seconds
        # So just extract this value, and convert it to minutes
        #
		elsif($array[$mascotname2col{"pep_scan_title"}] =~ /scanid=(\d+)/i)
		{
			$mrt = $1;
			$mrt = $mrt/1000; # Convert scanId to retention (scanId ~= 1000*RT(sec))
			$mrt = $mrt/60; # Convert seconds to minutes
			$has_rt = 1;
		}
        # ABSciex Triple-TOF have a pep_scan_title that looks like this:
        # Locus:1:1:1:1806:4
        # "1806 is the cycle number, which needs to be converted to retention time as per sequest/LTQ data
        # NOTE: This must appear before the LTQ .RAW section below, as the regular expression in that section will also match this format, and return a scan of "1" every time
        elsif($array[$mascotname2col{"pep_scan_title"}] =~ /Locus\:(\d+)\.(\d+)\.(\d+)\.(\d+)\.(\d+)/)
        {
            $mrt = $4;
            $has_rt = 0;
        }
        #
        # LTQ .RAW files converted to mgf in the Siuzdak lab have a pep_scan_title that looks like this:
        # 7OS.2828.2828.2.dta
        # The matching numbers are scan numbers that need to be converted to retention time in a manner similar to the sequest data
        # Converting LTQ .RAW to mgf using pwiz results in a similar format
        # Sometimes, in the 70S_mascot.xls data from sschen, the numbers do not match exactly, differing by 2 (the second value is always 2 higher than the first), hence the second check for equality below.
        elsif($array[$mascotname2col{"pep_scan_title"}] =~ /.*?\.(\d+)\.(\d+)/)
        {
            $mrt_a = $1;
            $mrt_b = $2;
            if($mrt_a == $mrt_b || ($mrt_a+2) == $mrt_b)
            {
                $mrt = $mrt_a;
                $has_rt = 0
            }
            else
            {
                print "ERROR: Could not detect RT from $array[$mascotname2col{'pep_scan_title'}]\n";
                STDOUT->autoflush(1);
                die;
            }
        }
		else
		{
			print "ERROR: Could not detect RT from $array[$mascotname2col{'pep_scan_title'}]\n";
			STDOUT->autoflush(1);
			die;
		}

		my $mpro2;

		# $msms_proteinkeep = 0 -> check for accession number in %acc2pro (ribosomal only)
		# $msms_proteinkeep = 1 -> check for accession number in %acc2pro_all (all annot.)
		# $msms_proteinkeep = 2 -> keep all regardless of accession
		if($msms_proteinkeep == 0)
		{
			if(defined($acc2pro{$mpro}))
			{
				$mpro2 = $acc2pro{$mpro};
			}
			else
			{
				$were_noacc2pro_peptides = 1;
				$mdesc_short = substr($mdesc, 0, 25);
				print "Skipping $mpro: $mdesc_short\n";
				STDOUT->autoflush(1);
				next;
			}
		}
		elsif($msms_proteinkeep == 1)
		{
			if(defined($acc2pro_all{$mpro}))
			{
				$mpro2 = $acc2pro_all{$mpro};
			}
			else
			{
				$were_noacc2pro_peptides = 1;
				$mdesc_short = substr($mdesc, 0, 25);
				print "Skipping $mpro: $mdesc_short\n";
				STDOUT->autoflush(1);
				next;
			}
		}
		elsif($msms_proteinkeep == 2)
		{
			if(defined($acc2pro_all{$mpro}))
			{
				$mpro2 = $acc2pro_all{$mpro};
			}
			else
			{
				$mpro2 = $mpro;
			}
		}

		my $key = join '_', $count, $mpro2, $mseq, $mcrg, $mmz, $mrt;

		$mascotdata[$count][0] = $key;
		$mascotdata[$count][1] = $mcrg;
		$mascotdata[$count][2] = $mmz;
		$mascotdata[$count][3] = $mrt;
		$mascotdata[$count][4] = $filearray[$ctr];

		if(defined($mascothash{$key}))
		{
			die "Error, $key already defined\n";
			STDOUT->autoflush(1);
		}
		$mascothash{$key} = $filearray[$ctr];

		print TEMP "$key $mcrg $mmz $mrt\n";

		$count++;
	}
	if($were_noacc2pro_peptides)
	{
		print "Warning: If you didn't want to skip some of those, adjust the Configuration for Step 0 and/or talk to Mike\n";
	}

	close TEMP;
}

sub parse_sequest
{
	my @filearray = readfile($mascotfile);

	open TEMP, ">$tempfile" or die "Can't open $tempfile\n";

	while(1)
	{
		my $line = shift(@filearray);
		my @array = split /\t/, $line;

		if($#filearray == 0 || $#filearray == -1)
		{
			print "ERROR: Does not appear to be a sequest export file\n";
			STDOUT->autoflush(1);
			die;
		}

		%sequestname2col1 = ();
		for(my $ctr = 0; $ctr <= $#array; $ctr++)
		{
			$sequestname2col1{lc($array[$ctr])} = $ctr;
		}

		if(defined($sequestname2col1{"locus"}) && defined($sequestname2col1{"sequence count"}))
		{
			$sequestheader1 = $line;
			$sequestheader2 = shift(@filearray);
			$sequestheader = join "\t", $sequestheader1, $sequestheader2;
			last;
		}
	}

	my @array = split /\t/, $sequestheader2;
	for(my $ctr = 0; $ctr <= $#array; $ctr++)
	{
		$sequestname2col2{lc($array[$ctr])} = $ctr;
	}

	$count = 0;

	my $mpro;

	my $num_seq_skip = 0;
	my $num_seq_keep = 0;

	$were_noacc2pro_peptides = 0;
	for(my $ctr = 0; $ctr <= $#filearray; $ctr++)
	{
		my @array = split /\t/, $filearray[$ctr];

		# Check for summary info indicating end of entries
        # This may occur at the end of a sequest data file
		if($array[0] eq "" && $array[1] eq "Proteins")
		{
			last;
		}
        # An asterix (*) or a blank entry indicate a continuation of a protein
        # Anything else is the protein locus, including accession number
		elsif($array[0] ne "*" && $array[0] ne "")
		{
			$sequestproteinline = $filearray[$ctr];
			my @array2 = split /\|/, $array[$sequestname2col1{"locus"}];

			# Looking for format crap|accession|crap
			# Sometimes get crap|accession
			# Check in case you get "stuff" (no |'s)
			if($#array2 >= 1)
			{
				$mpro = $array2[1];
			}
			else
			{
				$mpro = $array[$sequestname2col1{"locus"}];
			}

			if(substr($mpro, 0, 7) eq "REFSEQ:")
			{
				substr($mpro, 0, 7) = "";
			}

			# Replace underscores by dashes since underscores are used as delimiters elsewhere (that's probably not the best choice)
			$mpro =~ s/_/-/g;

			$mdesc = $array2[2];
			next;
		}
        # Here either an asterix or blank as detected, so this is assumed to be a continuation of the entry for a previously listed protein
		else
		{
            # @array2 splits the sequence entry, which looks like this: R.VTVQDAVEK.I
            # So $array2[1] is the actual sequence
			my @array2 = split /\./, $array[$sequestname2col2{"sequence"}];
			$mseq = $array2[1];

            # Sometimes the sequence looks like this: K.GAEAALTALEM(15.9949)INVLK.A
            # This implies a modification
            # These are skipped, but could be adjusted to be handled
            # The regular expression looks for any non A-Z characters
            # NOTE: Splitting on the periods actually makes $array2[1] be "GAEAALTALEM(15" so a better way of splitting is needed if modifications are going to be handled properly
			if($mseq =~ /\d/)
			{
				print "Skipping $filearray[$ctr] due to sequence |$mseq|\n";
				STDOUT->autoflush(1);
				$num_seq_skip++;
				next;
			}

			$num_seq_keep++;

            # @array3 splits the filename entry: long-04.6701.6701.1
            # It turns out the charge is embedded in the file name as the last entry
            # The [1] and [2] elements are the scan ID, and should be the same
            # If these are the same, they are saved as RT
            # This is converted to a proper RT by msc00_adjustmsmsrt_blu.pl (see the do_mzml_scan2rt subroutine for starters)
			my @array3 = split /\./, $array[$sequestname2col2{"filename"}];
			$mcrg = $array3[3];
			#$mmz = $array[$sequestname2col2{"m+h+"}] / $mcrg;
			$mmz = ( $array[$sequestname2col2{"calcm+h+"}] + ($mcrg-1) ) / $mcrg;
			if($array3[1] == $array3[2])
			{
				$mrt = $array3[1];
			}
			else
			{
				die "RT Detection Error: |$array3[1]| != |$array3[2]|\n";
			}
		}

		my $mpro2;

		# $msms_proteinkeep = 0 -> check for accession number in %acc2pro (ribosomal only)
		# $msms_proteinkeep = 1 -> check for accession number in %acc2pro_all (all annot.)
		# $msms_proteinkeep = 2 -> keep all regardless of accession
		if($msms_proteinkeep == 0)
		{
			if(defined($acc2pro{$mpro}))
			{
				$mpro2 = $acc2pro{$mpro};
			}
			else
			{
				$were_noacc2pro_peptides = 1;
				$mdesc_short = substr($mdesc, 0, 25);
				print "Skipping $mpro: $mdesc_short\n";
				STDOUT->autoflush(1);
				next;
			}
		}
		elsif($msms_proteinkeep == 1)
		{
			if(defined($acc2pro_all{$mpro}))
			{
				$mpro2 = $acc2pro_all{$mpro};
			}
			else
			{
				$were_noacc2pro_peptides = 1;
				$mdesc_short = substr($mdesc, 0, 25);
				print "Skipping $mpro: $mdesc_short\n";
				STDOUT->autoflush(1);
				next;
			}
		}
		elsif($msms_proteinkeep == 2)
		{
			if(defined($acc2pro_all{$mpro}))
			{
				$mpro2 = $acc2pro_all{$mpro};
			}
			else
			{
				$mpro2 = $mpro;
			}
		}

		my $key = join '_', $count, $mpro2, $mseq, $mcrg, $mmz, $mrt;

		if(defined($sequesthash{$key}))
		{
			die "Error, $key already defined\n";
			STDOUT->autoflush(1);
		}
		$sequesthash{$key} = join "\t", $sequestproteinline, $filearray[$ctr];

		print TEMP "$key $mcrg $mmz $mrt\n";

		$count++;
	}
	if($were_noacc2pro_peptides)
	{
		print "Warning: If you didn't want to skip some of those, adjust the Configuration for Step 0 and/or talk to Mike\n";
	}
	if($num_seq_skip)
	{
		my $num_seq_total = $num_seq_skip + $num_seq_keep;
		print "Skipped $num_seq_skip/$num_seq_total peptide entries due to weird sequence\n";
		STDOUT->autoflush(1);
	}

	close TEMP;
}

sub by_key_number
{
    @aa = split /\_/, $a;
    @bb = split /\_/, $b;

    $aa[0] <=> $bb[0];
}
