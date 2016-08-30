#!/usr/bin/perl

use Getopt::Long;
use IO::Handle;
use Text::ParseWords;
require "$ENV{'MASSACRE_PATH'}/msc/mikesubs.pl";

#
# "ail" line
# ID      multno  Z       Ion     RT      M/Z             Abundance
# 8       1       3       M+3H    22.139  371.8736        2.32
# [0]     [1]     [2]     [3]     [4]     [5]             [6]
#

#
# theoretical digest line
# digestionid digestpepid n14mass n15mass protein startres endres iontype charge missed seq mod n14prox n15prox
# 1 1 300.4888117489 304.8092962861 S2 11 25 M+6H 6 1 AGVHFGHQTRYWNPK --------------- 0.015398 0.021013
#

#############
# Variables #
#############
$ppm_offset = 0; # Machine-based offset for accuracy
$ppm_thresh = 50; # Final Threshold to keep or reject a match at this stage
$sample_id = "0"; # Sample ID -> goes in as "0min" or something like that

# mode = 0 -> Both N14 and N15 matches
# mode = 1 -> N14 Only
# mode = 2 -> N15 Only
$mode = 0;

$digestfile = "digestfile.txt";
$ailfile = "ailfile.xls";
$outputfile = "output01.txt";

$ailtype = 0; # 0 for original, 1 for MassHunter2 (Qualitative Analysis) Type

# This is the line where options are read from the command line
# This line can overwrite all the default values above
&GetOptions("ppm_offset=i" => \$ppm_offset, "ppm_thresh=i" => \$ppm_thresh, "id=s" => \$sample_id, "mode=i" => \$mode, "digest=s" => \$digestfile, "ail=s" => \$ailfile, "output=s" => \$outputfile, "ailtype=i" => \$ailtype);

#
# Read in the theoretical digest file
# Store the index of each column in %pep_header
# This allows us to operate with arbitrary array indices, allowing for easy insertion of additional columns into the data
#
print "Reading Theoretical Digest File\n";
STDOUT->autoflush(1);

@pep_array = readfile($digestfile);

print "Done Reading Theoretical Digest File\n";
STDOUT->autoflush(1);

@pep_array = stripcomments(@pep_array);
@pep_headerarray = split ' ', shift(@pep_array);
for($ctr = 0; $ctr <= $#pep_headerarray; $ctr++)
{
	#print "$pep_headerarray[$ctr]\n";
	$pep_name2col{ $pep_headerarray[$ctr] } = $ctr;
	$pep_col2name{ $ctr } = $pep_headerarray[$ctr];
}

#
# Set up ordered arrays (one for N14, one for N15) for the theoretical digest for easy searching
# Include an index based on integer mass
#
@pep_n14 = sort by_n14 @pep_array;
@pep_n15 = sort by_n15 @pep_array;

for($ctr = 0; $ctr <= $#pep_n14; $ctr++)
{
	@array_n14 = split ' ', $pep_n14[$ctr];
	@array_n15 = split ' ', $pep_n15[$ctr];
	
	for($ctr2 = 0; $ctr2 <= $#array_n14; $ctr2++)
	{
		$digest_n14[$ctr]{ $pep_col2name{$ctr2} } = $array_n14[$ctr2];
		$digest_n15[$ctr]{ $pep_col2name{$ctr2} } = $array_n15[$ctr2];
	}

	#
	# Get out the mass values as integers, and then use them to index the array for faster searching
	# $n14_index{$n14_mz_int} -> This contains the first array index of theoretical peptides with that N14 integer mass
	# So you would start searching at ($n14_mz_int - 1) to be safe
	#
	$n14_mz_int = sprintf("%d", $digest_n14[$ctr]{"n14mass"});
	$n15_mz_int = sprintf("%d", $digest_n15[$ctr]{"n15mass"});

	if(defined($n14_index{ $n14_mz_int }))
	{
	}
	else
	{
		$n14_index{ $n14_mz_int } = $ctr;
	}
	
	if(defined($n15_index{ $n15_mz_int }))
	{
	}
	else
	{
		$n15_index{ $n15_mz_int } = $ctr;
	}
}

#
# Determine the m/z bounds of the digest file
#
$smallest_n14 = $digest_n14[0]{"n14mass"};
$smallest_n15 = $digest_n15[0]{"n15mass"};
$largest_n14 = $digest_n14[ $#pep_n14 ]{"n14mass"};
$largest_n15 = $digest_n15[ $#pep_n15 ]{"n15mass"};

if($smallest_n14 <= $smallest_n15)
{
	$smallest_mz = $smallest_n14;
}
else
{
	$smallest_mz = $smallest_n15;
}

if($largest_n14 >= $largest_n15)
{
	$largest_mz = $largest_n14;
}
else
{
	$largest_mz = $largest_n15;
}

# Read in Peak List
#
# This is stupid - $ailtype is used by both msc01_findmatches.pl and msc010203_mascotmsmsmatches.pl
#
# 0, 1, 2, 3 -> Processed Here
# 4, 5 -> Mascot/Sequest, used by msc010203
# 6, 7 -> Processed Here
#
#
if($ailtype == 0)
{
	# Peak List Type 0: AIL File
	print "Reading ail File\n";
	STDOUT->autoflush(1);
	@peakfiledata = parse_ailfile($ailfile);
}
elsif($ailtype == 1)
{
	print "Reading Qualitative Analysis Peak List\n";
	STDOUT->autoflush(1);
	@peakfiledata = parse_mh2file($ailfile);
}
elsif($ailtype == 2)
{
	# Peak List Type 1: MassHunter 2 (Qualitative Analysis)
	print "Reading xcms monoisotopic peaks\n";
	STDOUT->autoflush(1);
	#@peakfiledata = parse_mh2file($ailfile);
	@peakfiledata = parse_camera_mono($ailfile);
}
elsif($ailtype == 3)
{
	# featureML output from OpenMS
	print "Reading OpenMS featureML peaks\n";
	STDOUT->autoflush(1);
	@peakfiledata = parse_featureml($ailfile);
}
elsif($ailtype == 6)
{
	# Peak List Type: AIL File from new Mass Profile
	print "Reading ail File (new Mass Profiler)\n";
	STDOUT->autoflush(1);
	@peakfiledata = parse_ailfile_newmassprofiler($ailfile);
}
elsif($ailtype == 7)
{
	# Peak List Type: Anna Popova Matlab Type for MALDI RNA Data
	print "Reading AMP Matlab MALDI Peaklist\n";
	STDOUT->autoflush(1);
	@peakfiledata = parse_ailfile_ampmatlabmaldi($ailfile);
}
else
{
	print "ERROR: Unknown peak list type |$ailtype|\n";
	STDOUT->autoflush(1);
}

for($ctr = 0; $ctr <= $#{$peakfiledata[0]}; $ctr++)
{
	$peakfile_headerarray[$ctr] = $peakfiledata[0][$ctr];
}

#
# Read in the ail file
# Store the header titles in an array - probably not necessary, as ail file is consistent in terms of format
#
#print "Reading ail File\n";
#STDOUT->autoflush(1);
#@ail_array = readfile($ailfile);

#
# The ail header line is not actually space delimited, so we will build our own version with new labels
#
#@ail_headerarray = split ' ', shift(@ail_array);
#@ail_headerarray = ();
#@ail_headerarray = ( "featureid", "multimer", "ailcharge", "ailiontype", "rt", "mz", "abundance" );
#for($ctr = 0; $ctr <= $#ail_headerarray; $ctr++)
#{
	#print "$ail_headerarray[$ctr]\n";
#	$ail_name2col{ $ail_headerarray[$ctr] } = $ctr;
#	$ail_col2name{ $ctr } = $ail_headerarray[$ctr];
#}

open OUTPUT, ">$outputfile" or die "Can't open $outputfile\n";

#
# Print out the complete header line, as we start printing matches below
#
#$ail_headerline = join ' ', @ail_headerarray;
#$pep_headerline = join ' ', @pep_headerarray;
#print OUTPUT "ailid $ail_headerline sampleid $pep_headerline matchtype ppm\n";
$peakfile_headerline = join ' ', @peakfile_headerarray;
$pep_headerline = join ' ', @pep_headerarray;
print OUTPUT "$peakfile_headerline sampleid $pep_headerline matchtype ppm\n";

print "Searching for Matches\n";
STDOUT->autoflush(1);
$num_match = 0;
# Start from 1, 0 index is heaader index
PEAK: for($ctr = 1; $ctr <= $#peakfiledata; $ctr++)
{
	$ail_mz = $peakfiledata[$ctr][2];
	$ail_mz_int = sprintf("%.0f", $ail_mz);
	$ail_charge = $peakfiledata[$ctr][4];
	
	if($ail_mz < ($smallest_mz-1) || $ail_mz > ($largest_mz+1))
	{
		next PEAK;
	}
	
	# N14 Matches
	if($mode == 0 || $mode == 1)
	{
		$incr = 2;
		N14START: while(1)
		{
			if(defined( $n14_index{ $ail_mz_int - $incr } ))
			{
				$start = $n14_index{ $ail_mz_int - $incr };
				last N14START;
			}
			elsif( ($ail_mz_int - $incr) < $smallest_mz )
			{
				$start = 0;
				last N14START;
			}			
			$incr++;
		}
		
		$incr = 2;
		N14END: while(1)
		{
			if(defined( $n14_index{ $ail_mz_int + $incr } ))
			{
				$end = $n14_index{ $ail_mz_int + $incr };
				last N14END;
			}
			elsif( ($ail_mz_int + $incr) > $largest_mz )
			{
				$end = $#pep_n14;
				last N14END;
			}			
			$incr++;
		}		
		
		if($start < 0)
		{
			$start = 0;
		}
		if($end > $#pep_n14)
		{
			$end = $#pep_n14;
		}
	
		N14SCAN: for($ctr2 = $start; $ctr2 <= $end; $ctr2++)
		{
			#if($ail_charge != $digest_n14[$ctr2]{ "charge" })
			if(abs($ail_charge) != abs($digest_n14[$ctr2]{ "charge" }))
			{
				next N14SCAN;
			}
		
			$mzdiff = $ail_mz - $digest_n14[$ctr2]{ "n14mass" };
			$ppm = 1000000*($mzdiff / $digest_n14[$ctr2]{ "n14mass" });

			#if(abs($ppm) > $ppm_thresh)
			#{
			#	next N14SCAN;
			#}
			if($ppm < ($ppm_offset - $ppm_thresh) || $ppm > ($ppm_offset + $ppm_thresh))
			{
				next N14SCAN;
			}

			$ppm_int = sprintf("%.1f", $ppm);
			
			$outstring = join ' ', @{$peakfiledata[$ctr]};
			
			print OUTPUT "$outstring $sample_id $pep_n14[$ctr2] N14 $ppm_int\n";
			$num_match++;
		}
	}

	# N15 Matches
	if($mode == 0 || $mode == 2)
	{
		$incr = 2;
		N15START: while(1)
		{
			if(defined( $n15_index{ $ail_mz_int - $incr } ))
			{
				$start = $n15_index{ $ail_mz_int - $incr };
				last N15START;
			}
			elsif( ($ail_mz_int - $incr) < $smallest_mz )
			{
				$start = 0;
				last N15START;
			}			
			$incr++;
		}
		$incr = 2;
		
		N15END: while(1)
		{
			if(defined( $n15_index{ $ail_mz_int + $incr } ))
			{
				$end = $n15_index{ $ail_mz_int + $incr };
				last N15END;
			}
			elsif( ($ail_mz_int + $incr) > $largest_mz )
			{
				$end = $#pep_n15;
				last N15END;
			}			
			$incr++;
		}	
		
		if($start < 0)
		{
			$start = 0;
		}
		if($end > $#pep_n15)
		{
			$end = $#pep_n15;
		}
		
		N15SCAN: for($ctr2 = $start; $ctr2 <= $end; $ctr2++)
		{
			#if($ail_charge != $digest_n15[$ctr2]{ "charge" })
			if(abs($ail_charge) != abs($digest_n15[$ctr2]{ "charge" }))
			{
				next N15SCAN;
			}
		
			$mzdiff = $ail_mz - $digest_n15[$ctr2]{ "n15mass" };
			$ppm = 1000000*($mzdiff / $digest_n15[$ctr2]{ "n15mass" });

			#if(abs($ppm) > $ppm_thresh)
			#{
			#	next N15SCAN;
			#}
			if($ppm < ($ppm_offset - $ppm_thresh) || $ppm > ($ppm_offset + $ppm_thresh))
			{
				next N15SCAN;
			}
						
			$ppm_int = sprintf("%.1f", $ppm);

			$outstring = join ' ', @{$peakfiledata[$ctr]};
			
			print OUTPUT "$outstring $sample_id $pep_n15[$ctr2] N15 $ppm_int\n";
			$num_match++;
		}
	}	
}



print "Found $num_match matches to within $ppm_thresh ppm\n";
print "Finished Searching for Matches\n";
STDOUT->autoflush(1);

################################################################################
################################################################################

####################
# Reading in Peak Files with Subroutines
#
# Subroutines must:
#	1) Accept the file name as an argument
#	2) Return the data formatted as a single 2D array
#		$array[0][j] is the header line containing (j+1) column titles
#		$array[1][j] -> $array[n][j] are n data lines
#		
#		$array[0][0] -> "ailid" this is generated by the subroutine as line number
#						Really just needs to be any unique number identifying the peak
#		$array[0][1] -> "rt" retention time of peak
#		$array[0][2] -> "mz" mass/charge of monoisotopic peak
#		$array[0][3] -> "abundance" magnitude of peak
#		$array[0][4] -> "ailcharge" charge state of peak
#
#	Note, you must preserve the header names "ailid" etc or they will be ignored.
#
#	These first five (j=0 -> j=4) elements are mandatory, as are the column names
#	All other information from the peak file can be included/discarded as desired
#
#	If other columns are included, their names must be prefaced with "pks_"
#   Otherwise, they will be discarded during consolidation in step 4
####################

# Agilent MassProfiler Format
# Acquired Included List
sub parse_ailfile
{
	my @params = @_;
	my $peakfilename = shift(@params);
	
	my @peakfilearray = readfile($peakfilename);
	
	my $peakfileheader = shift(@peakfilearray);

	my @peakdata = ();
	
	# feature ID	multimer no.	z	mono peak	rt	m/z	abundance
	# 0             1               2   3           4   5   6
	# We build this by hand since the AIL file isn't properly delimited
	# Also, it has lame column names
	my @peakfileheaderarray = ("ailid", "rt", "mz", "abundance", "ailcharge", "pks_featureid", "pks_multimer", "pks_ailiontype" );
	
	for(my $ctr = 0; $ctr <= $#peakfileheaderarray; $ctr++)
	{
		$peakdata[0][$ctr] = $peakfileheaderarray[$ctr];
	}
		
	for(my $ctr = 0; $ctr <= $#peakfilearray; $ctr++)
	{
		my @array = split ' ', $peakfilearray[$ctr];
		# Check in case ail file is in csv format for some strange reason
		if($#array == 0)
		{
			@array = split /\,/, $peakfilearray[$ctr];
		}
		
		my $peakid = $ctr + 1;
		$peakdata[$peakid][0] = $peakid; # "ailid" ID Generated from line number 
		$peakdata[$peakid][1] = $array[4]; # "rt"
		$peakdata[$peakid][2] = $array[5]; # "mz"
		$peakdata[$peakid][3] = $array[6]; # "abundance"
		$peakdata[$peakid][4] = $array[2]; # "ailcharge" Charge as defined in peak file
		$peakdata[$peakid][5] = $array[0]; # "featureid"
		$peakdata[$peakid][6] = $array[1]; # "multimer"
		$peakdata[$peakid][7] = $array[3]; # "ailiontype"
	}
	
	return @peakdata;
}

sub parse_ailfile_newmassprofiler
{
	my @params = @_;
	my $peakfilename = shift(@params);
	
	my @peakfilearray = readfile($peakfilename);
	
	my $peakfilejunk = shift(@peakfilearray);
	my $peakfileheader = shift(@peakfilearray);
	
	my @peakdata = ();
	
	# TargetedMSMSTable
	# On,Prec. m/z,Z,Ret. time (min),Delta ret. time (min)
	# 0  1         2 3               4
	# TRUE,185.943426565839,6,6.683803,0.1
	# TRUE,747.440993827851,3,64.03803,0.1
	# TRUE,1120.65785251593,2,64.03803,0.1
	# TRUE,269.066780599291,1,5.762038,0.1
	
	my @peakfileheaderarray = ("ailid", "rt", "mz", "abundance", "ailcharge" );

	for(my $ctr = 0; $ctr <= $#peakfileheaderarray; $ctr++)
	{
		$peakdata[0][$ctr] = $peakfileheaderarray[$ctr];
	}
	
	for(my $ctr = 0; $ctr <= $#peakfilearray; $ctr++)
	{
		my @array = split /\,/, $peakfilearray[$ctr];
		
		my $peakid = $ctr + 1;
		$peakdata[$peakid][0] = $peakid; # "ailid" ID Generated from line number 
		$peakdata[$peakid][1] = $array[3]; # "rt"
		$peakdata[$peakid][2] = $array[1]; # "mz"
		$peakdata[$peakid][3] = 1; # "abundance" NO ABUNDANCE IN THIS FORMAT
		$peakdata[$peakid][4] = $array[2]; # "ailcharge" Charge as defined in peak file
	}
	
	return @peakdata;
}

# Agilent Qualitative Analysis
sub parse_mh2file
{
	my @params = @_;
	my $peakfilename = shift(@params);
	
	my @peakfilearray = readfile($peakfilename);
	
	my %peakname2col = ();
	
	while(1)
	{
		$peakfileheader = shift(@peakfilearray);
		#my @temparray = split /\,/, $peakfileheader;
		my @temparray = quotewords(",", 0, $peakfileheader);

		
		%peakname2col = ();
		for(my $tempctr = 0; $tempctr <= $#temparray; $tempctr++)
		{
			$peakname2col{$temparray[$tempctr]} = $tempctr;
			#print "$temparray[$tempctr] -> $tempctr\n";
		}
		
		if(defined($peakname2col{"RT"}) && defined($peakname2col{"Base Peak"}))
		{
			#print "LAST\n";
			last;
		}
	
		if($#peakfilearray == 0 || $#peakfilearray == -1)
		{
			print "ERROR: Does not appear to be a Qualitative Analysis CSV file\n";
			STDOUT->autoflush(1);
			die;
		}
	}

	my @peakfileheaderarray = ("ailid", "rt", "mz", "abundance", "ailcharge", "pks_featureid" );

	for(my $ctr = 0; $ctr <= $#peakfileheaderarray; $ctr++)
	{
		$peakdata[0][$ctr] = $peakfileheaderarray[$ctr];
	}

	# Base charge on "Min Z"
	for(my $ctr = 0; $ctr <= $#peakfilearray; $ctr++)
	{
		@array = quotewords(",", 0, $peakfilearray[$ctr]);
		
		my $peakid = $ctr + 1;
		$peakdata[$peakid][0] = $peakid; # "ailid" ID Generated from line number 
		$peakdata[$peakid][1] = $array[$peakname2col{"RT"}]; # "rt"
		$peakdata[$peakid][2] = $array[$peakname2col{"Base Peak"}]; # "mz"
		$peakdata[$peakid][3] = $array[$peakname2col{"Vol"}]; # "abundance"
		$peakdata[$peakid][4] = $array[$peakname2col{"Min Z"}]; # "ailcharge" Charge as defined in peak file
		$peakdata[$peakid][5] = $array[$peakname2col{"Group"}]; # "featureid"
	}

	# Base charge on nearest integer to Mass/"Base Peak"
	for(my $ctr = 0; $ctr <= $#peakfilearray; $ctr++)
	{
		@array = quotewords(",", 0, $peakfilearray[$ctr]);
		
		my $charge_float = $array[$peakname2col{"Mass"}]/$array[$peakname2col{"Base Peak"}];
		my $charge_int = sprintf("%.0f", $charge_float);
		
		my $peakid = $ctr + 1;
		$peakdata[$peakid][0] = $peakid; # "ailid" ID Generated from line number 
		$peakdata[$peakid][1] = $array[$peakname2col{"RT"}]; # "rt"
		$peakdata[$peakid][2] = $array[$peakname2col{"Base Peak"}]; # "mz"
		$peakdata[$peakid][3] = $array[$peakname2col{"Vol"}]; # "abundance"
		$peakdata[$peakid][4] = $charge_int; # integer of Mass/"Base Peak"
		$peakdata[$peakid][5] = $array[$peakname2col{"Group"}]; # "featureid"
	}

	
	return @peakdata;
}

# xcms/CAMERA .csv file -> extract monoisotopic peaks from annotation
sub parse_camera_mono
{
	print "Parsing xcms/CAMERA peaklist\n";
	STDOUT->autoflush(1);

	my @params = @_;
	my $peakfilename = shift(@params);
	
	my @peakfilearray = readfile($peakfilename);
	
	my %peakname2col = ();
	
	my $header = shift(@peakfilearray);
	my @headerarray = quotewords(",", 0, $header);
	for(my $ctr = 0; $ctr <= $#headerarray; $ctr++)
	{
		$peakname2col{$headerarray[$ctr]} = $ctr;
	}

	# Build New Array for Actual Use
	my @peakfileheaderarray = ("ailid", "rt", "mz", "abundance", "ailcharge", "pks_sn");

	for(my $ctr = 0; $ctr <= $#peakfileheaderarray; $ctr++)
	{
		$peakdata[0][$ctr] = $peakfileheaderarray[$ctr];
	}
	
	my $mono_count = 0;
	for(my $ctr = 0; $ctr <= $#peakfilearray; $ctr++)
	{
		#my @array = split /\,/, $peakfilearray[$ctr];
		my @array = quotewords(",", 0, $peakfilearray[$ctr]);
		my $peakid = $ctr + 1;
		
		my $isotopes = $array[$peakname2col{"isotopes"}];
		
		if($isotopes !~ /\[M\]/)
		{
			next;
		}
		
        
		# Extract charge from the isotopes annotation
		#[435][M]+ -> BOO
        #[438][M]2+ -> 2
        #[5047][M]3+ -> 3
        # +1 charge not being detected properly (shows up as 0) (5/20/11) thanks sschen
        $isotopes =~ /\[M\](\d+)\+/;
		my $charge = $1;
        
        if($charge eq "")
        {
            $charge = 1;
        }
		
		# Maintain original $peakid so it corresponds to file line numbers
		# Generate this separate counter so we don't get a sparse array
		$mono_count++;
		
		$peakdata[$mono_count][0] = $peakid; # "Generated from line number"
		$peakdata[$mono_count][1] = $array[$peakname2col{"rt"}] / 60; # rt # convert s->min
		$peakdata[$mono_count][2] = $array[$peakname2col{"mz"}]; # m/z
		$peakdata[$mono_count][3] = $array[$peakname2col{"into"}]; # intensity(original)
		$peakdata[$mono_count][4] = $charge; # charge
		$peakdata[$mono_count][5] = $array[$peakname2col{"sn"}]; # signal/noise
	
		#$outline = join ' ', @{$peakdata[$mono_count]};
		#print "$outline\n";
		#STDOUT->autoflush(1);
	}
	
	print "$mono_count monoisotopic peaks in xcms peaks file\n";
	STDOUT->autoflush(1);
	
	return @peakdata;
}

sub parse_featureml
{
	print "Parsing featureML peaklist\n";
	STDOUT->autoflush(1);
	
	my @params = @_;
	my $peakfilename = shift(@params);
	
	my @peakfilearray = readfile($peakfilename);
	
	my @peakfileheaderarray = ("ailid", "rt", "mz", "abundance", "ailcharge");

	for(my $ctr = 0; $ctr <= $#peakfileheaderarray; $ctr++)
	{
		$peakdata[0][$ctr] = $peakfileheaderarray[$ctr];
	}

	my $feature_count = 0;
	for(my $ctr = 0; $ctr <= $#peakfilearray; $ctr++)
	{
		my $rt = -1;
		my $mz = -1;
		my $abundance = -1;
		my $charge = -1;
	
		if($peakfilearray[$ctr] =~ /feature id/)
		{
			$feature_count++;
		
			$ctr++;
			$peakfilearray[$ctr] =~ /\<.*\>(\d+)\.(\d+)\<.*\>/;
			$rt = join '.', $1, $2;
			
			$ctr++;
			$peakfilearray[$ctr] =~ /\<.*\>(\d+)\.(\d+)\<.*\>/;
			$mz = join '.', $1, $2;
			
			$ctr++;
			$peakfilearray[$ctr] =~ /\<.*\>(.*)\<.*\>/;
			$abundance = $1;
			
			$ctr += 4;
			$peakfilearray[$ctr] =~ /\<.*\>(.*)\<.*\>/;
			$charge = $1;
			
			$peakdata[$feature_count][0] = $feature_count - 1; #
			$peakdata[$feature_count][1] = $rt / 60;
			$peakdata[$feature_count][2] = $mz;
			$peakdata[$feature_count][3] = $abundance;
			$peakdata[$feature_count][4] = $charge;
		}
	}
	
	print "$feature_count features in the featureML file\n";
	STDOUT->autoflush(1);
	
	return @peakdata;
}

sub parse_ailfile_ampmatlabmaldi
{
	print "Parsing AMP Matlab MALDI peaklist\n";
	STDOUT->autoflush(1); # Just flushes the buffer to be output
	
	# Reading in the input parameters that were sent to this subroutine
	my @params = @_;
	# shift takes the first array element and assigns it to $peakfilename
	my $peakfilename = shift(@params);
	
	# Read the actual file from disk
	# Each line of the file gets stored in the @peakfilearray
	my @peakfilearray = readfile($peakfilename);
	
	# Now we parse the actual peak list file
	my @peakfileheaderarray = ("ailid", "rt", "mz", "abundance", "ailcharge");
		
	# Apply that peak list header to the actual data we will return
	for(my $ctr = 0; $ctr <= $#peakfileheaderarray; $ctr++)
	{
		$peakdata[0][$ctr] = $peakfileheaderarray[$ctr];
	}
	
	my $num_lines = 0;
	
	# Go through the peak list file line by line
	for(my $ctr = 0; $ctr <= $#peakfilearray; $ctr++)
	{
		$num_lines++;
		
		# Let's split the line into individual columns
	
		# This would split on commas
		#my @array = quotewords(",", 0, $peakfilearray[$ctr]);

		# This splits on spaces and places individual column elements into @array
		#my @array = quotewords(" ", 0, $peakfilearray[$ctr]);
		my @array = split ' ', $peakfilearray[$ctr];

		# Offset the input by 1 so we don't overwrite the header in the first row of @peakdata
		my $row_number = $ctr + 1;
		
		my $ailid = $array[0];
		# There is no RT!  Let's just set an arbitrary value
		# Let's totally not use 0 in case any issue arise from having a negative extraction range, since you extract +/- 0.1 min
		# Another fun detail - always use minutes!  You might have to convert from seconds sometimes
		my $rt = 1;
		my $mz = $array[1];
		my $abundance = $array[2];
		# No charge is given but Anna says they are mostly +1
		my $ailcharge = 1;

        $peakdata[$row_number][0] = $ailid;
        $peakdata[$row_number][1] = $rt; #
        $peakdata[$row_number][2] = $mz;
        $peakdata[$row_number][3] = $abundance; #
        $peakdata[$row_number][4] = $ailcharge; #
		
		if(0 == 1)
		{
			print "$peakdata[$row_number][0] | $peakdata[$row_number][1] | $peakdata[$row_number][2] | $peakdata[$row_number][3] | $peakdata[$row_number][4]\n";
		}
	}
	
	print "$num_lines lines found in this crazy peak list file\n";
	STDOUT->autoflush(1);
	
	# Return the actual data that we built
	return @peakdata;
}

sub by_n14
{
	my @array_a = split ' ', $a;
	my @array_b = split ' ', $b;
	$array_a[ $pep_name2col{"n14mass"} ] <=> $array_b[ $pep_name2col{"n14mass"} ];
}

sub by_n15
{
	my @array_a = split ' ', $a;
	my @array_b = split ' ', $b;
	$array_a[ $pep_name2col{"n15mass"} ] <=> $array_b[ $pep_name2col{"n15mass"} ];
}


