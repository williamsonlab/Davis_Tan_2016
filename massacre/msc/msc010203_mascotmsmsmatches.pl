#!/usr/bin/perl

# Perform all of Steps 1 2 and 3 for Mascot-based MS/MS IDs

use Getopt::Long;
use IO::Handle;
use Text::ParseWords;
require "$ENV{'MASSACRE_PATH'}/msc/mikesubs.pl";
require "$ENV{'MASSACRE_PATH'}/msc/mikevars.pl";

$pid = $$;

#$digestfile = "";
$mascotfile = "";
$outputfile = "";

#%digestname2col = ();
%digest = ();
#$digestheader = "";

%mascotname2col = ();
%mascotraw = ();
%mascotdata = ();

# ssc101911 changed mzrange_low from 400 to 200
$mzrange_low = 200;
$mzrange_high = 2000;

$ailtype = 0;

#&GetOptions("digest=s" => \$digestfile, "ail=s" => \$mascotfile, "output=s" => \$outputfile);
&GetOptions("ail=s" => \$mascotfile, "output=s" => \$outputfile, "mzrange_low=i" => \$mzrange_low, "mzrange_high=i" => \$mzrange_high, "ailtype=i" => \$ailtype, "rt_max_int_fraction=f" => \$rt_max_int_fraction);

#if($digestfile eq "" || $mascotfile eq "" || $outputfile eq "")
if($mascotfile eq "" || $outputfile eq "")
{
	#print "ERROR: Must define --digest=something --ail=something --output=something\n";
	print "ERROR: Must define --ail=something --output=something\n";
	STDOUT->autoflush(1);
	die;
}

#####
# No longer going to use a digest file
# Generate all digest-related parameters dynamically
#####
#parse_digest($digestfile);

#@temparray = keys (%digest);
#$num = $#temparray + 1;
#print "$num entries in theoretical digest\n";
#STDOUT->autoflush(1);
#####
# /
#####

# This is stupid - $ailtype is used by both msc01_findmatches.pl and msc010203_mascotmsmsmatches.pl
#
# 0, 1, 2, 3 -> msc01
# 4, 5 -> Mascot/Sequest, used by msc010203 (Here)
# 6 -> msc01
#
#
print "Sequencing scan intensity threshold: $rt_max_int_fraction\n";
STDOUT->autoflush(1);

if($ailtype == 4)
{
	parse_mascot($mascotfile);


}
elsif($ailtype == 5)
{
	parse_sequest($mascotfile);
}
else
{
	print "ERROR: ailtype |$ailtype| unknown\n";
	STDOUT->autoflush(1);
	die;
}

@temparray = keys (%mascotraw);
$num = $#temparray + 1;
print "$num peptides in mascot export file\n";
STDOUT->autoflush(1);

@temparray = keys (%mascotdata);
$num = $#temparray + 1;
print "$num peptides processed\n";
STDOUT->autoflush(1);

output_mascot($mascotfile);

##################
# parse_mascot() #
##################
#
# Need:
#	ailid
#	rt
#	mz
#	abundance
#	ailcharge
#
sub parse_mascot
{
	my @params = @_;
	
	my @filearray = readfile($mascotfile);
	
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
			if(defined($mascotname2col{"corr_rt"}))
			{
				print "Corrected RTs Detected\n";
				STDOUT->autoflush(1);
			}
			last;
		}
	}

	my $mpro;
	
	# Organize all entries into a hash (%mascotraw) sorted by ion (protein_peptide_charge)
	for(my $ctr = 0; $ctr <= $#filearray; $ctr++)
	{
		my $mseq;
		my $mcrg;

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
			$mdesc = $array[$mascotname2col{"prot_desc"}];
			
            # In the LTQ data that sschen got from the Siuzdak lab, gi| was appended to the start of accession codes
            # Strip gi| out
            $mpro =~ s/gi\|//;
		}
		
		$mseq = $array[$mascotname2col{"pep_seq"}];
		$mcrg = $array[$mascotname2col{"pep_exp_z"}];
		$mmiss = $array[$mascotname2col{"pep_miss"}];
		$mstart = $array[$mascotname2col{"pep_start"}];
		$mend = $array[$mascotname2col{"pep_end"}];

		$rtcorrfile = "";
		$rtcorrfile = $array[$mascotname2col{"rtcorrfile"}];

		my $mpro2;
		# msc00_adjustmsmsrt.pl does filtering based on msms_proteinkeep
		# Check for the accession number in %acc2pro_all (mikevars.pl)
		# %acc2pro_all is a superset of %acc2pro
		# If not present, just use boring accession number
		if(defined($acc2pro_all{$mpro}))
		{
			$mpro2 = $acc2pro_all{$mpro};
		}
		else
		{
			$mpro2 = $mpro;
		}
		my $key = join '_', $mpro2, $mseq, $mcrg;

		push(@{ $mascotraw{$key} }, $filearray[$ctr]);
		$mascotmiss{$key} = $mmiss;
		$mascotstart{$key} = $mstart;
		$mascotend{$key} = $mend;
		$mascotdesc{$key} = $mdesc;
		$mascotacc{$key} = $mpro; # Map translated protein back to accession
		
		push(@{ $mascotrtcorrfile{$key} }, $rtcorrfile);
	}
	
	# Process %mascot_raw, grouping and processing multiple entries from the same ion
	# Essentially the only interesting thing is RT - this needs to be averaged
	# Also worth tracking up to 2 possible pep_exp_mz, 14N and 15N versions
	my $ailid_n14 = 0;
	my $ailid_n15 = 1;
	foreach my $key (keys %mascotraw)
	{
		#@{$mascotraw{$key}} = sort by_score @{$mascotraw{$key}};

		#
		# Sort this my the maximum intensity observed in the RT correction plots
		# We can then filter based on a threshold of the maximum maximum intensity
		# This should eliminate long tails of elution, especially for really sensitive specs like the Triple-TOF
		#
		@sorted_raw = sort by_rtmaxint_mascot @{$mascotraw{$key}};
		
		#print "IN: $#{$mascotraw{$key}}\n";
		
		@max_rt_array = quotewords(",", 0, $sorted_raw[ $#sorted_raw ]);
		$max_rt_max_int = $max_rt_array[$mascotname2col{"rt_max_int"}];
		$rt_max_int_cutoff = $rt_max_int_fraction * $max_rt_max_int;
		
		#print "CUTOFF: $rt_max_int_cutoff\n";
		
		@new_array = ();
		@new_array_rtcorrfile = ();
		
		for(my $raw_ctr = 0; $raw_ctr <= $#{$mascotraw{$key}}; $raw_ctr++)
		{
			my @temp_array = quotewords(",", 0, $mascotraw{$key}[$raw_ctr]);
			my $rt_max_int = $temp_array[$mascotname2col{"rt_max_int"}];
			#print "\tRT_MAX_INT: $rt_max_int\n";
			if($rt_max_int >= $rt_max_int_cutoff)
			{
				#print "PUSH\n";
				push(@new_array, $mascotraw{$key}[$raw_ctr]);
				push(@new_array_rtcorrfile, $mascotrtcorrfile{$key}[$raw_ctr]);
			}
		}
		
		@{$mascotraw{$key}} = @new_array;
		
		#print "OUT: $#{$mascotraw{$key}} $#new_array\n";
		#print "OUT: $#{$mascotraw{$key}}\n";

		#
		# End filtering by $rt_max_int
		#
		
		my @array = split /\_/, $key;
		my $protein = $array[0];
		my $sequence = $array[1];
		my $charge = $array[2];
				
		my $rt = 0;
		my $orig_rt = 0;
		my $abundance_n14 = 0;
		my $abundance_n15 = 0;
		my $mz_n14 = 0;
		my $mz_n15 = 0;
		
		my $num_keep = 0;
		my $num_keep_n14 = 0;
		my $num_keep_n15 = 0;
		
		do_pseudo_digest($key);
		
		#my $rtcorrfile_line = join '#', @{$mascotrtcorrfile{$key}};
		my $rtcorrfile_line = join '#', @new_array_rtcorrfile;
		
		#####
		#####
		#####
		# Change mz < mz_low, mz > mz_high to only check for N14 OR N15, not both
		# Should there be averaging of mz_n14 or mz_n15?  Just overwrites with latest one
		#####
		#####
		#####
		
		MASCOTRAW: for(my $ctr = 0; $ctr <= $#{$mascotraw{$key}}; $ctr++)
		{
			my $n14n15 = -1; # 6/20/11 This doesn't actually seem to do anything
			my @array1 = quotewords(",", 0, $mascotraw{$key}[$ctr]);
			my @array2 = split ' ', $array1[$mascotname2col{"pep_scan_title"}];
			
			# Check the modifications column
			# Previously, anything here had caused things to be skipped, as the only thing ever present was Oxidation
			# Now, this column reports labeled amino acids.  These are generally reported with the prefix "Label"
			# So this skips peptides if they have a modification, but they don't have a modification of the type "Label"
			# In theory this could cause a problem if both an oxidation and a labeled amino acid occurred in the same peptide.  So far heavy modifications something more sophisticated needs to be done
			if($array1[$mascotname2col{'pep_var_mod'}] ne "" && !($array1[$mascotname2col{'pep_var_mod'}] =~ m/Label/))
			{
				print "Skipping $key due to mod |$array1[$mascotname2col{'pep_var_mod'}]|\n";
				next;
			}
			
			my $mz = $array1[$mascotname2col{"pep_exp_mz"}];
			# Consult digest to determine if 14N or 15N
            # Note (5/4/2011) Previously I had used a threshold of 0.1 to determine if a peptide matched the 14N or 15N version
            # The LTQ data from sschen has really bad mass accuracy, so I am increasing this to 1.0
            # Since 14N and 15N peptides are always more than 1 away from each other, this should not cause any problems
            # An alternative would be to use the pep_calc_mr (calculated mass of peptide, ignoring charge) but I like to keep the observed values in play
			#if(abs($mz - $digest{$key}{"n14mass"}) < 0.1)
			if(abs($mz - $digest{$key}{"n14mass"}) < 1.0)
			{
				$n14n15 = 0;
				$temp_mz_n14 = $mz;
				$abundance_n14 = 1;
				
				if($temp_mz_n14 < $mzrange_low || $temp_mz_n14 > $mzrange_high)
				{
					print "Skipping $key, m/z out of range ($temp_mz_n14, $mzrange_low <-> $mzrange_high)\n";
					STDOUT->autoflush(1);
					next MASCOTRAW;
				}
				
				$mz_n14 += $temp_mz_n14;
				$num_keep_n14++;
			}
			#elsif(abs($mz - $digest{$key}{"n15mass"}) < 0.1)
			elsif(abs($mz - $digest{$key}{"n15mass"}) < 1.0)
			{
				$n14n15 = 1;
				$temp_mz_n15 = $mz;
				$abundance_n15 = 1;

				if($temp_mz_n15 < $mzrange_low || $temp_mz_n15 > $mzrange_high)
				{
					print "Skipping $key, m/z out of range ($temp_mz_n15, $mzrange_low <-> $mzrange_high)\n";
					STDOUT->autoflush(1);
					next MASCOTRAW;
				}

				$mz_n15 += $temp_mz_n15;
				$num_keep_n15++;
			}
            # 6/20/11 - Amino acid labeling is causing more problems.  It can be properly sequenced using a SILAC quantitation method in mascot, but then the m/z value doesn't match either the n14 or n15 mass properly, so massacre was rejecting these
            # So now we just look for the "Label" entry inside the pep_var_mod field, and arbitrarily assign these values to the N15 category
            # These of course aren't really N15, but they kind of mimic N15 in the classic N14/N15 pairs scenario
            elsif($array1[$mascotname2col{'pep_var_mod'}] =~ m/Label/)
            {
				$n14n15 = 2; # Set this to 2 here to distinguish it from true N15, but this is completely ignored so this is really pointless
				$temp_mz_n15 = $mz;
				$abundance_n15 = 1;
                
				if($temp_mz_n15 < $mzrange_low || $temp_mz_n15 > $mzrange_high)
				{
					print "Skipping $key, m/z out of range ($temp_mz_n15, $mzrange_low <-> $mzrange_high)\n";
					STDOUT->autoflush(1);
					next MASCOTRAW;
				}
                
				$mz_n15 += $temp_mz_n15;
				$num_keep_n15++;
			}
			else
			{
				print "Warning: $key mz does not match either n14mass or n15mass ";
				print "|$array1[$mascotname2col{'pep_var_mod'}]|\n";
				print "\t$digest{$key}{'n14mass'} or $digest{$key}{'n15mass'} vs. $mz\n";
				STDOUT->autoflush(1);
				next MASCOTRAW;
			}

			# Corrected RT values as derived in step0
			if(defined($mascotname2col{"corr_rt"}))
			{
				$rt += $array1[$mascotname2col{"corr_rt"}];
				$orig_rt += $array1[$mascotname2col{"rt"}];
			}
			else
			{
				print "Warning: Corrected RT not found!\n";
				STDOUT->autoflush(1);
				# pep_scan_title
				# 1_30S_1to10dil_060409_cap_MSMS.d, MS/MS of 439.7410431 2+ at 33.7646166666667 mins 
				# So RT is the penultimate array element of @array2
				$rt += $array2[$#array2-1];
				$orig_rt += $array2[$#array2-1];
			}
			
			#if($mz_n14 < $mzrange_low)
			#{
			#	print "Skipping $key, m/z out of range ($mz_n14 < $mzrange_low)\n";
			#	STDOUT->autoflush(1);
			#	next MASCOTRAW;
			#}
			#if($mz_n15 > $mzrange_high)
			#{
			#	print "Skipping $key, m/z out of range ($mz_n15 > $mzrange_high)\n";
			#	STDOUT->autoflush(1);
			#	next MASCOTRAW;
			#}
			
			$num_keep++;
		}
		
		if($num_keep <= 0)
		{
			print "Warning: all matches eliminated for $key\n";
			next;
		}
		
		# Average RT and observed m/z values across number of peptides
		# m/z is averaged over N14 and N15 separately.  This means that if there are no N14 or N15 peptides, it will get an m/z values of 0.  This is okay because extraction (msc04) is done based on calculated mass values, not actual mass values
		# RT is averaged over N14 and N15 together.  This is done because having a 0 value for RT would completely throw off the extraction in msc04 (recall that MS/MS data always masquerades as peaks and so both an N14 and N15 entry are read in during step msc04)
		$rt = $rt / $num_keep;
		$orig_rt = $orig_rt / $num_keep;
		if($num_keep_n14)
		{
			$mz_n14 = $mz_n14 / $num_keep_n14;
		}
		if($num_keep_n15)
		{
			$mz_n15 = $mz_n15 / $num_keep_n15;
		}
		
		$mascotdata{$key}{"rt"} = $rt;
		$mascotdata{$key}{"orig_rt"} = $orig_rt;
		$mascotdata{$key}{"ailcharge"} = $charge;
		$mascotdata{$key}{"ailid_n14"} = $ailid_n14;
		$mascotdata{$key}{"ailid_n15"} = $ailid_n15;
		$mascotdata{$key}{"mz_n14"} = $mz_n14;
		$mascotdata{$key}{"mz_n15"} = $mz_n15;
		$mascotdata{$key}{"abundance_n14"} = $abundance_n14;
		$mascotdata{$key}{"abundance_n15"} = $abundance_n15;
		$mascotdata{$key}{"rtcorrfile"} = $rtcorrfile_line;
		
		$ailid_n14 += 2;
		$ailid_n15 += 2;
	}
}

sub parse_sequest
{
	my @params = @_;
	
	my @filearray = readfile($mascotfile);
	
	my $line = shift(@filearray);
	my @array = split /\t/, $line;
	
	%mascotname2col = ();
	for(my $ctr = 0; $ctr <= $#array; $ctr++)
	{
		my $key = lc($array[$ctr]);
		$mascotname2col{$key} = $ctr;
		#print "$key -> $ctr\n";
	}
	
	my $mpro;
	
	# Generate %mascotraw
	for(my $ctr = 0; $ctr <= $#filearray; $ctr++)
	{
		my $mseq;
		my $mcrg;
		
		my @array = split /\t/, $filearray[$ctr];
		
		my @seqarray = split /\./, $array[$mascotname2col{"sequence"}];
		$mseq = $seqarray[1];
		
		my @crgarray = split /\./, $array[$mascotname2col{"filename"}];
		$mcrg = $crgarray[3];
		
		$rtcorrfile = $array[$mascotname2col{"rtcorrfile"}];
		
		my @locusarray = split /\|/, $array[$mascotname2col{"locus"}];
		
		# Looking for format crap|accession|crap
		# Sometimes get crap|accession
		# Check in case you get "stuff" (no |'s)
		if($#locusarray >= 1)
		{
			$mpro = $locusarray[1];
		}
		else
		{
			$mpro = $array[$mascotname2col{"locus"}];
		}
		
		if(substr($mpro, 0, 7) eq "REFSEQ:")
		{
			substr($mpro, 0, 7) = "";
		}
		
		# Replace underscores by dashes since underscores are used as delimiters elsewhere (that's probably not the best choice)
		$mpro =~ s/_/-/g;
		
		$mdesc = $array[$mascotname2col{"descriptive name"}];
		
		if(defined($acc2pro_all{$mpro}))
		{
			$mpro2 = $acc2pro_all{$mpro};
		}
		else
		{
			$mpro2 = $mpro;
		}
		
		my $key = join '_', $mpro2, $mseq, $mcrg;
		
		push(@{ $mascotraw{$key} }, $filearray[$ctr]);
		
		$mascotdesc{$key} = $mdesc;
		$mascotacc{$key} = $mpro;
		
		# Values not available from Sequest output
		$mascotmiss{$key} = -1;
		$mascotstart{$key} = -1;
		$mascotend{$key} = -1;
		
		push(@{ $mascotrtcorrfile{$key} }, $rtcorrfile);
	}
	
	my $ailid_n14 = 0;
	my $ailid_n15 = 1;
	foreach my $key (keys %mascotraw)
	{		
		#
		# Sort this my the maximum intensity observed in the RT correction plots
		# We can then filter based on a threshold of the maximum maximum intensity
		# This should eliminate long tails of elution, especially for really sensitive specs like the Triple-TOF
		#
		@sorted_raw = sort by_rtmaxint_sequest @{$mascotraw{$key}};
		
		#print "IN: $#{$mascotraw{$key}}\n";
		
		#@max_rt_array = quotewords(",", 0, $sorted_raw[ $#sorted_raw ]);
		@max_rt_array = split /\t/, $sorted_raw[ $#sorted_raw ];		
		$max_rt_max_int = $max_rt_array[$mascotname2col{"rt_max_int"}];
		$rt_max_int_cutoff = $rt_max_int_fraction * $max_rt_max_int;
		
		#print "CUTOFF: $rt_max_int_cutoff\n";
		
		@new_array = ();
		@new_array_rtcorrfile = ();
		
		for(my $raw_ctr = 0; $raw_ctr <= $#{$mascotraw{$key}}; $raw_ctr++)
		{
			#my @temp_array = quotewords(",", 0, $mascotraw{$key}[$raw_ctr]);
			my @temp_array = split /\t/, $mascotraw{$key}[$raw_ctr];
			my $rt_max_int = $temp_array[$mascotname2col{"rt_max_int"}];
			#print "\tRT_MAX_INT: $rt_max_int\n";
			if($rt_max_int >= $rt_max_int_cutoff)
			{
				#print "PUSH\n";
				push(@new_array, $mascotraw{$key}[$raw_ctr]);
				push(@new_array_rtcorrfile, $mascotrtcorrfile{$key}[$raw_ctr]);
			}
		}
		
		@{$mascotraw{$key}} = @new_array;
		
		#print "OUT: $#{$mascotraw{$key}} $#new_array\n";
		#print "OUT: $#{$mascotraw{$key}}\n";
		
		#
		# End filtering by $rt_max_int
		#

		
		my @array = split /\_/, $key;
		my $protein = $array[0];
		my $sequence = $array[1];
		my $charge = $array[2];
		
		my $rt = 0;
		my $orig_rt = 0;
		my $abundance_n14 = 0;
		my $abundance_n15 = 0;
		my $mz_n14 = 0;
		my $mz_n15 = 0;
		
		my $num_keep = 0;
		my $num_keep_n14 = 0;
		my $num_keep_n15 = 0;
		
		do_pseudo_digest($key);
		
		#my $rtcorrfile_line = join '#', @{$mascotrtcorrfile{$key}};
		my $rtcorrfile_line = join '#', @new_array_rtcorrfile;

		SEQUESTRAW: for(my $ctr = 0; $ctr <= $#{$mascotraw{$key}}; $ctr++)
		{
			my $n14n15 = -1;
			my @array1 = split /\t/, $mascotraw{$key}[$ctr];
			
			#my $mz = $array1[$mascotname2col{"m+h+"}] / $charge;
			my $mz = ( $array1[$mascotname2col{"calcm+h+"}] + ($charge-1) ) / $charge;
			
			# Consult digest to determine if 14N or 15N
			if(abs($mz - $digest{$key}{"n14mass"}) < 0.1)
			{
				$n14n15 = 0;
				$temp_mz_n14 = $mz;
				$abundance_n14 = 1;
				
				if($temp_mz_n14 < $mzrange_low || $temp_mz_n14 > $mzrange_high)
				{
					print "Skipping $key, m/z out of range ($temp_mz_n14, $mzrange_low <-> $mzrange_high)\n";
					STDOUT->autoflush(1);
					next SEQUESTRAW;
				}
				
				$mz_n14 += $temp_mz_n14;
				$num_keep_n14++;
			}
			elsif(abs($mz - $digest{$key}{"n15mass"}) < 0.1)
			{
				$n14n15 = 1;
				$temp_mz_n15 = $mz;
				$abundance_n15 = 1;

				if($temp_mz_n15 < $mzrange_low || $temp_mz_n15 > $mzrange_high)
				{
					print "Skipping $key, m/z out of range ($temp_mz_n15, $mzrange_low <-> $mzrange_high)\n";
					STDOUT->autoflush(1);
					next SEQUESTRAW;
				}

				$mz_n15 += $temp_mz_n15;
				$num_keep_n15++;
			}
			else
			{
				print "Warning: $key mz does not match either n14mass or n15mass ";
				print "|$array1[$mascotname2col{'pep_var_mod'}]|\n";
				print "\t$digest{$key}{'n14mass'} or $digest{$key}{'n15mass'} vs. $mz\n";
				STDOUT->autoflush(1);
				next SEQUESTRAW;
			}
		
			# Corrected RT values as derived in step0
			if(defined($mascotname2col{"corr_rt"}))
			{
				$rt += $array1[$mascotname2col{"corr_rt"}];
				$orig_rt += $array1[$mascotname2col{"rt"}];
			}
			else
			{
				print "Warning: Corrected RT not found!\n";
				STDOUT->autoflush(1);
				# pep_scan_title
				# 1_30S_1to10dil_060409_cap_MSMS.d, MS/MS of 439.7410431 2+ at 33.7646166666667 mins 
				# So RT is the penultimate array element of @array2
				$rt += $array2[$#array2-1];
				$orig_rt += $array2[$#array2-1];
			}
			
			$num_keep++;
		}
	
		if($num_keep <= 0)
		{
			print "Warning: all matches eliminated for $key\n";
			next;
		}

		$rt = $rt / $num_keep;
		$orig_rt = $orig_rt / $num_keep;
		if($num_keep_n14)
		{
			$mz_n14 = $mz_n14 / $num_keep_n14;
		}
		if($num_keep_n15)
		{
			$mz_n15 = $mz_n15 / $num_keep_n15;
		}

		$mascotdata{$key}{"rt"} = $rt;
		$mascotdata{$key}{"orig_rt"} = $orig_rt;
		$mascotdata{$key}{"ailcharge"} = $charge;
		$mascotdata{$key}{"ailid_n14"} = $ailid_n14;
		$mascotdata{$key}{"ailid_n15"} = $ailid_n15;
		$mascotdata{$key}{"mz_n14"} = $mz_n14;
		$mascotdata{$key}{"mz_n15"} = $mz_n15;
		$mascotdata{$key}{"abundance_n14"} = $abundance_n14;
		$mascotdata{$key}{"abundance_n15"} = $abundance_n15;
		$mascotdata{$key}{"rtcorrfile"} = $rtcorrfile_line;
		
		$ailid_n14 += 2;
		$ailid_n15 += 2;	
	}
}

##################
# parse_digest() #
##################
sub parse_digest
{
	my @params = @_;
	
	my @filearray = readfile($digestfile);
	@filearray = stripcomments(@filearray);
	
	$digestheader = shift(@filearray);
	my @array = split ' ', $digestheader;
	
	for(my $ctr = 0; $ctr <= $#array; $ctr++)
	{
		$digestname2col{$array[$ctr]} = $ctr;
	}
	
	for(my $ctr = 0; $ctr <= $#filearray; $ctr++)
	{
		@array = split ' ', $filearray[$ctr];
		my $protein = $array[$digestname2col{"protein"}];
		my $sequence = $array[$digestname2col{"seq"}];
		my $charge = $array[$digestname2col{"charge"}];
		
		my $key = join '_', $protein, $sequence, $charge;
				
		$digest{$key}{"line"} = $filearray[$ctr];
		$digest{$key}{"n14mass"} = $array[$digestname2col{"n14mass"}];
		$digest{$key}{"n15mass"} = $array[$digestname2col{"n15mass"}];
	}
}

###################
# output_mascot() #
###################
sub output_mascot
{
	open OUTPUT, ">$outputfile" or die "ERROR: Can't open $outputfile\n";

	# Description is unique to MSMS
	# Need the pks_ prefix to avoid being filtered out in Step 4
	print OUTPUT "ailid rt mz abundance ailcharge pks_description pks_acc pks_orig_rt pks_rtcorrfile ";

	print OUTPUT "$digestheader\n";

	# Output two lines as this simulates matchmode=0, Unlabeled/Labeled pairs
	# Pairs such as this appear on two separate lines of .match3 output files
	foreach my $key (keys %mascotdata)
	{
		@output_description_array = split ' ', $mascotdesc{$key};
		$output_description = join '_', @output_description_array;
		$output_description =~ s/,/_/g;
	
		my $out1 = join ' ', $mascotdata{$key}{"ailid_n14"}, $mascotdata{$key}{"rt"}, $mascotdata{$key}{"mz_n14"}, $mascotdata{$key}{"abundance_n14"}, $mascotdata{$key}{"ailcharge"}, $output_description, $mascotacc{$key}, $mascotdata{$key}{"orig_rt"}, $mascotdata{$key}{"rtcorrfile"}, $digest{$key}{"line"};
	
		print OUTPUT "$out1\n";

		my $out2 = join ' ', $mascotdata{$key}{"ailid_n15"}, $mascotdata{$key}{"rt"}, $mascotdata{$key}{"mz_n15"}, $mascotdata{$key}{"abundance_n15"}, $mascotdata{$key}{"ailcharge"}, $output_description, $mascotacc{$key}, $mascotdata{$key}{"orig_rt"}, $mascotdata{$key}{"rtcorrfile"}, $digest{$key}{"line"};
	
		print OUTPUT "$out2\n";
	}
	
	close OUTPUT;
}

# sort by pep_score
sub by_score
{
	my @aa = quotewords(",", 0, $a);
	my @bb = quotewords(",", 0, $b);
	
	$aa[$mascotname2col{"pep_score"}] <=> $bb[$mascotname2col{"pep_score"}];
}

sub by_rtmaxint_mascot
{
	my @aa = quotewords(",", 0, $a);
	my @bb = quotewords(",", 0, $b);
	
	$aa[$mascotname2col{"rt_max_int"}] <=> $bb[$mascotname2col{"rt_max_int"}];
}

sub by_rtmaxint_sequest
{
	#my @aa = quotewords(",", 0, $a);
	#my @bb = quotewords(",", 0, $b);
	my @aa = split /\t/, $a;
	my @bb = split /\t/, $b;
	
	$aa[$mascotname2col{"rt_max_int"}] <=> $bb[$mascotname2col{"rt_max_int"}];
}

sub do_pseudo_digest
{
	my @params = @_;
	my $key = shift(@params);
	my @key_array = split /\_/, $key;

    # In the Williamson database, the protein names
    # are UNIPROT accession keys, such as P0CE47 
    # In the SwissProt database, the protein names
    # are EntryNames EFTU1_ECOLI. This additional underscore
    # screws things up
    my $protein = $key_array[0];
    my $sequence = $key_array[$#key_array - 1];
    my $charge = $key_array[$#key_array];

	my $miss = $mascotmiss{$key};
	my $start = $mascotstart{$key};
	my $end = $mascotend{$key};
	
	my @sequence_array = split '', $sequence;
	my @mod_array = ();
	my $seq_ctr;
	for($seq_ctr = 0; $seq_ctr <= $#sequence_array; $seq_ctr++)
	{
		if($sequence_array[$seq_ctr] eq "C")
		{
			$mod_array[$seq_ctr] = "c";
		}
		else
		{
			$mod_array[$seq_ctr] = "";
		}
	}
	my $mod_string = join '-', @mod_array;
	if($mod_string eq "")
	{
		$mod_string = "-";
	}
	
	open TEMP, ">pseudodigest.$pid" or die "Can't open pseudodigest.$pid\n";
	#print TEMP "seq missed num_copy loc\n";
	print TEMP "seq mod missed num_copy digest_loc protein startres endres digest_notes\n";
	#print TEMP "$sequence $miss 1 $protein,$start,$end\n";
	print TEMP "$sequence $mod_string $miss 1 $protein,$start,$end $protein $start $end 0_PSEUDO\n";
	close TEMP;
	
	# Specify carbamido=1 here as this is considered a fixed modification in mascot
	# At least in the default parameters generally used
	# If somehow it's not, it's okay because masses won't match up
	# 9/26/11 - New digest script doesn't handle modification (it's handled earlier in the digest pipeline). So we ignore "carbamido" and "modmethod" and instead for every "C" in the sequence we explicitly add a "c" modification (see $mod_string above)
	#my @pseudodigest = `$ENV{'MASSACRE_PATH'}/digest/03_calc_mass.pl --carbamido=1 --modmethod=0 --minz=$charge --z=$charge --mzlow=0 --mzhigh=10000 --pseudo=1 < pseudodigest.$pid`;
	my @pseudodigest = `$ENV{'MASSACRE_PATH'}/digest/03_calc_mass.pl --minz=$charge --z=$charge --mzlow=0 --mzhigh=10000 --pseudo=1 < pseudodigest.$pid`;
	# Actual digest is one line, rest is headers and crap, so just take last line
	my $digestline = pop(@pseudodigest);
	chomp($digestline);
	my @linearray = split ' ', $digestline;
	
	# Don't use the local variable here, as we need $digestheader for output
	# Yes, rewriting it every single time is kind of silly
	$digestheader = pop(@pseudodigest);
	chomp($digestheader);
	my @harray = split ' ', $digestheader;
	my $dctr;
	my %dname2col;
	for($dctr = 0; $dctr <= $#harray; $dctr++)
	{
		$dname2col{$harray[$dctr]} = $dctr;
	}
	
	$digest{$key}{"line"} = $digestline;
	$digest{$key}{"n14mass"} = $linearray[$dname2col{"n14mass"}];
	$digest{$key}{"n15mass"} = $linearray[$dname2col{"n15mass"}];
	
	`rm pseudodigest.$pid`;
}
