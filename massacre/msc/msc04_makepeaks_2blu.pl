#!/usr/bin/perl

use Getopt::Long;
use IO::Handle;
use MIME::Base64;
use XML::Twig;
use IO::Zlib;
use Compress::Zlib;

require "$ENV{'MASSACRE_PATH'}/msc/mikesubs.pl";

$t_start = time();

$pid = $$;

$pbstmpdir = `echo \$PBS_O_WORKDIR`;
chomp($pbstmpdir);
print "pbstmpdir -> $pbstmpdir\n";

# Note! This is the same mode as used in msc01 - it helps determine the approach used here (ie read in one or two lines at a time)
# mode = 0 -> Both N14 and N15 matches
# mode = 1 -> N14 Only
# mode = 2 -> N15 Only
$matchmode = 0;

# Extract this many units on either side of the N14/N15 monoisotopic peaks
$mz_expand = 4;
# Extract an additional this many units beyond the N15 monoisotopic peak
$mz_iso_expand = 0; # CURRENTLY UNUSED 2/12/2009
# Include this much in minutes above/below the average retention time
$rt_expand = 0.1;
$rt_expand_lo = 0.0;
$rt_expand_hi = 0.0;
# Threshold for summing points
$mz_bin_tol = 0.005;
# Do we filter based on...
#$filter_abundance = 1;
#$filter_rt = 0;
#$filter_mz = 0;
#$abundance_mincut = 2.5;
#$rt_mincut = 5;
#$rt_maxcut = 50;
#$mz_mincut = 300;
#$mz_maxcut = 1300;

# New Variables 8/8/2007
$sortstyle = 0;
$slice_per_pt = 5;
$mz_per_pt = 1;

# Stuff to deal with residue labeling
$special_range = 0;
$special_range_residues = "";
$special_range_values = "";

#$spectra3D = 0;
$extractstyle = 0;
$extractbin = 0.05;

$ailtype = 0; # This is used to help calculate mz_bin_tol

$gnuplot = "gnuplot"; # urocyclon + xtopher both work fine now

&GetOptions("input=s" => \$inputfile, "output=s" => \$outputfile, "csv=s" => \$csvfile, "dir=s" => \$dir, "matchmode=i" => \$matchmode, "filter_abundance=i" => \$filter_abundance, "filter_mz=i" => \$filter_mz, "filter_rt=i" => \$filter_rt, "abundance_min=f" => \$abundance_mincut, "rt_min=f" => \$rt_mincut, "rt_max=f" => \$rt_maxcut, "mz_min=f" => \$mz_mincut, "mz_max=f" => \$mz_maxcut, "id=s" => \$sample_id, "plot=s" => \$plotfile, "pid=i" => \$pid, "sortstyle=i" => \$sortstyle, "slice_per_pt=f" => \$slice_per_pt, "mz_per_pt=f" => \$mz_per_pt, "rt_expand=f" => \$rt_expand, "rt_expand_lo=f" => \$rt_expand_lo, "rt_expand_hi=f" => \$rt_expand_hi, "special_range=i" => \$special_range, "special_range_residues=s" => \$special_range_residues, "special_range_values=s" => \$special_range_values, "rofile=s" => \$rofile, "blueid=s" => \$blueid, "hostname=s" => \$hostname, "ailtype=i" => \$ailtype, "extractstyle=i" => \$extractstyle, "extractbin=f" => \$extractbin);

# Debug Variable Insertion
#$filter_rt = 1;
#$rt_mincut = 6;
#$rt_maxcut = 8;

#
# Add trailing slash to peaks directory if not present
# Not necessary here, but important later when specifying output files
#
if( substr($dir, -1, 1) ne "/" )
{
	$dir = join '', $dir, "/";
}
$rval = mkdir($dir);

#
# Convert $plotfile into a second variable for the raw text dump of the low-res contour data
# Sloppy conversion here just trims off the last 4 characters and replaces with ".txt"
# Should work fine as long as the initial file has a three character extension like ".png"
#
#$plotfile_text = $plotfile;
#substr($plotfile_text, -4, 4) = "";
#$plotfile_text = join '', $plotfile_text, ".txt";
# Not necessary right now - we read the .txt name in directly
$plotfile_text = $plotfile;

@input_array = readfile($inputfile);
@header_array = split ' ', shift(@input_array);
$header_line = join ' ', @header_array;
for($ctr = 0; $ctr <= $#header_array; $ctr++)
{
	$name2col{ $header_array[$ctr] } = $ctr;
	$col2name{ $ctr } = $header_array[$ctr];
}

@zero_array = ();
for($ctr = 0; $ctr <= $#header_array; $ctr++)
{
	push(@zero_array, 0);
}

if($matchmode == 0)
{
	$num_spectra = $#input_array + 1;
	if($num_spectra % 2 != 0)
	{
		roprint("Error: Odd number of lines but looking for feature *pairs*");
		die;
	}
	$num_spectra = $num_spectra / 2;
}
elsif($matchmode == 1 || $matchmode == 2)
{
	$num_spectra = $#input_array + 1;
}

if($filter_abundance == 1)
{
	roprint(sprintf("Abundance filtering enabled - %f -> inf", $abundance_mincut));
}
if($filter_mz == 1)
{
	roprint(sprintf("MZ filtering enabled - %f -> %f", $mz_mincut, $mz_maxcut));
}
if($filter_rt == 1)
{
	roprint(sprintf("RT filtering enabled - %f -> %f", $rt_mincut, $rt_maxcut));
}

STDOUT->autoflush(1);

#
# This is just parsing the .match3 input file
#
# For each spectra that we will excise from the csv (or mzML or whatever), extract all the necessary parameters
# Put these parameters into the @spectra array
#
$num_good_spectra = 0; # Spectra kept after various filtering steps
$num_skip = 0;
$num_skip_rt = 0;
$num_skip_mz = 0;
$num_skip_ab = 0;
$num_skip_sr = 0; # Skipped due to lack of residues from special_range_residues
SPECTRUM: for($ctr = 0; $ctr < $num_spectra; $ctr++)
{
	if($matchmode == 0)
	{
		$ctra = 2*$ctr;
		$ctrb = 2*$ctr + 1;
	}
	elsif($matchmode == 1 || $matchmode == 2)
	{
		$ctra = $ctr;
		$ctrb = $ctr;
	}

	@arraya = split ' ', $input_array[$ctra];
	@arrayb = split ' ', $input_array[$ctrb];
	
	$rt_avg = ( $arraya[ $name2col{"rt"} ] + $arrayb[ $name2col{"rt"} ] ) / 2;
	$rt_min = $rt_avg - $rt_expand - $rt_expand_lo;
	$rt_max = $rt_avg + $rt_expand + $rt_expand_hi;

	$charge = $arraya[ $name2col{"charge"} ];
	$seq = $arraya[ $name2col{"seq"} ];
	$mod = $arraya[ $name2col{"mod"} ];
			
	$mz_n14 = $arraya[ $name2col{"n14mass"} ];
	$mz_n15 = $arrayb[ $name2col{"n15mass"} ];
	#$mz_min = $mz_n14 - $mz_expand;
	#$mz_max = $mz_n15 + $mz_expand + $mz_iso_expand;
	
	# If amino acid labeling is being used, we will calculate a new mz_n15 value
	if($special_range == 0)
	{
		$mz_n15range = $mz_n15;
	}
	else
	{
		$mz_n15range = $mz_n14;
		
		@range_residues = split '_', $special_range_residues;
		@range_values = split '_', $special_range_values;

		# Scan sequence for relevant residues, and increment by appropriate amount
		for($seqctr1 = 0; $seqctr1 < length($seq); $seqctr1++)
		{
			$current_residue = substr($seq, $seqctr1, 1);
			
			for($seqctr2 = 0; $seqctr2 <= $#range_residues; $seqctr2++)
			{
				if($current_residue eq $range_residues[$seqctr2])
				{
					$mz_n15range += $range_values[$seqctr2] / $charge;
				}
			}
		}
	
		# means that no relevant amino acids were detected
		if($mz_n15range == $mz_n14)
		{
			# option to skip this line completely
			$num_skip++;
			$num_skip_sr++;
			next SPECTRUM;
		}
	}
	
	
	# Change these to be more like the old limits that were read in by isodist in read_dino
	# set mzlow [expr ($mw-1.5)/$z]
	# set mzhi [expr (($mw+$nn)*1.005+2.5)/$z]
	$mz_min = $mz_n14 - (1.5 / $charge);
	$mz_max = $mz_n15range*1.005 + (2.5 / $charge);

	# Additional Expansion for the plot
	$rt_min_plot = $rt_avg - 5*$rt_expand - $rt_expand_lo;
	$rt_max_plot = $rt_avg + 5*$rt_expand + $rt_expand_hi;
	#$mz_min_plot = $mz_n14 - 2*$mz_expand;
	#$mz_max_plot = $mz_n15range + 2*$mz_expand + $mz_iso_expand;
	$mz_min_plot = $mz_min - $mz_expand;
	$mz_max_plot = $mz_max + $mz_expand;
	
	$protein = $arraya[ $name2col{"protein"} ];
	$startres = $arraya[ $name2col{"startres"} ];
	$endres = $arraya[ $name2col{"endres"} ];

	$abundancea = $arraya[ $name2col{"abundance"} ];
	$abundanceb = $arrayb[ $name2col{"abundance"} ];

	# Join the sequence and modification fields, and then get rid of the non-modified flags (hyphens "-")
	$seqmod = join '', $seq, $mod;
	$seqmod =~ s/-//g;
	#print "$seq\n";

	# Skip if either peak has abundance too small
	if($filter_abundance == 1)
	{
		if($abundancea < $abundance_mincut || $abundanceb < $abundance_mincut)
		{
			$num_skip++;
			$num_skip_ab++;
			#print BATCH "SKIP abundance\n";
			next SPECTRUM;
		}
	}
	# Skip if ANY part of the RT region for this spectrum is outside the defined bounds
	if($filter_rt == 1)
	{
		if($rt_min < $rt_mincut || $rt_max > $rt_maxcut)
		{
			$num_skip++;
			$num_skip_rt++;
			#print BATCH "SKIP rt\n";
			next SPECTRUM;
		}
	}
	# Skip if n14mass or n15mass outside the bounds - expanded region is okay to be outside
	if($filter_mz == 1)
	{
		if($mz_n14 < $mz_mincut || $mz_n15 > $mz_maxcut)
		{
			$num_skip++;
			$num_skip_mz++;
			#print BATCH "SKIP mz\n";
			next SPECTRUM;
		}
	}

	$label = join '_', $protein, $startres, $endres, $sample_id, $mz_n14, $rt_avg;

	while(defined($labelhash{$label}))
	{
		$label = join '', $label, "a";
	}
	
	$labelhash{$label} = 1;

	#$spectra[$num_good_spectra]{"label"} = join '_', $protein, $startres, $endres, $sample_id, $mz_n14, $rt_avg;
	$spectra[$num_good_spectra]{"label"} = $label;
	$spectra[$num_good_spectra]{"rt_avg"} = $rt_avg;
	$spectra[$num_good_spectra]{"rt_min"} = $rt_min;
	$spectra[$num_good_spectra]{"rt_max"} = $rt_max;
	$spectra[$num_good_spectra]{"mz_min"} = $mz_min;
	$spectra[$num_good_spectra]{"mz_max"} = $mz_max;
	$spectra[$num_good_spectra]{"seqmod"} = $seqmod;
	$spectra[$num_good_spectra]{"charge"} = $charge;
	$spectra[$num_good_spectra]{"slices"} = 0;
	$spectra[$num_good_spectra]{"index"} = $ctr; # Reference the @spectra array to the original input file

	$spectra[$num_good_spectra]{"mz_min_plot"} = $mz_min_plot;
	$spectra[$num_good_spectra]{"mz_max_plot"} = $mz_max_plot;
	$spectra[$num_good_spectra]{"rt_min_plot"} = $rt_min_plot;
	$spectra[$num_good_spectra]{"rt_max_plot"} = $rt_max_plot;
	$spectra[$num_good_spectra]{"slices_plot"} = 0;
	
	# Insert new stuff keeping track of previous information

	#print "Spec: $spectra[$num_good_spectra]{\"rt_min\"} -> $spectra[$num_good_spectra]{\"rt_max\"}\n";
	#print "Plot: $spectra[$num_good_spectra]{\"rt_min_plot\"} -> $spectra[$num_good_spectra]{\"rt_max_plot\"}\n";
	$num_good_spectra++;
}

#print "Excising $num_good_spectra minispectra\n";
#print "Skipping $num_skip minispectra\n";
#print "  $num_skip_ab -> Abundance Filtering\n";
#print "  $num_skip_rt -> RT Filtering\n";
#print "  $num_skip_mz -> M/Z Filtering\n";
roprint(sprintf("Excising %d minispectra", $num_good_spectra));
roprint(sprintf("Skipping %d minispectra", $num_skip));
roprint(sprintf("  %d -> Abundance Filtering", $num_skip_ab));
roprint(sprintf("  %d -> RT Filtering", $num_skip_rt));
roprint(sprintf("  %d -> M/Z Filtering", $num_skip_mz));
roprint(sprintf("  %d -> Residue Filtering (no residues from Set Special Range found)", $num_skip_sr));
$t_elapsed = (time() - $t_start)/60;
roprint(sprintf("Done Initial Spectral Sweep - Elapsed Time: %.1f min", $t_elapsed));
STDOUT->autoflush(1);

#for($ctr = 0; $ctr < $num_good_spectra; $ctr++)
#{
#	print "$spectra[$ctr]{'label'}\n";
#}

#
# Open the text file for low resolution contour data
#
open IMAGETEXT, ">$plotfile_text" or die "Error: Can't open $plotfile_text\n";

# Split on periods to determine file type from file extension
@ext_array = split /\./, $csvfile;
$filetype1 = lc($ext_array[$#ext_array]); # this is the last array element, so the actual file extension, in lowercase
$filetype2 = lc(join '.', $ext_array[$#ext_array-1], $ext_array[$#ext_array]);

$gz = 0;
if(lc($ext_array[$#ext_array] eq "gz"))
{
    $gz = 1;
    pop(@ext_array);
    $csvfile_nogz = join '.', @ext_array;
}

#
# This is where the type of data to be processed is determined
#
if($filetype1 eq "csv")
{
	do_csv();
	$fileformat eq "csv";
}
elsif($filetype2 eq "mzdata.xml")
{
	do_mzdata();
	$fileformat = "mzdata";
}
elsif($filetype1 eq "mzml" || $filetype2 eq "mzml.gz")
{
	# The output of msconvert does not report these accurately
	# So just automatically peg these at 400 -> 2000 based on typical usage
	# ssc101911 changed to 200 -> 2000
	$mass_min = 200;
	$mass_max = 2000;

	do_mzml();
	$fileformat = "mzml";
}
elsif($filetype1 eq "prn")
{
	do_prn();
	$fileformat = "prn";	
}
else
{
	roprint(sprintf("Error, unrecognized file format for %s", $csvfile));
}

#if(substr($csvfile, -3, 3) eq "csv")
#{
#	$fileformat = "csv";
#}
#elsif($csvfile =~ /mzdata/)
#{
#	$fileformat = "mzdata";
#}
#else
#{
#	$fileformat = "";
#}

#if($fileformat eq "csv")

#elsif($fileformat eq "mzdata")

#else
#{
#	roprint(sprintf("Error, unrecognized file format for %s", $csvfile));
	#STDOUT->autoflush(1);
#}

close IMAGETEXT;

# 2/25/2011 - Update precalculation range for all types of data


# What is @tolarray?
# The problem is, from scan to scan the m/z values used to get intensity readings do not exactly line up
# *    *    *    *    *    *
#  *    *    *    *    *    *
#   *    *    *    *    *    *
# You can easily see that these cluster into little groups, but we need an algorith to determine this, and yield the following:
#  *    *    *    *    *    *
# So what we really do is...
# 1) Sort the points in terms of m/z, and plop them in an array
# ***  ***  ***  ***  ***  ***
# 2) Then we go through the array, and see how close sequential points are
#    -> If these points are within some value, they are lumped together as a single point
#    -> In the old days this value was called "mz_bin_tol"
#    -> Now it's called @tolarray, and has a different value depending on m/z (the spacing between readings varies with m/z in most spectrometers we use)
# 3) There are two possible problems:
#    a) mz_bin_tol is too large -> this means that too many points are lumped together into one.  This is easily detected because that means num_points > num_slices.  So we just reduce mz_bin_tol by a small amount, and try it all again.  This overestimation of mz_bin_tol is more-or-less desirable, and is an expected result
#    b) mz_bin_tol is too small -> this means that not enough points are lumped together into one.  This should be easily detectable because num_points < num_slices.  Unfortunately, it's not that simple because there could be some sparse missing readings (Anne used to see this in the ESI-TOF data), or there simply could be a bunch of values not present (pretty common in our orbitrap data which seems to omit near-baseline values entirely).  So the solution isn't a real solution, but in practice we want to OVERESTIMATE mz_bin_tol (see "a"), which can then be easily handled

#
# Right now @tolarray is calculated based on data type (.csv, mzML) because data type has generally corresponded to machine
# However, we are considering .wiff -> .d -> mzML, which would mean that mzML could correspond to both ESI-TOF and qTOF
# It might be better to check on the peak list type ($ailtype) as a safer check
# 

#
# Precalculate mz_bin_tol for the m/z range 0 -> 1600
# mz_bin_tol = 0.0007 * sqrt(m/z)
# When m/z = 100, mz_bin_tol = 0.007
# This was empirically determined by Anne and later confirmed by Sykes to be a better starting point than 0.005
# These values generally apply to the ESI-TOF instrument
#
for($ctr = 0; $ctr <= 2000; $ctr++)
{
	$tolval = 0.0007 * sqrt($ctr);
	$tolarray[$ctr] = $tolval;
}
	
# This mz_bin_tol doesn't work very well for the QTOF mzdata or mzml files
# SHIFT reduces mz_bin_tol
# SHIFT for one of Andrea's datasets: 2.48 X 10^-6 + 0.00163 (relative to mz_n15)
# SHIFT for Zahra 2.45 X 10^-6(mz) + 0.00164 (relative to mz_n15)
# Shift for Zahra 2.48 X 10^-6(mz) + 0.00164
# Because SHIFT is calculated relative to the csv parameters, this adjustment is also done relative to the csv parameters
if($fileformat eq "mzdata" || $fileformat eq "mzml")
{
	for($ctr = 0; $ctr <= 2000; $ctr++)
	{
		$tolarray[$ctr] -= 0.0000025*$ctr + 0.00164;
	}
}

# 2/25/11
# Orbitrap/Sequest Data
# Here we make the selection to calculate @tolarray based on the Sequest data type
# This isn't great, as really the issue is the particular orbitrap machine, and not the fact that we are using Sequest, but it works for now
# This overwrites @tolarray pre-calculated above
if($ailtype == 5)
{
	my $tol_A = 0.000107929;
	my $tol_B = 1.51567;
	my $tol_C = 3000000;
	
	for($ctr = 0; $ctr <= 2000; $ctr++)
	{
		$tolarray[$ctr] = $tol_A + (($ctr**$tol_B)/$tol_C);
	}
}

# Precalculation of mz_bin_tol for the Triple-TOF machine
# f(x) = (x**A)/B - 0.0001
if($blueid eq "sykes")
{
    my $tol_A = 0.497585;
    my $tol_B = 7023.89;
    
    for($ctr = 0; $ctr <= 3000; $ctr++)
    {
        $tolarray[$ctr] = (($ctr**$tol_A)/$tol_B) - 0.0001;
    }
}

#@batcharray = ();
#
# Go through each spectrum, write out the x/y values and generate the plot
#
open TXTLIST, ">txtlist.$pid" or die "Can't open txtlist.$pid\n";
open CONLIST, ">conlist.$pid" or die "Can't open conlist.$pid\n";
#open NEWCONLIST, ">newconlist.$pid" or die "Can't open newconlist.$pid\n";

#
# For each spectrum try to do the summing (via @tolarray or binning)
#
SPEC: for($spec_ctr = 0; $spec_ctr < $num_good_spectra; $spec_ctr++)
{
	#if($spectra3D) # Special case for Jamie's 3D fitting stuff
    if($extractstyle == 2) # Special case for Jamie's 3D fitting stuff
	{
		@sorted_values = sort bymz @{ $spectra_values[$spec_ctr] };

		if($#sorted_values < 0)
		{
			roprint(sprintf("Skipping %s%s", $dir, $spectra[$spec_ctr]{'label'}));
			#STDOUT->autoflush(1);
			#$batcharray[$spec_ctr] = "SKIP";
			$spectra[$spec_ctr]{"keep"} = 0;
			next SPEC;
		}
		$spectra[$spec_ctr]{"keep"} = 1;
	
		print TXTLIST "$dir$spectra[$spec_ctr]{'label'}.txt\n";
		open OUT, ">$dir$spectra[$spec_ctr]{'label'}.txt" or die "Can't open $dir$spectra[$spec_ctr]{'label'}.txt\n";
		for($ctr3D = 0; $ctr3D <= $#{$spectra_values[$spec_ctr]}; $ctr3D++)
		{
			print OUT "$spectra_values[$spec_ctr][$ctr3D]\n";
		}
		close OUT;
	}
	else # Just do a normal summed minispectra
	{
		@sorted_values = sort bymz @{ $spectra_values[$spec_ctr] };
		#print "Length: $#sorted_values\n";
		#print "Length: $#{ $spectra_values[$spec_ctr] }\n";
		#print "Min->Max: $spectra[$spec_ctr]{\"mz_min\"} $spectra[$spec_ctr]{\"mz_max\"}\n";
			
		#
		# Sometimes we end up with an empty spectrum (ie pair outside the range of the csv)
		# So skip it, and make a notation in the {"keep"} key
		#
		if($#sorted_values < 0)
		{
			roprint(sprintf("Skipping %s%s", $dir, $spectra[$spec_ctr]{'label'}));
			#STDOUT->autoflush(1);
			#$batcharray[$spec_ctr] = "SKIP";
			$spectra[$spec_ctr]{"keep"} = 0;
			next SPEC;
		}
		$spectra[$spec_ctr]{"keep"} = 1;
		
		for($ctr = 0; $ctr <= $#sorted_values; $ctr++)
		{
			@array = split ' ', $sorted_values[$ctr];
			$mz_sort[$ctr] = $array[0];
			$int_sort[$ctr] = $array[1];
		}

		print TXTLIST "$dir$spectra[$spec_ctr]{'label'}.txt\n";

        # 6/9/11 Triple-Tof data is causing problems with the mz_bin_tol method.  Basically, the data spacing is both erratic and very tight.  So rather than incrementing mz_bin_tol until you get num_points = num_slices, you end up with a variable num_points that is something < num_slices, and you get jaggies
        # On the other hand the data is quite high resolution and there are a lot of points.  Let's just do some simple binning of the data in arbitrary sized bins.  Based on experimentation with a few tests 0.05 m/z looks like a good bin size to start with
        if($extractstyle == 1)
        {
            $bin = $extractbin;
            $min_mz = $mz_sort[0];
            
            @bin_data = ();
            @bin_data_mz = ();
            @bin_data_count = ();
            
            for($bin_ctr = 0; $bin_ctr <= $#sorted_values; $bin_ctr++)
            {
                $bin_index1 = ($mz_sort[$bin_ctr] - $min_mz) / $bin;
                $bin_index2 = sprintf("%.0f", $bin_index1);
                $bin_data[$bin_index2] += $int_sort[$bin_ctr];
                $bin_data_mz[$bin_index2] += $mz_sort[$bin_ctr];
                $bin_data_count[$bin_index2] += 1;
            }
            
            open OUT, ">$dir$spectra[$spec_ctr]{'label'}.txt" or die "Can't open $dir$spectra[$spec_ctr]{'label'}.txt\n";
            for($bin_ctr = 0; $bin_ctr <= $#bin_data; $bin_ctr++)
            {
                if(!defined($bin_data[$bin_ctr]))
                {
                    #$bin_data[$bin_ctr] = 0;
                    next;
                }
                
                #$bin_mz = $min_mz + ($bin_ctr * $bin);
                #$bin_int = $bin_data[$bin_ctr] / $spectra[$spec_ctr]{"slices"};
                
                $bin_int = $bin_data[$bin_ctr] / $bin_data_count[$bin_ctr];
                $bin_mz = $bin_data_mz[$bin_ctr] / $bin_data_count[$bin_ctr];
                
                print OUT "$bin_mz $bin_int $bin_data_count[$bin_ctr]\n";
            }
            close OUT;
        }
        # Otherwise revert back to the mz_bin_tol method ($extractstyle == 0)
        elsif($extractstyle == 0)
        {
            $shift = 0.00000;
            
            SHIFT: while(1)
            {
                #if($shift == 0.00000)
                #{
                    #print "Writing $spectra[$spec_ctr]{'label'}\n";
                    #STDOUT->autoflush(1);
                #	print TXTLIST "$dir$spectra[$spec_ctr]{'label'}.txt\n";
                #}
                
                #$batcharray[$spec_ctr] = "$spectra[$spec_ctr]{'seq'} $spectra[$spect_ctr]{'charge'} $dir$spectra[$spec_ctr]{'label'}.txt";
                open OUT, ">$dir$spectra[$spec_ctr]{'label'}.txt" or die "Can't open $dir$spectra[$spec_ctr]{'label'}.txt\n";
                
                #
                # Group points that are <= $mz_bin_tol away from each other
                #
                $flag = 0;
                $num = 0;
                $int_sum = 0;
                $mz_sum = 0;
                for($ctr2 = 0; $ctr2 <= $#sorted_values; $ctr2++)
                {
                    if($ctr2 == 0)
                    {
                        $int_sum += $int_sort[$ctr2];
                        $mz_sum += $mz_sort[$ctr2];
                        $num = 1;
                    }
                    else
                    {
                        $diff = $mz_sort[$ctr2] - $last_mz;
                        if($diff <= ($tolarray[$mz_sort[$ctr2]] - $shift))
                        {
                            $int_sum += $int_sort[$ctr2];
                            $mz_sum += $mz_sort[$ctr2];
                            $num++;
                        }
                        else
                        {
                            $mz_out = $mz_sum / $num;
                            $int_out = $int_sum / $num;
                            $int_out = $int_sum / $spectra[$spec_ctr]{"slices"};
                            #if($blueid eq "sykes")
                            #{
							#   print OUT "$mz_out $int_out $num $spectra[$spec_ctr]{'slices'}\n";
                            #}
                            #else
                            #{
                                print OUT "$mz_out $int_out\n";
                            #}

                            #
                            # This checks to see if we have accumulated more data points than spectra
                            # If so, we have lumped too many data points together, so increment $shift and try again
                            #
                            if($num > $spectra[$spec_ctr]{"slices"})
                            {
                                #roprint(sprintf("ID of Slice Adjustment: %d", $ctr2));
                                if($shift == 0.00000)
                                {
                                    roprint("Adjusting mz_bin_tol");
                                }
                                $shift = $shift + 0.00002;
                                #printf("%.5f\n", $shift);
                                next SHIFT;
                            }
                            
                            $mz_sum = $mz_sort[$ctr2];
                            $int_sum = $int_sort[$ctr2];
                            $num = 1;
                        }

                    }
                    
                    $last_mz = $mz_sort[$ctr2];
                }
                
                # Do it one last time to handle the final cluster
                $mz_out = $mz_sum / $num;
                $int_out = $int_sum / $num;
                $int_out = $int_sum / $spectra[$spec_ctr]{"slices"};
                #if($blueid eq "sykes")
                #{
                #    print OUT "$mz_out $int_out $num $spectra[$spec_ctr]{'slices'}\n";
                #}
                #else
                #{
                    print OUT "$mz_out $int_out\n";
                #}
                
                if($num > $spectra[$spec_ctr]{"slices"})
                {
                    if($shift == 0.00000)
                    {
                        roprint("Adjusting mz_bin_tol");
                        #STDOUT->autoflush(1);
                    }
                    $shift = $shift + 0.00002;
                    #printf("%.5f\n", $shift);
                    next SHIFT;
                }
                
                close OUT;
                
                last SHIFT;
            }
            
            if($shift > 0.0)
            #if($shift > 0.0 || $blueid eq "sykes")
            {
                roprint(sprintf("SHIFT: %.5f", $shift));
            }
        
        }
        else
        {
            roprint("ERROR: extractstyle |$extractstyle| unrecognized");
            die;
        }
	}
	
	#
	# Write out the high-res contour plot for the immediate area
	# 
	#
	open CONTOUR, ">$dir$spectra[$spec_ctr]{'label'}.con" or die "Can't open $dir$spectra[$spec_ctr]{'label'}.con\n";

	print CONLIST "$dir$spectra[$spec_ctr]{'label'}.con\n";

	$zthresh_mini = 25000;
	$ysize = $spectra[$spec_ctr]{"slices_plot"};
	$xsize = 100*($spectra[$spec_ctr]{"mz_max_plot"} - $spectra[$spec_ctr]{"mz_min_plot"});
	#@plot_values = @{ $spectra_values_plot[$spec_ctr] };
	#@plot_values = split ' ', $spectra_values_plot[$spec_ctr];
	#$num_plot_val = ($#plot_values + 1) / 3;
	#$num_plot_val = $#plot_values + 1;

	$out_mz_min_plot = $spectra[$spec_ctr]{'mz_min_plot'} * $spectra[$spec_ctr]{'charge'};
	$out_mz_max_plot = $spectra[$spec_ctr]{'mz_max_plot'} * $spectra[$spec_ctr]{'charge'};

	print CONTOUR "$xsize,$ysize,$out_mz_min_plot,$out_mz_max_plot\n";
	
	# Scan the entire thing, to scale $zthresh_mini dynamically
	for($slice_ctr = 0; $slice_ctr < $spectra[$spec_ctr]{"slices_plot"}; $slice_ctr++)
	{
		@plot_values = split ' ', $spectra_values_plot[$spec_ctr][$slice_ctr];
		# Searching only over spectra_values will scale to within only the actual spectrum
		# This might make contour plots a little better overall (2/8/2009)
		# Ummm, nope.  Instead of all white ones, we get all black ones
		#@plot_values = split ' ', $spectra_values[$spec_ctr][$slice_ctr];
		$num_plot_val = ($#plot_values + 1) / 2;

		for($plot_ctr = 0; $plot_ctr < $num_plot_val; $plot_ctr++)
		{
			$z = $plot_values[2*$plot_ctr + 1];

			if($plot_ctr == 0)
			{
				$zthresh_mini = $z;
			}
			elsif($z > $zthresh_mini)
			{
				$zthresh_mini = $z;
			}
		}
	}

	# 12/10/2010 Sequest data, some yield zthresh_mini = 0
	# Empty spectrum?  Set this to avoid divide by zero errors
	if($zthresh_mini == 0)
	{
		$zthresh_mini = 1;
	}

	for($slice_ctr = 0; $slice_ctr < $spectra[$spec_ctr]{"slices_plot"}; $slice_ctr++)
	{
		$y = $slice_ctr;
		$y = ($ysize - 1) - $y; # Use $ysize - 1 as the values are indexed from 0

		@plot_values = split ' ', $spectra_values_plot[$spec_ctr][$slice_ctr];
		$num_plot_val = ($#plot_values + 1) / 2;

		#print "num_plot_val: $num_plot_val\n";

		for($plot_ctr = 0; $plot_ctr < $num_plot_val; $plot_ctr++)
		{
			$x = $plot_values[2*$plot_ctr];
			$z = $plot_values[2*$plot_ctr + 1];

			$x = 100*($x - $spectra[$spec_ctr]{"mz_min_plot"});

			#print "(x,y,z) $x $y $z\n";
			
			$color_fraction = $z / $zthresh_mini;
			if($color_fraction > 1)
			{
				$color_fraction = 1;
			}
			
			print CONTOUR "$x,$y,$color_fraction\n";
		}
	
	}
	
	$rec_x1 = 100*($spectra[$spec_ctr]{"mz_min"} - $spectra[$spec_ctr]{"mz_min_plot"});
	$rec_x2 = 100*($spectra[$spec_ctr]{"mz_max"} - $spectra[$spec_ctr]{"mz_min_plot"});

	$rec_y1 = ($ysize - 1) - $spectra[$spec_ctr]{"slices_before"}; # ($ysize-1) is the bottom of the plot - work up by slices_before
	$rec_y2 = $spectra[$spec_ctr]{"slices_after"};	

	print CONTOUR "$rec_x1,$rec_y1,$rec_x2,$rec_y2\n";	

	close CONTOUR;
	
	# New style contour data for on the fly contour plot generation and exploration
	
	open NEWCONTOUR, ">$dir$spectra[$spec_ctr]{'label'}.newcon" or die "Can't open $dir$spectra[$spec_ctr]{'label'}.newcon\n";

	#print NEWCONLIST "$dir$spectra[$spec_ctr]{'label'}.newcon\n";

	$nc_mz_min = $spectra[$spec_ctr]{"mz_min"};
	$nc_mz_max = $spectra[$spec_ctr]{"mz_max"};
	$nc_rt_min = $spectra[$spec_ctr]{"rt_min"};
	$nc_rt_max = $spectra[$spec_ctr]{"rt_max"};
	$nc_mz_min_plot = $spectra[$spec_ctr]{"mz_min_plot"};
	$nc_mz_max_plot = $spectra[$spec_ctr]{"mz_max_plot"};
	$nc_rt_min_plot = $spectra[$spec_ctr]{"rt_min_plot"};
	$nc_rt_max_plot = $spectra[$spec_ctr]{"rt_max_plot"};

	$mzshift = sprintf("%.1f", $nc_mz_min_plot);
	$mzshift = 10*$mzshift;
	$mzshift_max = sprintf("%.1f", $nc_mz_max_plot);
	$mzshift_max = 10*$mzshift_max;

	print NEWCONTOUR "$nc_mz_min $nc_mz_max $nc_rt_min $nc_rt_max\n";
	print NEWCONTOUR "$nc_mz_min_plot $nc_mz_max_plot $nc_rt_min_plot $nc_rt_max_plot\n";

	@ncdata = ();

	for($slice_ctr = 0; $slice_ctr < $spectra[$spec_ctr]{"slices_plot"}; $slice_ctr++)
	{
		$y = $slice_ctr;

		@plot_values = split ' ', $spectra_values_plot[$spec_ctr][$slice_ctr];
		$num_plot_val = ($#plot_values + 1) / 2;

		for($plot_ctr = 0; $plot_ctr < $num_plot_val; $plot_ctr++)
		{
			$x = sprintf("%.1f", $plot_values[2*$plot_ctr]);
			$x = 10*$x - $mzshift;
			
			$z = $plot_values[2*$plot_ctr + 1];

			$ncdata[$x][$y] += $z;
		}
	
	}
	
	$nc_mzcount = 0;
	@mzvals = ();
	for($nc_ctr = $mzshift; $nc_ctr <= $mzshift_max; $nc_ctr++)
	{
		push(@mzvals, $nc_ctr/10);
		$nc_mzcount++;
	}
	
	printf(NEWCONTOUR "%d\n", $nc_mzcount); # Number of m/z values
	printf(NEWCONTOUR "%s\n", join ' ', @mzvals); # Actual m/z values
	
	printf(NEWCONTOUR "%d\n", $spectra[$spec_ctr]{"slices_plot"}); # Number of rt values
	printf(NEWCONTOUR "%s\n", join ' ', @{$spectra[$spec_ctr]{"rt_vals"}});
	
	for($nc_ctr1 = 0; $nc_ctr1 < $spectra[$spec_ctr]{"slices_plot"}; $nc_ctr1++)
	{
		for($nc_ctr2 = 0; $nc_ctr2 < $nc_mzcount; $nc_ctr2++)
		{
			if(defined($ncdata[$nc_ctr2][$nc_ctr1]))
			{
				print NEWCONTOUR "$ncdata[$nc_ctr2][$nc_ctr1] ";
			}
			else
			{
				print NEWCONTOUR "0 ";
			}
		}
		print NEWCONTOUR "\n";
	}
	
	close NEWCONTOUR;

	$t_elapsed = (time() - $t_start)/60;
	#roprint(sprintf("Spectrum: %d Time: %.1f min", $spec_ctr, $t_elapsed));
}

close TXTLIST;
close CONLIST;
#close NEWCONLIST;

#
# NOTE: The following method of collecting the two input lines into one consolidated line is pretty hideous
#

#
# Generate a new array of hashes with consolidated information
#
# ailid featureid multimer ailcharge ailiontype rt mz abundance sampleid digestionid digestpepid n14mass n15mass protein startres endres iontype charge missed seq mod n14prox n15prox matchtype ppm
#
# In the N14/N15 case, we need to keep both copies of:
# ailid featureid multimer ailiontype rt mz abundance ppm
#
# Things that should be the same for both N14 and N15 lines, so only need to keep one
# ailcharge sampleid digestionid digestpepid n14mass n15mass protein startres endres iontype charge missed seq mod 
#
# Other
# n14prox n15prox -> keep only the relevant one in each case
# matchtype -> information that is kind of encoded in the $matchmode variable, and is constant for all spectra
#

# 6/11/09 New (Semi-)Automatic Generation of @output_labels
@output_labels = ();
# First, add in the mandatory 
@output_labels1 = ("ailid_n14", "ailid_n15", "rt_n14", "rt_n15", "mz_n14", "mz_n15", "abundance_n14", "abundance_n15", "ailcharge");
# Now go through and fetch out all the "pks_" stuff
@output_labels2 = ();
@output_labels2a = ();
foreach $name (keys %name2col)
{
	if(substr($name, 0, 4) eq "pks_")
	{
		$pks_n14 = join '', $name, "_n14";
		$pks_n15 = join '', $name, "_n15";
		push(@output_labels2, $pks_n14);
		push(@output_labels2, $pks_n15);
		push(@output_labels2a, $name);
	}
}
# Now, add on the rest
@output_labels3 = ( "matchmode", "ppm_n14", "ppm_n15", "sampleid", "digestionid", "digestpepid", "n14mass", "n15mass", "protein", "startres", "endres", "iontype", "charge", "missed", "seq", "mod", "n14prox", "n15prox", "seqmod", "file");

@output_labels = (@output_labels1, @output_labels2, @output_labels3);

$keep_ctr = 0;
for($ctr = 0; $ctr < $num_good_spectra; $ctr++)
{
	if($spectra[$ctr]{"keep"} == 0)
	{
		next;
	}
		
	if($matchmode == 0)
	{
		$index_n14 = 2*$spectra[$ctr]{"index"};
		$index_n15 = 2*$spectra[$ctr]{"index"} + 1;
		
		@n14array = split ' ', $input_array[$index_n14];
		@n15array = split ' ', $input_array[$index_n15];
		@goodarray = @n14array;
	}
	elsif($matchmode == 1)
	{
		$index_n14 = $spectra[$ctr]{"index"};
		
		@n14array = split ' ', $input_array[$index_n14];
		@n15array = @zero_array;
		@goodarray = @n14array;
	}
	elsif($matchmode == 2)
	{
		$index_n15 = $spectra[$ctr]{"index"};

		@n14array = @zero_array;
		@n15array = split ' ', $input_array[$index_n15];
		@goodarray = @n15array;
	}


	##########
	# These items originate from the peak file (Mandatory)
	# ailid is actually generated from the peak file line numbers in step 1
	##########
	$output_data[$keep_ctr]{"ailid_n14"} = $n14array[ $name2col{"ailid"} ];
	$output_data[$keep_ctr]{"ailid_n15"} = $n15array[ $name2col{"ailid"} ];

	$output_data[$keep_ctr]{"rt_n14"} = $n14array[ $name2col{"rt"} ];
	$output_data[$keep_ctr]{"rt_n15"} = $n15array[ $name2col{"rt"} ];

	$output_data[$keep_ctr]{"mz_n14"} = $n14array[ $name2col{"mz"} ];
	$output_data[$keep_ctr]{"mz_n15"} = $n15array[ $name2col{"mz"} ];

	$output_data[$keep_ctr]{"abundance_n14"} = $n14array[ $name2col{"abundance"} ];
	$output_data[$keep_ctr]{"abundance_n15"} = $n15array[ $name2col{"abundance"} ];

	$output_data[$keep_ctr]{"ailcharge"} = $goodarray[ $name2col{"ailcharge"} ];
	##########
	##########

	##########
	# These items originate from the peak file (Optional)
	# featureid/multimer/ailiontype are part of the AIL file format
	# Need a way to automatically detect and incorporate these
	##########
	#$output_data[$keep_ctr]{"featureid_n14"} = $n14array[ $name2col{"featureid"} ];
	#$output_data[$keep_ctr]{"featureid_n15"} = $n15array[ $name2col{"featureid"} ];
	
	#$output_data[$keep_ctr]{"multimer_n14"} = $n14array[ $name2col{"multimer"} ];
	#$output_data[$keep_ctr]{"multimer_n15"} = $n15array[ $name2col{"multimer"} ];
	
	#$output_data[$keep_ctr]{"ailiontype_n14"} = $n14array[ $name2col{"ailiontype"} ];
	#$output_data[$keep_ctr]{"ailiontype_n15"} = $n15array[ $name2col{"ailiontype"} ];

	for($olctr = 0; $olctr <= $#output_labels2a; $olctr++)
	{
		$label_n14 = join '', $output_labels2a[$olctr], "_n14";
		$label_n15 = join '', $output_labels2a[$olctr], "_n15";
		$output_data[$keep_ctr]{$label_n14} = $n14array[$name2col{$output_labels2a[$olctr]}];
		$output_data[$keep_ctr]{$label_n15} = $n15array[$name2col{$output_labels2a[$olctr]}];
	}
	
	##########
	# Internally Generated
	##########
	$output_data[$keep_ctr]{"matchmode"} = $matchmode;
	
	$output_data[$keep_ctr]{"ppm_n14"} = $n14array[ $name2col{"ppm"} ];
	$output_data[$keep_ctr]{"ppm_n15"} = $n15array[ $name2col{"ppm"} ];

	$output_data[$keep_ctr]{"sampleid"} = $goodarray[ $name2col{"sampleid"} ];
	$output_data[$keep_ctr]{"digestionid"} = $goodarray[ $name2col{"digestionid"} ];
	$output_data[$keep_ctr]{"digestpepid"} = $goodarray[ $name2col{"digestpepid"} ];
	$output_data[$keep_ctr]{"n14mass"} = $goodarray[ $name2col{"n14mass"} ];
	$output_data[$keep_ctr]{"n15mass"} = $goodarray[ $name2col{"n15mass"} ];
	$output_data[$keep_ctr]{"protein"} = $goodarray[ $name2col{"protein"} ];
	$output_data[$keep_ctr]{"startres"} = $goodarray[ $name2col{"startres"} ];
	$output_data[$keep_ctr]{"endres"} = $goodarray[ $name2col{"endres"} ];
	$output_data[$keep_ctr]{"iontype"} = $goodarray[ $name2col{"iontype"} ];
	$output_data[$keep_ctr]{"charge"} = $goodarray[ $name2col{"charge"} ];
	$output_data[$keep_ctr]{"missed"} = $goodarray[ $name2col{"missed"} ];
	$output_data[$keep_ctr]{"seq"} = $goodarray[ $name2col{"seq"} ];
	$output_data[$keep_ctr]{"mod"} = $goodarray[ $name2col{"mod"} ];

	$output_data[$keep_ctr]{"n14prox"} = $n14array[ $name2col{"n14prox"} ];
	$output_data[$keep_ctr]{"n15prox"} = $n15array[ $name2col{"n15prox"} ];

	$output_data[$keep_ctr]{"seqmod"} = $spectra[$ctr]{"seqmod"};
	$output_data[$keep_ctr]{"file"} = join '', $dir, $spectra[$ctr]{"label"}, ".txt";

	$keep_ctr++;
}

#
# Note: multimer, ailiontype and ailcharge are probably all useless
#
# Note: This is now defined above
#@output_labels = ( "matchmode", "ailid_n14", "ailid_n15", "featureid_n14", "featureid_n15", "multimer_n14", "multimer_n15", "ailiontype_n14", "ailiontype_n15", "rt_n14", "rt_n15", "mz_n14", "mz_n15", "abundance_n14", "abundance_n15", "ppm_n14", "ppm_n15", "ailcharge", "sampleid", "digestionid", "digestpepid", "n14mass", "n15mass", "protein", "startres", "endres", "iontype", "charge", "missed", "seq", "mod", "n14prox", "n15prox", "seqmod", "file");

if($sortstyle == 0)
{
	#@sorted_data = sort bymz2 @output_data;
	@output_data = sort bymz2 @output_data;
}
elsif($sortstyle == 1)
{
	#@sorted_data = sort byrt @output_data;
	@output_data = sort byrt @output_data;
}
elsif($sortstyle == 2)
{
	#@sorted_data = sort bypep @output_data;
	@output_data = sort bypep @output_data;
}

open BATCH, ">$outputfile" or die "Error: Can't open $outputfile\n";
print BATCH "isoid ";
for($ctr2 = 0; $ctr2 <= $#output_labels; $ctr2++)
{
	print BATCH "$output_labels[$ctr2] ";
	$output_name2col{ $output_labels[$ctr2] } = $ctr2;
}
print BATCH "\n";

for($ctr = 0; $ctr <= $#output_data; $ctr++)
{
	#print "Go: $ctr\n";
	print BATCH "$ctr ";
	for($ctr2 = 0; $ctr2 <= $#output_labels; $ctr2++)
	{
		print BATCH "$output_data[$ctr]{ $output_labels[$ctr2] } ";
	}
	print BATCH "\n";
}
close BATCH;


# Processing Data Files
# 1) Need to retrieve retention time in minutes ($rt)
# 2) Need to establish two arrays:
#      a) @mz_array (m/z values)
#      b) @int_array (intensity values)
# 3) Need to allow for retention time filtering
# 4) Need to define $mass_max and $mass_min if not done elsewhere
#
# Once $rt, @mz_array and @int_array are setup, process_slice() is called

sub do_csv
{
	open CSV, "$csvfile" or die "Error: Can't open $csvfile\n";
	roprint("Reading csv file");
	#STDOUT->autoflush(1);
	$slice_ctr = 0;
	$flag = 0;
	CSVLINE: while(defined($input=<CSV>))
	{
		@array = split /\,/, $input;
		
		#
		# Skip a bunch of the header information
		#
		if($array[0] eq "mass range")
		{
			$mass_min = $array[1];
			$mass_max = $array[2];
		}
		if($array[0] eq "time range")
		{
		
		}
		if($array[0] eq "number of spectra")
		{
			$nspectra = $array[1];
		}
		if(substr($input, 0, 9) eq "[spectra]")
		{
			#print "MATCH\nMATCH\nMATCH\nMATCH\n";
			$flag = 1;
			next;
		}
		if(substr($input, 0, 10) eq "[\spectra]")
		{
			$flag = 0;
			next;
		}
		
		$plot_xsize = ($mass_max - $mass_min) / 0.2;
		$plot_ysize = $nspectra;
		
		if($flag == 0)
		{
			next CSVLINE;
		}

		$rt = shift(@array);
		
		#
		# Filter by RT
		# Not strictly necessary as we have already filtered peak spectra by RT, but this will speed things up
		#
		if($filter_rt == 1)
		{
			if($rt < $rt_mincut)
			{
				if($slice_ctr % 100 == 0)
				{
					#print "Slice: $slice_ctr  RT: $rt\n";
					$t_elapsed = (time() - $t_start)/60;
					roprint(sprintf("Slice: %d RT: %f Elapsed Time: %.1f min", $slice_ctr, $rt, $t_elapsed));
					#STDOUT->autoflush(1);
				}

				$slice_ctr++;
				next CSVLINE;
			}
			elsif($rt > $rt_maxcut)
			{
				if($slice_ctr % 100 == 0)
				{
					#print "Slice: $slice_ctr  RT: $rt\n";
					$t_elapsed = (time() - $t_start)/60;
					roprint(sprintf("Slice: %d RT: %f Elapsed Time: %.1f min", $slice_ctr, $rt, $t_elapsed));
					#STDOUT->autoflush(1);
				}
				
				$slice_ctr++;
				last CSVLINE;
			}
		}
		
		#
		# Get rid of some miscellaneous padding on the spectrum
		#
		for($ctr = 0; $ctr < 6; $ctr++)
		{
			shift(@array);
		}
		$last = pop(@array);
		
		if( ($#array + 1) % 2 != 0)
		{
			roprint(sprintf("Warning: Unmatched x,y values in Slice %d", $slice_ctr));
			#STDOUT->autoflush(1);
		}

		$num_xy = ($#array + 1) / 2;

		# Convert @array into @mz_array and @int_array
		# This better matches the format of mzdata and mzML
		for($num_xy_ctr = 0; $num_xy_ctr < $num_xy; $num_xy_ctr++)
		{
			$mz_array[$num_xy_ctr] = $array[2*$num_xy_ctr];
			$int_array[$num_xy_ctr] = $array[2*$num_xy_ctr + 1];;
		}
		
		process_slice();
		
		$slice_ctr++;
	}
	close CSV;
	roprint("Done reading csv file");
	#STDOUT->autoflush(1);
}

# We need these three things
# $rt
# @mz_array
# @int_array
sub do_prn
{
	roprint("Start reading PRN file");

	@prn_array = readfile($csvfile);

	$slice_ctr = 0;
	
	# We arbitrarily set this as 1 because that's what we did for the AMP Matlab MALDI format in Step 1
	$rt = 1;

	@mz_array = ();
	@int_array = ();
	
	# Go through the @prn_array line by line
	for(my $ctr = 0; $ctr <= $#prn_array; $ctr++)
	{
		# Split the array into hopefully (mz, intensity) pairs
		my @array = split ' ', $prn_array[$ctr];
		
		# This checks and sees if we have two array elements
		# If not, exit the loop
		if($#array != 1)
		{
			last;
		}
		
		my $mz = $array[0];
		my $int = $array[1];
		
		# Push just takes the value and tacks it on the end of an array
		push(@mz_array, $mz);
		push(@int_array, $int);
	}
	
	$mass_min = $mz_array[0]; # Just take the first m/z value
	$mass_max = $mz_array[$#mz_array]; # Just take the last m/z value
	
	process_slice();

	roprint("Done reading PRN file");
}

sub do_mzdata
{
	open MZDATA, "$csvfile" or die "Error: Can't open $csvfile\n";
	roprint("Reading mzdata file");
	#STDOUT->autoflush(1);

	$slice_ctr = 0;
	$flag1 = 0;
	$flag2 = 0;

	$mass_min = -1;
	$mass_max = -1;

	MZQUICKSCAN: while(defined($input=<MZDATA>))
	{
		if($input =~ /<spectrumInstrument msLevel=\"1\" mzRangeStart=\"(\d+\.\d+)\" mzRangeStop=\"(\d+\.\d+)\">/) #"
		{
			$temp_mass_min = $1;
			$temp_mass_max = $2;
			
			if($mass_min == -1)
			{
				$mass_min = $temp_mass_min;
				$mass_max = $temp_mass_max;
			}
			else
			{
				if($temp_mass_min < $mass_min)
				{
					$mass_min = $temp_mass_min;
					#print "MIN: $mass_min\n";
				}
				if($temp_mass_max > $mass_max)
				{
					$mass_max = $temp_mass_max;
					#print "MAX: $mass_max\n";
				}
			}
		}
	}

	close MZDATA;
	open MZDATA, "$csvfile" or die "Error: Can't open $csvfile\n";

	MZDATALINE: while(defined($input=<MZDATA>))
	{
		if($input =~ /<spectrum id=\"\d+\">/) #"
		{
			$flag1 = 1;
			#print "$input";
			next MZDATALINE;
		}
		if($input =~ /msLevel=\"1\"/) #"
		{
			$flag2 = 1;
			#print "$input";
			# Extract mass_min and mass_max
			
			#$input =~ /<spectrumInstrument msLevel=\"1\" mzRangeStart=\"(\d+\.\d+)\" mzRangeStop=\"(\d+\.\d+)\">/; #"
			#$temp_mass_min = $1;
			#$temp_mass_max = $2;
			#print "$temp_mass_min -> $temp_mass_max\n";

			#if($slice_ctr == 0)
			#{
			#	$mass_min = $temp_mass_min;
			#	$mass_max = $temp_mass_max;
			#}
			#else
			#{
			#	if($temp_mass_min < $mass_min)
			#	{
			#		$mass_min = $temp_mass_min;
			#	}
			#	if($temp_mass_max > $mass_max)
			#	{
			#		$mass_max = $temp_mass_max;
			#	}
			#}
			
			$mzline = "";
			$intline = "";
			
			# We have located an appropriate spectrum entry, so now go through and process it
			MZDSPECTRUM: while(defined($input2=<MZDATA>))
			{
				if($input2 =~ /TimeInMinutes/)
				{
					#print "$input2";
					$input2 =~ /.*TimeInMinutes\" value=\"(\d+\.\d+)\".*/; #"
					$rt = $1;
					#print "$rt\n";
					
					if($filter_rt == 1)
					{
						if($rt < $rt_mincut)
						{
							last MZDSPECTRUM;
						}
						elsif($rt > $rt_maxcut)
						{
							last MZDSPECTRUM;
						}
					}				
				}
				elsif($input2 =~ /<mzArrayBinary>/)
				{
					$mzline = <MZDATA>;
				}
				elsif($input2 =~ /<intenArrayBinary>/)
				{
					$intline = <MZDATA>;
				}
				# Terminate the spectrum
				elsif($input2 =~ /<\/spectrum>/ && $flag1 == 1 && $flag2 == 1)
				{
					# Process $mzline and $intline
					if($mzline eq "" || $intline eq "")
					{
						roprint("ERROR: Missing mz or int data");
						#STDOUT->autoflush(1);
					}
				
					$mzline =~ /.*precision=\"(\d+)\".*length=\"(\d+)\".*/; #"
					$precision = $1;
					$num_xy = $2;
					#print "$precision $num_xy\n";
					
					$mzline =~ /<.*>(.*)<.*>/;
					$b64_mz = $1;
					
					if($precision == 32)
					{
						@mz_array = unpack("f$num_xy", decode_base64($b64_mz));
					}
					elsif($precision == 64)
					{
						@mz_array = unpack("d$num_xy", decode_base64($b64_mz));
					}
					else
					{
						roprint(sprintf("ERROR: Unrecognized precision value |%d|", $precision));
						#STDOUT->autoflush(1);
					}

					$intline =~ /.*precision=\"(\d+)\".*length=\"(\d+)\".*/; #"
					$precision = $1;
					$num_xy2 = $2;
					#print "$precision $num_xy\n";
					if($num_xy != $num_xy2)
					{
						roprint(sprintf("ERROR: num_mz (%d) != num_int (%d)", $num_xy, $num_xy2));
						#STDOUT->autoflush(1);
					}
					
					$intline =~ /<.*>(.*)<.*>/;
					$b64_int = $1;
					
					if($precision == 32)
					{
						@int_array = unpack("f$num_xy2", decode_base64($b64_int));
					}
					elsif($precision == 64)
					{
						@int_array = unpack("d$num_xy2", decode_base64($b64_int));
					}
					else
					{
						roprint(sprintf("ERROR: Unrecognized precision value |%d|", $precision));
						#STDOUT->autoflush(1);
					}
					
					process_slice();
					
					$slice_ctr++;
					last MZDSPECTRUM;
				}
			}
			
			next MZDATALINE;
		}
		elsif($input =~ /msLevel=\"2\"/) #"
		{
			$flag1 = 0;
			$flag2 = 0;
			next MZDATALINE;
		}
	}
}

sub do_mzml_regexp
{
	#@ext_array = split /\./, $csvfile;
	#if($ext_array[$#ext_array] eq "gz")
	#{
	#	open MZMLF, "gunzip -c $csvfile |" or die "Error, Can't open $csvfile\n";
	#}
	#else
	#{
	#	open MZMLF, "$csvfile" or die "Error: Can't open $csvfile\n";
	#}
    if($gz)
    {
		open MZMLF, "gunzip -c $csvfile |" or die "Error, Can't open $csvfile\n";
    }
    else
    {
		open MZMLF, "$csvfile" or die "Error: Can't open $csvfile\n";
    }
    
	MZML: while(defined($input=<MZMLF>))
	{
		chomp($input);

		# Start processing a new spectrum entry
		if($input =~ /<spectrum index/)
		{
			#print "START $input\n";
			
			$mslevel = 0;
			
			MZML_SPEC: while(defined($input2=<MZMLF>))
			{
				chomp($input2);
			
				$cvparam_name = "";
				$cvparam_value = "";
				
				if($input2 =~ /<cvParam.*name=\"(.*?)\".*.*/)
				{
					$cvparam_name = $1;
					
					# Skip MS2 scans if they are present
					if(lc($cvparam_name) eq "ms level")
					{
						$input2 =~ /<cvParam.*name=\"(.*?)\".*value=\"(.*?)\".*/;
						$mslevel = $2;
					
						# Search for end of this <spectrum> entry
						if($mslevel > 1)
						{
							while(defined($input3=<MZMLF>))
							{
								chomp($input3);
								
								if($input3 =~ /<\/spectrum>/)
								{
									next MZML;
								}
							}
						}
					}
					
					# Get out RT in minutes
					if(lc($cvparam_name) eq "scan start time")
					{
						$input2 =~ /<cvParam.*name=\"(.*?)\".*value=\"(.*?)\".*unitName=\"(.*?)\".*/;
					
						$rt = $2;
						$rt_unit = $3;
						
						if(lc($rt_unit) ne "minute")
						{
							roprint("Warning, RT in units of |$rt_unit|.");
						}
					}
					
					#print "|$cvparam_name| -> |$cvparam_value|\n";
				}
				elsif($input2 =~ /<binaryDataArrayList/)
				{
					MZML_BLIST: while(defined($input3=<MZMLF>))
					{
						chomp($input3);
						
						$cvparam_mz = 0;
						$cvparam_int = 0;
						$cvparam_precision = 0;
						
						if($input3 =~ /<binaryDataArray encodedLength/)
						{
							MZML_BDATA: while(defined($input4=<MZMLF>))
							{
								chomp($input4);
							
								if($input4 =~ /<\/binaryDataArray>/)
								{
									last MZML_BDATA;
								}
								elsif($input4 =~ /<cvParam.*name=\"(.*?)\".*.*/)
								{
									$cvparam_name = $1;
									
									if($cvparam_name eq "m/z array")
									{
										#print "MZ\n";
										$cvparam_mz = 1;
										$cvparam_int = 0;
									}
									elsif($cvparam_name eq "intensity array")
									{
										#print "INT\n";
										$cvparam_int = 1;
										$cvparam_mz = 0;
									}
									elsif($cvparam_name eq "32-bit float")
									{
										$cvparam_precision = 32;
									}
									elsif($cvparam_name eq "64-bit float")
									{
										$cvparam_precision = 64;
									}
								}
								elsif($input4 =~ /<binary>/)
								{
									#print "BINARY\n";
								
									$input4 =~ /<.*>(.*)<.*>/;
									$b64_line = $1;
									
									if($cvparam_precision == 32)
									{
										# The < (endian flag) caused errors on garibaldi
										#@array = unpack("f<*", decode_base64($b64_line));
										@array = unpack("f*", decode_base64($b64_line));
									}
									elsif($cvparam_precision == 64)
									{
										#@array = unpack("d<*", decode_base64($b64_line));
										@array = unpack("d*", decode_base64($b64_line));
									}
									else
									{
										roprint("Error, precision not detected.");
										die;
									}
									
									if($cvparam_mz == 1)
									{
										#print "MZ2\n";
										@mz_array = @array;
									}
									elsif($cvparam_int == 1)
									{
										#print "INT2\n";
										@int_array = @array;
									}
								}
							}
						}
						
						if($input3 =~ /<\/binaryDataArrayList>/)
						{
							last MZML_BLIST;
						}
					}
				}
			
				if($input2 =~ /<\/spectrum>/)
				{
					#print "($mz_array[0],$int_array[0]) ($mz_array[1],$int_array[1]) ($mz_array[2],$int_array[2])\n";
				
					#print "END $input2\n----------------------------------\n";
					
					if($filter_rt == 1)
					{
						if($rt < $rt_mincut)
						{
							last MZML_SPEC;
						}
						elsif($rt > $rt_maxcut)
						{
							last MZML_SPEC;
						}
					}
					
					process_slice();
					
					$slice_ctr++;
					last MZML_SPEC;
				}
			}
		}
	}
	
	close MZMLF;
}

sub do_mzml
{
    roprint("Start do_mzml");
    
    $spectrum_ctr = 0;
    
	$xml = XML::Twig->new( twig_handlers =>
    {
        spectrum => \&do_mzml_spectrum,
    },
    );
    
	if($gz)
	{
		#$gzfile = new IO::Zlib;
		#$gzfile->open($csvfile, "rb");
		#$xml2->parse($gzfile);
        
        # Parsing dies when using a gzipped file (no idea why)
        # So just unzip it here and proceed
        # Make sure it hasn't already been gunzip'd
        if(!(-e $csvfile_nogz))
        {
            `gunzip -c $csvfile > $csvfile_nogz`;            
        }
        $xml->parsefile($csvfile_nogz);
    }
	else
	{
		$xml->parsefile($csvfile);
	}
    
	$xml->purge;
    
    roprint("End do_mzml");
}

sub do_mzml_spectrum
{
    #roprint("do_mzml_spectrum");
    
    my($xml, $spectrum) = @_;
    
	my @cvparam = $spectrum->children('cvParam');
    
    # For mzML files generated by proteowizard, this will eliminate MS2 scans
	for(my $ctr = 0; $ctr <= $#cvparam; $ctr++)
	{
		my $name = $cvparam[$ctr]->att('name');
		my $value = $cvparam[$ctr]->att('value');
		if($name eq "ms level" && $value > 1)
		{
			$spectrum->purge;
			return
		}
	}
    
    # For the ABSciex mzML files, there is no "ms level" parameter as per above
    # For these, parse out the experiment, and eliminate those > 1
    # This is the same technique used in do_mzml_scan2rt_spectrum()
	my $id = $spectrum->att('id');
    
	my $scan_num;
    my $experiment = 0;
    
	#if($id =~ /controllerType=(\d*) controllerNumber=(\d*) scan=(\d*)/)
	if($id =~ /cycle=(\d+) experiment=(\d+)/)
    {
        $scan_num = $1;
        $experiment = $2;
    }
    # This ignores > MS1 scans (see notes above)
    if($experiment > 1)
    {
        $spectrum->purge;
        return;
    }
    
    
    # Get the RT value from "scan start time"
	my @cvparam2 = $spectrum->get_xpath('scanList/scan/cvParam');
    
    # Actually can't have this as a local variable as process_slice() uses it
	#my $rt;
    $rt = -1;
    
	my $did_find_rt = 0;
    
	for(my $ctr = 0; $ctr <= $#cvparam2; $ctr++)
	{
		my $name = $cvparam2[$ctr]->att('name');
		my $value;
		my $unit;
		
		if($name eq "scan start time")
		{
			$value = $cvparam2[$ctr]->att('value');
			$unit = $cvparam2[$ctr]->att('unitName');
			
			#roprint("value -> unit $value -> $unit"); 
			
			$did_find_rt = 1;
            
			if($unit eq "second")
			{
				$rt = $value/60;
			}
			else
			{
				$rt = $value;
				
				#if($spectrum_ctr % 100 == 0)
				#{
				#	roprint("RT1: $rt");
				#}
			}		
		}
	}
    
    if($did_find_rt == 0)
    {
		roprint("ERROR: Could not find RT value do_mzml_spectrum");
		die;
    }
    
    # Build @mz_array and @int_array
    # Handle 32 or 64 bit, zlib/no compression, and either order
    my @cvparam3 = $spectrum->get_xpath('binaryDataArrayList/binaryDataArray/cvParam');
    
    my $which_mz = -1;
    my $which_int = -1;
    my @precision = ();
    my @compression = ();
    
    for(my $ctr = 0; $ctr <= $#cvparam3; $ctr++)
    {
        my $name = $cvparam3[$ctr]->att('name');
        #print "$name\n";
        if(lc($name) eq "m/z array")
        {
            # m/z array is first (index 0)
            if($which_mz == -1 && $which_int == -1)
            {
                $which_mz = 0;
            }
            # m/z array is second (index 1)
            elsif($which_mz == -1)
            {
                $which_mz = 1;
            }
            # Should not happen - this means which_mz has already been found, and this is a second one, so die
            elsif($which_int == -1)
            {
                roprint("Error, second instance of m/z array in this spectrum");
                die;
            }
        }
        elsif(lc($name) eq "intensity array")
        {
            # m/z array is first (index 0)
            if($which_mz == -1 && $which_int == -1)
            {
                $which_int = 0;
            }
            # m/z array is second (index 1)
            elsif($which_int == -1)
            {
                $which_int = 1;
            }
            # Should not happen - this means which_int has already been found, and this is a second one, so die
            elsif($which_mz == -1)
            {
                roprint("Error, second instance of intensity array in this spectrum");
                die;
            }
        }
        elsif(lc($name) eq "32-bit float")
        {
            push(@precision, 32);
        }
        elsif(lc($name) eq "64-bit float")
        {
            push(@precision, 64);
        }
        elsif(lc($name) eq "no compression")
        {
            push(@compression, 0);
        }
        elsif(lc($name) eq "zlib compression")
        {
            push(@compression, 1);
        }
    }
    
    #print "MZ: $which_mz INT: $which_int\n";
    
    if($which_mz == -1)
    {
        roprint("Error, m/z array not detected");
        die;
    }
    if($which_int == -1)
    {
        roprint("Error, intensity array not detected");
        die;
    }
    
	my @binary = $spectrum->get_xpath('binaryDataArrayList/binaryDataArray/binary');
	
	#print "$#binary\n";
	
	my $mz_data = $binary[$which_mz]->field;
	my $int_data = $binary[$which_int]->field;
	
	@mz_array = ();
    if($precision[$which_mz] == 32)
    {
        if($compression[$which_mz] == 0)
        {
            #$test_t1 = time();
            @mz_array = unpack("f*", decode_base64($mz_data));
            #$test_t2 = time();
            #$test_diff = $test_t2 - $test_t1;
            #roprint("Elapsed: $test_diff");
        }
        elsif($compression[$which_mz] == 1)
        {
            #@mz_array = unpack("f*", uncompress(decode_base64($mz_data)));
            my $uncompress_data = uncompress(decode_base64($mz_data));
            @mz_array = unpack("f*", $uncompress_data);
            $uncompress_data = "";
        }
    }
    elsif($precision[$which_mz] == 64)
    {
        if($compression[$which_mz] == 0)
        {
            @mz_array = unpack("d*", decode_base64($mz_data));
        }
        elsif($compression[$which_mz] == 1)
        {
            @mz_array = unpack("d*", uncompress(decode_base64($mz_data)));
        }
    }
    else
    {
        roprint("Error, unrecognized precision: |$precision[$which_mz]|");
        die;
    }
    $mz_data = "";
    
	@int_array = ();
    if($precision[$which_int] == 32)
    {
        if($compression[$which_int] == 0)
        {
            @int_array = unpack("f*", decode_base64($int_data));
        }
        elsif($compression[$which_int] == 1)
        {
            #@int_array = unpack("f*", uncompress(decode_base64($int_data)));
            my $uncompress_data = uncompress(decode_base64($int_data));
            @int_array = unpack("f*", $uncompress_data);
            $uncompress_data = "";
        }
    }
    elsif($precision[$which_int] == 64)
    {
        if($compression[$which_int] == 0)
        {
            @int_array = unpack("d*", decode_base64($int_data));
        }
        elsif($compression[$which_int] == 1)
        {
            @int_array = unpack("d*", uncompress(decode_base64($int_data)));
        }
    }
    else
    {
        roprint("Error, unrecognized precision: |$precision[$which_int]|");
        die;
    }
    $int_data = "";
    
    #$test_t1 = time();
    process_slice();
    #$test_t2 = time();
    #$test_diff = $test_t2 - $test_t1;
    
    #$diff1 = $ps_t2 - $ps_t1;
    #$diff2 = $ps_t3 - $ps_t2;
    #roprint("Elapsed: $test_diff $diff1 $diff2 $#mz_array $#int_array $num_good_spectra");
    
    $slice_ctr++;
    
    $spectrum->purge;
}

#
# This takes the data for a particular slice, and checks it against each minispectrum to be extracted as defined in @spectra
# If the $rt is in the proper range, it puts data into @spectra at the appropriate location
#
# Does the RT range of the minispectrum as previously defined match the RT of the slice?
#
sub process_slice
{
	if($slice_ctr %100 == 0)
	{
		#print "Slice: $slice_ctr RT: $rt NumXY: $num_xy\n";
		$t_elapsed = (time() - $t_start)/60;
		roprint(sprintf("Slice: %d RT: %f Elapsed Time: %.1f min", $slice_ctr, $rt, $t_elapsed));
		#STDOUT->autoflush(1);
	}
					
	# Zero this array every $slice_per_pt times
	if($slice_ctr % $slice_per_pt == 0)
	{
		$contour_rt_sum = 0;
		@contour_array = ();
	}
	
	$contour_rt_sum += $rt;

	# 12/16/10 Need to reset this every slice due to orbitrap data
	# orbitrap data is sparse, and so many indices are missing
	# Without resetting, was just reading from a previous slice's values
	%array_index = ();

	#for($pt_ctr = 0; $pt_ctr < $num_xy; $pt_ctr++)
	for($pt_ctr = 0; $pt_ctr <= $#mz_array; $pt_ctr++)
	{
		$mz_val = $mz_array[$pt_ctr];

		$x = ($mz_val - $mass_min) / $mz_per_pt;
		$contour_array[$x] += $int_array[$pt_ctr];
	
		$mz_int = sprintf("%d", $mz_val);
		$array_index{$mz_int} = $pt_ctr;
	}

	# One-time only printing of mz labels to IMAGETEXT
	if($slice_ctr == 0)
	{
		print IMAGETEXT "X ";
		for($ctr = 0; $ctr <= ($mass_max - $mass_min) / $mz_per_pt; $ctr++)
		{
			$contour_mz = $mass_min + ($ctr / ( ($mass_max - $mass_min) / $mz_per_pt) ) * ($mass_max - $mass_min);
			printf(IMAGETEXT "%.3f ", $contour_mz);
		}
		print IMAGETEXT "\n";
	}

	if( ($slice_ctr + 1) % $slice_per_pt == 0)
	{
		$contour_rt = $contour_rt_sum / $slice_per_pt;
		printf(IMAGETEXT "%.2f ", $contour_rt);
		for($ctr = 0; $ctr <= ($mass_max - $mass_min) / $mz_per_pt; $ctr++)
		{
			if(defined($contour_array[$ctr]))
			{
				print IMAGETEXT "$contour_array[$ctr] ";
			}
			else
			{
				print IMAGETEXT "0 ";
			}
		}
		print IMAGETEXT "\n";
	}

	#
	# Main loop checking against each minispectrum
	#
	for($spec_ctr = 0; $spec_ctr < $num_good_spectra; $spec_ctr++)
	{
		if($rt >= $spectra[$spec_ctr]{"rt_min_plot"} && $rt < $spectra[$spec_ctr]{"rt_min"})
		{
			#print "RTYES1\n";
			
			#newcontour
			push(@{$spectra[$spec_ctr]{"rt_vals"}}, $rt);
		
			if($rt < $spectra[$spec_ctr]{"rt_min"})
			{
				$spectra[$spec_ctr]{"slices_before"} += 1;
			}
			elsif($rt > $spectra[$spec_ctr]{"rt_max"})
			{
				$spectra[$spec_ctr]{"slices_after"} += 1;
			}
		
			$spectra[$spec_ctr]{"slices_plot"} += 1;
			
			$startloc = $spectra[$spec_ctr]{"mz_min_plot"} - 2;
			$startloc = sprintf("%d", $startloc);
			
			# 12/16/10 loop to scan for defined $array_index{} locations
			# This is necessitated by the sparse orbitrap data
			# If nothing is defined, $startloc increments smaller, $endloc larger
			while(1)
			{
				if(defined($array_index{$startloc}))
				{
					$startloc = $array_index{$startloc};
					last;
				}
				elsif($startloc < $mz_array[0])
				{
					$startloc = 0;
					last;
				}
				else
				{
					$startloc -= 1;
				}
			}

			$endloc = $spectra[$spec_ctr]{"mz_max_plot"} + 2;
			$endloc = sprintf("%d", $endloc);
			
			while(1)
			{
				if(defined($array_index{$endloc}))
				{
					$endloc = $array_index{$endloc};				
					last;
				}
				elsif($endloc > $mz_array[$#mz_array])
				{
					$endloc = $#mz_array;
					last;
				}
				else
				{
					$endloc += 1;
				}
			}
		
			#if($blueid eq "sykes")
			#{
			#	#print "1s->e $startloc $endloc\n";
			#	roprint("$spec_ctr 1s->e $startloc $endloc");
			#}
		
			for($xy_ctr = $startloc; $xy_ctr < $endloc; $xy_ctr++)
			#for($xy_ctr = 0; $xy_ctr < $num_xy; $xy_ctr++)
			{
				#$mz_val = sprintf("%.4f", $array[2*$xy_ctr]);
				$mz_val = $mz_array[$xy_ctr];
				$int_val = $int_array[$xy_ctr];
				
				if($mz_val > $spectra[$spec_ctr]{"mz_min_plot"} && $mz_val < $spectra[$spec_ctr]{"mz_max_plot"})
				{
					#print "MZ_YES\n";
					$slice_num = $spectra[$spec_ctr]{"slices_plot"} - 1;
					#$xy_val = join ' ', $slice_num, $mz_val, $int_val;
					$xy_val = join ' ', $mz_val, $int_val;
					#push @{ $spectra_values_plot[$spec_ctr] }, $xy_val;
					
					$spectra_values_plot[$spec_ctr][$slice_num] = join ' ', $spectra_values_plot[$spec_ctr][$slice_num], $xy_val;
					
					#push @{ $spectra_values_plot[$spec_ctr] }, $slice_num;
					#push @{ $spectra_values_plot[$spec_ctr] }, $mz_val;
					#push @{ $spectra_values_plot[$spec_ctr] }, $int_val;
				}
			}			
		}
		elsif($rt >= $spectra[$spec_ctr]{"rt_min"} && $rt <= $spectra[$spec_ctr]{"rt_max"})
		{
			#if($spectra3D)
            if($extractstyle == 2)
            {
				push(@{$spectra_values[$spec_ctr]}, "START");
			}
		
			#print "RTYES2\n";
		
			#newcontour
			push(@{$spectra[$spec_ctr]{"rt_vals"}}, $rt);

			# Consolidate the Two Regions
			$spectra[$spec_ctr]{"slices"} += 1;

			if($rt < $spectra[$spec_ctr]{"rt_min"})
			{
				$spectra[$spec_ctr]{"slices_before"} += 1;
			}
			elsif($rt > $spectra[$spec_ctr]{"rt_max"})
			{
				$spectra[$spec_ctr]{"slices_after"} += 1;
			}
		
			$spectra[$spec_ctr]{"slices_plot"} += 1;
			
			# Use the wider range from mz_min_plot and mz_max_plot
			$startloc = $spectra[$spec_ctr]{"mz_min_plot"} - 2;
			$startloc = sprintf("%d", $startloc);
			#$startloc = $array_index{$startloc};

			while(1)
			{
				if(defined($array_index{$startloc}))
				{
					$startloc = $array_index{$startloc};
					last;
				}
				elsif($startloc < $mz_array[0])
				{
					$startloc = 0;
					last;
				}
				else
				{
					$startloc -= 1;
				}
			}


			$endloc = $spectra[$spec_ctr]{"mz_max_plot"} + 2;
			$endloc = sprintf("%d", $endloc);
			#$endloc = $array_index{$endloc};				

			while(1)
			{
				if(defined($array_index{$endloc}))
				{
					$endloc = $array_index{$endloc};				
					last;
				}
				elsif($endloc > $mz_array[$#mz_array])
				{
					$endloc = $#mz_array;
					last;
				}
				else
				{
					$endloc += 1;
				}
			}

			#if($blueid eq "sykes")
			#{
			#	$mz_print_label = $spectra[$spec_ctr]{"label"};
			#	$mz_print_low1 = $spectra[$spec_ctr]{"mz_min_plot"} - 2;
			#	$mz_print_high1 = $spectra[$spec_ctr]{"mz_max_plot"} + 2;
			#	$mz_print_low2 = $mz_array[$startloc];
			#	$mz_print_high2 = $mz_array[$endloc];
			#	#print "2s->e $startloc $endloc\n";
			#	roprint("$spec_ctr 2s->e |$startloc| |$endloc| |$mz_print_low1| |$mz_print_high1| |$mz_print_low2| |$mz_print_high2| |$mz_print_label|");
			#}
			
			for($xy_ctr = $startloc; $xy_ctr < $endloc; $xy_ctr++)
			{
				#$mz_val = sprintf("%.4f", $array[2*$xy_ctr]);
				$mz_val = $mz_array[$xy_ctr];
				$int_val = $int_array[$xy_ctr];
				
				if($mz_val > $spectra[$spec_ctr]{"mz_min_plot"} && $mz_val < $spectra[$spec_ctr]{"mz_max_plot"})
				{
					#print "MZ_YES\n";
					$slice_num = $spectra[$spec_ctr]{"slices_plot"} - 1;
					#$xy_val = join ' ', $slice_num, $mz_val, $int_val;
					$xy_val = join ' ', $mz_val, $int_val;
					#push @{ $spectra_values_plot[$spec_ctr] }, $xy_val;

					$spectra_values_plot[$spec_ctr][$slice_num] = join ' ', $spectra_values_plot[$spec_ctr][$slice_num], $xy_val;
					
					#push @{ $spectra_values_plot[$spec_ctr] }, $slice_num;
					#push @{ $spectra_values_plot[$spec_ctr] }, $mz_val;
					#push @{ $spectra_values_plot[$spec_ctr] }, $int_val;
				}				
			
				if($mz_val > $spectra[$spec_ctr]{"mz_min"} && $mz_val < $spectra[$spec_ctr]{"mz_max"})
				{
					#print "MZ_YES\n";
					$xy_val = join ' ', $mz_val, $int_val;
					push @{ $spectra_values[$spec_ctr] }, $xy_val;
				}
			}
			
			#if($spectra3D)
            if($extractstyle == 2)
			{
				push(@{$spectra_values[$spec_ctr]}, "END");
			}			
		
		}
		elsif($rt > $spectra[$spec_ctr]{"rt_max"} && $rt <= $spectra[$spec_ctr]{"rt_max_plot"})
		{
			#print "RTYES3\n";

			#newcontour
			push(@{$spectra[$spec_ctr]{"rt_vals"}}, $rt);
		
			if($rt < $spectra[$spec_ctr]{"rt_min"})
			{
				$spectra[$spec_ctr]{"slices_before"} += 1;
			}
			elsif($rt > $spectra[$spec_ctr]{"rt_max"})
			{
				$spectra[$spec_ctr]{"slices_after"} += 1;
			}
		
			$spectra[$spec_ctr]{"slices_plot"} += 1;
			
			$startloc = $spectra[$spec_ctr]{"mz_min_plot"} - 2;
			$startloc = sprintf("%d", $startloc);
			#$startloc = $array_index{$startloc};

			while(1)
			{
				if(defined($array_index{$startloc}))
				{
					$startloc = $array_index{$startloc};
					last;
				}
				elsif($startloc < $mz_array[0])
				{
					$startloc = 0;
					last;
				}
				else
				{
					$startloc -= 1;
				}
			}

			$endloc = $spectra[$spec_ctr]{"mz_max_plot"} + 2;
			$endloc = sprintf("%d", $endloc);
			#$endloc = $array_index{$endloc};				

			while(1)
			{
				if(defined($array_index{$endloc}))
				{
					$endloc = $array_index{$endloc};				
					last;
				}
				elsif($endloc > $mz_array[$#mz_array])
				{
					$endloc = $#mz_array;
					last;
				}
				else
				{
					$endloc += 1;
				}
			}

			#if($blueid eq "sykes")
			#{
			#	#print "3s->e $startloc $endloc\n";
			#	roprint("$spec_ctr 3s->e $startloc $endloc");
			#}
		
			for($xy_ctr = $startloc; $xy_ctr < $endloc; $xy_ctr++)
			#for($xy_ctr = 0; $xy_ctr < $num_xy; $xy_ctr++)
			{
				#$mz_val = sprintf("%.4f", $array[2*$xy_ctr]);
				$mz_val = $mz_array[$xy_ctr];
				$int_val = $int_array[$xy_ctr];
				
				if($mz_val > $spectra[$spec_ctr]{"mz_min_plot"} && $mz_val < $spectra[$spec_ctr]{"mz_max_plot"})
				{
					#print "MZ_YES\n";
					$slice_num = $spectra[$spec_ctr]{"slices_plot"} - 1;
					#$xy_val = join ' ', $slice_num, $mz_val, $int_val;
					$xy_val = join ' ', $mz_val, $int_val;
					#push @{ $spectra_values_plot[$spec_ctr] }, $xy_val;

					$spectra_values_plot[$spec_ctr][$slice_num] = join ' ', $spectra_values_plot[$spec_ctr][$slice_num], $xy_val;
					
					#push @{ $spectra_values_plot[$spec_ctr] }, $slice_num;
					#push @{ $spectra_values_plot[$spec_ctr] }, $mz_val;
					#push @{ $spectra_values_plot[$spec_ctr] }, $int_val;
				}
			}			
		}
	}
}

# Are these working?!  10/10/2007 - Sykes
# All this sorting is now done in massive 10/22/2010 - Sykes

sub bymz
{
	@array_a = split ' ', $a;
	@array_b = split ' ', $b;
	$array_a[0] <=> $array_b[0];
}

sub bymz2
{
	@array_a = split ' ', $a;
	@array_b = split ' ', $b;
	$array_a[ $output_name2col{"n14mass"} ] <=> $array_b[ $output_name2col{"n14mass"} ];
}

sub byrt
{
	@array_a = split ' ', $a;
	@array_b = split ' ', $b;
	$array_a[ $output_name2col{"rt_n14"} ] <=> $array_b[ $output_name2col{"rt_n14"} ];
}

sub bypep
{
	@array_a = split ' ', $a;
	@array_b = split ' ', $b;
	$array_a[ $output_name2col{"digestpepid"} ] <=> $array_b[ $output_name2col{"digestpepid"} ]
		or
	$array_a[ $output_name2col{"n14mass"} ] <=> $array_b[ $output_name2col{"n14mass"} ];
}

sub roprint
{
	my @vars = @_;
	my $outstring = join "\n", @vars;
	
	`echo "$outstring" | cat - >> $rofile`;
}
