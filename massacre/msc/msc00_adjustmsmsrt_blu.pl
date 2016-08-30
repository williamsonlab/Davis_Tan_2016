#!/usr/bin/perl

use Getopt::Long;
use IO::Handle;
use Text::ParseWords;
use MIME::Base64;
use XML::Twig;
use IO::Zlib;
use Compress::Zlib;

require "$ENV{'MASSACRE_PATH'}/msc/mikesubs.pl";

$t_start = time();
$pid = $$;

$whoami = `whoami`;
chomp($whoami);

&GetOptions("input=s" => \$inputfile, "csv=s" => \$csvfile, "dir=s" => \$dir, "blueid=s" => \$blueid, "hostname=s" => \$hostname, "rofile=s" => \$rofile, "ailtype=i" => \$ailtype, "has_rt=i" => \$has_rt);

$outputfile = join '.', $inputfile, "2";

if( substr($dir, -1, 1) ne "/" )
{
	$dir = join '', $dir, "/";
}
$rval = mkdir($dir);

# Check for file type and process based on result
$gz = 0;
@csvfilearray = split /\./, lc($csvfile);
@csvfilearray_nolc = split /\./, $csvfile;

if($csvfilearray[$#csvfilearray] eq "gz")
{
	$gz = 1;
	pop(@csvfilearray);
    pop(@csvfilearray_nolc);
    $csvfile_nogz = join '.', @csvfilearray_nolc; # Use this later when parsing
}

$filetype1 = lc($csvfilearray[$#csvfilearray]);
$filetype2 = lc(join '.', $csvfilearray[$#csvfilearray-1], $csvfilearray[$#csvfilearray]);

%scan2rt = ();

# Convert scan numbers to actual RT for sequest data or mascot data that only has scan numbers
if($ailtype == 5 || $has_rt <= 0)
{    
	if($filetype1 eq "mzml")
	{
		do_mzml_scan2rt();
    }
	else
	{
		roprint("ERROR: Filetype |$filetype1| |$filetype2| unrecognized.");
		die;
	}
}

#foreach $scan_key (keys %scan2rt)
#{
#	roprint("$scan_key -> $scan2rt{$scan_key}");
#}

# Process the input file
roprint("Processing input file.");
@input_array = readfile($inputfile);
for($ctr = 0; $ctr <= $#input_array; $ctr++)
{
	@array = split ' ', $input_array[$ctr];
	#roprint("$array[0].rt");
	#`touch $dir/$array[0].rt`;

	$key = shift(@array);
	$z = shift(@array);
	$mz = shift(@array);
	$rt = shift(@array);

	#roprint("KEY Z MZ RT |$key|$z|$mz|$rt|");

    # Convert the scan number into an actual retention time
    # The tricky thing is that the scan numbers given are usually MS2 scans, but these are not kept in the conversion (file size reasons)
    # So we need to find the closest MS1 match based on the scan number
	if($ailtype == 5 || $has_rt <= 0)
	{
		$scan2rt_offset = 0;
		while(1)
		{
            # First check and see if the scan number is defined
			if(defined($scan2rt{$rt - $scan2rt_offset}))
			{
                # If it is, keep the rt value that is returned
				$rt = $scan2rt{$rt - $scan2rt_offset};
				last;
			}
            # If not, increment the offset until it matches
            # This is fine as we do an RT correction step anyways and the actual time difference is < 0.01 minutes in most cases
			else
			{
				$scan2rt_offset++;

				if($scan2rt_offset > $rt)
				{
					roprint("ERROR: scan2rt_offset ($scan2rt_offset) > scan ($rt)");
					die;
				}
			}
		}
	}

	#if($ctr % 10 == 0)
	#{
	#	roprint("KEY Z MZ RT |$key|$z|$mz|$rt|");
	#}

	# z +/-
	# 1 0.50
	# 2 0.25
	# 3 0.167

	$mz_expand = 1/(2*$z);
	$rt_expand = 0.5;
	
	$spectra[$ctr]{"rt"} = $rt;
	$spectra[$ctr]{"mz_min"} = $mz - $mz_expand;
	$spectra[$ctr]{"mz_max"} = $mz + $mz_expand;
	$spectra[$ctr]{"rt_min"} = $rt - $rt_expand;
	$spectra[$ctr]{"rt_max"} = $rt + $rt_expand;
	$spectra[$ctr]{"key"} = $key;
    $spectra[$ctr]{"output"} = 0; # Use this to keep track of if spectra has been written out
	
	$spectra[$ctr]{"mz_min_plot"} = $spectra[$ctr]{"mz_min"} - 5;
	$spectra[$ctr]{"mz_max_plot"} = $spectra[$ctr]{"mz_max"} + 5;
	$spectra[$ctr]{"rt_min_plot"} = $spectra[$ctr]{"rt_min"} - 2;
	$spectra[$ctr]{"rt_max_plot"} = $spectra[$ctr]{"rt_max"} + 2;
    
	#if($whoami eq "sykes")
	#{
	#	roprint("$spectra[$ctr]{'rt_min'} - $spectra[$ctr]{'rt_max'} * $spectra[$ctr]{'mz_min'} - $spectra[$ctr]{'mz_max'}");
	#}

}
roprint("Finished processing input file.");

$num_good_spectra = $ctr;


# Sort by RT.
# We can now write out spectra and purge them from memory as we go through the spectrum.  This is important for huge datasets, which take up a ton of memory.  It should also make the search process quicker.
#for(my $sortctr = 0; $sortctr <= $#spectra; $sortctr++)
#{
#    print "A $sortctr $spectra[$sortctr]{'rt'}\n";
#}
@spectra = sort by_rt @spectra;
#for(my $sortctr = 0; $sortctr <= $#spectra; $sortctr++)
#{
#    print "B $sortctr $spectra[$sortctr]{'rt'}\n";
#}



# Open this here, so that output can be written bit by bit as slices are processed
open OUTFILE, ">$outputfile" or die "Can't open $outputfile\n";




# Process the data file
# This could be .mzdata.xml or .mzML
$slice_ctr = 0;
$flag1 = 0;
$flag2 = 0;
$mass_min = -1;
$mass_max = -1;


if($filetype2 eq "mzdata.xml")
{
	roprint("Processing mzdata file.");
	do_mzdata();
}
elsif($filetype1 eq "mzml")
{
	# The output of msconvert does not report these accurately
	# <cvParam cvRef="MS" accession="MS:1000501" name="scan window lower limit" value="80.001051846596482" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
	# <cvParam cvRef="MS" accession="MS:1000500" name="scan window upper limit" value="2499.3999360293742" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
	# So just automatically peg these at 400 -> 2000 based on typical usage
	# ssc101911 changed from 400 -> 2000 to 200 -> 2000
	$mass_min = 200;
	$mass_max = 2000;
	
	roprint ("Processing mzML file.");
	do_mzml();
}
elsif($filetype1 eq "mzxml")
{
	$mass_min = 200;
	$mass_max = 2000;
	
	roprint ("Processing mzxml file.");
	do_mzxml();
}
else
{
	roprint("Error, file type |$filetype| not recognized.");
	die;
}


# We now have two important arrays
# 1) spectra_values[$spec_ctr][xy_vals]
# These are already properly averaged - just need to write them out
# This is in contrast to step msc04 where binning is needed
# 2) spectra_values_plot[$spec_ctr]
# Use the step4 style output for newcon plots

roprint("Generating output file.");


#for($spec_output_ctr = 0; $spec_output_ctr < $num_good_spectra; $spec_output_ctr++)
#{
#    output_spectrum();
#}

# 7/26/11 - Spectra are output during process_slice() when $rt > $rt_max_plot
# However, that condition is not always satisfied for all spectra, especially those very close to the end of the .mzML (or whatever) file
# Here we will not output the remaining spectra to tidy this up

for($spec_output_ctr = 0; $spec_output_ctr <= $#spectra; $spec_output_ctr++)
{
	output_spectrum();
}

close OUTFILE;

roprint("Finished generating output file.");


####################
####################
####################
####################
####################
####################
####################
####################



# Processing Data Files
# 1) Need to retrieve retention time in minutes ($rt)
# 2) Need to establish two arrays:
#      a) @mz_array (m/z values)
#      b) @int_array (intensity values)

#############
# do_mzdata #
#############
sub do_mzdata
{
	roprint("Processing mzdata file.");

	open MZDATA, "$csvfile" or die "Error: Can't open $csvfile\n";
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
				}
				if($temp_mass_max > $mass_max)
				{
					$mass_max = $temp_mass_max;
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
			next MZDATALINE;
		}
		if($input =~ /msLevel=\"1\"/) #"
		{
			$flag2 = 1;
			
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
					}

					$intline =~ /.*precision=\"(\d+)\".*length=\"(\d+)\".*/; #"
					$precision = $1;
					$num_xy2 = $2;
					#print "$precision $num_xy\n";
					if($num_xy != $num_xy2)
					{
						roprint(sprintf("ERROR: num_mz (%d) != num_int (%d)", $num_xy, $num_xy2));
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
					}
									
					# Go through all the spectra and see if you need to accumulate intensity
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
	close MZDATA;
}

# Old way to parse mzml files based on regular expressions
# This worked fine only if mzml files were generated in a particular way
# This is being replaced by a proper XML based parsing which should work with all proper mzML files
sub do_mzml_regexp
{
	#@ext_array = split /\./, $csvfile;
	#if($ext_array[$#ext_array] eq "gz")
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
						
						if(lc($rt_unit) eq "second")
						{
							$rt = $rt/60;
						}
						elsif(lc($rt_unit) ne "minute")
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
					#roprint("($mz_array[0],$int_array[0]) ($mz_array[1],$int_array[1]) ($mz_array[2],$int_array[2])");
				
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

# Virtually identical to do_mzml_scan2rt except it calls a different subroutine to handle the "spectrum" XML element
sub do_mzml
{
	roprint("Start do_mzml");
    
	$spectrum_ctr = 0;
    
	$xml2 = XML::Twig->new( twig_handlers =>
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
           #`gunzip $csvfile `;            
        }
        $xml2->parsefile($csvfile_nogz);
    }
	else
	{
		$xml2->parsefile($csvfile);
	}
    
	$xml2->purge;
    
    roprint("End do_mzml");
}

sub do_mzml_scan2rt
{
	roprint("Scanning for scan numbers and retention times.");

	$spectrum_ctr = 0;

	$xml = XML::Twig->new( twig_handlers =>
		{
			spectrum => \&do_mzml_scan2rt_spectrum,
		},
		);

	if($gz)
	{
		#$gzfile = new IO::Zlib;
		#$gzfile->open($csvfile, "rb");
		#$xml->parse($gzfile);

        # Parsing dies when using a gzipped file (no idea why)
        # So just unzip it here and proceed
        # Make sure it hasn't already been gunzip'd
        if(!(-e $csvfile_nogz))
        {
            `gunzip -c $csvfile > $csvfile_nogz`;            
            #`gunzip -c $csvfile `;            
        }
        $xml->parsefile($csvfile_nogz);
    }
	else
	{
		$xml->parsefile($csvfile);
	}

	$xml->purge;

    roprint("Finished scanning for scan numbers and retention times.");
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

sub do_mzml_scan2rt_spectrum
{
	my($xml, $spectrum) = @_;
	
	my $id = $spectrum->att('id');

	my $scan_num;
    my $experiment = 0;

	#if($id =~ /controllerType=(\d*) controllerNumber=(\d*) scan=(\d*)/)
	if($id =~ /scan=(\d+)/)
	{
		$scan_num = $1;
	}
	elsif($id =~ /scanId=(\d+)/)
	{
		$scan_num = $1;
	}
    # Triple-TOF Data from ABSciex looks like this:
    # <spectrum id="sample=1 period=1 cycle=2470 experiment=1" defaultArrayLength="58757" index="23413">
    # <spectrum id="sample=1 period=1 cycle=2470 experiment=9" defaultArrayLength="509" index="23421">
    # "cycle" appears to take the place of scanId/scan in this case
    # It's not clear, but based on the defaultArrayLength, "experiment" seems related to the MS level, which is not reported in these mzML files
    # So experiment=1 means MS1 scan, and >1 means MS2+ scan
    # The mascot output contains a pep_scan_title that looks like this:
    # Locus:1:1:1:1806:4
    # This seems to map to sample/period/cycle/experiment, except the last number often exceeds that of any experiment reported in the mzML file
    elsif($id =~ /cycle=(\d+) experiment=(\d+)/)
    {
        $scan_num = $1;
        $experiment = $2;
    }
	else
	{
		roprint("ERROR: Could not detect scan number |$id|.");
		die;
	}
	
    # This ignores > MS1 scans (see notes above)
    if($experiment > 1)
    {
        $spectrum->purge;
        return;
    }
    
	my @cvparam = $spectrum->children('cvParam');

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
	
	#print "Scan Num: $scan_num\n";
	
	for(my $ctr = 0; $ctr <= $#cvparam; $ctr++)
	{
		my $name = $cvparam[$ctr]->att('name');
		my $value = $cvparam[$ctr]->att('value');
		#print "$name -> $value\n";
	}

	my @cvparam2 = $spectrum->get_xpath('scanList/scan/cvParam');

	my $rt;

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
		roprint("ERROR: Could not find RT value do_mzml_scan2rt_spectrum");
		die;
	}
	
	#if($spectrum_ctr % 100 == 0)
	#{
	#	roprint("RT2: $rt");
	#}

	#print "----------\n----------\n";
	
	$scan2rt{$scan_num} = $rt;
	
	if($spectrum_ctr % 100 == 0)
	{
		$t_elapsed = (time() - $t_start)/60;
		$t_elapsed = sprintf("%.1f", $t_elapsed);
		roprint("Spectrum: $spectrum_ctr ($scan_num -> $rt) Time: $t_elapsed");
	}
	
	$spectrum->purge;
	
	$spectrum_ctr++;
}





#################
# process_slice #
#################
# Requires $rt, @mz_array, @int_array, $slice_ctr
sub process_slice
{
	#if($whoami eq "sykes")
	#{
	#	roprint("$mz_array[0] $int_array[0] $mz_array[1] $int_array[1] $mz_array[2] $int_array[2]");
	#}

	# Output info
	if($slice_ctr %100 == 0)
	{
		$t_elapsed = (time() - $t_start)/60;
		roprint(sprintf("Slice: %d RT: %f Time: %.1f min", $slice_ctr, $rt, $t_elapsed));
	}
	
	%array_index = ();
	
    $ps_t1 = time();
    
	# Index array for faster searching
	#for($pt_ctr = 0; $pt_ctr < $num_xy; $pt_ctr++)
	for($pt_ctr = 0; $pt_ctr <= $#mz_array; $pt_ctr++)
	{
		$mz_val = $mz_array[$pt_ctr];

		$mz_int = sprintf("%d", $mz_val);
		$array_index{$mz_int} = $pt_ctr;
	}

    $ps_t2 = time();
    
    # Use this to keep track of which spectra have been output in this round
    @spec_output_array = ();
    
	#for($spec_ctr = 0; $spec_ctr < $num_good_spectra; $spec_ctr++)
    # Change to <= $#spectra as we are going to start punting these
	for($spec_ctr = 0; $spec_ctr <= $#spectra; $spec_ctr++)
	{
		#roprint("$spectra[$spec_ctr]{'mz_min'} -> $spectra[$spec_ctr]{'mz_max'}");
		#roprint("$spectra[$spec_ctr]{'rt_min'} -> $spectra[$spec_ctr]{'rt_max'}");
	
        # RT has exceeded spectrum, output and purge it
        if($rt > $spectra[$spec_ctr]{"rt_max_plot"})
        {
            $spec_output_ctr = $spec_ctr;
            output_spectrum();
            $spectra[$spec_ctr]{"output"} = 1;
            push(@spec_output_array, $spec_ctr);
			next;
        }
        
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
		
			#print "1s->e $startloc $endloc\n";
		
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
			if($spectra3D)
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

			#print "2s->e $startloc $endloc\n";

			$num_spectra_mz = 0;
			$sum_spectra_mz = 0;

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
					$num_spectra_mz++;
					$sum_spectra_mz += $int_val;
					#roprint("MZ_YES");
					#$xy_val = join ' ', $mz_val, $int_val;
					#push @{ $spectra_values[$spec_ctr] }, $xy_val;
				}
			}
			
			# RAW files tend to simply not output a lot of m/z values when int=0
			# Since we are really summing intensity (don't let the variable names fool you, see above $sum_spectra_mz += $int_val), just set $avg_spectra_mz to 0
			if($num_spectra_mz == 0)
			{
				#roprint("$spectra[$spec_ctr]{'rt_min'} - $spectra[$spec_ctr]{'rt_max'} * $spectra[$spec_ctr]{'mz_min'} - $spectra[$spec_ctr]{'mz_max'}");
				#roprint("RT: $rt");
				#roprint("MZ($#mz_array): $mz_array[0] -> $mz_array[$#mz_array]");
				#roprint("INT($#int_array): $int_array[0] -> $int_array[$#int_array]");
				#roprint("Start/End $startloc $endloc");

				#for($xy_ctr = $startloc; $xy_ctr < $endloc; $xy_ctr++)
				#{
				#	roprint("MZ,INT: $mz_array[$xy_ctr] $int_array[$xy_ctr]");
				#}
				
				$avg_spectra_mz = 0;
			}
			# For RAW files, we can't actually average, since some RT slices may report 0s, some may report nothing, others may report numbers
			# Let's just take the sum, rather than the average, since the actual intensity is arbitrary
			else
			{
				#if($whoami eq "nobody")
				#{
					$avg_spectra_mz = $sum_spectra_mz;
				#}
				#else
				#{
				#	$avg_spectra_mz = $sum_spectra_mz / $num_spectra_mz;
				#}
			}
			
			$xy_val = join ' ', $rt, $avg_spectra_mz;
			push(@{$spectra_values[$spec_ctr]}, $xy_val);
			
			if($spectra3D)
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

			#print "3s->e $startloc $endloc\n";
		
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
    
    # Eliminate the spectra that were output
    # Because the spectra are sorted on RT, this array will always be a contiguous list of number starting at 0.  However, just in case let's remove the spliced elements one by one starting from the end.
    # It's also possible to be an empty array;
    # Also need to splice @spectra_values_plot, which is created in parallel
    while($#spec_output_array >= 0)
    {
        $splice_offset = pop(@spec_output_array);
        splice(@spectra, $splice_offset, 1);
        splice(@spectra_values, $splice_offset, 1);
        splice(@spectra_values_plot, $splice_offset, 1);
    }
    
    $ps_t3 = time();
}

###########
# roprint #
###########
sub roprint
{
	my @vars = @_;
	my $outstring = join "\n", @vars;
	
    #print "$outstring\n";
	`echo "$outstring" | cat - >> $rofile`;
}

###################
# output_spectrum #
###################
# This used to be part of a loop on $spec_ctr outputting everything at the end
# for($spec_ctr = 0; $spec_ctr < $num_good_spectra; $spec_ctr++)
# Now factor this code out so that individual spectra can be output at any stage, freeing up memory and eliminating them from processing for individual slices
# This used to require the variable $spec_ctr to be set properly
# Replace this with $spec_output_ctr
sub output_spectrum
{
	$spectra_filename = join '.', $spectra[$spec_output_ctr]{"key"}, "rt";
	$newcon_filename = join '.', $spectra[$spec_output_ctr]{"key"}, "newcon";
	open OUT, ">$dir$spectra_filename" or die "Can't open $dir$spectra_filename\n";

	# 7/22/11
    # Maximum observed intensity for a single point in the summed RT spectrum.
	# This is now going to be used to filter out a bunch of the crappy RT fits
	# Basically get $max_int for all sequencing attempts
	# Find max_max_int (max of all max_int values for that ion), and if max_int < (1/2)max_max_int, reject it
	$max_int = 0; 
    
	for($ctr = 0; $ctr <= $#{$spectra_values[$spec_output_ctr]}; $ctr++)
	{
		@array = split ' ', $spectra_values[$spec_output_ctr][$ctr];
		if($array[1] > $max_int)
		{
			$max_int = $array[1];
		}
		print OUT "$spectra_values[$spec_output_ctr][$ctr]\n";
	}
	
	close OUT;
	
	print OUTFILE "$spectra[$spec_output_ctr]{'key'} $spectra[$spec_output_ctr]{'rt'} $max_int\n";
	
	open NEWCONTOUR, ">$dir$newcon_filename" or die "Can't open $dir$newcon_filename\n";
    
	$nc_mz_min = $spectra[$spec_output_ctr]{"mz_min"};
	$nc_mz_max = $spectra[$spec_output_ctr]{"mz_max"};
	$nc_rt_min = $spectra[$spec_output_ctr]{"rt_min"};
	$nc_rt_max = $spectra[$spec_output_ctr]{"rt_max"};
	$nc_mz_min_plot = $spectra[$spec_output_ctr]{"mz_min_plot"};
	$nc_mz_max_plot = $spectra[$spec_output_ctr]{"mz_max_plot"};
	$nc_rt_min_plot = $spectra[$spec_output_ctr]{"rt_min_plot"};
	$nc_rt_max_plot = $spectra[$spec_output_ctr]{"rt_max_plot"};
    
	$mzshift = sprintf("%.1f", $nc_mz_min_plot);
	$mzshift = 10*$mzshift;
	$mzshift_max = sprintf("%.1f", $nc_mz_max_plot);
	$mzshift_max = 10*$mzshift_max;
    
	print NEWCONTOUR "$nc_mz_min $nc_mz_max $nc_rt_min $nc_rt_max\n";
	print NEWCONTOUR "$nc_mz_min_plot $nc_mz_max_plot $nc_rt_min_plot $nc_rt_max_plot\n";
    
	@ncdata = ();
    
	for($local_slice_ctr = 0; $local_slice_ctr < $spectra[$spec_output_ctr]{"slices_plot"}; $local_slice_ctr++)
	{
		$y = $local_slice_ctr;
        
		@plot_values = split ' ', $spectra_values_plot[$spec_output_ctr][$local_slice_ctr];
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
	
	printf(NEWCONTOUR "%d\n", $spectra[$spec_output_ctr]{"slices_plot"}); # Number of rt values
	printf(NEWCONTOUR "%s\n", join ' ', @{$spectra[$spec_output_ctr]{"rt_vals"}});
	
	for($nc_ctr1 = 0; $nc_ctr1 < $spectra[$spec_output_ctr]{"slices_plot"}; $nc_ctr1++)
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
}

sub by_rt
{
    $a->{"rt"} <=> $b->{"rt"};
}
