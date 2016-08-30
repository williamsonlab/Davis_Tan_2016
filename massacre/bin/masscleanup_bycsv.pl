#!/usr/bin/perl

#
# 1) Takes a .csv file as input
# 2) From the .csv file, identifies the _peaks (input) and _fits (output) directories
#		-probably the same, but no guarantee
# 3) Keep:
#		a) _peaks/basename.txt
#		b) _fits/basename.fit.png
#		c) _fits/basename.contour.png
#

$pid = $$;

$header = <STDIN>;
chomp($header);
@array = split /\,/, $header;
for($ctr = 0; $ctr <= $#array; $ctr++)
{
	$name2col{$array[$ctr]} = $ctr;
}

$flag1 = 0;
$flag2 = 0;

# Used to only have a "file" record
# This has been renamed "isofile" -> base name, with no extension
# "file" now specifies the actual ".txt" file, including extension
if(defined($name2col{"isofile"}))
{
	$flag1 = 1;
}
if(defined($name2col{"file"}))
{
	$flag2 = 1;
}

%dirs = ();
%newdirs = ();

print "Reading .csv file\n";
while(defined($input=<STDIN>))
{
	chomp($input);
	@array = split /\,/, $input;
	
	$isofile = "";
	$file = "";

	if($flag1 == 1)
	{
		$isofile = $array[ $name2col{"isofile"} ];
		@isofile_array = split /\//, $isofile;
		$isofile_dir = $isofile_array[0]; # possible additional base directory for fits
		if($isofile_dir eq "")
		{
			die "Error, isofile entry broken\n";
		}
		$dirs{$isofile_dir} = 1;
		$new_isofile_dir = join '.', $isofile_dir, "tempclean", $pid;
		$newdirs{$new_isofile_dir} = 1;
		
		if(!(-e $new_isofile_dir))
		{
			`mkdir $new_isofile_dir`;
		}
	}
	
	if($flag2 == 1)
	{
		$file = $array[ $name2col{"file"} ];
		@file_array = split /\//, $file;
		$file_dir = $file_array[0]; # base directory for peaks
		if($file_dir eq "")
		{
			die "Error, file entry broken\n";
		}
		$dirs{$file_dir} = 1;
		$new_file_dir = join '.', $file_dir, "tempclean", $pid;
		$newdirs{$new_file_dir} = 1;
		
		if(!(-e $new_file_dir))
		{
			`mkdir $new_file_dir`;
		}
	} 
	
	if($flag1 == 0 && $flag2 == 0)
	{
		die "Error, no isofile or file entry detected in your .csv\n";
	}
	
	
	if($flag1 == 1 && $flag2 == 1)
	{
		if(substr($file_array[1], -4, 4) eq ".txt")
		{
			substr($file_array[1], -4, 4) = "";
		}
		if(substr($isofile_array[1], -4, 4) eq ".txt")
		{
			substr($isofile_array[1], -4, 4) = "";
		}
	
		$txtfile = join '', $file_array[1], ".txt";
		$confile = join '', $file_array[1], ".contour.png";
		$fitfile = join '', $isofile_array[1], ".fit.png";
		
		`mv $file_dir/$txtfile $new_file_dir/$txtfile`;
		`mv $file_dir/$confile $new_file_dir/$confile`;
		`mv $isofile_dir/$fitfile $new_isofile_dir/$fitfile`;
	}
	elsif($flag1 == 1)
	{
		if(substr($isofile_array[1], -4, 4) eq ".txt")
		{
			substr($isofile_array[1], -4, 4) = "";
		}
	
		$txtfile = join '', $isofile_array[1], ".txt";
		$confile = join '', $isofile_array[1], ".contour.png";
		$fitfile = join '', $isofile_array[1], ".fit.png";
		
		`mv $isofile_dir/$txtfile $new_isofile_dir/$txtfile`;
		`mv $isofile_dir/$confile $new_isofile_dir/$confile`;
		`mv $isofile_dir/$fitfile $new_isofile_dir/$fitfile`;
	}
	elsif($flag2 == 1)
	{
		if(substr($file_array[1], -4, 4) eq ".txt")
		{
			substr($file_array[1], -4, 4) = "";
		}
		
		$txtfile = join '', $file_array[1], ".txt";
		$confile = join '', $file_array[1], ".contour.png";
		$fitfile = join '', $file_array[1], ".fit.png";

		`mv $file_dir/$txtfile $new_file_dir/$txtfile`;
		`mv $file_dir/$confile $new_file_dir/$confile`;
		`mv $file_dir/$fitfile $new_file_dir/$fitfile`;
	}
}

foreach $dir (keys %dirs)
{
	$deldir = join '', $dir, ".DELETEME";
	print "Moving discarded files to $deldir\n";
	`mv $dir $deldir`;
}

foreach $dir (keys %newdirs)
{
	@array = split /\./, $dir;
	pop(@array); # pop off $pid
	pop(@array); # pop off tempclean
	
	$freshdir = join '.', @array;
	print "Keeping good files in $freshdir\n";
	`mv $dir $freshdir`;
}

print "Done\n";