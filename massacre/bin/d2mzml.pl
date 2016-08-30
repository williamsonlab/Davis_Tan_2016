#!/usr/bin/perl

$pid = $$;
$who = `whoami`;
chomp($who);

$dir = join '_', $who, $pid;

@ls_d = `ls -d *.d`;
chomp(@ls_d);

@ls_raw = `ls *.raw *.RAW`;
chomp(@ls_raw);

#@ls_mzml = `ls *.mzml *.mzML`;
#push(@ls_d, @ls_raw, @ls_mzml);

push(@ls_d, @ls_raw);

for($ctr_d = 0; $ctr_d <= $#ls_d; $ctr_d++)
{
	print "$ctr_d: $ls_d[$ctr_d]\n";
	
	@base_array = split /\./, $ls_d[$ctr_d];
	pop(@base_array);
	$base = join '.', @base_array;	
	#$base = $base_array[0];
	
	$mzml = join '.', $base, "mzML", "gz";
	
	if(!(-e $mzml))
	{
		#
		print "Making $dir directory\n";
		$out1 = `ssh mass\@vaginicola.scripps.edu "mkdir vbox_ms/$dir"`;		
		print "$out1\n";
		
		#
		print "Copying $ls_d[$ctr_d]\n";
		$out2 = `scp -r $ls_d[$ctr_d] mass\@vaginicola.scripps.edu:vbox_ms/$dir/`;
		print "$out2\n";
		print "Done copying\n";
		
		#
		print "Generating mzML file\n";
		$out3 = `ssh -p 2222 mass\@vaginicola.scripps.edu "pwiz/msconvert.exe //VBOXSVR/vbox_ms/$dir/$ls_d[$ctr_d] --filter 'msLevel 1' --mz32 --inten32 --filter 'mzWindow [400,2000]' --o //VBOXSVR/vbox_ms/$dir --gzip"`;
		print "$out3\n";
		print "Done generating mzML file\n";
		
		#
		print "Copying mzML file back\n";
		$out4 = `scp mass\@vaginicola.scripps.edu:vbox_ms/$dir/$mzml .`;
		print "$out4\n";
		print "Done Copying\n";
		
		#
		print "Cleaning up\n";
		$out5 = `ssh mass\@vaginicola.scripps.edu "rm -r vbox_ms/$dir/"`;
		print "$out5\n";
		print "Done cleaning up\n";
	}
	else
	{
		print "Skipping, $mzml already exists\n";
	}
}

