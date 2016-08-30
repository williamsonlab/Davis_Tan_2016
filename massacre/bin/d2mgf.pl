#!/usr/bin/perl

$pid = $$;
$who = `whoami`;
chomp($who);

$dir = join '_', $who, $pid;

@ls_d = `ls -d *.d`;
chomp(@ls_d);

@ls_raw = `ls *.raw *.RAW`;
chomp(@ls_raw);

push(@ls_d, @ls_raw);

for($ctr_d = 0; $ctr_d <= $#ls_d; $ctr_d++)
{
	print "$ctr_d: $ls_d[$ctr_d]\n";
	
	@base_array = split /\./, $ls_d[$ctr_d];
	pop(@base_array);
	$base = join '.', @base_array;	
	#$base = $base_array[0];
	
	$mgf = join '.', $base, "mgf";
	
	if(!(-e $mgf))
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
		print "Generating MGF file\n";
		#$out3 = `ssh mass\@vaginicola.scripps.edu "wine /usr/local/bin/pwiz/msconvert.exe $dir/$ls_d[$ctr_d] --filter 'peakPicking true 1-' --mgf --o $dir"`;
		$out3 = `ssh -p 2222 mass\@vaginicola.scripps.edu "pwiz/msconvert.exe //VBOXSVR/vbox_ms/$dir/$ls_d[$ctr_d] --filter 'peakPicking true 1-' --mgf --o //VBOXSVR/vbox_ms/$dir"`;
		# 9/26/11 - Updated to include new settings via sschen
		# deisotope, counts > 25, top 100 most intense peaks
		$out3 = `ssh -p 2222 mass\@vaginicola.scripps.edu "pwiz/msconvert.exe //VBOXSVR/vbox_ms/$dir/$ls_d[$ctr_d] --filter 'peakPicking true 1-' --mgf --filter MS2Deisotope --filter 'threshold absolute 25 most-intense' --filter 'threshold count 100 most-intense' --o //VBOXSVR/vbox_ms/$dir"`;
		print "$out3\n";	
		print "Done generating MGF file\n";
		
		#
		print "Copying MGF file back\n";
		$out4 = `scp mass\@vaginicola.scripps.edu:vbox_ms/$dir/$mgf .`;
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
		print "Skipping, $mgf already exists\n";
	}
}

