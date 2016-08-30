###############
# SUBROUTINES #
###############

#
# readfile($rf_input_file)
# reads files in both Unix or old-style Mac format
# returns the file as an array, with each line as a different array element
# newlines at end of lines are removed
#
sub readfile
{
	my @rf_parameters = @_;
	my @rf_array1;
	my @rf_array2;
	my $rf_ctr1;
	my $rf_ctr2;
	my $rf_input_file;
	my @rf_output_array = ();	
	my $rf_inputline;
	
	if($#rf_parameters != 0)
	{
		die "Error: Subroutine \"readfile()\" accepts only one parameter\n";
	}

	$rf_input_file = $rf_parameters[0];
	
	open INPUT, "$rf_input_file" or die "Error: Can't open $rf_input_file\n";
	while(defined($rf_inputline=<INPUT>))
	{
		push(@rf_array1, $rf_inputline);
	}
	#@rf_array1 = <INPUT>;
	close INPUT;

	for ($rf_ctr1 = 0; $rf_ctr1 <= $#rf_array1; $rf_ctr1++)
	{
		chomp($rf_array1[$rf_ctr1]);
		@rf_array2 = split /\r/, $rf_array1[$rf_ctr1];
		for($rf_ctr2 = 0; $rf_ctr2 <= $#rf_array2; $rf_ctr2++)
		{
			push(@rf_output_array, $rf_array2[$rf_ctr2]);
		}
	}
	
	return @rf_output_array;
}

sub stripcomments
{
	my @oldarray = @_;
	my @newarray = ();
	my $ctr;
	
	for($ctr = 0; $ctr <= $#oldarray; $ctr++)
	{
		if(substr($oldarray[$ctr], 0, 1) eq "#")
		{
			next;
		}
		else
		{
			push(@newarray, $oldarray[$ctr]);
		}
	}

	return @newarray;
}

1;