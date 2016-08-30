#!/usr/bin/perl

use Getopt::Long;
use Text::ParseWords;

require "$ENV{'MASSACRE_PATH'}/msc/mikevars.pl";
require "$ENV{'MASSACRE_PATH'}/msc/mikesubs.pl";

# %acc2pro_all - use this to map accession numbers to protein names at the end

&GetOptions("plist=s" => \$plist, "mlist=s" => \$mlistfile);

print "plist = $plist, mlist = $mlistfile\n";

#####
# Process list of proteins and get out an array of accession numbers @proteins
#####
if($plist eq "")
{
	die "Error, must define --plist=something\n";
}
elsif(lc($plist) eq "human")
{
	$plistfile = "$ENV{'MASSACRE_PATH'}/msc/plist/human_ribosome.txt";
}
elsif(lc($plist) eq "ecoli")
{
    $plistfile = "$ENV{'MASSACRE_PATH'}/msc/plist/ecoli_ribosome.txt";
}
elsif(lc($plist) eq "yeast")
{
    $plistfile = "$ENV{'MASSACRE_PATH'}/msc/plist/yeast_ribosome.txt";
}
else
{
	print "Do not recognize plist \"$plist\", setting as file path\n";
	$plistfile = $plist;
}

#####
# Process list of proteins into an array
#####
@plistlines = readfile($plistfile);
@plistlines = stripcomments(@plistlines);

@proteins = ();

for($ctr = 0; $ctr <= $#plistlines; $ctr++)
{
	@array = quotewords(" ", 0, $plistlines[$ctr]);
	push(@proteins, $array[0]);
    
    if($array[1] eq "=>")
    {
        $prot_name = $array[2];
    }
    else
    {
        $prot_name = $array[1];
    }
    
    # Build up a local mapping of accession -> protein
    # This way you can build a custom list of proteins that includes some not present in mikevars.pl
    $local_acc2pro{$array[0]} = $prot_name;
}

#####
# Process list of mascot files and get out an array of file names @mfiles
# Also map file names to real names %mfile2name
# Also store real names in @mnames
#####
if($mlistfile eq "")
{
	die "Error, must define --mlist=something\n";
}
else
{
	@mlistarray = readfile($mlistfile);
}

@mfiles = ();
@mnames = ();

for($ctr = 0; $ctr <= $#mlistarray; $ctr++)
{
	@array = split ' ', $mlistarray[$ctr];
	if($#array < 0)
	{
		die "Error, blank line in list of mascot files? |$mlistarray[$ctr]|\n";
	}
	
	push(@mfiles, $array[0]);
	
	if($#array == 0)
	{
		$mfile2name{$array[0]} = $array[0];
		push(@mnames, $array[0]);
	}
	else
	{
		$mfile2name{$array[0]} = $array[1];
		push(@mnames, $array[1]);
	}
}

#####
# Process each mascot file
#####

for($mctr = 0; $mctr <= $#mfiles; $mctr++)
{
	@filearray = readfile($mfiles[$mctr]);

	while(1)
	{
		my $line = shift(@filearray);
		my @array = quotewords(",", 0, $line);
		
		if($#filearray == 0 || $#filearray == -1)
		{
			print "ERROR: |$mfiles[$ctr]| does not appear to be a mascot export file\n";
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
	
	for($ctr = 0; $ctr <= $#filearray; $ctr++)
	{
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
			
			$mascot_count[$mctr]{$mpro} = $array[$mascotname2col{"prot_matches"}]
		}
	}
}

#####
# Print header line
#####

push(@header, "X");
push(@header, @mnames);

$header_line = join ',', @header;
print "$header_line\n";

#####
# Print individual protein lines
#####
for($ctr = 0; $ctr <= $#proteins; $ctr++)
{
    # Get name for current protein
    # Give priority to the global acc2pro from mikevars.pl, then try local_acc2pro, then just use the accession number
    if(defined($acc2pro_all{$proteins[$ctr]}))
    {
        $current_protein = $acc2pro_all{$proteins[$ctr]};
    }
    elsif(defined($local_acc2pro{$proteins[$ctr]}))
    {
        $current_protein = $local_acc2pro{$proteins[$ctr]};
    }
    else
    {
        $current_protein = $proteins[$ctr];
    }
    
    # Get name for previous protein
    if(defined($acc2pro_all{$proteins[$ctr-1]}))
    {
        $previous_protein = $acc2pro_all{$proteins[$ctr-1]};
    }
    elsif(defined($local_acc2pro{$proteins[$ctr-1]}))
    {
        $previous_protein = $local_acc2pro{$proteins[$ctr-1]};
    }
    else
    {
        $previous_protein = $proteins[$ctr-1];
    }
    
    # If the current line points to the same protein as the previous line, add in the values rather than starting a new entry
    # This is useful for something like yeast, where isoforms are treated as a single protein (L1A and L1B -> L1X)
    $same_as_previous = 0;
    #if(defined($acc2pro_all{$proteins[$ctr]}) && defined($acc2pro_all{$proteins[$ctr-1]}))
    #{
    #    if($acc2pro_all{$proteins[$ctr]} eq $acc2pro_all{$proteins[$ctr-1]})
    #    {
    #        $same_as_previous = 1;
    #    }
    #}
    
    if($current_protein eq $previous_protein)
    {
        $same_as_previous = 1;
    }
    
    if($same_as_previous)
    {
        for($mctr = 0; $mctr <= $#mfiles; $mctr++)
        {
            # Increment $line[$mctr+1] by the appropriate amount
            # use $mctr+1 since the first array element is the protein name
            if(defined($mascot_count[$mctr]{$proteins[$ctr]}))
            {
                $line[$mctr+1] += $mascot_count[$mctr]{$proteins[$ctr]};
            }
        }
    }
    else
    {
        # Print the line from the previous entry since we are moving onto a new protein
        if($ctr > 0)
        {
            $lineout = join ',', @line;
            print "$lineout\n";
        }
            
        @line = ();
        #if(defined($acc2pro_all{$proteins[$ctr]}))
        #{
        #    push(@line, $acc2pro_all{$proteins[$ctr]});
        #}
        #else
        #{
        #    push(@line, $proteins[$ctr]);
        #}
        push(@line, $current_protein);
        
        for($mctr = 0; $mctr <= $#mfiles; $mctr++)
        {
            if(defined($mascot_count[$mctr]{$proteins[$ctr]}))
            {
                push(@line, $mascot_count[$mctr]{$proteins[$ctr]});
            }
            else
            {
                push(@line, 0);
            }
        }
	}
}

# Print the line from the final entry
$lineout = join ',', @line;
print "$lineout\n";
