#!/usr/bin/perl

use Getopt::Long;
use Text::ParseWords;

&GetOptions("a=s" => \$filea, "b=s" => \$fileb);

@array = split /\./, $filea;
$suba = $array[0];
@array = split /\./, $fileb;
$subb = $array[0];

open FILEA, "$filea" or die "Can't open $filea\n";
while(defined($input=<FILEA>))
{
	chomp($input);
	@array = quotewords(",", 0, $input);
	
	%name2col = ();
	for($ctr = 0; $ctr <= $#array; $ctr++)
	{
		$name2col{$array[$ctr]} = $ctr;
	}
	
	if(defined($name2col{"prot_acc"}) && defined($name2col{"prot_desc"}))
	{
		while(defined($line=<FILEA>))
		{
			chomp($line);
			@array2 = quotewords(",", 0, $line);
			
			$acc = $array2[$name2col{"prot_acc"}];
			$desc = $array2[$name2col{"prot_desc"}];
			
			if($acc ne "" && $desc ne "")
			{
				$data_a{$acc} = join ' ', $acc, $desc;
				$data{$acc} = 1;
			}
		}
	}
}
close FILEA;

open FILEB, "$fileb" or die "Can't open $fileb\n";
while(defined($input=<FILEB>))
{
	chomp($input);
	@array = quotewords(",", 0, $input);
	
	%name2col = ();
	for($ctr = 0; $ctr <= $#array; $ctr++)
	{
		$name2col{$array[$ctr]} = $ctr;
	}
	
	if(defined($name2col{"prot_acc"}) && defined($name2col{"prot_desc"}))
	{
		while(defined($line=<FILEB>))
		{
			chomp($line);
			@array2 = quotewords(",", 0, $line);
			
			$acc = $array2[$name2col{"prot_acc"}];
			$desc = $array2[$name2col{"prot_desc"}];
			
			if($acc ne "" && $desc ne "")
			{
				$data_b{$acc} = join ' ', $acc, $desc;
				$data{$acc} = 1;
			}
		}
	}
}
close FILEA;

$count = 0;

foreach $acc (keys %data)
{
	if(defined($data_a{$acc}) && defined($data_b{$acc}))
	{
	
	}
	elsif(defined($data_a{$acc}))
	{
		$line = join ' ', "Only", $suba, $data_a{$acc};
		push(@{$output[$count]}, $line);
		#print "Only $suba $data_a{$acc}\n";

		$string = $data_a{$acc};
		$search_string1 = get_search_string($string);
		
		foreach $acc_b (keys %data_b)
		{
			$string = $data_b{$acc_b};
			$search_string2 = get_search_string($string);;
					
			if($search_string1 eq $search_string2)
			{
				$line = join ' ', "\t ->", $subb, $data_b{$acc_b};
				push(@{$output[$count]}, $line);
				#print "    -> $subb $data_b{$acc_b}\n";
			}
		}

		$count++;
	}
	elsif(defined($data_b{$acc}))
	{
		$line = join ' ', "Only", $subb, $data_b{$acc};
		push(@{$output[$count]}, $line);
		#print "Only $subb $data_b{$acc}\n";

		$string = $data_b{$acc};
		$search_string1 = get_search_string($string);
		
		foreach $acc_a (keys %data_a)
		{
			$string = $data_a{$acc_a};
			$search_string2 = get_search_string($string);
			
			if($search_string1 eq $search_string2)
			{
				$line = join ' ', "\t ->", $suba, $data_a{$acc_a};
				push(@{$output[$count]}, $line);
				#print "    -> $suba $data_a{$acc_a}\n";
			}
		}

		$count++;
	}
	else
	{
		print "WTF?!\n";
	}
}

# Sort @output by number of entries
for($ctr = 0; $ctr <= $#output; $ctr++)
{
	$output2[$ctr] = join "\n", @{$output[$ctr]};
}

@sort = sort by_newline @output2;

sub by_newline
{
	@aa = split /\n/, $a;
	@bb = split /\n/, $b;
	
	$#aa <=> $#bb;
}

for($ctr = 0; $ctr <= $#sort; $ctr++)
{
	print "$sort[$ctr]\n";
	print "----------\n";
}

sub get_search_string
{
	my @vars = @_;
	my $thestring = shift(@vars);
	my @array = split ' ', $thestring;
	my $junk = shift(@array);
	
	$new_string = "";
	
	while($#array >= 0)
	{
		if(substr($array[0], 0, 3) eq "OS=")
		{
			last;
		}
		else
		{
			$new_string = join ' ', $new_string, shift(@array);
		}
	}
	
	my @array2 = split ' ', $new_string;
		
	if($array2[$#array2] =~ m/^\w$/)
	{
		pop(@array2);
	}
		
	my $return_string = join ' ', @array2;
	
	return $return_string;
}
