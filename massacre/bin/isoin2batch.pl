#!/usr/bin/perl

$header = <STDIN>;
chomp($header);
@array = split ' ', $header;
for($ctr = 0; $ctr <= $#array; $ctr++)
{
	$name2col{$array[$ctr]} = $ctr;
}

while(defined($input=<STDIN>))
{
	chomp($input);
	@array = split ' ', $input;
	$seq = $array[$name2col{"seqmod"}];
	$z = $array[$name2col{"charge"}];
	$file = $array[$name2col{"file"}];
	
	print "$seq $z $file\n";
}