#!/usr/bin/perl

#
# Take aa_elem.txt as an input, and from this calculate the masses of each amino acid residue
# elem2mass.pl < aa_elem.txt > aa_mass.txt
#
# see 03_calc_mass.pl for masses

$c_mass = 12.0000000;
$h_mass = 1.00782503207;
$n_mass = 14.0030740048;
$o_mass = 15.99491461956;
$s_mass = 31.97207100; 
$p_mass = 30.97376163;

$n15_mass = 15.0001088982;

@c_m = (12, 13.0033548378);
@h_m = (1.0078250320735, 2.0141017778);
@n_m = (14.0030740048, 15.0001088982);
@o_m = (15.99491461956, 16.99913170, 17.9991610);
@s_m = (31.97207100, 32.97145876, 33.96786690, 35.96708076);
@p_m = (30.97376163);

@c_a = (0.9893, 0.0107);
@h_a = (0.999885, 0.000115);
@n_a = (0.99636, 0.00364);
@o_a = (0.99757, 0.00038, 0.00205);
@s_a = (0.9499, 0.0075, 0.0425, 0.0001);
@p_a = (1.0);

@n993_a = (0.007, 0.993);

$c_avg = 0;
for($ctr = 0; $ctr <= $#c_m; $ctr++)
{
	$c_avg += $c_m[$ctr]*$c_a[$ctr];
}
#print "C $c_avg\n";

$h_avg = 0;
for($ctr = 0; $ctr <= $#h_m; $ctr++)
{
	$h_avg += $h_m[$ctr]*$h_a[$ctr];
}
#print "H $h_avg\n";

$n_avg = 0;
for($ctr = 0; $ctr <= $#n_m; $ctr++)
{
	$n_avg += $n_m[$ctr]*$n_a[$ctr];
}
#print "N $n_avg\n";

$o_avg = 0;
for($ctr = 0; $ctr <= $#o_m; $ctr++)
{
	$o_avg += $o_m[$ctr]*$o_a[$ctr];
}
#print "O $o_avg\n";

$s_avg = 0;
for($ctr = 0; $ctr <= $#s_m; $ctr++)
{
	$s_avg += $s_m[$ctr]*$s_a[$ctr];
}
#print "S $s_avg\n";

$p_avg = 0;
for($ctr = 0; $ctr <= $#p_m; $ctr++)
{
	$p_avg += $p_m[$ctr]*$p_a[$ctr];
}
#print "P $p_avg\n";

$n993_avg = 0;
for($ctr = 0; $ctr <= $#n_m; $ctr++)
{
	$n993_avg += $n_m[$ctr]*$n993_a[$ctr];
}
#print "N993 $n993_avg\n";

$headerline = <STDIN>;

print "id n14mass n15mass n14avgmass n993avgmass\n";

while(defined($input=<STDIN>))
{
	chomp($input);
	@array = split ' ', $input;

	$mass = $array[1]*$c_mass + $array[2]*$h_mass + $array[3]*$n_mass + $array[4]*$o_mass + $array[5]*$s_mass + $array[6]*$p_mass;

	$mass_n15 = $array[1]*$c_mass + $array[2]*$h_mass + $array[3]*$n15_mass + $array[4]*$o_mass + $array[5]*$s_mass + $array[6]*$p_mass;

	$mass_avg = $array[1]*$c_avg + $array[2]*$h_avg + $array[3]*$n_avg + $array[4]*$o_avg + $array[5]*$s_avg + $array[6]*$p_avg;

	$mass_n993 = $array[1]*$c_avg + $array[2]*$h_avg + $array[3]*$n993_avg + $array[4]*$o_avg + $array[5]*$s_avg + $array[6]*$p_avg;

	# The iodoacetamide based modification of cysteine happens after the fact, and thus does not include labeled nitrogen but rather 14N nitrogen
	# So we have a special exception here in the calculation of the "c" modification
	if($array[0] eq "c")
	{
		print "$array[0] $mass $mass $mass_avg $mass_n993\n";		
	}
	else
	{
		print "$array[0] $mass $mass_n15 $mass_avg $mass_n993\n";		
	}
}
