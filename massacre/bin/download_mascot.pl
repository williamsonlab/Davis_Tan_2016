#!/usr/bin/perl

#http://snowwhite.scripps.edu/mascot/cgi/export_dat_2.pl?file=..%2Fdata%2F20101119%2FF014264.dat&do_export=1&prot_hit_num=1&prot_acc=1&pep_query=1&pep_rank=1&pep_isbold=1&pep_exp_mz=1&_showallfromerrortolerant=0&_onlyerrortolerant=0&_noerrortolerant=0&_show_decoy_report=0&export_format=CSV&_sigthreshold=0.05&REPORT=AUTO&_server_mudpit_switch=0.000000001&_ignoreionsscorebelow=35&_showsubsets=0&_requireboldred=1&search_master=1&show_header=1&show_decoy=1&show_mods=1&show_params=1&show_format=1&protein_master=1&prot_score=1&prot_desc=1&prot_mass=1&prot_matches=1&peptide_master=1&pep_exp_mr=1&pep_exp_z=1&pep_calc_mr=1&pep_delta=1&pep_start=1&pep_end=1&pep_miss=1&pep_score=1&pep_expect=1&pep_seq=1&pep_var_mod=1&pep_scan_title=1

use LWP::Simple;
use Getopt::Long;

#@mascot_list = "http://snowwhite.scripps.edu/mascot/cgi/master_results.pl?file=../data/20101119/F014264.dat";

$query_base = "http://snowwhite.scripps.edu/mascot/cgi/export_dat_2.pl?";

#$query_param = "do_export=1&prot_hit_num=1&prot_acc=1&pep_query=1&pep_rank=1&pep_isbold=1&pep_exp_mz=1&_showallfromerrortolerant=0&_onlyerrortolerant=0&_noerrortolerant=0&_show_decoy_report=0&export_format=CSV&_sigthreshold=0.05&REPORT=AUTO&_server_mudpit_switch=0.000000001&_ignoreionsscorebelow=35&_showsubsets=0&_requireboldred=1&search_master=1&show_header=1&show_decoy=1&show_mods=1&show_params=1&show_format=1&protein_master=1&prot_score=1&prot_desc=1&prot_mass=1&prot_matches=1&peptide_master=1&pep_exp_mr=1&pep_exp_z=1&pep_calc_mr=1&pep_delta=1&pep_start=1&pep_end=1&pep_miss=1&pep_score=1&pep_expect=1&pep_seq=1&pep_var_mod=1&pep_scan_title=1";

$query_param = "do_export=1&prot_hit_num=1&prot_acc=1&pep_query=1&pep_rank=1&pep_isbold=1&pep_exp_mz=1&_showallfromerrortolerant=0&_onlyerrortolerant=0&_noerrortolerant=0&_show_decoy_report=0&export_format=CSV&_sigthreshold=0.05&REPORT=AUTO&_server_mudpit_switch=0.000000001&_showsubsets=0&_requireboldred=1&search_master=1&show_header=1&show_decoy=1&show_mods=1&show_params=1&show_format=1&protein_master=1&prot_score=1&prot_desc=1&prot_mass=1&prot_matches=1&peptide_master=1&pep_exp_mr=1&pep_exp_z=1&pep_calc_mr=1&pep_delta=1&pep_start=1&pep_end=1&pep_miss=1&pep_score=1&pep_expect=1&pep_seq=1&pep_var_mod=1&pep_scan_title=1";

# _ignoreionsscorebelow=35&

$threshold = 20;

&GetOptions("threshold=i" => \$threshold);

print "Downloading mascot results with a threshold of $threshold\n";

while(defined($input=<STDIN>))
{
	chomp($input);
	#print "$input\n";;
	@array = split ' ', $input;
	$mascot_path = shift(@array);
	if($#array >= 0)
	{
		$output_file = shift(@array);
	}
	else
	{
		@array2 = split /\//, $mascot_path;
		$output_file = $array2[$#array2];
		$output_file = join '', $output_file, ".csv";
	}

	print "$mascot_path -> $output_file\n";

	if($mascot_path =~ /file=(.*)/)
	{
		$file = $1;
		#print "\t$file\n";
	}
	else
	{
		print "$mascot_path does not appear to be a valid path\n";
	}
	
	$query = join '', $query_base, $query_param, "&_ignoreionsscorebelow=", $threshold, "&file=", $file;
	#print "$query\n";
	
	$result = get($query);
	open OUT, ">$output_file" or die "Can't open $output_file\n";
	print OUT "$result";
	close OUT;
}