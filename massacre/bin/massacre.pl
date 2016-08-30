#!/usr/bin/perl

use strict;
require Tk;
require Tk::ROText;
require Tk::X11Font;
use Tk;
use Getopt::Long;

#######################
# Operation Variables #
#######################
my $program_status = 0;

##########################
# Command Line Variables #
##########################
my $go = -1;
my $paramfile_input = "";
my $output_root_input = "";
# Other variables which can be altered on the command line but are defined elsewhere
# $input_digest


my @datearray = ();
@datearray = `ls -ltT $ENV{'MASSACRE_PATH'}/msc/msc*.pl $ENV{'MASSACRE_PATH'}/bin/massacre.pl`;
my @datearray2 = split ' ', $datearray[0];
my $date_display = join ' ', $datearray2[5], $datearray2[6], $datearray2[7], $datearray2[8];
my $title_display = join '', "massacre - ", $date_display;

##############################
# Filename Related Variables #
##############################
my @log;

#my $input_csv = "somefile.csv";
#my $input_ail = "somefile_ail.xls";
my $input_digest = "somedigestfile.txt";
my $digestselect_text;
my $digestselect_val;
#my $output_root = "outputroot";

my $output_root_line = `pwd`;
chomp($output_root_line);
my @output_root_array = split /\//, $output_root_line;
my $output_root = $output_root_array[$#output_root_array];

my $input_csv = join '', $output_root, ".csv";
my $input_ail = join '', $output_root, "_ail.xls";

my $output_match1 = join '', $output_root, ".match1"; # Output of Step 1 - Input of Step 2
my $output_match2 = join '', $output_root, ".match2"; # Output of Step 2
my $output_match3 = join '', $output_root, ".match3"; # Output of Step 3
my $output_specdir = join '', $output_root, "_peaks"; # Output of Step 4
my $output_batch = join '', $output_root, ".isoin"; # Output of Step 4
my $output_plot = join '', $output_root, "_plot.txt"; # Output of Step 4
#my $output_fitsdir = join '', $output_root, "_fits"; # Output of Step 5
my $output_fitsdir = $output_specdir; # Output of Step 5
my $output_fits = join '', $output_root, "_iso.csv"; # Output of Step 5

#my $blueid = "sykes";
my $blueid = `whoami`;
chomp($blueid);

my $logfilename = join '', $output_root, ".log"; # Logfile

my $notifymessage = "Ready";

my $time_start;
my $time_end;
my $elapsed_minutes;

my $paramfile = "file.massparam";

###################################
# Configuration Related Variables #
###################################
# All(?) of these variables need to be stored/loaded when desired
#

# Step 1
# 0 = N14 and N15
# 1 = N14 only
# 2 = N15 only
my $matchmode = 0;
my $match_ppm_thresh = 50;
my $sample_id = "0";

# Step 2
my $match_rt_filter = 1;
my $match_rt_minval = 10;
my $match_rt_maxval = 50;

# Step 3
my $pairmode = 0;
my $pair_rt_thresh = 0.1;

# Step 4
my $filter_abundance = 1;
my $filter_rt = 0;
my $filter_mz = 0;
my $abundance_mincut = 2.5;
my $rt_mincut = 10;
my $rt_maxcut = 50;
my $mz_mincut = 300;
my $mz_maxcut = 1300;
my $sortstyle = 0; # 0 for n14mz, 1 for n14rt
my $slice_per_pt = 5;
my $mz_per_pt = 1;

# Step 5
my $fit_model = 3;
my $last_model = 3; # Doesn't need to be output, read in as a duplicate of $fit_model
my $niter = 4;
my $sig_global = 100;
my $baseline = 1.0;
my $ul_amp = 1.2;
my $l_amp = 0.5;
my $offset = 0.01;
my $gw = 0.0003;
my $frac_lab = 0.993;
my $nround = 5;
my $nfitpar = 6;
#my $fitschedule = "1 0 0 1 1 0 1 0 0 0 0 1 1 0 1 1 1 1 1 1 1 0 0 1 1 1 1 1 1 1"; # From isodist_batch.in
#my $fitschedule = "0 1 1 0 0 0 0 0 0 0 0 1 0 1 1 0 1 1 0 0 0 1 1 1 1 1 1 1 1 1"; # From isodist_marq.tcl
my $fitschedule = "0 0 0 1 1 0 0 0 0 0 0 1 0 0 1 1 1 1 0 1 1 0 0 1 1 1 1 1 1 1";
my @fitarray_1d;
my @fitarray_2d;

# This is not really a variable per se - we do not need to save/load
my %nfitparhash; # links $nfitpar with $fit_model
$nfitparhash{1} = 4;
$nfitparhash{2} = 5;
$nfitparhash{3} = 6;

# This just stores some built-in defaults - we do not need to save/load 
my %defaultschedule;
$defaultschedule{1} = "0 1 0 0 0 1 0 1 0 1 1 1 1 1 1 1";
$defaultschedule{2} = "0 1 0 0 0 0 0 0 0 1 0 1 0 1 1 0 0 1 1 1 1 1 1 1 1";
#$defaultschedule{3} = "1 0 0 1 1 0 1 0 0 0 0 1 1 0 1 1 1 1 1 1 1 0 0 1 1 1 1 1 1 1";
#$defaultschedule{3} = "0 1 1 0 0 0 0 0 0 0 0 1 0 1 1 0 1 1 0 0 0 1 1 1 1 1 1 1 1 1";
$defaultschedule{3} = "0 0 0 1 1 0 0 0 0 0 0 1 0 0 1 1 1 1 0 1 1 0 0 1 1 1 1 1 1 1";

# Upon the first call of the config window, the appropriate key/value pair will be populated with $fitschedule
# Therefor we do not need to save/load
my %lastschedule;
$lastschedule{1} = "0 1 0 0 0 1 0 1 0 1 1 1 1 1 1 1";
$lastschedule{2} = "0 1 0 0 0 0 0 0 0 1 0 1 0 1 1 0 0 1 1 1 1 1 1 1 1";
#$lastschedule{3} = "1 0 0 1 1 0 1 0 0 0 0 1 1 0 1 1 1 1 1 1 1 0 0 1 1 1 1 1 1 1";
#$lastschedule{3} = "0 1 1 0 0 0 0 0 0 0 0 1 0 1 1 0 1 1 0 0 0 1 1 1 1 1 1 1 1 1";
$lastschedule{3} = "0 0 0 1 1 0 0 0 0 0 0 1 0 0 1 1 1 1 0 1 1 0 0 1 1 1 1 1 1 1";

# Again, not really variables per se - we do not need to save/load
my %fitvars;
@{ $fitvars{1} } = ( "b", "ul_amp", "off", "gw" );
@{ $fitvars{2} } = ( "b", "l_amp", "off", "gw", "frac" );
#@{ $fitvars{3} } = ( "b", "ul_amp", "l_amp", "off", "gw", "frac" );
@{ $fitvars{3} } = ( "b", "off", "gw", "ul_amp", "amp", "frac" );

#"0 0 0 1 1 0 0 0 0 0 0 1 0 0 1 1 1 1 0 1 1 0 0 1 1 1 1 1 1 1";

#############
# GUI Stuff #
#############

# The main window
my $mw = MainWindow->new();
$mw->configure(-title=>$title_display, -background=>'white');
$mw->minsize( qw(1024 768));

# The configuration window (will be called/setup later)
my $cw;
# Is there a better way to allow access to these than putting them as global variables?
# Can you pass these Tk-related variables to a subroutine?
my $cw_s5_fitvar1entry;
my $cw_s5_fitvar2entry;
my $cw_s5_fitvar3entry;
my $cw_s5_fitvar4entry;
my $cw_s5_fitvar5entry;
my $cw_s5_fitvar6entry;

# The Fit Schedule Window
my $sw;

#
# The Main Window is Divided into Two Frames
# The Left Frame will contain all the buttons, text inputs etc for the user to interact with
# The Right Frame will contain only the contents of the logs and output
#
my $frame_left = $mw->Frame(-border=>'5', -background=>'white')->pack(-fill=>'both', -expand=>'1', -anchor=>'nw', -side=>'left');
my $frame_right = $mw->Frame(-border=>'5', -background=>'white')->pack(-fill=>'y', -expand=>'0', -anchor=>'ne', -side=>'right');

#
# Right Frame Log Stuff
#
my $frame_logview_border = $frame_right->Frame(-background=>'black', -borderwidth=>'3')->pack(-fill=>'both', -expand=>'0', -anchor=>'ne', -side=>'right');
my $logview = $frame_logview_border->ROText(-state=>'normal', -width=>'80', -spacing1=>'1', -spacing2=>'0', -spacing3=>'2', -wrap=>'word')->pack(-side=>'left', -fill=>'both', -expand=>'1');
my $logscrolly = $frame_logview_border->Scrollbar(-orient=>'vertical')->pack(-side=>'right', -fill=>'both');
$logview->configure(-yscrollcommand=>['set'=>$logscrolly]);
$logscrolly->configure(-command=>['yview'=>$logview]);

#
# Left Frame Buttons Text etc
#
my $frame_go_border = $frame_left->Frame(-background=>'black', -borderwidth=>'3')->pack(-fill=>'both', -expand=>'0', -anchor=>'nw', -side=>'top');
my $frame_go = $frame_go_border->Frame(-borderwidth=>'2')->pack(-fill=>'both', -expand=>'0', -anchor=>'nw', -side=>'top');

my $frame_go_top = $frame_go->Frame()->pack(-side=>'top', -fill=>'both');
my $frame_go_bottom = $frame_go->Frame()->pack(-side=>'bottom', -fill=>'both');
my $frame_go_bottom_left = $frame_go_bottom->Frame()->pack(-side=>'left', -fill=>'both');
my $frame_go_bottom_right = $frame_go_bottom->Frame()->pack(-side=>'right', -fill=>'both');

	my $go_label = $frame_go_top->Label(-text=>'Start Here to Run the Complete Analysis', -font=>'*-helvetica-bold-r-normal-18-*')->pack(-side=>'left');

	my $inputlabel = $frame_go_bottom_left->Label(-text=>'Input Files:')->pack(-side=>'top', -anchor=>'nw');

	my $csvframe = $frame_go_bottom_left->Frame()->pack(-side=>'top', -fill=>'x');
	my $csventry = $csvframe->Entry(-textvariable=>\$input_csv)->pack(-side=>'right');
	my $csvlabel = $csvframe->Label(-text=>'csv',)->pack(-side=>'right');

	my $ailframe = $frame_go_bottom_left->Frame()->pack(-side=>'top', -fill=>'x');
	my $ailentry = $ailframe->Entry(-textvariable=>\$input_ail)->pack(-side=>'right');
	my $aillabel = $ailframe->Label(-text=>'ail')->pack(-side=>'right');

	my $digestframe = $frame_go_bottom_left->Frame()->pack(-side=>'top', -fill=>'x');
	my $digestentry = $digestframe->Entry(-textvariable=>\$input_digest)->pack(-side=>'right');
	my $digestlabel = $digestframe->Label(-text=>'Digest')->pack(-side=>'right');
	my $digestframe2 = $frame_go_bottom_left->Frame(-width=>'4')->pack(-side=>'top', -expand=>'0', -fill=>'none');
	my $digestselect1 = $digestframe2->Optionmenu(-options=>
	[
	["Select a Digest"=>0],
	["30S ms2 z6 c1 m0 3-13"=>1],
	["30S ms4 z6 c1 m0 3-13"=>2],
	["50S ms2 z6 c1 m0 3-13"=>3],
	["50S ms4 z6 c1 m0 3-13"=>4],
	["70S ms2 z6 c1 m0 3-13"=>5],
	["70S ms4 z6 c1 m0 3-13"=>6],
	["30S ms4 z6 c1 m0 2-13"=>7]
	], 
	-command=>\&set_digest, -textvariable=>\$digestselect_text, -variable=>\$digestselect_val)->pack(-side=>'right');

	my $rootlabel = $frame_go_bottom_right->Label(-text=>'Root Name of Output Files:')->pack(-side=>'top', -anchor=>'nw');
	my $rootframe = $frame_go_bottom_right->Frame()->pack(-side=>'top', -fill=>'x');
		my $rootentry = $rootframe->Entry(-textvariable=>\$output_root)->pack(-side=>'left', -anchor=>'ne');
		my $popbutton = $rootframe->Button(-text=>'pop', -command=>\&populate)->pack(-side=>'right');

	my $idframe = $frame_go_bottom_right->Frame()->pack(-side=>'top', -fill=>'x');
		my $idspacer = $idframe->Label(-text=>'          ')->pack(-side=>'right');
		my $identry = $idframe->Entry(-textvariable=>\$blueid, -width=>'8')->pack(-side=>'right', -anchor=>'ne');
		my $idlabel = $idframe->Label(-text=>'blue ID:')->pack(-side=>'right');

	my $gobutton = $frame_go_bottom_right->Button(-text=>'Go!', -command=>\&go)->pack(-side=>'bottom', -anchor=>'se');

my $frame_spacer1a = $frame_left->Frame(-height=>'3')->pack(-side=>'top');

my $frame_config_border = $frame_left->Frame(-background=>'black', -borderwidth=>'3')->pack(-fill=>'x', -expand=>'0', -anchor=>'nw', -side=>'top');
my $frame_config = $frame_config_border->Frame(-borderwidth=>'2')->pack(-fill=>'both', -expand=>'0', -anchor=>'nw', -side=>'top');
	my $paramlabel = $frame_config->Label(-text=>'Parameters:')->pack(-side=>'left');
	my $configbutton = $frame_config->Button(-text=>'Configure', -command=>\&config_vars)->pack(-side=>'left', -anchor=>'nw');
	my $paramsavebutton = $frame_config->Button(-text=>'Save', -command=>\&save_param)->pack(-side=>'right');
	my $paramloadbutton = $frame_config->Button(-text=>'Load', -command=>\&load_param)->pack(-side=>'right');
	my $paramentry = $frame_config->Entry(-textvariable=>\$paramfile, -width=>'12')->pack(-side=>'right');

my $frame_spacer1b = $frame_left->Frame(-height=>'3')->pack(-side=>'top');

my $frame_steps_border = $frame_left->Frame(-background=>'black', -borderwidth=>'3')->pack(-fill=>'x', -expand=>'0', -anchor=>'nw', -side=>'top');
my $frame_steps = $frame_steps_border->Frame(-borderwidth=>'2')->pack(-fill=>'both', -expand=>'0', -anchor=>'nw', -side=>'top');

my $frame_steps_top = $frame_steps->Frame()->pack(-side=>'top', -fill=>'both');
	my $steps_label = $frame_steps_top->Label(-text=>'Start Here to Run Individual Steps', -font=>'*-helvetica-bold-r-normal-18-*')->pack(-side=>'left');

my $step1_border = $frame_steps->Frame(-background=>'black', -borderwidth=>'2')->pack(-side=>'top', -fill=>'both');
my $step1 = $step1_border->Frame()->pack(-side=>'top', -fill=>'both');

	my $step1_top = $step1->Frame()->pack(-fill=>'both', -side=>'top');
	my $step1_label = $step1_top->Label(-text=>'1) Compare ail with Theoretical Digest')->pack(-side=>'left');

	my $step1_bottom = $step1->Frame()->pack(-fill=>'both', -side=>'bottom');
	my $step1_bottom_left = $step1_bottom->Frame()->pack(-side=>'left', -fill=>'both');
	my $step1_bottom_right = $step1_bottom->Frame()->pack(-side=>'right', -fill=>'both');

	my $s1_ailframe = $step1_bottom_left->Frame()->pack(-side=>'top', -fill=>'x');
	my $s1_ailentry = $s1_ailframe->Entry(-textvariable=>\$input_ail)->pack(-side=>'right');
	my $s1_aillabel = $s1_ailframe->Label(-text=>'ail')->pack(-side=>'right');	

	my $s1_digestframe = $step1_bottom_left->Frame()->pack(-side=>'top', -fill=>'x');
	my $s1_digestentry = $s1_digestframe->Entry(-textvariable=>\$input_digest)->pack(-side=>'right');
	my $s1_digestlabel = $s1_digestframe->Label(-text=>'Digest')->pack(-side=>'right');
	
	my $s1_outframe = $step1_bottom_right->Frame()->pack(-side=>'top', -fill=>'x');
	my $s1_outentry = $s1_outframe->Entry(-textvariable=>\$output_match1)->pack(-side=>'right');
	my $s1_outlabel = $s1_outframe->Label(-text=>'Output')->pack(-side=>'right');
	
	my $button1 = $step1_bottom_right->Button(-text=>'Go (1)', -command=>\&step1)->pack(-side=>'top', -anchor=>'se');

my $step1_spacer = $frame_steps->Frame(-height=>'2')->pack(-side=>'top');

my $step2_border = $frame_steps->Frame(-background=>'black', -borderwidth=>'2')->pack(-side=>'top', -fill=>'both');
my $step2 = $step2_border->Frame()->pack(-side=>'top', -fill=>'both');

	my $step2_top = $step2->Frame()->pack(-fill=>'both', -side=>'top');
	my $step2_label = $step2_top->Label(-text=>'2) Filter Matches')->pack(-side=>'left');

	my $step2_bottom = $step2->Frame()->pack(-fill=>'both', -side=>'bottom');
	my $step2_bottom_left = $step2_bottom->Frame()->pack(-side=>'left', -fill=>'both');
	my $step2_bottom_right = $step2_bottom->Frame()->pack(-side=>'right', -fill=>'both');

	my $s2_inframe = $step2_bottom_left->Frame()->pack(-side=>'top', -fill=>'both', -anchor=>'nw');
	my $s2_inentry = $s2_inframe->Entry(-textvariable=>\$output_match1)->pack(-side=>'right');
	my $s2_inlabel = $s2_inframe->Label(-text=>'Input')->pack(-side=>'right');	
	
	my $s2_outframe = $step2_bottom_right->Frame()->pack(-side=>'top', -fill=>'x');
	my $s2_outentry = $s2_outframe->Entry(-textvariable=>\$output_match2)->pack(-side=>'right');
	my $s2_outlabel = $s2_outframe->Label(-text=>'Output')->pack(-side=>'right');
	
	my $button2 = $step2_bottom_right->Button(-text=>'Go (2)', -command=>\&step2)->pack(-side=>'top', -anchor=>'se');

my $step2_spacer = $frame_steps->Frame(-height=>'2')->pack(-side=>'top');

my $step3_border = $frame_steps->Frame(-background=>'black', -borderwidth=>'2')->pack(-side=>'top', -fill=>'both');
my $step3 = $step3_border->Frame()->pack(-side=>'top', -fill=>'both');

	my $step3_top = $step3->Frame()->pack(-fill=>'both', -side=>'top');
	my $step3_label = $step3_top->Label(-text=>'3) Find Good Features or Feature Pairs')->pack(-side=>'left');

	my $step3_bottom = $step3->Frame()->pack(-fill=>'both', -side=>'bottom');
	my $step3_bottom_left = $step3_bottom->Frame()->pack(-side=>'left', -fill=>'both');
	my $step3_bottom_right = $step3_bottom->Frame()->pack(-side=>'right', -fill=>'both');

	my $s3_inframe = $step3_bottom_left->Frame()->pack(-side=>'top', -fill=>'both', -anchor=>'nw');
	my $s3_inentry = $s3_inframe->Entry(-textvariable=>\$output_match2)->pack(-side=>'right');
	my $s3_inlabel = $s3_inframe->Label(-text=>'Input')->pack(-side=>'right');	
	
	my $s3_outframe = $step3_bottom_right->Frame()->pack(-side=>'top', -fill=>'x');
	my $s3_outentry = $s3_outframe->Entry(-textvariable=>\$output_match3)->pack(-side=>'right');
	my $s3_outlabel = $s3_outframe->Label(-text=>'Output')->pack(-side=>'right');
	
	my $button3 = $step3_bottom_right->Button(-text=>'Go (3)', -command=>\&step3)->pack(-side=>'top', -anchor=>'se');

my $step3_spacer = $frame_steps->Frame(-height=>'2')->pack(-side=>'top');

my $step4_border = $frame_steps->Frame(-background=>'black', -borderwidth=>'2')->pack(-side=>'top', -fill=>'both');
my $step4 = $step4_border->Frame()->pack(-side=>'top', -fill=>'both');

	my $step4_top = $step4->Frame()->pack(-fill=>'both', -side=>'top');
	my $step4_label = $step4_top->Label(-text=>'4) Extract the Minispectra')->pack(-side=>'left');

	my $step4_bottom = $step4->Frame()->pack(-fill=>'both', -side=>'bottom');
	my $step4_bottom_left = $step4_bottom->Frame()->pack(-side=>'left', -fill=>'both');
	my $step4_bottom_right = $step4_bottom->Frame()->pack(-side=>'right', -fill=>'both');

	my $s4_inframe = $step4_bottom_left->Frame()->pack(-side=>'top', -fill=>'x');
	my $s4_inentry = $s4_inframe->Entry(-textvariable=>\$output_match3)->pack(-side=>'right');
	my $s4_inlabel = $s4_inframe->Label(-text=>'Input')->pack(-side=>'right');	

	my $s4_csvframe = $step4_bottom_left->Frame()->pack(-side=>'top', -fill=>'x');
	my $s4_csventry = $s4_csvframe->Entry(-textvariable=>\$input_csv)->pack(-side=>'right');
	my $s4_csvlabel = $s4_csvframe->Label(-text=>'csv')->pack(-side=>'right');

	my $s4_plotframe = $step4_bottom_left->Frame()->pack(-side=>'top', -fill=>'x');
	my $s4_plotentry = $s4_plotframe->Entry(-textvariable=>\$output_plot)->pack(-side=>'right');
	my $s4_plotlabel = $s4_plotframe->Label(-text=>'plot')->pack(-side=>'right');
	
	my $s4_outframe = $step4_bottom_right->Frame()->pack(-side=>'top', -fill=>'x');
	my $s4_outentry = $s4_outframe->Entry(-textvariable=>\$output_batch)->pack(-side=>'right');
	my $s4_outlabel = $s4_outframe->Label(-text=>'Output')->pack(-side=>'right');

	my $s4_specframe = $step4_bottom_right->Frame()->pack(-side=>'top', -fill=>'x');
	my $s4_specentry = $s4_specframe->Entry(-textvariable=>\$output_specdir)->pack(-side=>'right');
	my $s4_speclabel = $s4_specframe->Label(-text=>'Peaks Dir')->pack(-side=>'right');
	
	my $button4 = $step4_bottom_right->Button(-text=>'Go (4)', -command=>\&step4)->pack(-side=>'top', -anchor=>'se');

my $step4_spacer = $frame_steps->Frame(-height=>'2')->pack(-side=>'top');

my $step5_border = $frame_steps->Frame(-background=>'black', -borderwidth=>'2')->pack(-side=>'top', -fill=>'both');
my $step5 = $step5_border->Frame()->pack(-side=>'top', -fill=>'both');

	my $step5_top = $step5->Frame()->pack(-fill=>'both', -side=>'top');
	my $step5_label = $step5_top->Label(-text=>'5) Fit Spectra')->pack(-side=>'left');

	my $step5_bottom = $step5->Frame()->pack(-fill=>'both', -side=>'bottom');
	my $step5_bottom_left = $step5_bottom->Frame()->pack(-side=>'left', -fill=>'both');
	my $step5_bottom_right = $step5_bottom->Frame()->pack(-side=>'right', -fill=>'both');

	my $s5_inframe = $step5_bottom_left->Frame()->pack(-side=>'top', -fill=>'x');
	my $s5_inentry = $s5_inframe->Entry(-textvariable=>\$output_batch)->pack(-side=>'right');
	my $s5_inlabel = $s5_inframe->Label(-text=>'Input')->pack(-side=>'right');	

#	my $s5_csvframe = $step5_bottom_left->Frame()->pack(-side=>'top', -fill=>'x');
#	my $s5_csventry = $s5_csvframe->Entry(-textvariable=>\$input_batch)->pack(-side=>'right');
#	my $s5_csvlabel = $s5_csvframe->Label(-text=>'Batchfile')->pack(-side=>'right');

	my $s5_idframe = $step5_bottom_left->Frame()->pack(-side=>'top', -fill=>'x');
	my $s5_identry = $s5_idframe->Entry(-textvariable=>\$blueid)->pack(-side=>'right');
	my $s5_idlabel = $s5_idframe->Label(-text=>'blue ID')->pack(-side=>'right');

	my $s5_peaksframe = $step5_bottom_left->Frame()->pack(-side=>'top', -fill=>'x');
	my $s5_peaksentry = $s5_peaksframe->Entry(-textvariable=>\$output_specdir)->pack(-side=>'right');
	my $s5_peakslabel = $s5_peaksframe->Label(-text=>'Peaks Dir')->pack(-side=>'right');
	
	my $s5_outframe = $step5_bottom_right->Frame()->pack(-side=>'top', -fill=>'x');
	my $s5_outentry = $s5_outframe->Entry(-textvariable=>\$output_fits)->pack(-side=>'right');
	my $s5_outlabel = $s5_outframe->Label(-text=>'Fitsfile')->pack(-side=>'right');

	my $s5_specframe = $step5_bottom_right->Frame()->pack(-side=>'top', -fill=>'x');
	my $s5_specentry = $s5_specframe->Entry(-textvariable=>\$output_fitsdir)->pack(-side=>'right');
	my $s5_speclabel = $s5_specframe->Label(-text=>'Fits Dir')->pack(-side=>'right');
	
	my $button5 = $step5_bottom_right->Button(-text=>'Go (5)', -command=>\&step5)->pack(-side=>'top', -anchor=>'se');

my $frame_spacer2 = $frame_left->Frame(-height=>'3')->pack(-side=>'top');

my $frame_log_border = $frame_left->Frame(-background=>'black', -borderwidth=>'3')->pack(-fill=>'x', -expand=>'0', -anchor=>'nw', -side=>'top');
my $frame_log = $frame_log_border->Frame(-borderwidth=>'2')->pack(-fill=>'both', -expand=>'0', -anchor=>'nw', -side=>'top');
	my $logentry = $frame_log->Entry(-textvariable=>\$logfilename, -width=>'15')->pack(-side=>'left');
	my $logclearbutton = $frame_log->Button(-text=>'Clear Log', -command=>\&clear_log)->pack(-side=>'left');
	my $logsavebutton = $frame_log->Button(-text=>'Save Log', -command=>\&save_log)->pack(-side=>'left');
	my $quitbutton = $frame_log->Button(-text=>'Quit', -command=>\&quit)->pack(-side=>'right', -anchor=>'nw');

my $frame_spacer3 = $frame_left->Frame(-height=>'3')->pack(-side=>'top');

my $frame_notify_border = $frame_left->Frame(-background=>'black', -borderwidth=>'3')->pack(-fill=>'x', -expand=>'0', -anchor=>'nw', -side=>'top');
my $frame_notify = $frame_notify_border->Frame(-background=>'green', -borderwidth=>'5')->pack(-fill=>'both', -expand=>'1', -anchor=>'nw', -side=>'top');
my $frame_notify_label = $frame_notify->Label(-background=>'green', -textvariable=>\$notifymessage)->pack(-fill=>'both', -expand=>'1');

################################
# Misc GUI Configuration Stuff #
################################
$logview->tagConfigure('red', -foreground=>'#CC0000');
$logview->tagConfigure('green', -foreground=>'#008800');
$logview->tagConfigure('blue', -foreground=>'blue');


###############
# Subroutines #
###############

#
# Update both the $logview widget and the @log array
#
sub ulog
{
	my $ctr;
	my $color;
	my @logentries = @_;
	for($ctr = 0; $ctr <= $#logentries; $ctr++)
	{
		chomp($logentries[$ctr]);
		my $entry = join '', $logentries[$ctr], "\n";
		
		if( index( lc($entry), "error") >= 0)
		{
			$color = "red";
			$program_status = -1;
		}
		elsif( index( lc($entry), "can't") >= 0)
		{
			$color = "red";
			$program_status = -1;
		}
		elsif( index( lc($entry), "terminated") >= 0)
		{
			$color = "red";
			$program_status = -1;
		}
		elsif( index( lc($entry), "finished step") >= 0)
		{
			$color = "green";
		}
		
		push(@log, $entry);
		$logview->insert('end', $entry, $color);
	}
	$logview->yview('end');
	$logview->update;
}

sub save_log
{
	my $ctr;	
	
	ulog("Writing logfile to $logfilename");
	
	open LOGFILE, ">$logfilename" or die "Can't open $logfilename\n";
	for($ctr = 0; $ctr <= $#log; $ctr++)
	{
		print LOGFILE "$log[$ctr]";
	}
	close LOGFILE;
}

sub clear_log
{
	ulog("Clearing log (does not affect display or saved files)");
	@log = ();
}

sub set_digest
{
	if($digestselect_val == 0)
	{
		return;
	}

	if($digestselect_val == 1)
	{
		$input_digest = "$ENV{'MASSACRE_PATH'}/digest/files/30S_miss2_cons_z6_1_0_300_1300.txt";
	}
	elsif($digestselect_val == 2)
	{
		$input_digest = "$ENV{'MASSACRE_PATH'}/digest/files/30S_miss4_cons_z6_1_0_300_1300.txt";
	}
	elsif($digestselect_val == 3)
	{
		$input_digest = "$ENV{'MASSACRE_PATH'}/digest/files/50S_miss2_cons_z6_1_0_300_1300.txt";
	}
	elsif($digestselect_val == 4)
	{
		$input_digest = "$ENV{'MASSACRE_PATH'}/digest/files/50S_miss4_cons_z6_1_0_300_1300.txt";
	}
	elsif($digestselect_val == 5)
	{
		$input_digest = "$ENV{'MASSACRE_PATH'}/digest/files/70S_miss2_cons_z6_1_0_300_1300.txt";
	}
	elsif($digestselect_val == 6)
	{
		$input_digest = "$ENV{'MASSACRE_PATH'}/digest/files/70S_miss4_cons_z6_1_0_300_1300.txt";
	}
	elsif($digestselect_val == 7)
	{
		$input_digest = "$ENV{'MASSACRE_PATH'}/digest/files/30S_miss4_cons_z6_1_0_200_1300.txt";
	}
		
	ulog("Setting Digest File:");
	#ulog("$digestselect_text");
	ulog("$input_digest");


}

sub populate
{
	$output_match1 = join '', $output_root, ".match1"; # Output of Step 1 - Input of Step 2
	$output_match2 = join '', $output_root, ".match2"; # Output of Step 2
	$output_match3 = join '', $output_root, ".match3"; # Output of Step 3
	$output_specdir = join '', $output_root, "_peaks"; # Output of Step 4
	$output_batch = join '', $output_root, ".isoin"; # Output of Step 4
	$output_plot = join '', $output_root, "_plot.txt"; # Output of Step 4
	#$output_fitsdir = join '', $output_root, "_fits"; # Output of Step 5
	$output_fitsdir = $output_specdir; # Output of Step 5
	$output_fits = join '', $output_root, "_iso.csv"; # Output of Step 5

	$logfilename = join '', $output_root, ".log"; # Logfile
}

sub quit
{
	print "Goodbye\n";
	exit(0);
}

#
# Run all steps from one command
#
sub go
{
	$time_start = time();
	$program_status = 0;

	populate();
	
	step1();
	
	if($program_status == 0)
	{
		 step2();
	}
	else
	{
		ulog("Skipping Step 2");
	}
	
	if($program_status == 0)
	{
		 step3();
	}
	else
	{
		ulog("Skipping Step 3");
	}
		
	if($program_status == 0)
	{
		 step4();
	}
	else
	{
		ulog("Skipping Step 4");
	}
		
	if($program_status == 0)
	{
		 step5();
	}
	else
	{
		ulog("Skipping Step 5");
	}
		
	save_log();
	
	if($program_status == -1)
	{
		ulog("Error: Program failed to complete.  Please check log for details.");
	}
	
	$time_end = time();
	$elapsed_minutes = ($time_end - $time_start) / 60;
	ulog("Total Time Elapsed: $elapsed_minutes minutes");
}

sub step1
{
	my $line;

	notify("massacre is running - please wait");
	ulog("Starting Step 1");
	ulog("$ENV{'MASSACRE_PATH'}msc/msc01_findmatches.pl --mode=$matchmode --sample_id=$sample_id --ppm_thresh=$match_ppm_thresh --digest=$input_digest --ail=$input_ail --output=$output_match1");

	open LOG, "$ENV{'MASSACRE_PATH'}/msc/msc01_findmatches.pl --mode=$matchmode --id=$sample_id --ppm_thresh=$match_ppm_thresh --digest=$input_digest --ail=$input_ail --output=$output_match1 2>&1 |" or die "Can't Fork: $!\n";
	while(defined($line=<LOG>))
	{
		ulog($line);
	}
	close LOG;
	
	ulog("Finished Step 1");
	ulog("----------");
	stop_notify();
}

#
# Note - if adding additional processing steps, need to add intermediate input/output files
# First Input: $output_match1
# Last Output: $output_match2
#
sub step2
{
	my $line;
	
	notify("massacre is running - please wait");
	ulog("Starting Step 2");
	
	if($match_rt_filter == 1)
	{
		ulog("$ENV{'MASSACRE_PATH'}msc/msc02a_rtfilter.pl --input=$output_match1 --output=$output_match2 --rt_min=$match_rt_minval --rt_max=$match_rt_maxval");
		
		open LOG, "$ENV{'MASSACRE_PATH'}/msc/msc02a_rtfilter.pl --input=$output_match1 --output=$output_match2 --rt_min=$match_rt_minval --rt_max=$match_rt_maxval 2>&1 |" or die "Can't Fork: $!\n";
		while(defined($line=<LOG>))
		{
			ulog($line);
		}
		close LOG;
	}
	else
	{
		`cp $output_match1 $output_match2`;
	}

	ulog("Finished Step 2");
	ulog("----------");
	stop_notify();
}

sub step3
{
	my $line;

	notify("massacre is running - please wait");
	ulog("Starting Step 3");
	ulog("$ENV{'MASSACRE_PATH'}/msc/msc03_getsinglesorpairs.pl --input=$output_match2 --output=$output_match3 --matchmode=$matchmode --pairmode=$pairmode --rt_thresh=$pair_rt_thresh");
	
	open LOG, "$ENV{'MASSACRE_PATH'}/msc/msc03_getsinglesorpairs.pl --input=$output_match2 --output=$output_match3 --matchmode=$matchmode --pairmode=$pairmode --rt_thresh=$pair_rt_thresh 2>&1 |" or die "Can't Fork: $!\n";
	while(defined($line=<LOG>))
	{
		ulog($line);
	}
	close LOG;
	
	ulog("Finished Step 3");
	ulog("----------");
	stop_notify();
}

sub step4
{
	my $line;

	notify("massacre is running - please wait");
	ulog("Starting Step 4");
	ulog("$ENV{'MASSACRE_PATH'}/msc/msc04_makepeaks_1local.pl --input=$output_match3 --output=$output_batch --dir=$output_specdir --csv=$input_csv --matchmode=$matchmode --filter_abundance=$filter_abundance --filter_mz=$filter_mz --filter_rt=$filter_rt --abundance_min=$abundance_mincut --mz_min=$mz_mincut --mz_max=$mz_maxcut --rt_min=$rt_mincut --rt_max=$rt_maxcut --id=$sample_id --plot=$output_plot --blueid=$blueid");

	open LOG, "$ENV{'MASSACRE_PATH'}/msc/msc04_makepeaks_1local.pl --input=$output_match3 --output=$output_batch --dir=$output_specdir --csv=$input_csv --matchmode=$matchmode --filter_abundance=$filter_abundance --filter_mz=$filter_mz --filter_rt=$filter_rt --abundance_min=$abundance_mincut --mz_min=$mz_mincut --mz_max=$mz_maxcut --rt_min=$rt_mincut --rt_max=$rt_maxcut --id=$sample_id --plot=$output_plot --blueid=$blueid --sortstyle=$sortstyle --slice_per_pt=$slice_per_pt --mz_per_pt=$mz_per_pt 2>&1 |" or die "Can't Fork: $!\n";
	while(defined($line=<LOG>))
	{
		ulog($line);
	}
	close LOG;

	ulog("Finished Step 4");
	ulog("----------");
	stop_notify();
}

sub step5
{
	my $line;

	notify("massacre is running - please wait");
	ulog("Starting Step 5");
	ulog("$ENV{'MASSACRE_PATH'}/msc/msc05_runisodist.pl --input=$output_batch --output=$output_fits --dir=$output_specdir --fit_model=$fit_model --baseline=$baseline --ul_amp=$ul_amp --l_amp=$l_amp --offset=$offset --gw=$gw --frac_lab=$frac_lab --nround=$nround --nfitpar=$nfitpar --fitschedule=\"$fitschedule\" --blueid=$blueid");

	# Note, we pass $output_specdir here as the new isodist_batch uses the same output directory

	open LOG, "$ENV{'MASSACRE_PATH'}/msc/msc05_runisodist.pl --input=$output_batch --output=$output_fits --dir=$output_specdir --fitsdir=$output_fitsdir --fit_model=$fit_model --baseline=$baseline --ul_amp=$ul_amp --l_amp=$l_amp --offset=$offset --gw=$gw --frac_lab=$frac_lab --nround=$nround --nfitpar=$nfitpar --fitschedule=\"$fitschedule\" --blueid=$blueid 2>&1 |" or die "Can't Fork: $!\n";
	
	while(defined($line=<LOG>))
	{
		ulog($line);
	}
	close LOG;
	
	ulog("Finished Step 5");
	ulog("----------");
	stop_notify();
}

###################################################
# Subroutines related to the configuration window #
###################################################

sub config_vars
{
	$cw = MainWindow->new();
	$cw->configure(-title=>'Configure Variables', -background=>'white', -border=>'5');
	#$cw->minsize( qw(800 600));
	
	my $cw_s1_border = $cw->Frame(-border=>'3', -background=>'black')->pack(-side=>'top', -fill=>'both');
	my $cw_s1 = $cw_s1_border->Frame(-border=>'2')->pack(-side=>'top', -fill=>'both');
	my $cw_s1_labelframe = $cw_s1->Frame()->pack(-side=>'top', -fill=>'both');
	my $cw_s1_label = $cw_s1_labelframe->Label(-text=>'1) Compare ail with Theoretical Digest')->pack(-side=>'left');

	my $cw_s1_varframe1 = $cw_s1->Frame(-pady=>'2')->pack(-side=>'top', -fill=>'both');
		my $cw_s1_matchbutton0 = $cw_s1_varframe1->Radiobutton(-text=>'N14 + N15', -variable=>\$matchmode, -value=>'0')->pack(-side=>'left');
		my $cw_s1_matchbutton1 = $cw_s1_varframe1->Radiobutton(-text=>'N14 Only', -variable=>\$matchmode, -value=>'1')->pack(-side=>'left');
		my $cw_s1_matchbutton2 = $cw_s1_varframe1->Radiobutton(-text=>'N15 Only', -variable=>\$matchmode, -value=>'2')->pack(-side=>'left');

	my $cw_s1_varframe2 = $cw_s1->Frame(-pady=>'2')->pack(-side=>'top', -fill=>'both');
		my $cw_s1_var2label = $cw_s1_varframe2->Label(-text=>'PPM Threshold for Matches:')->pack(-side=>'left');
		my $cw_s2_var2val = $cw_s1_varframe2->Entry(-textvariable=>\$match_ppm_thresh)->pack(-side=>'left');

	my $cw_s1_varframe3 = $cw_s1->Frame(-pady=>'2')->pack(-side=>'top', -fill=>'both');
		my $cw_s1_var3label = $cw_s1_varframe3->Label(-text=>'Sample ID:')->pack(-side=>'left');
		my $cw_s2_var3val = $cw_s1_varframe3->Entry(-textvariable=>\$sample_id)->pack(-side=>'left');

	my $cw_spacer1 = $cw->Frame(-height=>'3')->pack(-side=>'top');

	my $cw_s2_border = $cw->Frame(-border=>'3', -background=>'black')->pack(-side=>'top', -fill=>'both');
	my $cw_s2 = $cw_s2_border->Frame(-border=>'2')->pack(-side=>'top', -fill=>'both');
	my $cw_s2_labelframe = $cw_s2->Frame()->pack(-side=>'top', -fill=>'both');
	my $cw_s2_label = $cw_s2_labelframe->Label(-text=>'2) Filter Matches')->pack(-side=>'left');

	my $cw_s2_varframe1 = $cw_s2->Frame(-pady=>'2')->pack(-side=>'top', -fill=>'both');
		my $cw_s2_var1label = $cw_s2_varframe1->Label(-text=>'Filter by RT?')->pack(-side=>'left');
		my $cw_s2_var1button1 = $cw_s2_varframe1->Radiobutton(-text=>'Yes', -variable=>\$match_rt_filter, -value=>'1')->pack(-side=>'left');
		my $cw_s2_var1button0 = $cw_s2_varframe1->Radiobutton(-text=>'No', -variable=>\$match_rt_filter, -value=>'0')->pack(-side=>'left');

		my $cw_s2_var1label2 = $cw_s2_varframe1->Label(-text=>' -> Min RT:')->pack(-side=>'left');
		my $cw_s2_var1entry2 = $cw_s2_varframe1->Entry(-textvariable=>\$match_rt_minval, -width=>'5')->pack(-side=>'left');
		my $cw_s2_var1label3 = $cw_s2_varframe1->Label(-text=>'Max RT:')->pack(-side=>'left');
		my $cw_s2_var1entry3 = $cw_s2_varframe1->Entry(-textvariable=>\$match_rt_maxval, -width=>'5')->pack(-side=>'left');

	my $cw_spacer2 = $cw->Frame(-height=>'3')->pack(-side=>'top');

	my $cw_s3_border = $cw->Frame(-border=>'3', -background=>'black')->pack(-side=>'top', -fill=>'both');
	my $cw_s3 = $cw_s3_border->Frame(-border=>'2')->pack(-side=>'top', -fill=>'both');
	my $cw_s3_labelframe = $cw_s3->Frame()->pack(-side=>'top', -fill=>'both');
	my $cw_s3_label = $cw_s3_labelframe->Label(-text=>'3) Find Good Features or Feature Pairs')->pack(-side=>'left');

	my $cw_s3_varframe1 = $cw_s3->Frame(-pady=>'2')->pack(-side=>'top', -anchor=>'nw');
		my $cw_s3_matchbutton0 = $cw_s3_varframe1->Radiobutton(-text=>'N14 + N15', -variable=>\$matchmode, -value=>'0')->pack(-side=>'left');
		my $cw_s3_matchbutton1 = $cw_s3_varframe1->Radiobutton(-text=>'N14 Only', -variable=>\$matchmode, -value=>'1')->pack(-side=>'left');
		my $cw_s3_matchbutton2 = $cw_s3_varframe1->Radiobutton(-text=>'N15 Only', -variable=>\$matchmode, -value=>'2')->pack(-side=>'left');

	my $cw_s3_varframe2 = $cw_s3->Frame(-pady=>'2')->pack(-side=>'top', -anchor=>'nw');
		my $cw_s3_pairtypebutton0 = $cw_s3_varframe2->Radiobutton(-text=>'Non-Redundant m/z Only', -variable=>\$pairmode, -value=>'0')->pack(-side=>'top', -anchor=>'nw');
		my $cw_s3_pairtypebutton1 = $cw_s3_varframe2->Radiobutton(-text=>'Redundant OK (Same Protein)', -variable=>\$pairmode, -value=>'1', -state=>'disabled')->pack(-side=>'top', -anchor=>'nw');
		my $cw_s3_pairtypebutton2 = $cw_s3_varframe2->Radiobutton(-text=>'Redundant OK (All)', -variable=>\$pairmode, -value=>'2')->pack(-side=>'top', -anchor=>'nw');

	my $cw_s3_varframe3 = $cw_s3->Frame(-pady=>'2')->pack(-side=>'top', -anchor=>'nw');
		my $cw_s3_var3label1 = $cw_s3_varframe3->Label(-text=>'Pairs Must be Within')->pack(-side=>'left');
		my $cw_s3_var3entry = $cw_s3_varframe3->Entry(-textvariable=>\$pair_rt_thresh, -width=>'5')->pack(-side=>'left');
		my $cw_s3_var3label2 = $cw_s3_varframe3->Label(-text=>'Minutes')->pack(-side=>'left');

	my $cw_spacer3 = $cw->Frame(-height=>'3')->pack(-side=>'top');

	my $cw_s4_border = $cw->Frame(-border=>'3', -background=>'black')->pack(-side=>'top', -fill=>'both');
	my $cw_s4 = $cw_s4_border->Frame(-border=>'2')->pack(-side=>'top', -fill=>'both');
	my $cw_s4_labelframe = $cw_s4->Frame()->pack(-side=>'top', -fill=>'both');
	my $cw_s4_label = $cw_s4_labelframe->Label(-text=>'4) Extract the Minispectra')->pack(-side=>'left');

	my $cw_s4_varframe1 = $cw_s4->Frame(-pady=>'2')->pack(-side=>'top', -fill=>'both');
		my $cw_s4_matchbutton0 = $cw_s4_varframe1->Radiobutton(-text=>'N14 + N15', -variable=>\$matchmode, -value=>'0')->pack(-side=>'left');
		my $cw_s4_matchbutton1 = $cw_s4_varframe1->Radiobutton(-text=>'N14 Only', -variable=>\$matchmode, -value=>'1')->pack(-side=>'left');
		my $cw_s4_matchbutton2 = $cw_s4_varframe1->Radiobutton(-text=>'N15 Only', -variable=>\$matchmode, -value=>'2')->pack(-side=>'left');

	my $cw_s4_varframe2 = $cw_s4->Frame(-pady=>'2')->pack(-side=>'top', -fill=>'both');
		my $cw_s4_var2label = $cw_s4_varframe2->Label(-text=>'Filter by RT?')->pack(-side=>'left');
		my $cw_s4_var2button1 = $cw_s4_varframe2->Radiobutton(-text=>'Yes', -variable=>\$filter_rt, -value=>'1')->pack(-side=>'left');
		my $cw_s4_var2button0 = $cw_s4_varframe2->Radiobutton(-text=>'No', -variable=>\$filter_rt, -value=>'0')->pack(-side=>'left');

		my $cw_s4_var2label2 = $cw_s4_varframe2->Label(-text=>' -> Min RT:')->pack(-side=>'left');
		my $cw_s4_var2entry2 = $cw_s4_varframe2->Entry(-textvariable=>\$rt_mincut, -width=>'5')->pack(-side=>'left');
		my $cw_s4_var2label3 = $cw_s4_varframe2->Label(-text=>'Max RT:')->pack(-side=>'left');
		my $cw_s4_var2entry3 = $cw_s4_varframe2->Entry(-textvariable=>\$rt_maxcut, -width=>'5')->pack(-side=>'left');

	my $cw_s4_varframe3 = $cw_s4->Frame(-pady=>'2')->pack(-side=>'top', -fill=>'both');
		my $cw_s4_var3label = $cw_s4_varframe3->Label(-text=>'Filter by MZ?')->pack(-side=>'left');
		my $cw_s4_var3button1 = $cw_s4_varframe3->Radiobutton(-text=>'Yes', -variable=>\$filter_mz, -value=>'1')->pack(-side=>'left');
		my $cw_s4_var3button0 = $cw_s4_varframe3->Radiobutton(-text=>'No', -variable=>\$filter_mz, -value=>'0')->pack(-side=>'left');

		my $cw_s4_var3label2 = $cw_s4_varframe3->Label(-text=>' -> Min MZ:')->pack(-side=>'left');
		my $cw_s4_var3entry2 = $cw_s4_varframe3->Entry(-textvariable=>\$mz_mincut, -width=>'5')->pack(-side=>'left');
		my $cw_s4_var3label3 = $cw_s4_varframe3->Label(-text=>'Max MZ:')->pack(-side=>'left');
		my $cw_s4_var3entry3 = $cw_s4_varframe3->Entry(-textvariable=>\$mz_maxcut, -width=>'5')->pack(-side=>'left');

	my $cw_s4_varframe4 = $cw_s4->Frame(-pady=>'2')->pack(-side=>'top', -fill=>'both');
		my $cw_s4_var4label = $cw_s4_varframe4->Label(-text=>'Filter by AB?')->pack(-side=>'left');
		my $cw_s4_var4button1 = $cw_s4_varframe4->Radiobutton(-text=>'Yes', -variable=>\$filter_abundance, -value=>'1')->pack(-side=>'left');
		my $cw_s4_var4button0 = $cw_s4_varframe4->Radiobutton(-text=>'No', -variable=>\$filter_abundance, -value=>'0')->pack(-side=>'left');

		my $cw_s4_var4label2 = $cw_s4_varframe4->Label(-text=>' -> Min Abundance:')->pack(-side=>'left');
		my $cw_s4_var4entry2 = $cw_s4_varframe4->Entry(-textvariable=>\$abundance_mincut, -width=>'5')->pack(-side=>'left');

	my $cw_s4_varframe5 = $cw_s4->Frame(-pady=>'4')->pack(-side=>'top', -fill=>'both');
		my $cw_s4_var5label = $cw_s4_varframe5->Label(-text=>'Sort Output By:')->pack(-side=>'left');
		my $cw_s4_var5button0 = $cw_s4_varframe5->Radiobutton(-text=>'MZ', -variable=>\$sortstyle, -value=>'0')->pack(-side=>'left');
		my $cw_s4_var5button1 = $cw_s4_varframe5->Radiobutton(-text=>'RT', -variable=>\$sortstyle, -value=>'1')->pack(-side=>'left');
		my $cw_s4_var5button1 = $cw_s4_varframe5->Radiobutton(-text=>'Peptide', -variable=>\$sortstyle, -value=>'2')->pack(-side=>'left');

	my $cw_s4_varframe6 = $cw_s4->Frame(-pady=>'4')->pack(-side=>'top', -fill=>'both');
		my $cw_s4_var6label1 = $cw_s4_varframe6->Label(-text=>'Contour Plot: RTSlices/pt')->pack(-side=>'left');
		my $cw_s4_var6entry1 = $cw_s4_varframe6->Entry(-textvariable=>\$slice_per_pt, -width=>'5')->pack(-side=>'left');
		my $cw_s4_var6label2 = $cw_s4_varframe6->Label(-text=>'MZ/pt')->pack(-side=>'left');
		my $cw_s4_var6entry2 = $cw_s4_varframe6->Entry(-textvariable=>\$mz_per_pt, -width=>'5')->pack(-side=>'left');

	my $cw_spacer4 = $cw->Frame(-height=>'3')->pack(-side=>'top');

	my $cw_s5_border = $cw->Frame(-border=>'3', -background=>'black')->pack(-side=>'top', -fill=>'both');
	my $cw_s5 = $cw_s5_border->Frame(-border=>'2')->pack(-side=>'top', -fill=>'both');
	my $cw_s5_labelframe = $cw_s5->Frame()->pack(-side=>'top', -fill=>'both');
	my $cw_s5_label = $cw_s5_labelframe->Label(-text=>'5) Fit Spectra')->pack(-side=>'left');

	my $cw_s5_varframe1a = $cw_s5->Frame(-pady=>'0')->pack(-side=>'top', -fill=>'both');
		my $cw_s5_fitbutton1 = $cw_s5_varframe1a->Radiobutton(-text=>'Natural Abundance (b,ul_amp,gw,off)', -variable=>\$fit_model, -value=>'1', -command=>\&set_fitvars, -state=>'disabled')->pack(-side=>'left');
	my $cw_s5_varframe1b = $cw_s5->Frame(-pady=>'0')->pack(-side=>'top', -fill=>'both');
		my $cw_s5_fitbutton2 = $cw_s5_varframe1b->Radiobutton(-text=>'Fractional Label (b,l_amp,gw,off,frac)', -variable=>\$fit_model, -value=>'2', -command=>\&set_fitvars, -state=>'disabled')->pack(-side=>'left');
	my $cw_s5_varframe1c = $cw_s5->Frame(-pady=>'0')->pack(-side=>'top', -fill=>'both');
		my $cw_s5_fitbutton3 = $cw_s5_varframe1c->Radiobutton(-text=>'Unlabeled + Labeled (b,ul_amp,l_amp,gw,off,frac)', -variable=>\$fit_model, -value=>'3', -command=>\&set_fitvars)->pack(-side=>'left');

	my $cw_s5_varframe2 = $cw_s5->Frame(-pady=>'2')->pack(-side=>'top', -fill=>'both');
	my $cw_s5_varframe2left = $cw_s5_varframe2->Frame(-background=>'red')->pack(-side=>'left');

		my $cw_s5_fitvar1frame = $cw_s5_varframe2left->Frame()->pack(-side=>'top', -fill=>'both');
		 $cw_s5_fitvar1entry = $cw_s5_fitvar1frame->Entry(-textvariable=>\$baseline, -width=>'6')->pack(-side=>'right');
		my $cw_s5_fitvar1label = $cw_s5_fitvar1frame->Label(-text=>'Baseline:')->pack(-side=>'right');

		my $cw_s5_fitvar2frame = $cw_s5_varframe2left->Frame()->pack(-side=>'top', -fill=>'both');
		 $cw_s5_fitvar2entry = $cw_s5_fitvar2frame->Entry(-textvariable=>\$ul_amp, -width=>'6')->pack(-side=>'right');
		my $cw_s5_fitvar2label = $cw_s5_fitvar2frame->Label(-text=>'Unlabeled Amp:')->pack(-side=>'right');

		my $cw_s5_fitvar3frame = $cw_s5_varframe2left->Frame()->pack(-side=>'top', -fill=>'both');
		 $cw_s5_fitvar3entry = $cw_s5_fitvar3frame->Entry(-textvariable=>\$l_amp, -width=>'6')->pack(-side=>'right');
		my $cw_s5_fitvar3label = $cw_s5_fitvar3frame->Label(-text=>'Labeled Amp:')->pack(-side=>'right');

	my $cw_s5_varframe2center = $cw_s5_varframe2->Frame()->pack(-side=>'left');

		my $cw_s5_fitvar4frame = $cw_s5_varframe2center->Frame()->pack(-side=>'top', -fill=>'both');
		 $cw_s5_fitvar4entry = $cw_s5_fitvar4frame->Entry(-textvariable=>\$offset, -width=>'6')->pack(-side=>'right');
		my $cw_s5_fitvar4label = $cw_s5_fitvar4frame->Label(-text=>'Offset:')->pack(-side=>'right');

		my $cw_s5_fitvar5frame = $cw_s5_varframe2center->Frame()->pack(-side=>'top', -fill=>'both');
		 $cw_s5_fitvar5entry = $cw_s5_fitvar5frame->Entry(-textvariable=>\$gw, -width=>'6')->pack(-side=>'right');
		my $cw_s5_fitvar5label = $cw_s5_fitvar5frame->Label(-text=>'Gaussian Width:')->pack(-side=>'right');

		my $cw_s5_fitvar6frame = $cw_s5_varframe2center->Frame()->pack(-side=>'top', -fill=>'both');
		 $cw_s5_fitvar6entry = $cw_s5_fitvar6frame->Entry(-textvariable=>\$frac_lab, -width=>'6')->pack(-side=>'right');
		my $cw_s5_fitvar6label = $cw_s5_fitvar6frame->Label(-text=>'Frac N:')->pack(-side=>'right');
	
	my $cw_s5_varframe2right = $cw_s5_varframe2->Frame(-padx=>'5')->pack(-side=>'right', -fill=>'both'); 
		my $cw_s5_schedulebutton = $cw_s5_varframe2right->Button(-text=>'Fit Schedule', -command=>\&fit_schedule)->pack(-side=>'bottom', -anchor=>'se');

	my $cw_s5_varframe3 = $cw_s5->Frame(-pady=>'2')->pack(-side=>'top', -fill=>'both');

	my $cw_spacer5 = $cw->Frame(-height=>'3')->pack(-side=>'top');

	my $cw_close_border = $cw->Frame(-border=>'3', -background=>'black')->pack(-side=>'top', -fill=>'both');
	my $cw_close = $cw_close_border->Frame(-border=>'2')->pack(-side=>'top', -fill=>'both');
	my $cw_closebutton = $cw_close->Button(-text=>'Save/Close', -command=>\&close_config)->pack(-side=>'right', -anchor=>'e');

	set_fitvars();
}

sub close_config
{
	$cw->destroy();
}

# Check the status of $fit_model and color variables appropriately
sub set_fitvars
{
	my @local_array;

	$lastschedule{$last_model} = $fitschedule;

	if($fit_model == 1)
	{
		$cw_s5_fitvar1entry->configure(-background=>'lightgrey');
		$cw_s5_fitvar2entry->configure(-background=>'lightgrey');
		$cw_s5_fitvar3entry->configure(-background=>'red');
		$cw_s5_fitvar4entry->configure(-background=>'lightgrey');
		$cw_s5_fitvar5entry->configure(-background=>'lightgrey');
		$cw_s5_fitvar6entry->configure(-background=>'red');
		
		$nfitpar = $nfitparhash{$fit_model};
		$fitschedule = $lastschedule{$fit_model};
		@local_array = split ' ', $fitschedule;
		$nround = ($#local_array + 1) / $nfitpar;
	}
	elsif($fit_model == 2)
	{
		$cw_s5_fitvar1entry->configure(-background=>'lightgrey');
		$cw_s5_fitvar2entry->configure(-background=>'red');
		$cw_s5_fitvar3entry->configure(-background=>'lightgrey');
		$cw_s5_fitvar4entry->configure(-background=>'lightgrey');
		$cw_s5_fitvar5entry->configure(-background=>'lightgrey');
		$cw_s5_fitvar6entry->configure(-background=>'lightgrey');	
		
		$nfitpar = $nfitparhash{$fit_model};
		$fitschedule = $lastschedule{$fit_model};
		@local_array = split ' ', $fitschedule;
		$nround = ($#local_array + 1) / $nfitpar;
	}
	elsif($fit_model == 3)
	{
		$cw_s5_fitvar1entry->configure(-background=>'lightgrey');
		$cw_s5_fitvar2entry->configure(-background=>'lightgrey');
		$cw_s5_fitvar3entry->configure(-background=>'lightgrey');
		$cw_s5_fitvar4entry->configure(-background=>'lightgrey');
		$cw_s5_fitvar5entry->configure(-background=>'lightgrey');
		$cw_s5_fitvar6entry->configure(-background=>'lightgrey');

		
		$nfitpar = $nfitparhash{$fit_model};
		$fitschedule = $lastschedule{$fit_model};
		@local_array = split ' ', $fitschedule;
		$nround = ($#local_array + 1) / $nfitpar;
	}
}

sub fit_schedule
{
	my $ctr1;
	my $ctr2;

	$sw = MainWindow->new();
	$sw->configure(-title=>'Fit Schedule', -background=>'white', -border=>'5');

	my $sw_borderframe = $sw->Frame(-background=>'black', -border=>'3')->pack();

	my $sw_buttonframe = $sw_borderframe->Frame()->pack(-side=>'top', -fill=>'both');
	my $sw_addbutton = $sw_buttonframe->Button(-text=>'Add Round', -command=>\&add_round)->pack(-side=>'left');
	my $sw_defaultbutton = $sw_buttonframe->Button(-text=>'Default Schedule')->pack(-side=>'left');
	my $sw_savebutton = $sw_buttonframe->Button(-text=>'Save/Close', -command=>\&close_schedule)->pack(-side=>'right');

	my $sw_scheduleframe = $sw_borderframe->Frame()->pack(-side=>'top', -fill=>'both');
	
	@fitarray_1d = split ' ', $fitschedule;
	$nround = ($#fitarray_1d + 1) / $nfitpar;
	
	my @sw_frames;
	my @sw_labels;
	my @sw_text;
	my @sw_checks;
	my @sw_buttons;
	
	for($ctr1 = 0; $ctr1 < $nround; $ctr1++)
	{
		$sw_frames[$ctr1] = $sw_scheduleframe->Frame()->pack(-side=>'top', -fill=>'both');
		
		$sw_text[$ctr1] = join '', "Round ", $ctr1 + 1, ": ";
		$sw_labels[$ctr1] = $sw_frames[$ctr1]->Label(-text=>$sw_text[$ctr1])->pack(-side=>'left');
		
		for($ctr2 = 0; $ctr2 < $nfitpar; $ctr2++)
		{
			$fitarray_2d[$ctr1][$ctr2] = $fitarray_1d[ $nfitpar*$ctr1 + $ctr2 ];

			$sw_checks[$ctr1][$ctr2] = $sw_frames[$ctr1]->Checkbutton(-text=>$fitvars{$fit_model}[$ctr2], -offvalue=>'0', -onvalue=>'1', -variable=>\$fitarray_2d[$ctr1][$ctr2])->pack(-side=>'left');
		}
	
		$sw_buttons[$ctr1] = $sw_frames[$ctr1]->Button(-text=>'Delete', -command=>[\&delete_round, $ctr1])->pack(-side=>'left');
	}
	
}

sub delete_round
{
	my $ctr1;
	my $ctr2;
	
	my @parameters = @_;
	my $which_round = $parameters[0];

	#ulog($param);

	#save_schedule(); # Not necessary -> Current state is already in @fitarray_2d;

	@fitarray_1d = ();
	
	for($ctr1 = 0; $ctr1 < $nround; $ctr1++)
	{
		if($ctr1 == $which_round)
		{
			next;
		}
	
		for($ctr2 = 0; $ctr2 < $nfitpar; $ctr2++)
		{
			push(@fitarray_1d, $fitarray_2d[$ctr1][$ctr2]);
		}
	}

	$fitschedule = join ' ', @fitarray_1d;

	$sw->destroy();
	
	fit_schedule();
}

sub add_round
{
	my $ctr;
	
	save_schedule();
	
	for($ctr = 0; $ctr < $nfitpar; $ctr++)
	{
		$fitschedule = join ' ', $fitschedule, "0";
	}

	$sw->destroy();
	
	fit_schedule();
}

sub save_schedule
{
	# Read the current state into $fitschedule
	my $ctr1;
	my $ctr2;
	
	@fitarray_1d = ();
	
	for($ctr1 = 0; $ctr1 < $nround; $ctr1++)
	{
		for($ctr2 = 0; $ctr2 < $nfitpar; $ctr2++)
		{
			push(@fitarray_1d, $fitarray_2d[$ctr1][$ctr2]);
		}
	}
	
	$fitschedule = join ' ', @fitarray_1d;
	
	ulog($fitschedule);
}

sub close_schedule
{
	save_schedule();

	$sw->destroy();
}

sub notify
{
	my @messages = @_;
	$notifymessage = $messages[0];
	
	$frame_notify->configure(-background=>'red');
	$frame_notify_label->configure(-background=>'red');
}

sub stop_notify
{
	$notifymessage = "Ready";
	
	$frame_notify->configure(-background=>'green');
	$frame_notify_label->configure(-background=>'green');
}

#
# load_param() and save_param() are really ugly
#
sub load_param
{
	my $ctr;
	my @parray;
	my $param;
	my $value;
	my %phash;
	my @local_array = ();

	ulog("Reading parameter file");

	my @param_array = readfile($paramfile);

	if($#param_array < 0)
	{
		ulog("Error loading parameters");
		return;
	}

	PLINE: for($ctr = 0; $ctr <= $#param_array; $ctr++)
	{
		if(substr($param_array[$ctr], 0, 1) eq "#")
		{
			next PLINE;
		}

		# Due to this split operation, the first thing is the parameter name, the second is the value
		# Third etc elements are ignored, regardless of if they are commented out in some way
		@parray = split ' ', $param_array[$ctr];
		$param = $parray[0];
		$value = $parray[1];
		#ulog("$param -> $value");
		$phash{$param} = $value;
	}
	
	foreach $param (keys %phash)
	{
		$value = $phash{$param};
		if($value eq "YES")
		{
			$value = 1;
		}
		elsif($value eq "NO")
		{
			$value = 0;
		}
	
		# Step 1
		if($param eq "matchmode")
		{
			$matchmode = $value;
		}
		elsif($param eq "match_ppm_thresh")
		{
			$match_ppm_thresh = $value;
		}
		elsif($param eq "sample_id")
		{
			$sample_id = $value;
		}

		# Step 2
		elsif($param eq "match_rt_filter")
		{
			$match_rt_filter = $value;
		}
		elsif($param eq "match_rt_minval")
		{
			$match_rt_minval = $value;
		}
		elsif($param eq "match_rt_maxval")
		{
			$match_rt_maxval = $value;
		}
		
		# Step 3
		elsif($param eq "pairmode")
		{
			$pairmode = $value;
		}		
		elsif($param eq "pair_rt_thresh")
		{
			$pair_rt_thresh = $value;
		}		
		
		# Step 4
		elsif($param eq "filter_abundance")
		{
			$filter_abundance = $value;
		}		
		elsif($param eq "filter_rt")
		{
			$filter_rt = $value;
		}
		elsif($param eq "filter_mz")
		{
			$filter_mz = $value;
		}
		elsif($param eq "abundance_mincut")
		{
			$abundance_mincut = $value;
		}
		elsif($param eq "rt_mincut")
		{
			$rt_mincut = $value;
		}
		elsif($param eq "rt_maxcut")
		{
			$rt_maxcut = $value;
		}
		elsif($param eq "mz_mincut")
		{
			$mz_mincut = $value;
		}
		elsif($param eq "mz_maxcut")
		{
			$mz_maxcut = $value;
		}
		elsif($param eq "sortstyle")
		{
			$sortstyle = $value;
		}
		elsif($param eq "slice_per_pt")
		{
			$slice_per_pt = $value;
		}
		elsif($param eq "mz_per_pt")
		{
			$mz_per_pt = $value;
		}
		
		# Step 5
		elsif($param eq "fit_model")
		{
			$fit_model = $value;
			$last_model = $value;
		}
		elsif($param eq "niter")
		{
			$niter = $value;
		}
		elsif($param eq "sig_global")
		{
			$sig_global = $value;
		}
		elsif($param eq "baseline")
		{
			$baseline = $value;
		}
		elsif($param eq "ul_amp")
		{
			$ul_amp = $value;
		}
		elsif($param eq "l_amp")
		{
			$l_amp = $value;
		}
		elsif($param eq "offset")
		{
			$offset = $value;
		}
		elsif($param eq "gw")
		{
			$gw = $value;
		}
		elsif($param eq "frac_lab")
		{
			$frac_lab = $value;
		}
		elsif($param eq "nround")
		{
			$nround = $value;
		}		
		elsif($param eq "nfitpar")
		{
			$nfitpar = $value;
		}				
		elsif($param eq "fitschedule")
		{
			# Fit Schedule is saved with commas, so deal with that
			@local_array = split /\,/, $value;
			$value = join ' ', @local_array;
			$fitschedule = $value;
		}						
	}
	
	ulog("Finished reading parameter file");
}

sub save_param
{
	my $error = 0;
	my @local_array;
	my $fitout;

	ulog("Writing parameter file");
	
	open PARAM, ">$paramfile" or $error = 1;
	if($error == 1)
	{
		ulog("Error: Can't open $paramfile");
		die;
	}
	
	# Step 1
	print PARAM "matchmode $matchmode\n";
	print PARAM "match_ppm_thresh $match_ppm_thresh\n";
	print PARAM "sample_id $sample_id\n";

	# Step 2
	print PARAM "match_rt_filter $match_rt_filter\n";
	print PARAM "match_rt_minval $match_rt_minval\n";
	print PARAM "match_rt_maxval $match_rt_maxval\n";
		
	# Step 3
	print PARAM "pairmode $pairmode\n";
	print PARAM "pair_rt_thresh $pair_rt_thresh\n";
	
	# Step 4
	print PARAM "filter_abundance $filter_abundance\n";
	print PARAM "filter_rt $filter_rt\n";
	print PARAM "filter_mz $filter_mz\n";
	print PARAM "abundance_mincut $abundance_mincut\n";
	print PARAM "rt_mincut $rt_mincut\n";
	print PARAM "rt_maxcut $rt_maxcut\n";
	print PARAM "mz_mincut $mz_mincut\n";
	print PARAM "mz_maxcut $mz_maxcut\n";
	print PARAM "sortstyle $sortstyle\n";
	print PARAM "slice_per_pt $slice_per_pt\n";
	print PARAM "mz_per_pt $mz_per_pt\n";
	
	# Step 5
	print PARAM "fit_model $fit_model\n";
	print PARAM "niter $niter\n";
	print PARAM "sig_global $sig_global\n";
	print PARAM "baseline $baseline\n";
	print PARAM "ul_amp $ul_amp\n";
	print PARAM "l_amp $l_amp\n";
	print PARAM "offset $offset\n";
	print PARAM "gw $gw\n";
	print PARAM "frac_lab $frac_lab\n";
	print PARAM "nround $nround\n";
	print PARAM "nfitpar $nfitpar\n";
	@local_array = split ' ', $fitschedule;
	$fitout = join ',', @local_array;
	print PARAM "fitschedule $fitout\n";
	
	close PARAM;
	
	ulog("Finished writing $paramfile");
}


sub readfile
{
	my @rf_parameters = @_;
	my @rf_array1;
	my @rf_array2;
	my $rf_ctr1;
	my $rf_ctr2;
	my $rf_input_file;
	my @rf_output_array = ();	
	my $rf_error = 0;
		
	if($#rf_parameters != 0)
	{
		ulog("Error: Subroutine \"readfile()\" accepts only one parameter");
		return;
	}

	$rf_input_file = $rf_parameters[0];
	open INPUT, "$rf_input_file" or $rf_error = 1;
	if($rf_error == 1)
	{
		ulog("Error: Can't open $rf_input_file");
		return;
	}
	@rf_array1 = <INPUT>;
	for ($rf_ctr1 = 0; $rf_ctr1 <= $#rf_array1; $rf_ctr1++)
	{
		chomp($rf_array1[$rf_ctr1]);
		@rf_array2 = split /\r/, $rf_array1[$rf_ctr1];
		for($rf_ctr2 = 0; $rf_ctr2 <= $#rf_array2; $rf_ctr2++)
		{
			push(@rf_output_array, $rf_array2[$rf_ctr2]);
		}
	}
	close INPUT;
	
	return @rf_output_array;
}

&GetOptions("param=s" => \$paramfile_input, "root=s" => \$output_root_input, "go=i" => \$go, "digest=s" => \$input_digest);

#
# Here we check on some possible command line input, and react accordingly
#
if($paramfile_input ne "")
{
	$paramfile = $paramfile_input;
	load_param();
}
if($output_root_input ne "")
{
	$output_root = $output_root_input;
	populate();
}
if($go >= 0)
{
	if($go == 0)
	{
		go();
	}
	elsif($go == 1)
	{
		step1();
	}
	elsif($go == 2)
	{
		step2();
	}
	elsif($go == 3)
	{
		step3();
	}
	elsif($go == 4)
	{
		step4();
	}
	elsif($go == 5)
	{
		step5();
	}
}

MainLoop();

