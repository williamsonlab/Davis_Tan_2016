#!/usr/bin/perl

use strict;
require Tk;
require Tk::ROText;
require Tk::X11Font;
use Tk;
use Tk::FileSelect;
use Getopt::Long;
use IPC::Open3;

#######################
# Operation Variables #
#######################
my $program_status = 0;

##########################
# Command Line Variables #
##########################
my $massacre_path = $ENV{'MASSACRE_PATH'};
my $go = 0;
my $paramfile_input = "";
my $output_root_input = "";
# Other variables which can be altered on the command line but are defined elsewhere
# $input_digest

my @datearray = ();
@datearray = `/bin/ls -lt $massacre_path/msc/msc*.pl $massacre_path/bin/massacre.pl`;
#@datearray = `/bin/ls -ltT $massacre_path/msc/msc*.pl $massacre_path/bin/massacre.pl`;
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
my $ailtype = 0;
my $ailtypetext;

my $output_match0 = join '', $output_root, ".match0";
my $output_rtdir = join '', $output_root, "_rtdir";
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
#my $machine = "bluefish";
#my $machine = "garibaldi";
my $machine = "garibaldi";

my $logfilename = join '', $output_root, ".log"; # Logfile

my $notifymessage = "Ready";

my $time_start;
my $time_end;
my $elapsed_minutes;

my $paramfile = "file.massparam";

my $datatype_flag = 0;

#if($blueid ne "sykes")
#{
#	print "12/2/2010 9:41PM - massacre is down for maintenance.  You will receive an email from Sykes when it is fully operational again\n";
#	die;
#}

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
my $match_ppm_offset = 0;
my $match_ppm_thresh = 50;
my $sample_id = "0";

# Step 2
my $match_rt_filter = 1;
my $match_rt_minval = 10;
my $match_rt_maxval = 50;

# Step 3
my $pairmode = 0; # Keep Completely Non-Redundant Only By Default
my $pair_rt_thresh = 0.1;

# Step 4
my $filter_abundance = 0;
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
my $rt_expand = 0.1;
my $rt_expand_lo = 0.0;
my $rt_expand_hi = 0.0;
my $special_range = 0; # by default, use standard N14/N15 range
my $special_range_residues = "";
my $special_range_values = "";
my $set_special_range = 0;
my $keepcon = 0;
#my $spectra3D = 0; # set to 1 to extract 3D (ie not summed) minispectra
# $extractstyle now replaces and extends $spectra3D
# 0 for the standard mz_bin_tol method
# 1 for simple binning in the m/z domain
# 2 for 3D minispectra
my $extractstyle = 0;
my $extractbin = 0.05; # The m/z range to sum intensity over for $extractstyle = 1

# Step 5
my $fit_model = 3;
my $last_model = 3; # Doesn't need to be output, read in as a duplicate of $fit_model
my $niter = 4;
my $sig_global = 100;
my $baseline = 1.0;
my $offset = 0.01;
my $gw = 0.003;
my $baseline_fit_type = 0;

my $input_model = "somemodelfile.txt";
my $modelselect_text;
my $modelselect_val;

my $input_atoms = "$massacre_path/msc/files/exp_atom_defs.txt";

my $isodist_binary = "isodist3_rsd2"; # not saved to parameter file

my $keepfit = 0; # By default, do not keep the .fit files

# Step 0 MSMS
my $msms_proteinkeep = 2;

# Step 1,2,3 MSMS
# ssc101911 changed mzrange_low from 400 to 200
my $mzrange_low = 200;
my $mzrange_high = 2000;
my $rt_max_int_fraction = 0.5;

#############
# GUI Stuff #
#############

# The main window
my $mw = MainWindow->new();
$mw->configure(-title=>$title_display, -background=>'white');
$mw->minsize( qw(1024 768));

#
# gnuplot stuff
#
my $gnupid1 = open3( \*GIN1, \*GOUT1, \*GERR1, "gnuplot" ) || die "Can't open gnuplot\n";
$mw->fileevent( \*GOUT1, readable => \&read_gout1 );
$mw->fileevent( \*GERR1, readable => \&read_gerr1 );
print GIN1 "set term X11 noraise nopersist\n";


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
my $logview = $frame_logview_border->ROText(-state=>'normal', -width=>'60', -spacing1=>'1', -spacing2=>'0', -spacing3=>'2', -wrap=>'word')->pack(-side=>'left', -fill=>'both', -expand=>'1');
my $logscrolly = $frame_logview_border->Scrollbar(-orient=>'vertical')->pack(-side=>'right', -fill=>'both');
$logview->configure(-yscrollcommand=>['set'=>$logscrolly]);
$logscrolly->configure(-command=>['yview'=>$logview]);

#
# Left Frame Buttons Text etc
#
my $frame_header = $frame_left->Frame(-background=>'white')->pack(-side=>'top');
my $frame_header_left = $frame_header->Frame(-background=>'white')->pack(-side=>'left');
my $frame_header_spacer = $frame_header->Frame(-width=>'5', -background=>'white')->pack(-side=>'left');
my $frame_header_right = $frame_header->Frame(-background=>'white')->pack(-side=>'right');

my $frame_go_border = $frame_header_left->Frame(-background=>'black', -borderwidth=>'3')->pack(-side=>'left', -anchor=>'nw');
my $frame_go = $frame_go_border->Frame(-borderwidth=>'2')->pack(-fill=>'both', -expand=>'0', -anchor=>'nw', -side=>'top');

my $frame_go_top = $frame_go->Frame()->pack(-side=>'top', -fill=>'both');
my $frame_go_bottom = $frame_go->Frame()->pack(-side=>'bottom', -fill=>'both');
my $frame_go_bottom_left = $frame_go_bottom->Frame()->pack(-side=>'left', -fill=>'both');
my $frame_go_bottom_right = $frame_go_bottom->Frame()->pack(-side=>'right', -fill=>'both');

	my $go_label = $frame_go_top->Label(-text=>'Start Here to Run the Complete Analysis', -font=>'*-helvetica-bold-r-normal-18-*')->pack(-side=>'left');

	my $inputlabel = $frame_go_bottom_left->Label(-text=>'Input Files:')->pack(-side=>'top', -anchor=>'nw');

	my $csvframe = $frame_go_bottom_left->Frame()->pack(-side=>'top', -fill=>'x');
	my $csventry = $csvframe->Entry(-textvariable=>\$input_csv)->pack(-side=>'right');
	my $csvlabel = $csvframe->Label(-text=>'MS Data',)->pack(-side=>'right');

	my $ailframe = $frame_go_bottom_left->Frame()->pack(-side=>'top', -fill=>'x');
	my $ailframe2 = $frame_go_bottom_left->Frame()->pack(-side=>'top', -fill=>'x');
	my $ailentry = $ailframe->Entry(-textvariable=>\$input_ail)->pack(-side=>'right');
	my $aillabel = $ailframe->Label(-text=>'Peaks')->pack(-side=>'right');
	#my $ailtypebutton1 = $ailframe2->Radiobutton(-text=>'MH2', -variable=>\$ailtype, -value=>'1')->pack(-side=>'right');
	#my $ailtypebutton0 = $ailframe2->Radiobutton(-text=>'MH1', -variable=>\$ailtype, -value=>'0')->pack(-side=>'right');
	my $ailtypemenu = $ailframe2->Optionmenu(-options=>
	[
	["MassProfiler AIL"=>0],
	["New MassProfiler AIL"=>6],
	["Qualitative Analysis CSV"=>1],
	["xcms Annotated Peak List"=>2],
	["OpenMS FeatureML"=>3],
	["AMP Matlab MALDI"=>7],
	["Mascot MS/MS IDs"=>4],
	["Sequest MS/MS IDs"=>5]
	],
	-variable=>\$ailtype, -textvariable=>\$ailtypetext, -command=>\&newailtype)->pack(-side=>'right');

	my $digestframe = $frame_go_bottom_left->Frame()->pack(-side=>'top', -fill=>'x');
	my $digestentry = $digestframe->Entry(-textvariable=>\$input_digest)->pack(-side=>'right');
	my $digestlabel = $digestframe->Label(-text=>'Digest')->pack(-side=>'right');
	my $digestframe2 = $frame_go_bottom_left->Frame(-width=>'4')->pack(-side=>'top', -expand=>'0', -fill=>'none');
	my $digestselectbutton = $digestframe2->Button(-text=>'Select Digest', -command=>\&select_digest)->pack(-side=>'right');

	my $rootlabel = $frame_go_bottom_right->Label(-text=>'Root Name of Output:')->pack(-side=>'top', -anchor=>'ne');
	my $rootframe = $frame_go_bottom_right->Frame()->pack(-side=>'top', -fill=>'x');
		my $popbutton = $rootframe->Button(-text=>'pop', -command=>\&populate)->pack(-side=>'right');
		my $rootentry = $rootframe->Entry(-textvariable=>\$output_root, -width=>
		'12')->pack(-side=>'right', -anchor=>'ne');

	my $idframe = $frame_go_bottom_right->Frame()->pack(-side=>'top', -fill=>'x');
		my $idspacer = $idframe->Label(-text=>'     ')->pack(-side=>'right');
		my $identry = $idframe->Entry(-textvariable=>\$blueid, -width=>'8')->pack(-side=>'right', -anchor=>'ne');
		my $idlabel = $idframe->Label(-text=>'Remote ID:')->pack(-side=>'right');

	my $machineframe = $frame_go_bottom_right->Frame()->pack(-side=>'top', -fill=>'x');
		my $machinespacer = $machineframe->Label(-text=>'     ')->pack(-side=>'right');
		my $machineentry = $machineframe->Entry(-textvariable=>\$machine, -width=>'8')->pack(-side=>'right', -anchor=>'ne');
		my $machinelabel = $machineframe->Label(-text=>'Remote CPU:')->pack(-side=>'right');

	my $goframe = $frame_go_bottom_right->Frame()->pack(-side=>'bottom', -anchor=>'se');

	my $gobutton = $goframe->Button(-text=>'Go!', -command=>\&go)->pack(-side=>'right');
	my $gobutton_matchhisto = $goframe->Button(-text=>'Get Offset', -command=>\&go_match_histo)->pack(-side=>'right');


my $frame_config_border = $frame_header_right->Frame(-background=>'black', -borderwidth=>'3')->pack(-side=>'top');

my $frame_config_top = $frame_config_border->Frame(-borderwidth=>'2')->pack(-fill=>'both', -expand=>'0', -anchor=>'nw', -side=>'top');
my $frame_config_middle = $frame_config_border->Frame(-borderwidth=>'2')->pack(-fill=>'both', -expand=>'0', -anchor=>'nw', -side=>'top');
my $frame_config_bottom = $frame_config_border->Frame(-borderwidth=>'2')->pack(-fill=>'both', -expand=>'0', -anchor=>'nw', -side=>'top');
	my $paramlabel = $frame_config_top->Label(-text=>'Parameters:')->pack(-side=>'left');
	my $configbutton = $frame_config_middle->Button(-text=>'Configure', -command=>\&config_vars)->pack(-side=>'left', -anchor=>'nw');
	my $paramsavebutton = $frame_config_bottom->Button(-text=>'Save', -command=>\&save_param, -width=>'4')->pack(-side=>'right');
	my $paramloadbutton = $frame_config_bottom->Button(-text=>'Load', -command=>\&load_param, -width=>'4')->pack(-side=>'right');
	my $paramentry = $frame_config_bottom->Entry(-textvariable=>\$paramfile, -width=>'12')->pack(-side=>'right');

my $frame_spacer3 = $frame_header_right->Frame(-height=>'3')->pack(-side=>'top');

my $frame_log_border = $frame_header_right->Frame(-background=>'black', -borderwidth=>'3')->pack(-fill=>'x', -expand=>'0', -side=>'top');
my $frame_log_top = $frame_log_border->Frame(-borderwidth=>'2')->pack(-fill=>'both', -expand=>'0', -anchor=>'nw', -side=>'top');
my $frame_log_bottom = $frame_log_border->Frame(-borderwidth=>'2')->pack(-fill=>'both', -expand=>'0', -anchor=>'nw', -side=>'top');

	my $loglabel = $frame_log_top->Label(-text=>'Logfile:')->pack(-side=>'left');
	my $logentry = $frame_log_bottom->Entry(-textvariable=>\$logfilename, -width=>'12')->pack(-side=>'left');
	my $logclearbutton = $frame_log_bottom->Button(-text=>'Clear', -command=>\&clear_log, -width=>'4')->pack(-side=>'left');
	my $logsavebutton = $frame_log_bottom->Button(-text=>'Save', -command=>\&save_log, -width=>'4')->pack(-side=>'left');

my $frame_spacer1b = $frame_left->Frame(-height=>'3')->pack(-side=>'top');

my $frame_steps_border = $frame_left->Frame(-background=>'black', -borderwidth=>'3')->pack(-fill=>'x', -expand=>'0', -anchor=>'nw', -side=>'top');
my $frame_steps = $frame_steps_border->Frame(-borderwidth=>'2')->pack(-fill=>'both', -expand=>'0', -anchor=>'nw', -side=>'top');

my $frame_steps_top = $frame_steps->Frame()->pack(-side=>'top', -fill=>'both');
	my $steps_label = $frame_steps_top->Label(-text=>'Start Here to Run Individual Steps', -font=>'*-helvetica-bold-r-normal-18-*')->pack(-side=>'left');

my $frame_steps_middle = $frame_steps->Frame()->pack(-side=>'top', -fill=>'both');
	my $frame_steps_left = $frame_steps_middle->Frame()->pack(-side=>'left', -fill=>'both');
	my $frame_steps_middle_space = $frame_steps_middle->Frame(-width=>'5')->pack(-side=>'left');
	my $frame_steps_right  = $frame_steps_middle->Frame()->pack(-side=>'right', -fill=>'both');

my $frame_steps_bottom = $frame_steps->Frame()->pack(-side=>'top', -fill=>'both');

my $step1_border = $frame_steps_left->Frame(-background=>'black', -borderwidth=>'2')->pack(-side=>'top', -fill=>'both');
my $step1 = $step1_border->Frame()->pack(-side=>'left', -fill=>'both', -expand=>'1');
my $step1r = $step1_border->Frame()->pack(-side=>'right', -fill=>'both');

	my $step1_top = $step1->Frame()->pack(-fill=>'both', -side=>'top');
	my $step1_label = $step1_top->Label(-text=>'1) Compare ail with Theoretical Digest')->pack(-side=>'left');

	my $step1_bottom = $step1->Frame()->pack(-fill=>'both', -side=>'bottom');
	my $step1_bottom_left = $step1_bottom->Frame()->pack(-side=>'left', -fill=>'both');
	my $step1_bottom_right = $step1_bottom->Frame()->pack(-side=>'right', -fill=>'both');

	my $s1_ailframe = $step1_bottom_left->Frame()->pack(-side=>'top', -fill=>'x');
	my $s1_aillabel = $s1_ailframe->Label(-text=>'ail')->pack(-side=>'top', -anchor=>'w');	
	my $s1_ailentry = $s1_ailframe->Entry(-textvariable=>\$input_ail)->pack(-side=>'top', -anchor=>'w');

	my $s1_digestframe = $step1_bottom_left->Frame()->pack(-side=>'top', -fill=>'x');
	my $s1_digestlabel = $s1_digestframe->Label(-text=>'Digest')->pack(-side=>'top', -anchor=>'w');
	my $s1_digestentry = $s1_digestframe->Entry(-textvariable=>\$input_digest)->pack(-side=>'top', -anchor=>'w');
	
	my $s1_outframe = $step1_bottom_right->Frame()->pack(-side=>'top', -fill=>'x');
	my $s1_outlabel = $s1_outframe->Label(-text=>'Output')->pack(-side=>'top', -anchor=>'w');
	my $s1_outentry = $s1_outframe->Entry(-textvariable=>\$output_match1)->pack(-side=>'top', -anchor=>'w');
	
	my $button1 = $step1r->Button(-text=>'Go (1)', -command=>\&step1)->pack(-side=>'bottom', -anchor=>'se');

my $step1_spacer = $frame_steps_left->Frame(-height=>'2')->pack(-side=>'top');

my $step2_border = $frame_steps_left->Frame(-background=>'black', -borderwidth=>'2')->pack(-side=>'top', -fill=>'both');
my $step2 = $step2_border->Frame()->pack(-side=>'left', -fill=>'both', -expand=>'1');
my $step2r = $step2_border->Frame()->pack(-side=>'right', -fill=>'both');

	my $step2_top = $step2->Frame()->pack(-fill=>'both', -side=>'top');
	my $step2_label = $step2_top->Label(-text=>'2) Filter Matches')->pack(-side=>'left');

	my $step2_bottom = $step2->Frame()->pack(-fill=>'both', -side=>'bottom');
	my $step2_bottom_left = $step2_bottom->Frame()->pack(-side=>'left', -fill=>'both');
	my $step2_bottom_right = $step2_bottom->Frame()->pack(-side=>'right', -fill=>'both');

	my $s2_inframe = $step2_bottom_left->Frame()->pack(-side=>'top', -fill=>'both', -anchor=>'nw');
	my $s2_inlabel = $s2_inframe->Label(-text=>'Input')->pack(-side=>'top', -anchor=>'w');	
	my $s2_inentry = $s2_inframe->Entry(-textvariable=>\$output_match1)->pack(-side=>'top', -anchor=>'w');
	
	my $s2_outframe = $step2_bottom_right->Frame()->pack(-side=>'top', -fill=>'x');
	my $s2_outlabel = $s2_outframe->Label(-text=>'Output')->pack(-side=>'top', -anchor=>'w');
	my $s2_outentry = $s2_outframe->Entry(-textvariable=>\$output_match2)->pack(-side=>'top', -anchor=>'w');
	
	my $button2 = $step2r->Button(-text=>'Go (2)', -command=>\&step2)->pack(-side=>'bottom', -anchor=>'se');

my $step2_spacer = $frame_steps_left->Frame(-height=>'2')->pack(-side=>'top');

my $step3_border = $frame_steps_left->Frame(-background=>'black', -borderwidth=>'2')->pack(-side=>'top', -fill=>'both');
my $step3 = $step3_border->Frame()->pack(-side=>'left', -fill=>'both', -expand=>'1');
my $step3r = $step3_border->Frame()->pack(-side=>'right', -fill=>'both');

	my $step3_top = $step3->Frame()->pack(-fill=>'both', -side=>'top');
	my $step3_label = $step3_top->Label(-text=>'3) Find Good Features or Feature Pairs')->pack(-side=>'left');

	my $step3_bottom = $step3->Frame()->pack(-fill=>'both', -side=>'bottom');
	my $step3_bottom_left = $step3_bottom->Frame()->pack(-side=>'left', -fill=>'both');
	my $step3_bottom_right = $step3_bottom->Frame()->pack(-side=>'right', -fill=>'both');

	my $s3_inframe = $step3_bottom_left->Frame()->pack(-side=>'top', -fill=>'both', -anchor=>'nw');
	my $s3_inlabel = $s3_inframe->Label(-text=>'Input')->pack(-side=>'top', -anchor=>'w');	
	my $s3_inentry = $s3_inframe->Entry(-textvariable=>\$output_match2)->pack(-side=>'top', -anchor=>'w');
	
	my $s3_outframe = $step3_bottom_right->Frame()->pack(-side=>'top', -fill=>'x');
	my $s3_outlabel = $s3_outframe->Label(-text=>'Output')->pack(-side=>'top', -anchor=>'w');
	my $s3_outentry = $s3_outframe->Entry(-textvariable=>\$output_match3)->pack(-side=>'top', -anchor=>'w');
	
	my $button3 = $step3r->Button(-text=>'Go (3)', -command=>\&step3)->pack(-side=>'bottom', -anchor=>'se');

my $step3_spacer = $frame_steps_left->Frame(-height=>'2')->pack(-side=>'top');

my $step0ms2_border = $frame_steps_right->Frame(-background=>'black', -borderwidth=>'2')->pack(-side=>'top', -fill=>'both', -anchor=>'ne');
my $step0ms2 = $step0ms2_border->Frame()->pack(-side=>'left', -fill=>'both', -expand=>'1');
my $step0ms2r = $step0ms2_border->Frame()->pack(-side=>'right', -fill=>'both');

	my $step0ms2_top = $step0ms2->Frame()->pack(-fill=>'both', -side=>'top');
	my $step0ms2_label = $step0ms2_top->Label(-text=>'0) Perform Retention Time Correction')->pack(-side=>'left');

	my $step0ms2_bottom = $step0ms2->Frame()->pack(-fill=>'both', -side=>'bottom');
	my $step0ms2_bottom_left = $step0ms2_bottom->Frame()->pack(-side=>'left', -fill=>'both');
	my $step0ms2_bottom_right = $step0ms2_bottom->Frame()->pack(-side=>'right', -fill=>'both');

	my $step0ms2_csvframe = $step0ms2_bottom_left->Frame()->pack(-side=>'top', -fill=>'x');
	my $step0ms2_csvlabel = $step0ms2_csvframe->Label(-text=>'csv')->pack(-side=>'top', -anchor=>'w');
	my $step0ms2_csventry = $step0ms2_csvframe->Entry(-textvariable=>\$input_ail)->pack(-side=>'top', -anchor=>'w');

	my $step0ms2_mzdataframe = $step0ms2_bottom_left->Frame()->pack(-side=>'top', -fill=>'x');
	my $step0ms2_mzdatalabel = $step0ms2_mzdataframe->Label(-text=>'MS Data')->pack(-side=>'top', -anchor=>'w');
	my $step0ms2_mzdataentry = $step0ms2_mzdataframe->Entry(-textvariable=>\$input_csv)->pack(-side=>'top', -anchor=>'w');
	
	my $step0ms2_outframe = $step0ms2_bottom_right->Frame()->pack(-side=>'top', -fill=>'x');
	my $step0ms2_outlabel = $step0ms2_outframe->Label(-text=>'Output')->pack(-side=>'top', -anchor=>'w');
	my $step0ms2_outentry = $step0ms2_outframe->Entry(-textvariable=>\$output_match0)->pack(-side=>'top', -anchor=>'w');

	my $step0ms2_rtdirframe = $step0ms2_bottom_left->Frame()->pack(-side=>'top', -fill=>'x');
	my $step0ms2_rtdirlabel = $step0ms2_outframe->Label(-text=>'RT Plot Dir')->pack(-side=>'top', -anchor=>'w');
	my $step0ms2_rtdirentry = $step0ms2_outframe->Entry(-textvariable=>\$output_rtdir)->pack(-side=>'top', -anchor=>'w');
	
	my $button0ms2 = $step0ms2r->Button(-text=>'Go (0)', -command=>\&step0ms2)->pack(-side=>'bottom', -anchor=>'se');

my $step0ms2_spacer = $frame_steps_right->Frame(-height=>'2')->pack(-side=>'top');

my $step123ms2_border = $frame_steps_right->Frame(-background=>'black', -borderwidth=>'2')->pack(-side=>'top', -fill=>'both', -anchor=>'ne');
my $step123ms2 = $step123ms2_border->Frame()->pack(-side=>'left', -fill=>'both', -expand=>'1');
my $step123ms2r = $step123ms2_border->Frame()->pack(-side=>'right', -fill=>'both');

	my $step123ms2_top = $step123ms2->Frame()->pack(-fill=>'both', -side=>'top');
	my $step123ms2_label = $step123ms2_top->Label(-text=>'1,2,3) Get Peptide Data from Digest and Consolidate')->pack(-side=>'left');

	my $step123ms2_bottom = $step123ms2->Frame()->pack(-fill=>'both', -side=>'bottom');
	my $step123ms2_bottom_left = $step123ms2_bottom->Frame()->pack(-side=>'left', -fill=>'both');
	my $step123ms2_bottom_right = $step123ms2_bottom->Frame()->pack(-side=>'right', -fill=>'both');

	#my $step123ms2_digestframe = $step123ms2_bottom_left->Frame()->pack(-side=>'top', -fill=>'x');
	#my $step123ms2_digestlabel = $step123ms2_digestframe->Label(-text=>'Digest')->pack(-side=>'top', -anchor=>'w');
	#my $step123ms2_digestentry = $step123ms2_digestframe->Entry(-textvariable=>\$input_digest)->pack(-side=>'top', -anchor=>'w');

	my $step123ms2_inputframe = $step123ms2_bottom_left->Frame()->pack(-side=>'top', -fill=>'x');
	my $step123ms2_inputlabel = $step123ms2_inputframe->Label(-text=>'Input')->pack(-side=>'top', -anchor=>'w');
	my $step123ms2_inputentry = $step123ms2_inputframe->Entry(-textvariable=>\$output_match0)->pack(-side=>'top', -anchor=>'w');
	
	my $step123ms2_outframe = $step123ms2_bottom_right->Frame()->pack(-side=>'top', -fill=>'x');
	my $step123ms2_outlabel = $step123ms2_outframe->Label(-text=>'Output')->pack(-side=>'top', -anchor=>'w');
	my $step123ms2_outentry = $step123ms2_outframe->Entry(-textvariable=>\$output_match3)->pack(-side=>'top', -anchor=>'w');
	
	my $button123ms2 = $step123ms2r->Button(-text=>'Go (1)', -command=>\&step123ms2)->pack(-side=>'bottom', -anchor=>'se');


my $step123ms2_spacer = $frame_steps_right->Frame(-height=>'2')->pack(-side=>'top');

my $step4_border = $frame_steps_bottom->Frame(-background=>'black', -borderwidth=>'2')->pack(-side=>'top');
my $step4 = $step4_border->Frame()->pack(-side=>'left', -fill=>'both', -expand=>'1');
my $step4r = $step4_border->Frame()->pack(-side=>'right', -fill=>'both');

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
	
	my $s4_outframe = $step4_bottom_right->Frame()->pack(-side=>'top', -fill=>'x');
	my $s4_outentry = $s4_outframe->Entry(-textvariable=>\$output_batch)->pack(-side=>'right');
	my $s4_outlabel = $s4_outframe->Label(-text=>'Output')->pack(-side=>'right');

	my $s4_specframe = $step4_bottom_right->Frame()->pack(-side=>'top', -fill=>'x');
	my $s4_specentry = $s4_specframe->Entry(-textvariable=>\$output_specdir)->pack(-side=>'right');
	my $s4_speclabel = $s4_specframe->Label(-text=>'Peaks Dir')->pack(-side=>'right');

	my $s4_plotframe = $step4_bottom_right->Frame()->pack(-side=>'top', -fill=>'x');
	my $s4_plotentry = $s4_plotframe->Entry(-textvariable=>\$output_plot)->pack(-side=>'right');
	my $s4_plotlabel = $s4_plotframe->Label(-text=>'Contour')->pack(-side=>'right');
	
	my $button4 = $step4r->Button(-text=>'Go (4)', -command=>\&step4)->pack(-side=>'bottom', -anchor=>'se');

my $step4_spacer = $frame_steps_bottom->Frame(-height=>'2')->pack(-side=>'top');

my $step5_border = $frame_steps_bottom->Frame(-background=>'black', -borderwidth=>'2')->pack(-side=>'top');
my $step5 = $step5_border->Frame()->pack(-side=>'left', -fill=>'both', -expand=>'1');
my $step5r = $step5_border->Frame()->pack(-side=>'right', -fill=>'both');

	my $step5_top = $step5->Frame()->pack(-fill=>'both', -side=>'top');
	my $step5_label = $step5_top->Label(-text=>'5) Fit Spectra')->pack(-side=>'left');

	my $step5_bottom = $step5->Frame()->pack(-fill=>'both', -side=>'bottom');
	my $step5_bottom_left = $step5_bottom->Frame()->pack(-side=>'left', -fill=>'both');
	my $step5_bottom_right = $step5_bottom->Frame()->pack(-side=>'right', -fill=>'both');

	my $s5_inframe = $step5_bottom_left->Frame()->pack(-side=>'top', -fill=>'x');
	my $s5_inentry = $s5_inframe->Entry(-textvariable=>\$output_batch)->pack(-side=>'right');
	my $s5_inlabel = $s5_inframe->Label(-text=>'Input')->pack(-side=>'right');	

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
	
	my $button5 = $step5r->Button(-text=>'Go (5)', -command=>\&step5)->pack(-side=>'bottom', -anchor=>'se');

my $frame_spacer2 = $frame_left->Frame(-height=>'3')->pack(-side=>'top');

my $frame_footer = $frame_left->Frame(-background=>'white')->pack(-side=>'top', -fill=>'both', -expand=>'1');

my $frame_notify_border = $frame_footer->Frame(-background=>'black', -borderwidth=>'3')->pack(-fill=>'x', -expand=>'1', -anchor=>'nw', -side=>'left');
my $frame_notify = $frame_notify_border->Frame(-background=>'green', -borderwidth=>'5')->pack(-fill=>'both', -expand=>'1', -anchor=>'nw', -side=>'top');
my $frame_notify_label = $frame_notify->Label(-background=>'green', -textvariable=>\$notifymessage)->pack(-fill=>'both', -expand=>'1');

my $quitbutton = $frame_footer->Button(-text=>'Quit', -command=>\&quit)->pack(-side=>'right', -anchor=>'nw');


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
		$input_digest = "$massacre_path/digest/files/30S_miss2_cons_z6_1_0_300_1300.txt";
	}
	elsif($digestselect_val == 2)
	{
		$input_digest = "$massacre_path/digest/files/30S_miss4_cons_z6_1_0_300_1300.txt";
	}
	elsif($digestselect_val == 3)
	{
		$input_digest = "$massacre_path/digest/files/50S_miss2_cons_z6_1_0_300_1300.txt";
	}
	elsif($digestselect_val == 4)
	{
		$input_digest = "$massacre_path/digest/files/50S_miss4_cons_z6_1_0_300_1300.txt";
	}
	elsif($digestselect_val == 5)
	{
		$input_digest = "$massacre_path/digest/files/70S_miss2_cons_z6_1_0_300_1300.txt";
	}
	elsif($digestselect_val == 6)
	{
		$input_digest = "$massacre_path/digest/files/70S_miss4_cons_z6_1_0_300_1300.txt";
	}
	elsif($digestselect_val == 7)
	{
		$input_digest = "$massacre_path/digest/files/30S_miss4_cons_z6_1_0_200_1300.txt";
	}
	elsif($digestselect_val == 8)
	{
		$input_digest = "$massacre_path/digest/files/16S.digest";
	}		
	
	ulog("Setting Digest File:");
	#ulog("$digestselect_text");
	ulog("$input_digest");
}

sub set_model
{
	if($modelselect_val == 0)
	{
		return;
	}
	
	if($modelselect_val == 1)
	{
		$input_model = "model_fix0.993N15.txt";
		$special_range_residues = "";
		$special_range_values = "";
	}
	elsif($modelselect_val == 2)
	{
		$input_model = "model_var0.8N15.txt";
		$special_range_residues = "";
		$special_range_values = "";
	}
	elsif($modelselect_val == 3)
	{
		$input_model = "model_var0.8N15_fix0.993N15.txt";	
		$special_range_residues = "";
		$special_range_values = "";
	}
	elsif($modelselect_val == 4)
	{
		$input_model = "model_fix0.5N15_fix0.993N15.txt	";	
		$special_range_residues = "";
		$special_range_values = "";
	}
	elsif($modelselect_val == 5)
	{
		$input_model = "model_fix0.2N15_fix0.5N15_fix0.8N15.txt	";	
		$special_range_residues = "";
		$special_range_values = "";
	}
	elsif($modelselect_val == 6)
	{
		$input_model = "model_fix0.2N15_fix0.4N15_fix0.6N15_fix0.8N15.txt";	
		$special_range_residues = "";
		$special_range_values = "";
	}
	elsif($modelselect_val == 7)
	{
		$input_model = "model_var0.5H2.txt";	
		$special_range_residues = "";
		$special_range_values = "";
	}
	elsif($modelselect_val == 8)
	{
		$input_model = "model_fix0.98H2var0.5Leu.txt";	
		$special_range_residues = "L";
		$special_range_values = "10";
	}
	elsif($modelselect_val == 9)
	{
		$input_model = "model_fix0.98H2var1.0Leu.txt";	
		$special_range_residues = "L";
		$special_range_values = "10";
	}
	elsif($modelselect_val == 10)
	{
		$input_model = "model_fix0.98H2var0.5Lys.txt";	
		$special_range_residues = "K";
		$special_range_values = "9";
	}
	elsif($modelselect_val == 11)
	{
		$input_model = "model_fix0.98H2var1.0Lys.txt";	
		$special_range_residues = "K";
		$special_range_values = "9";
	}
	elsif($modelselect_val == 12)
	{
		$input_model = "model_fix0.98H2var0.5Val.txt";	
		$special_range_residues = "V";
		$special_range_values = "8";
	}
	elsif($modelselect_val == 13)
	{
		$input_model = "model_fix0.98H2var1.0Val.txt";	
		$special_range_residues = "V";
		$special_range_values = "8";
	}
	elsif($modelselect_val == 14)
	{
		$input_model = "model_fix0.98C13fix0.98H2fix0.98N15var0.5Glu.txt";	
		$special_range_residues = "E";
		$special_range_values = "11";
	}
	elsif($modelselect_val == 15)
	{
		$input_model = "model_fix0.98C13fix0.98H2fix0.98N15var1.0Glu.txt";	
		$special_range_residues = "E";
		$special_range_values = "11";
	}
	elsif($modelselect_val == 16)
	{
		$input_model = "model_fix0.98H2var0.5Leu-fix0.98H2var0.5Val.txt";	
		$special_range_residues = "L_V";
		$special_range_values = "10_8";
	}
	elsif($modelselect_val == 17)
	{
		$input_model = "model_fix0.98H2var1.0Leu-fix0.98H2var1.0Val.txt";	
		$special_range_residues = "L_V";
		$special_range_values = "10_8";
	}
	elsif($modelselect_val == 18)
	{
		$input_model = "model_fix0.98H2var0.5Lys-fix0.98H2var0.5Val.txt";	
		$special_range_residues = "K_V";
		$special_range_values = "9_8";
	}
	elsif($modelselect_val == 19)
	{
		$input_model = "model_fix0.98H2var1.0Lys-fix0.98H2var1.0Val.txt";	
		$special_range_residues = "K_V";
		$special_range_values = "9_8";
	}
	elsif($modelselect_val == 20)
	{
		$input_model = "model_fix0.98H2var0.5Leu-fix0.98H2var0.5Lys-fix0.98H2var0.5Val.txt";	
		$special_range_residues = "L_K_V";
		$special_range_values = "10_9_8";
	}
	elsif($modelselect_val == 21)
	{
		$input_model = "model_fix0.98H2var1.0Leu-fix0.98H2var1.0Lys-fix0.98H2var1.0Val.txt";	
		$special_range_residues = "L_K_V";
		$special_range_values = "10_9_8";
	}
	elsif($modelselect_val == 22)
	{
		$input_model = "model_rna_fix0.993N15.txt";
		$special_range_residues = "";
		$special_range_values = "";
	}
	elsif($modelselect_val == 23)
	{
		$input_model = "model_rna_fix0.5N15.txt";
		$special_range_residues = "";
		$special_range_values = "";
	}
	elsif($modelselect_val == 24)
	{
		$input_model = "model_rna_var0.5N15.txt";
		$special_range_residues = "";
		$special_range_values = "";
	}

	$input_model = join '', "$massacre_path/msc/files/", $input_model;
	
	ulog("Setting Model:");
	ulog("$input_model");
}

sub populate
{
	$output_match0 = join '', $output_root, ".match0";
	$output_rtdir = join '', $output_root, "_rtdir";
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
	# If it detects a MS2 peak list format, run go_msms() instead\
	# 4 = MASCOT
	# 5 = SEQUEST
	if($ailtype == 4 || $ailtype == 5)
	{
		go_msms();
		return;
	}

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

sub go_msms
{
	$time_start = time();
	$program_status = 0;

	populate();
	
	step0ms2();
	
	if($program_status == 0)
	{
		step123ms2();
	}
	else
	{
		ulog("Skipping Step 1,2,3");
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

sub step0ms2
{
	my @ext_array = split /\./, $input_csv;
	my $local_input_csv;
	if($ext_array[$#ext_array] eq "d")
	{
		pop(@ext_array);
		$local_input_csv = join '.', @ext_array;
		$local_input_csv = join '', $local_input_csv, ".mzML.gz";

		if(!(-e $local_input_csv))
		{
			convert_d2mzmlgz();
		}
	}
	elsif($ext_array[$#ext_array] eq "RAW" || $ext_array[$#ext_array] eq "raw")
	{
		pop(@ext_array);
		$local_input_csv = join '.', @ext_array;
		$local_input_csv = join '', $local_input_csv, ".mzML.gz";

		if(!(-e $local_input_csv))
		{
			convert_d2mzmlgz();
		}
	}
	else
	{
		$local_input_csv = $input_csv;
	}

	my $line;
	#my $rtdir = join '_', $output_root, "rtdir";
	
	$program_status = 0;
	if(!(-e $input_ail))
	{
		ulog("Error: No $input_ail");
		$program_status = -1;
		return;
	}

	notify("massacre is running - please wait");
	ulog("Adjusting MSMS RTs");
	
	# This just writes to the output
	ulog("$massacre_path/msc/msc00_adjustmsmsrt.pl --ail=$input_ail --output=$output_match0 --blueid=$blueid --csvfile=$local_input_csv --dir=$output_rtdir --machine=$machine --msms_proteinkeep=$msms_proteinkeep --ailtype=$ailtype");
	
	# This actually calls the script
	open LOG, "$massacre_path/msc/msc00_adjustmsmsrt.pl --ail=$input_ail --output=$output_match0 --blueid=$blueid --csvfile=$local_input_csv --dir=$output_rtdir --machine=$machine --msms_proteinkeep=$msms_proteinkeep --ailtype=$ailtype 2>&1 |" or die "Can't Fork: $!\n";
	
	while(defined($line=<LOG>))
	{
		ulog($line);
	}
	close LOG;

	ulog("Finished Adjusting MSMS RTs");
	stop_notify();
}

sub step123ms2
{
	my $line;

	$program_status = 0;
	#if(!(-e $input_digest)) 
	#{
	#	ulog("Error: No $input_digest");
	#	$program_status = -1;
	#}
	#if(!(-e $input_ail))
	#{
	#	ulog("Error: No $input_ail");
	#	$program_status = -1;
	#}
	
	if(!(-e $output_match0))
	{
		ulog("Error: No $output_match0");
		$program_status = -1;
	}
	
	if($program_status == -1)
	{
		return;
	}

	notify("massacre is running - please wait");
	ulog("Starting Step 123_msms");
	
	#ulog("$ENV{'MASSACRE_PATH'}/msc/msc010203_mascotmsmsmatches.pl --digest=$input_digest --ail=$output_match0 --output=$output_match3");
	ulog("$massacre_path/msc/msc010203_mascotmsmsmatches.pl --ail=$output_match0 --output=$output_match3 --mzrange_low=$mzrange_low --mzrange_high=$mzrange_high --ailtype=$ailtype --rt_max_int_fraction=$rt_max_int_fraction");
	
	#open LOG, "$ENV{'MASSACRE_PATH'}/msc/msc010203_mascotmsmsmatches.pl --digest=$input_digest --ail=$output_match0 --output=$output_match3 2>&1 |" or die "Can't Fork: $!\n";
	open LOG, "$massacre_path/msc/msc010203_mascotmsmsmatches.pl --ail=$output_match0 --output=$output_match3 --mzrange_low=$mzrange_low --mzrange_high=$mzrange_high --ailtype=$ailtype --rt_max_int_fraction=$rt_max_int_fraction 2>&1 |" or die "Can't Fork: $!\n";
	
	while(defined($line=<LOG>))
	{
		ulog($line);
	}
	close LOG;
	
	ulog("Finished Step 123_msms");
	ulog("----------");
	stop_notify();
}

sub step1
{
	my $line;

	$program_status = 0;
	if(!(-e $input_digest)) 
	{
		ulog("Error: No $input_digest");
		$program_status = -1;
	}
	if(!(-e $input_ail))
	{
		ulog("Error: No $input_ail");
		$program_status = -1;
	}
	
	if($program_status == -1)
	{
		return;
	}

	notify("massacre is running - please wait");
	ulog("Starting Step 1");
	ulog("$massacre_path/msc/msc01_findmatches.pl --mode=$matchmode --sample_id=$sample_id --ppm_offset=$match_ppm_offset --ppm_thresh=$match_ppm_thresh --digest=$input_digest --ail=$input_ail --output=$output_match1 --ailtype=$ailtype");

	open LOG, "$massacre_path/msc/msc01_findmatches.pl --mode=$matchmode --id=$sample_id --ppm_offset=$match_ppm_offset --ppm_thresh=$match_ppm_thresh --digest=$input_digest --ail=$input_ail --output=$output_match1 --ailtype=$ailtype 2>&1 |" or die "Can't Fork: $!\n";
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

	$program_status = 0;
	if(!(-e $output_match1)) 
	{
		ulog("Error: No $output_match1");
		$program_status = -1;
	}
	
	if($program_status == -1)
	{
		return;
	}
	
	notify("massacre is running - please wait");
	ulog("Starting Step 2");
	
	if($match_rt_filter == 1)
	{
		ulog("$massacre_path/msc/msc02a_rtfilter.pl --input=$output_match1 --output=$output_match2 --rt_min=$match_rt_minval --rt_max=$match_rt_maxval");
		
		open LOG, "$massacre_path/msc/msc02a_rtfilter.pl --input=$output_match1 --output=$output_match2 --rt_min=$match_rt_minval --rt_max=$match_rt_maxval 2>&1 |" or die "Can't Fork: $!\n";
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

	$program_status = 0;
	if(!(-e $output_match2)) 
	{
		ulog("Error: No $output_match2");
		$program_status = -1;
	}
	
	if($program_status == -1)
	{
		return;
	}

	notify("massacre is running - please wait");
	ulog("Starting Step 3");
	ulog("$massacre_path/msc/msc03_getsinglesorpairs.pl --input=$output_match2 --output=$output_match3 --matchmode=$matchmode --pairmode=$pairmode --rt_thresh=$pair_rt_thresh");
	
	open LOG, "$massacre_path/msc/msc03_getsinglesorpairs.pl --input=$output_match2 --output=$output_match3 --matchmode=$matchmode --pairmode=$pairmode --rt_thresh=$pair_rt_thresh 2>&1 |" or die "Can't Fork: $!\n";
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
	# Splits the $input_csv (raw MS data file) filename on periods
	my @ext_array = split /\./, $input_csv;
	my $local_input_csv;
	
	# If it's a .d file, convert it to .mzML.gz
	if($ext_array[$#ext_array] eq "d")
	{
		pop(@ext_array);
		$local_input_csv = join '.', @ext_array;
		$local_input_csv = join '', $local_input_csv, ".mzML.gz";

		if(!(-e $local_input_csv))
		{
			convert_d2mzmlgz();
		}
	}
	# If it's a .raw file, convert to .mzML.gz
	elsif($ext_array[$#ext_array] eq "RAW" || $ext_array[$#ext_array] eq "raw")
	{
		pop(@ext_array);
		$local_input_csv = join '.', @ext_array;
		$local_input_csv = join '', $local_input_csv, ".mzML.gz";

		if(!(-e $local_input_csv))
		{
			convert_d2mzmlgz();
		}
	}	
	# $local_input_csv is the new name used here after potential conversion
	else
	{
		$local_input_csv = $input_csv;
	}

	my $line;

	# Check to make sure Step 3 happened
	$program_status = 0;
	if(!(-e $output_match3)) 
	{
		ulog("Error: No $output_match3");
		$program_status = -1;
	}
	# Make sure the file actually exists
	if(!(-e $input_csv)) 
	{
		ulog("Error: No $input_csv");
		$program_status = -1;
	}
		
	if($program_status == -1)
	{
		return;
	}

	if($special_range_residues ne "" && $special_range_values ne "" && $set_special_range == 1)
	{
		$special_range = 1;
		ulog("Setting Special Range");
	}
	else
	{
		$special_range = 0;
	}

	notify("massacre is running - please wait");
	ulog("Starting Step 4");
	ulog("$massacre_path/msc/msc04_makepeaks_1local.pl --input=$output_match3 --output=$output_batch --dir=$output_specdir --csv=$local_input_csv --matchmode=$matchmode --filter_abundance=$filter_abundance --filter_mz=$filter_mz --filter_rt=$filter_rt --abundance_min=$abundance_mincut --mz_min=$mz_mincut --mz_max=$mz_maxcut --rt_min=$rt_mincut --rt_max=$rt_maxcut --id=$sample_id --plot=$output_plot --blueid=$blueid --sortstyle=$sortstyle --slice_per_pt=$slice_per_pt --mz_per_pt=$mz_per_pt --rt_expand=$rt_expand --rt_expand_lo=$rt_expand_lo --rt_expand_hi=$rt_expand_hi --special_range=$special_range --special_range_residues=$special_range_residues --special_range_values=$special_range_values --machine=$machine --keepcon=$keepcon --ailtype=$ailtype --extractstyle=$extractstyle --extractbin=$extractbin");

	open LOG, "$massacre_path/msc/msc04_makepeaks_1local.pl --input=$output_match3 --output=$output_batch --dir=$output_specdir --csv=$local_input_csv --matchmode=$matchmode --filter_abundance=$filter_abundance --filter_mz=$filter_mz --filter_rt=$filter_rt --abundance_min=$abundance_mincut --mz_min=$mz_mincut --mz_max=$mz_maxcut --rt_min=$rt_mincut --rt_max=$rt_maxcut --id=$sample_id --plot=$output_plot --blueid=$blueid --sortstyle=$sortstyle --slice_per_pt=$slice_per_pt --mz_per_pt=$mz_per_pt --rt_expand=$rt_expand --rt_expand_lo=$rt_expand_lo --rt_expand_hi=$rt_expand_hi --special_range=$special_range --special_range_residues=$special_range_residues --special_range_values=$special_range_values --machine=$machine --keepcon=$keepcon --ailtype=$ailtype --extractstyle=$extractstyle --extractbin=$extractbin 2>&1 |" or die "Can't Fork: $!\n";
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

	$program_status = 0;
	if(!(-e $output_batch)) 
	{
		ulog("Error: No output_batch: |$output_batch|");
		$program_status = -1;
	}
	if(!(-e $output_specdir)) 
	{
		ulog("Error: No output_specdir: |$output_specdir|");
		$program_status = -1;
	}
	if(!(-e $input_model)) 
	{
		ulog("Error: No input_model: |$input_model|");
		$program_status = -1;
	}
	if(!(-e $input_atoms)) 
	{
		ulog("Error: No input_atoms: |$input_atoms|");
		$program_status = -1;
	}		
	
	if($program_status == -1)
	{
		return;
	}


	notify("massacre is running - please wait");
	ulog("Starting Step 5");
	ulog("$massacre_path/msc/msc05_runisodist_iso3.pl --input=$output_batch --output=$output_fits --dir=$output_specdir --fitsdir=$output_fitsdir --baseline=$baseline --offset=$offset --gw=$gw --blueid=$blueid --model=$input_model --atoms=$input_atoms --isodist_binary=$isodist_binary --machine=$machine --keepfit=$keepfit --baseline_fit_type=$baseline_fit_type --ailtype=$ailtype");

	# Note, we pass $output_specdir here as the new isodist_batch uses the same output directory

	open LOG, "$massacre_path/msc/msc05_runisodist_iso3.pl --input=$output_batch --output=$output_fits --dir=$output_specdir --fitsdir=$output_fitsdir --baseline=$baseline --offset=$offset --gw=$gw --blueid=$blueid --model=$input_model --atoms=$input_atoms --isodist_binary=$isodist_binary --machine=$machine --keepfit=$keepfit --baseline_fit_type=$baseline_fit_type --ailtype=$ailtype 2>&1 |" or die "Can't Fork: $!\n";
	
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
	
	my $cw_top = $cw->Frame()->pack(-side=>'top', -fill=>'both');
	my $cw_tbspacer = $cw->Frame(-height=>'3')->pack(-side=>'top');
	my $cw_bot = $cw->Frame()->pack(-side=>'top', -fill=>'both');
	my $cw_topl = $cw_top->Frame()->pack(-side=>'left', -anchor=>'w', -fill=>'both');
	my $cw_topr = $cw_top->Frame()->pack(-side=>'right', -anchor=>'e', -fill=>'both');
	
	my $cw_s1_border = $cw_topl->Frame(-border=>'3', -background=>'black')->pack(-side=>'top', -fill=>'both', -anchor=>'w');
	my $cw_s1 = $cw_s1_border->Frame(-border=>'2')->pack(-side=>'top', -fill=>'both');
	my $cw_s1_labelframe = $cw_s1->Frame()->pack(-side=>'top', -fill=>'both');
	my $cw_s1_label = $cw_s1_labelframe->Label(-text=>'1) Compare ail with Theoretical Digest')->pack(-side=>'left');

	my $cw_s1_varframe1 = $cw_s1->Frame(-pady=>'2')->pack(-side=>'top', -fill=>'both');
		my $cw_s1_matchbutton0 = $cw_s1_varframe1->Radiobutton(-text=>'N14 + N15', -variable=>\$matchmode, -value=>'0')->pack(-side=>'left');
		my $cw_s1_matchbutton1 = $cw_s1_varframe1->Radiobutton(-text=>'N14 Only', -variable=>\$matchmode, -value=>'1')->pack(-side=>'left');
		my $cw_s1_matchbutton2 = $cw_s1_varframe1->Radiobutton(-text=>'N15 Only', -variable=>\$matchmode, -value=>'2')->pack(-side=>'left');

	my $cw_s1_varframe2a = $cw_s1->Frame(-pady=>'2')->pack(-side=>'top', -fill=>'both');
		my $cw_s2_var2alabel = $cw_s1_varframe2a->Label(-text=>'PPM Offset for Matches:')->pack(-side=>'left');
		my $cw_s2_var2aval = $cw_s1_varframe2a->Entry(-textvariable=>\$match_ppm_offset)->pack(-side=>'left');

	my $cw_s1_varframe2 = $cw_s1->Frame(-pady=>'2')->pack(-side=>'top', -fill=>'both');
		my $cw_s1_var2label = $cw_s1_varframe2->Label(-text=>'PPM Threshold for Matches:')->pack(-side=>'left');
		my $cw_s2_var2val = $cw_s1_varframe2->Entry(-textvariable=>\$match_ppm_thresh)->pack(-side=>'left');

	my $cw_s1_varframe3 = $cw_s1->Frame(-pady=>'2')->pack(-side=>'top', -fill=>'both');
		my $cw_s1_var3label = $cw_s1_varframe3->Label(-text=>'Sample ID:')->pack(-side=>'left');
		my $cw_s2_var3val = $cw_s1_varframe3->Entry(-textvariable=>\$sample_id)->pack(-side=>'left');

	my $cw_spacer1 = $cw_topl->Frame(-height=>'3')->pack(-side=>'top');

	my $cw_s2_border = $cw_topl->Frame(-border=>'3', -background=>'black')->pack(-side=>'top', -fill=>'both', -anchor=>'w');
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

	my $cw_spacer2 = $cw_topl->Frame(-height=>'3')->pack(-side=>'top');

	my $cw_s3_border = $cw_topl->Frame(-border=>'3', -background=>'black')->pack(-side=>'top', -fill=>'both', -anchor=>'w');
	my $cw_s3 = $cw_s3_border->Frame(-border=>'2')->pack(-side=>'top', -fill=>'both');
	my $cw_s3_labelframe = $cw_s3->Frame()->pack(-side=>'top', -fill=>'both');
	my $cw_s3_label = $cw_s3_labelframe->Label(-text=>'3) Find Good Features or Feature Pairs')->pack(-side=>'left');

	my $cw_s3_varframe1 = $cw_s3->Frame(-pady=>'2')->pack(-side=>'top', -anchor=>'nw');
		my $cw_s3_matchbutton0 = $cw_s3_varframe1->Radiobutton(-text=>'N14 + N15', -variable=>\$matchmode, -value=>'0')->pack(-side=>'left');
		my $cw_s3_matchbutton1 = $cw_s3_varframe1->Radiobutton(-text=>'N14 Only', -variable=>\$matchmode, -value=>'1')->pack(-side=>'left');
		my $cw_s3_matchbutton2 = $cw_s3_varframe1->Radiobutton(-text=>'N15 Only', -variable=>\$matchmode, -value=>'2')->pack(-side=>'left');

	my $cw_s3_varframe2 = $cw_s3->Frame(-pady=>'2')->pack(-side=>'top', -anchor=>'nw');
		my $cw_s3_pairtypebutton0 = $cw_s3_varframe2->Radiobutton(-text=>'Non-Redundant Only', -variable=>\$pairmode, -value=>'0')->pack(-side=>'top', -anchor=>'nw');
		my $cw_s3_pairtypebutton1 = $cw_s3_varframe2->Radiobutton(-text=>'Same Protein Only Redundant OK', -variable=>\$pairmode, -value=>'1')->pack(-side=>'top', -anchor=>'nw');
		my $cw_s3_pairtypebutton2 = $cw_s3_varframe2->Radiobutton(-text=>'Redundant OK (Keep All -> Filter in Massive)', -variable=>\$pairmode, -value=>'2')->pack(-side=>'top', -anchor=>'nw');

	my $cw_s3_varframe3 = $cw_s3->Frame(-pady=>'2')->pack(-side=>'top', -anchor=>'nw');
		my $cw_s3_var3label1 = $cw_s3_varframe3->Label(-text=>'Pairs Must be Within')->pack(-side=>'left');
		my $cw_s3_var3entry = $cw_s3_varframe3->Entry(-textvariable=>\$pair_rt_thresh, -width=>'5')->pack(-side=>'left');
		my $cw_s3_var3label2 = $cw_s3_varframe3->Label(-text=>'Minutes')->pack(-side=>'left');
	
	my $cw_msms_s0_border = $cw_topr->Frame(-border=>'3', -background=>'black')->pack(-side=>'top', -fill=>'both', -anchor=>'ne');
	my $cw_msms_s0 = $cw_msms_s0_border->Frame(-border=>'2')->pack(-side=>'top', -fill=>'both');
	my $cw_msms_s0_labelframe = $cw_msms_s0->Frame()->pack(-side=>'top', -fill=>'both');
	my $cw_msms_s0_label = $cw_msms_s0_labelframe->Label(-text=>'0) Perform Retention Time Correction')->pack(-side=>'left');
	
	my $cw_msms_s0_varframe1 = $cw_msms_s0->Frame(-pady=>'2')->pack(-side=>'top', -anchor=>'nw');
		my $cw_msms_s0_keepbutton0 = $cw_msms_s0_varframe1->Radiobutton(-text=>'Only Keep Verified Ribosomal Proteins', -variable=>\$msms_proteinkeep, -value=>'0')->pack(-side=>'top', -anchor=>'nw');
		my $cw_msms_s0_keepbutton1 = $cw_msms_s0_varframe1->Radiobutton(-text=>'Only Keep Verified Proteins', -variable=>\$msms_proteinkeep, -value=>'1')->pack(-side=>'top', -anchor=>'nw');
		my $cw_msms_s0_keepbutton2 = $cw_msms_s0_varframe1->Radiobutton(-text=>'Keep All Proteins', -variable=>\$msms_proteinkeep, -value=>'2')->pack(-side=>'top', -anchor=>'nw');

	my $cw_spacer3 = $cw_topr->Frame(-height=>'3')->pack(-side=>'top');

	my $cw_msms_s1_border = $cw_topr->Frame(-border=>'3', -background=>'black')->pack(-side=>'top', -fill=>'both', -anchor=>'ne');
	my $cw_msms_s1 = $cw_msms_s1_border->Frame(-border=>'2')->pack(-side=>'top', -fill=>'both');
	my $cw_msms_s1_labelframe = $cw_msms_s1->Frame()->pack(-side=>'top', -fill=>'both');
	my $cw_msms_s1_label = $cw_msms_s1_labelframe->Label(-text=>'1,2,3) Get Peptide Data from Digest and Consolidate')->pack(-side=>'left');
	
	my $cw_msms_s1_varframe1 = $cw_msms_s1->Frame(-pady=>'2')->pack(-side=>'top', -anchor=>'nw', -fill=>'both');
		my $cw_msms_s1_var1entry2 = $cw_msms_s1_varframe1->Entry(-textvariable=>\$mzrange_high, -width=>'5')->pack(-side=>'right');
		my $cw_msms_s1_var1label2 = $cw_msms_s1_varframe1->Label(-text=>'M/Z Range High:')->pack(-side=>'right');
		my $cw_msms_s1_var1label1 = $cw_msms_s1_varframe1->Label(-text=>'M/Z Range Low:')->pack(-side=>'left');
		my $cw_msms_s1_var1entry1 = $cw_msms_s1_varframe1->Entry(-textvariable=>\$mzrange_low, -width=>'5')->pack(-side=>'left');

	my $cw_msms_s1_varframe2 = $cw_msms_s1->Frame(-pady=>'2')->pack(-side=>'top', -anchor=>'nw', -fill=>'both');
		my $cw_msms_s1_var2entry = $cw_msms_s1_varframe2->Entry(-textvariable=>\$rt_max_int_fraction, -width=>'5')->pack(-side=>'right');
		my $cw_msms_s1_var2label = $cw_msms_s1_varframe2->Label(-text=>'Sequencing Scan Intensity Threshold:')->pack(-side=>'right');
	
	#####
	
	my $cw_s4_border = $cw->Frame(-border=>'3', -background=>'black')->pack(-side=>'top', -fill=>'both');
	my $cw_s4 = $cw_s4_border->Frame(-border=>'2')->pack(-side=>'top', -fill=>'both');
	my $cw_s4_labelframe = $cw_s4->Frame()->pack(-side=>'top', -fill=>'both');
	my $cw_s4_label = $cw_s4_labelframe->Label(-text=>'4) Extract the Minispectra')->pack(-side=>'left');

	my $cw_s4_varframe1 = $cw_s4->Frame(-pady=>'2')->pack(-side=>'top', -fill=>'both');
		my $cw_s4_matchbutton0 = $cw_s4_varframe1->Radiobutton(-text=>'N14 + N15', -variable=>\$matchmode, -value=>'0')->pack(-side=>'left');
		my $cw_s4_matchbutton1 = $cw_s4_varframe1->Radiobutton(-text=>'N14 Only', -variable=>\$matchmode, -value=>'1')->pack(-side=>'left');
		my $cw_s4_matchbutton2 = $cw_s4_varframe1->Radiobutton(-text=>'N15 Only', -variable=>\$matchmode, -value=>'2')->pack(-side=>'left');
    
    my $cw_s4_varframe1b = $cw_s4->Frame(-pady=>'2')->pack(-side=>'top', -fill=>'both');
    my $cw_s4_extractlabel = $cw_s4_varframe1b->Label(-text=>'Extraction Style:')->pack(-side=>'left');
        my $cw_s4_extractbutton0 = $cw_s4_varframe1b->Radiobutton(-text=>'Standard', -variable=>\$extractstyle, -value=>'0')->pack(-side=>'left');
        my $cw_s4_extractbutton0 = $cw_s4_varframe1b->Radiobutton(-text=>'Bin', -variable=>\$extractstyle, -value=>'1')->pack(-side=>'left');
        my $cw_s4_extractbin = $cw_s4_varframe1b->Entry(-textvariable=>\$extractbin, -width=>'5')->pack(-side=>'left');
        my $cw_s4_extractbutton0 = $cw_s4_varframe1b->Radiobutton(-text=>'3D', -variable=>\$extractstyle, -value=>'2')->pack(-side=>'left');

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
		my $cw_s4_var5label1 = $cw_s4_varframe5->Label(-text=>'RT Extraction Radius:')->pack(-side=>'left');
		my $cw_s4_var5entry1 = $cw_s4_varframe5->Entry(-textvariable=>\$rt_expand, -width=>'5')->pack(-side=>'left');
		my $cw_s4_var5label2 = $cw_s4_varframe5->Label(-text=>'+ Lo RT:')->pack(-side=>'left');
		my $cw_s4_var5entry2 = $cw_s4_varframe5->Entry(-textvariable=>\$rt_expand_lo, -width=>'5')->pack(-side=>'left');
		my $cw_s4_var5label3 = $cw_s4_varframe5->Label(-text=>'+ Hi RT:')->pack(-side=>'left');
		my $cw_s4_var5entry3 = $cw_s4_varframe5->Entry(-textvariable=>\$rt_expand_hi, -width=>'5')->pack(-side=>'left');
		
		my $cw_s4_var5flag = $cw_s4_varframe5->Checkbutton(-text=>"Set Special Range:", -variable=>\$set_special_range)->pack(-side=>'left');
		my $cw_s4_var5label4 = $cw_s4_varframe5->Label(-text=>'M/Z Range Residues')->pack(-side=>'left');
		my $cw_s4_var5entry4 = $cw_s4_varframe5->Entry(-textvariable=>\$special_range_residues)->pack(-side=>'left');
		my $cw_s4_var5label5 = $cw_s4_varframe5->Label(-text=>'Values')->pack(-side=>'left');
		my $cw_s4_var5entry5 = $cw_s4_varframe5->Entry(-textvariable=>\$special_range_values)->pack(-side=>'left');

		
	my $cw_s4_varframe6 = $cw_s4->Frame(-pady=>'4')->pack(-side=>'top', -fill=>'both');
		my $cw_s4_var6label1 = $cw_s4_varframe6->Label(-text=>'Contour Plot: RTSlices/pt')->pack(-side=>'left');
		my $cw_s4_var6entry1 = $cw_s4_varframe6->Entry(-textvariable=>\$slice_per_pt, -width=>'5')->pack(-side=>'left');
		my $cw_s4_var6label2 = $cw_s4_varframe6->Label(-text=>'MZ/pt')->pack(-side=>'left');
		my $cw_s4_var6entry2 = $cw_s4_varframe6->Entry(-textvariable=>\$mz_per_pt, -width=>'5')->pack(-side=>'left');
		my $cw_s4_var6checkbutton3 = $cw_s4_varframe6->Checkbutton(-text=>'Keep Raw Minicontour Data', -variable=>\$keepcon)->pack(-side=>'left');
    #my $cw_s4_var6checkbutton4 = $cw_s4_varframe6->Checkbutton(-text=>'Extract 3D Minispectra', -variable=>\$spectra3D)->pack(-side=>'left');


	my $cw_spacer4 = $cw->Frame(-height=>'3')->pack(-side=>'top');

	my $cw_s5_border = $cw->Frame(-border=>'3', -background=>'black')->pack(-side=>'top', -fill=>'both');
	my $cw_s5 = $cw_s5_border->Frame(-border=>'2')->pack(-side=>'top', -fill=>'both');
	my $cw_s5_labelframe = $cw_s5->Frame()->pack(-side=>'top', -fill=>'both');
	my $cw_s5_label = $cw_s5_labelframe->Label(-text=>'5) Fit Spectra')->pack(-side=>'left');

	my $cw_s5_modelframe = $cw_s5->Frame()->pack(-side=>'top', -fill=>'x');


	my $modellabel = $cw_s5_modelframe->Label(-text=>'Model:', -width=>'6')->pack(-side=>'left');
	my $modelentry = $cw_s5_modelframe->Entry(-textvariable=>\$input_model)->pack(-side=>'left');

	my $modelselectbutton = $cw_s5_modelframe->Button(-text=>'Select Model', -command=>\&select_model)->pack(-side=>'left');
	
#	my $modelselect = $cw_s5_modelframe->Optionmenu(-width=>'15', -options=>
#	[
#	["Select a Model"=>0],
#	["U + fix 0.993 N15"=>1],
#	["U + var 0.8 N15"=>2],
#	["U + var 0.8 N15 + fix 0.993 N15"=>3],
#	["U + fix 0.5 N15 + fix 0.993 N15"=>4],
#	["U + fix 0.2 N15 + fix 0.5 N15 + fix 0.8 N15"=>5],
#	["U + fix 0.2 N15 + fix 0.4 N15 + fix 0.6 N15 + fix 0.8 N15"=>6],
#	["U + var 0.5 H2"=>7],
#	["U + fix 0.98 H2 var 0.5 Leu"=>8],
#	["U + fix 0.98 H2 var 1.0 Leu"=>9],
#	["U + fix 0.98 H2 var 0.5 Lys"=>10],
#	["U + fix 0.98 H2 var 1.0 Lys"=>11],
#	["U + fix 0.98 H2 var 0.5 Val"=>12],
#	["U + fix 0.98 H2 var 1.0 Val"=>13],
#	["U + fix 0.98 C13 fix 0.98 H2 fix 0.98 N15 var 0.5 Glu"=>14],
#	["U + fix 0.98 C13 fix 0.98 H2 fix 0.98 N15 var 1.0 Glu"=>15],
#	["U + (fix 0.98 H2 var 0.5 Leu, fix 0.98 H2 var 0.5 Val)"=>16],
#	["U + (fix 0.98 H2 var 1.0 Leu, fix 0.98 H2 var 1.0 Val)"=>17],
#	["U + (fix 0.98 H2 var 0.5 Lys, fix 0.98 H2 var 0.5 Val)"=>18],
#	["U + (fix 0.98 H2 var 1.0 Lys, fix 0.98 H2 var 1.0 Val)"=>19],
#	["U + (fix 0.98 H2 var 0.5 Leu, fix 0.98 H2 var 0.5 Lys, fix 0.98 H2 var 0.5 Val)"=>20],
#	["U + (fix 0.98 H2 var 1.0 Leu, fix 0.98 H2 var 1.0 Lys, fix 0.98 H2 var 1.0 Val)"=>21],
#	["RNA U + fix 0.993 N15"=>22],
#	["RNA U + fix 0.5 N15"=>23],
#	["RNA U + var 0.5 N15"=>24]	
#	], 
#	-command=>\&set_model, -textvariable=>\$modelselect_text, -variable=>\$modelselect_val)->pack(-side=>'left');

	my $cw_s5_atomframe = $cw_s5->Frame()->pack(-side=>'top', -fill=>'x');
	my $atomslabel = $cw_s5_atomframe->Label(-text=>'Atoms:', -width=>'6')->pack(-side=>'left');
	my $atomsentry = $cw_s5_atomframe->Entry(-textvariable=>\$input_atoms, -width=>'40')->pack(-side=>'left');
	my $atomsnote = $cw_s5_atomframe->Label(-text=>'OK 4 All')->pack(-side=>'left');

	my $isodist_binary_entry = $cw_s5_atomframe->Entry(-textvariable=>\$isodist_binary)->pack(-side=>'right');

	my $cw_s5_varframe2 = $cw_s5->Frame(-pady=>'2')->pack(-side=>'top', -fill=>'both');
	my $cw_s5_varframe2left = $cw_s5_varframe2->Frame(-background=>'red')->pack(-side=>'left');

		my $cw_s5_fitvar1frame = $cw_s5_varframe2left->Frame()->pack(-side=>'top', -fill=>'both');
		 $cw_s5_fitvar1entry = $cw_s5_fitvar1frame->Entry(-textvariable=>\$baseline, -width=>'6')->pack(-side=>'right');
		my $cw_s5_fitvar1label = $cw_s5_fitvar1frame->Label(-text=>'Baseline')->pack(-side=>'right');
		
		

		my $cw_s5_fitvar2frame = $cw_s5_varframe2left->Frame()->pack(-side=>'top', -fill=>'both');
		 $cw_s5_fitvar2entry = $cw_s5_fitvar2frame->Entry(-textvariable=>\$offset, -width=>'6')->pack(-side=>'right');
		my $cw_s5_fitvar2label = $cw_s5_fitvar2frame->Label(-text=>'Offset (var)')->pack(-side=>'right');

		my $cw_s5_fitvar3frame = $cw_s5_varframe2left->Frame()->pack(-side=>'top', -fill=>'both');
		 $cw_s5_fitvar3entry = $cw_s5_fitvar3frame->Entry(-textvariable=>\$gw, -width=>'6')->pack(-side=>'right');
		my $cw_s5_fitvar3label = $cw_s5_fitvar3frame->Label(-text=>'Gaussian Width (var)')->pack(-side=>'right');

	my $cw_s5_varframe2center = $cw_s5_varframe2->Frame()->pack(-side=>'left', -fill=>'both');


		my $cw_s5_fitvar4frame = $cw_s5_varframe2center->Frame(-background=>'red')->pack(-side=>'top', -fill=>'both');

		my $cw_s5_baselinefittype0 = $cw_s5_fitvar4frame->Radiobutton(-text=>'Auto', -variable=>\$baseline_fit_type, -value=>'0')->pack(-side=>'left', -fill=>'both');
		my $cw_s5_baselinefittype1 = $cw_s5_fitvar4frame->Radiobutton(-text=>'Fixed', -variable=>\$baseline_fit_type, -value=>'1')->pack(-side=>'left', -fill=>'both');


		#my $cw_s5_fitvar5frame = $cw_s5_varframe2center->Frame()->pack(-side=>'top', -fill=>'both');

		#my $cw_s5_fitvar6frame = $cw_s5_varframe2center->Frame()->pack(-side=>'top', -fill=>'both');

	my $cw_s5_varframe2right = $cw_s5_varframe2->Frame()->pack(-side=>'right');
		my $cw_s5_keepfit = $cw_s5_varframe2right->Checkbutton(-text=>'Keep .fit/.dat Files', -variable=>\$keepfit)->pack(-side=>'top');

	my $cw_s5_varframe3 = $cw_s5->Frame(-pady=>'2')->pack(-side=>'top', -fill=>'both');

	my $cw_spacer5 = $cw->Frame(-height=>'3')->pack(-side=>'top');

	my $cw_close_border = $cw->Frame(-border=>'3', -background=>'black')->pack(-side=>'top', -fill=>'both');
	my $cw_close = $cw_close_border->Frame(-border=>'2')->pack(-side=>'top', -fill=>'both');
	my $cw_closebutton = $cw_close->Button(-text=>'Save/Close', -command=>\&close_config)->pack(-side=>'right', -anchor=>'e');
}

sub select_digest
{
	my $digest_fs = $mw->FileSelect(-directory=>"$massacre_path/digest/files");
	$digest_fs->minsize( qw(800 480));
	$input_digest = $digest_fs->Show;
}

sub select_model
{
	my $model_fs = $mw->FileSelect(-directory=>"$massacre_path/msc/files");
	$model_fs->minsize( qw(800 480));
	$input_model = $model_fs->Show;
}

sub close_config
{
	$cw->destroy();
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
		elsif($param eq "match_ppm_offset")
		{
			$match_ppm_offset = $value;
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
		elsif($param eq "slice_per_pt")
		{
			$slice_per_pt = $value;
		}
		elsif($param eq "mz_per_pt")
		{
			$mz_per_pt = $value;
		}
		elsif($param eq "rt_expand")
		{
			$rt_expand = $value;
		}		
		elsif($param eq "rt_expand_lo")
		{
			$rt_expand_lo = $value;
		}
		elsif($param eq "rt_expand_hi")
		{
			$rt_expand_hi = $value;
		}
		elsif($param eq "special_range")
		{
			$special_range = $value;
		}
		elsif($param eq "special_range_residues")
		{
			$special_range_residues = $value;
		}
		elsif($param eq "special_range_values")
		{
			$special_range_values = $value;
		}
		elsif($param eq "set_special_range")
		{
			$set_special_range = $value;
		}
		elsif($param eq "keepcon")
		{
			$keepcon = $value;
		}
		# If we comment out the entire elsif statement we throw an error
        # We want to ignore this because as of 6/20/11 nobody has actually used spectra3D, and that's why it was removed and replaced by extractstyle
        elsif($param eq "spectra3D")
		{
		#	$spectra3D = $value;
		}
        elsif($param eq "extractstyle")
        {
            $extractstyle = $value;
        }
        elsif($param eq "extractbin")
        {
            $extractbin = $value;
        }
		
		# Step 5
		elsif($param eq "baseline")
		{
			$baseline = $value;
		}
		elsif($param eq "baseline_fit_type")
		{
			$baseline_fit_type = $value;
		}
		elsif($param eq "offset")
		{
			$offset = $value;
		}
		elsif($param eq "gw")
		{
			$gw = $value;
		}
		elsif($param eq "model")
		{
			$input_model = $value;
		}
		elsif($param eq "atoms")
		{
			$input_atoms = $value;
		}
		elsif($param eq "keepfit")
		{
			$keepfit = $value;
		}
		
		# Step 0 MSMS
		elsif($param eq "msms_proteinkeep")
		{
			$msms_proteinkeep = $value;
		}

		# Step 1,2,3 MSMS
		elsif($param eq "mzrange_low")
		{
			$mzrange_low = $value;
		}
		elsif($param eq "mzrange_high")
		{
			$mzrange_high = $value;
		}
		elsif($param eq "rt_max_int_fraction")
		{
			$rt_max_int_fraction = $value;			
		}
		
		
		else
		{
			ulog("Error: Can't find parameter $param");
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
	print PARAM "match_ppm_offset $match_ppm_offset\n";
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
	print PARAM "slice_per_pt $slice_per_pt\n";
	print PARAM "mz_per_pt $mz_per_pt\n";
	print PARAM "rt_expand $rt_expand\n";
	print PARAM "rt_expand_lo $rt_expand_lo\n";
	print PARAM "rt_expand_hi $rt_expand_hi\n";
	print PARAM "special_range $special_range\n";
	print PARAM "special_range_residues $special_range_residues\n";
	print PARAM "special_range_values $special_range_values\n";
	print PARAM "set_special_range $set_special_range\n";
	print PARAM "keepcon $keepcon\n";
	#print PARAM "spectra3D $spectra3D\n";
    print PARAM "extractstyle $extractstyle\n";
    print PARAM "extractbin $extractbin\n";
	
	# Step 5
	print PARAM "baseline $baseline\n";
	print PARAM "baseline_fit_type $baseline_fit_type\n";
	print PARAM "offset $offset\n";
	print PARAM "gw $gw\n";
	print PARAM "model $input_model\n";
	print PARAM "atoms $input_atoms\n";
	print PARAM "keepfit $keepfit\n";
	
	# Step 0 MSMS
	print PARAM "msms_proteinkeep $msms_proteinkeep\n";
	
	# Step 1,2,3 MSMS
	print PARAM "mzrange_low $mzrange_low\n";
	print PARAM "mzrange_high $mzrange_high\n";
	print PARAM "rt_max_int_fraction $rt_max_int_fraction\n";
	
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

# Read in .match3 file and generate histogram of ppm values
# Does Steps 1-3, and then calls calc_match_histo
sub go_match_histo
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
	
	calc_match_histo();
	
	if($program_status == -1)
	{
		ulog("Error: Program failed to complete.  Please check log for details.");
	}
	
	$time_end = time();
	$elapsed_minutes = ($time_end - $time_start) / 60;
	ulog("Total Time Elapsed: $elapsed_minutes minutes");
}

sub calc_match_histo
{
	open MATCH3, "$output_match3" or ulog("Error, could not open $output_match3.");
	my $match3_header = <MATCH3>;
	my @m3_harray = split ' ', $match3_header;
	my %m3_name2col = ();
	for(my $m3_ctr = 0; $m3_ctr <= $#m3_harray; $m3_ctr++)
	{
		$m3_name2col{$m3_harray[$m3_ctr]} = $m3_ctr;
	}
	close MATCH3;
	
	my $histo_offset = $match_ppm_offset - $match_ppm_thresh;
	if($histo_offset > 0)
	{
		$histo_offset = 0;
	}
	$histo_offset = abs($histo_offset);
	
	my $histo_numbins = ($match_ppm_offset + $match_ppm_thresh) - ($match_ppm_offset - $match_ppm_thresh) + 2;
		
	#`$ENV{'MASSACRE_PATH'}/bin/histo.pl --numbins=101 --binsize=1 --offset=50 --header=1 --colname=ppm < $output_match3 > $output_match3.histoa`;
	`$massacre_path/bin/histo.pl --numbins=$histo_numbins --binsize=1 --offset=$histo_offset --header=1 --colname=ppm --xval=1 < $output_match3 > $output_match3.histo`;
	print GIN1 "set style data histogram\n";
	print GIN1 "set style histogram gap 0\n";
	print GIN1 "set xtic rotate by 90\n";

	#print GIN1 "g(x) = a*exp( (-1*(x-m)**2) / (2*w**2) ) / (w*sqrt(2*pi))\n";
	#print GIN1 "m=5\n";
	#print GIN1 "a=100\n";
	#print GIN1 "w=5\n";
	#print GIN1 "fit g(x) \"$output_match3.histoa\" via m,a,w\n";

	print GIN1 "plot \"$output_match3.histo\" using 1:xtic(2)\n";
}

sub read_gout1
{
	my $buffer = <GOUT1>;
	print "$buffer\n";
}

sub read_gerr1
{
	my $buffer = <GERR1>;
	print "$buffer\n";
}

# More file naming stuff where massacre tries to be clever and help the user by automatically naming things.
sub newailtype
{
	# For MassProfiler AIL, use .csv
	if($ailtype == 0)
	{
		$input_csv = join '', $output_root, ".csv";
	}
	# For all others, use .mzdata.xml
	elsif($ailtype > 0)
	{
		#my @dir_ls = `ls | egrep -i "mzdata.xml|mzxml|mzml|mzml.gz|\.d"`;
		#chomp(@dir_ls);
		#$input_csv = $dir_ls[$#dir_ls];
		autodetect_msdata();
	}
}

sub convert_d2mzmlgz
{
	my $line;

	$program_status = 0;
	if(!(-e $input_csv))
	{
		ulog("Error: No $input_csv");
		$program_status = -1;
		return;
	}

	notify("massacre is running - please wait");
	ulog("Converting .d/.raw to .mzML.gz");
	
	ulog("$massacre_path/msc/mscM1_d2mzmlgz.pl --input=$input_csv");
	
	open LOG, "$massacre_path/msc/mscM1_d2mzmlgz.pl --input=$input_csv 2>&1 |" or die "Can't Fork: $!\n";
	
	while(defined($line=<LOG>))
	{
		ulog($line);
	}
	close LOG;

	ulog("Finished Converting .d/.raw to .mzML.gz");
	stop_notify();
}

sub autodetect_msdata
{
	my @dir_ls;
	my @dir_lsB;

	if(@dir_ls = `ls | egrep -i .mzdata.xml`)
	{
		chomp(@dir_ls);
		$datatype_flag = 1;
		$input_csv = shift(@dir_ls);
		print "Detected .mzdata.xml $input_csv\n";	
	}
	elsif(@dir_ls = `ls | egrep -i .mzml.gz`)
	{
		chomp(@dir_ls);
		$datatype_flag = 1;
		$input_csv = shift(@dir_ls);
		print "Detected .mzML.gz $input_csv\n";
		
		if(@dir_lsB = `ls *.sequest`)
		{
			$datatype_flag = 2;
		}
	}
	elsif(@dir_ls = `ls | egrep -i .mzml`)
	{
		chomp(@dir_ls);
		$datatype_flag = 1;
		$input_csv = shift(@dir_ls);
		print "Detected .mzML $input_csv\n";
	}
	#elsif(@dir_ls = `ls | egrep -i '\.d'`)
	elsif(@dir_ls = `ls -d *.d`)
	{
		chomp(@dir_ls);
		$datatype_flag = 1;
		$input_csv = shift(@dir_ls);
		print "Detected .d $input_csv\n";
	}
	elsif(@dir_ls = `ls *.raw *.RAW`)
	{
		chomp(@dir_ls);
		$datatype_flag = 2;
		$input_csv = shift(@dir_ls);
		print "Detected .raw $input_csv\n";
	}

}

################################################################################
################################################################################

&GetOptions("param=s" => \$paramfile_input, "root=s" => \$output_root_input, "go=i" => \$go, "digest=s" => \$input_digest);

# 1) Check for a .mzdata.xml file
# 2) Check for a .mzML.gz file
# 3) Check for a .mzML file
# 4) Check for a .d file

autodetect_msdata();

# This is just a bunch of code to try to automatically detect the right file formats and data types and name things properly

# Check for a .mzdata.xml file
# If present, change ailtype
#my @dir_ls = `ls *.mzdata.xml`;
#my @dir_ls = `ls | egrep -i "mzdata.xml|mzxml"`;
#chomp(@dir_ls);
#if($#dir_ls >= 0)
if($datatype_flag == 1) # We think we have MASCOT, let's name stuff accordingly
{
	#print "mzXML file detected ($dir_ls[$#dir_ls]), adjusting input\n";

	$ailtype = 4; # Sets the MASCOT $ailtype
	$ailtypetext = "Mascot MS/MS IDs";
	#newailtype();
	#$input_csv = $dir_ls[$#dir_ls];
	
	my @dir_ls2 = `ls *.xls *.mascot`;
	chomp(@dir_ls2);
	if($#dir_ls2 >= 0)
	{
		$input_ail = $dir_ls2[$#dir_ls2];
	}
}
elsif($datatype_flag == 2)
{
	$ailtype = 5;
	$ailtypetext = "Sequest MS/MS IDs";
	
	my @dir_ls2 = `ls *.sequest`;
	chomp(@dir_ls2);
	if($#dir_ls2 >= 0)
	{
		$input_ail = $dir_ls2[$#dir_ls2];
	}
}

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
if($go > 0)
{
	if($go == 1)
	{
		go();
	}
	elsif($go == 11)
	{
		step1();
	}
	elsif($go == 12)
	{
		step2();
	}
	elsif($go == 13)
	{
		step3();
	}
	elsif($go == 14)
	{
		step4();
	}
	elsif($go == 15)
	{
		step5();
	}
	elsif($go == 20)
	{
		step0ms2();
	}
	elsif($go = 21)
	{
		step123ms2();
	}
}

MainLoop();

