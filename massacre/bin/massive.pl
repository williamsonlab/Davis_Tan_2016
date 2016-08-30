#!/usr/bin/perl

#isofile,isoprotein,isopep,isomf,mw,isoz_charge,tim,chisq,b,b_err,ul_amp,ul_amp_err,l_amp,l_amp_err,off,off_err,gw,gw_err,frac_n,frac_n_err,frac_lab,frac_lab_err,symb,mz,isoid,matchmode,ailid_n14,ailid_n15,featureid_n14,featureid_n15,multimer_n14,multimer_n15,ailiontype_n14,ailiontype_n15,rt_n14,rt_n15,mz_n14,mz_n15,abundance_n14,abundance_n15,ppm_n14,ppm_n15,ailcharge,sampleid,digestionid,digestpepid,n14mass,n15mass,protein,startres,endres,iontype,charge,missed,seq,mod,n14prox,n15prox,seqmod,file,comment

use Getopt::Long;
require Tk;
use Tk;
require Tk::Image;
use Tk::Image;
use Tk::PNG;
use IPC::Open3;

require "$ENV{'MASSACRE_PATH'}/msc/mikesubs.pl";
require "$ENV{'MASSACRE_PATH'}/msc/mikevars.pl";

$dateline = `ls -ltT $ENV{'MASSACRE_PATH'}/bin/massive.pl`;
@datearray = split ' ', $dateline;
$date_display = join ' ', $datearray[5], $datearray[6], $datearray[7], $datearray[8];
$title_display = join '', "massive: Mass Spectrometry: Interactive Viewing Environment - ", $date_display;

#############
# Variables #
#############
$sortstyle = 0; # No sorting on launch
$infile = "";
$outfile = "massfilt.csv";
$plotframe = 2; # Defaults to 70S
$last_isoid = -1; # isoid of the last thing selected
$pid = $$;
#$pngtemp = join '', "pngtemp.", $pid, ".png";
$display_abx = 0; # Set to 1 for variable A/B/X protein display style

$bg = '#FFAAAA';
$bg1 = "#FFFFFF";
$bg2 = "#FFDDFF";
$bg3 = "#DDFFDD";
$bg4 = "#DDDDFF";

$default_fitschedule = "0 0 0 1 1 0 0 0 0 0 0 1 0 0 1 1 1 1 0 1 1 0 0 1 1 1 1 1 1 1";
$fitschedule = $default_fitschedule;
$nvar = 6;
@fitvars = ( "b", "off", "gw", "ul_amp", "l_amp", "frac_n" );
$nround;
$tryitfitit = 0; # Defaults to "Fit It"

$datplotstyle = 0; # Default to points
$fitplotstyle = 1; # Default to lines
$xrangemin = 0; # Default to auto
$xrangemax = 0; # Default to auto
$yrangemin = 0; # Default to auto
$yrangemin = 0; # Default to auto
$global_frac_n = 0; # Stephen Variable

$nvlp_user = "user";
$nvlp_machine = "machine";

$nc_user = `whoami`;
chomp($nc_user);
$nc_machine = "";

$obey_display_flag = 0;
$only_save_displayed = 0;

$iso3_modelfile = "model_or_massparam";
$iso3_atomsfile = "$ENV{'MASSACRE_PATH'}/msc/files/exp_atom_defs.txt";

%name2col = ();

# Automatic Filtering Variables
$auto_type = -1;
$auto_max = 1;

# Redundancy Filtering Variables
$show_redunpeaks = 1;

# Specifies that the load window exists or not
# Used to toggle destroying it for command-line based file loading
$lw_exists = 0;

# Isodist3 does not calculate frac_lab by default
# This flag describes what to do in case f_lab is not found in the data read in
# Note that column headers are all transformed to lowercase
# 0: amp_l / (amp_l + amp_u)
# 1: (amp_u + amp_l) / (amp_u + amp_l + amp_f) (Unlabeled Labeled Fully labeled)
$frac_lab_type = 0;

&GetOptions("flab=i" => \$frac_lab_type, "file=s" => \$infile);

#########################################
# %pro2loc variables now in mikevars.pl #
#########################################

$mw = MainWindow->new();
$mw->configure(-title=>$title_display);
$mw->minsize( qw(1280 800) );

$ftop = $mw->Frame()->pack(-side=>'top', -fill=>'both');
	$ftopl = $ftop->Frame()->pack(-side=>'left', -fill=>'both');
	$ftopr = $ftop->Frame()->pack(-side=>'right', -fill=>'both', -expand=>'1');

$fbot = $mw->Frame()->pack(-side=>'top', -fill=>'both', -expand=>'1');
	$fbot1 = $fbot->Frame()->pack(-side=>'left', -fill=>'both', -expand=>'1');
		$fbot1a = $fbot1->Frame()->pack(-side=>'top', -fill=>'both', -expand=>'1');
		$fbot1b = $fbot1->Frame()->pack(-side=>'top', -fill=>'both');
	$fbot2 = $fbot->Frame()->pack(-side=>'left', -fill=>'both', -expand=>'1');
		$fbot2a = $fbot2->Frame()->pack(-side=>'top', -fill=>'both', -fill=>'both', -expand=>'1');
		$fbot2b = $fbot2->Frame()->pack(-side=>'top', -fill=>'both', -fill=>'both', -expand=>'1');
		$fbot2c = $fbot2->Frame()->pack(-side=>'top', -fill=>'both', -fill=>'both', -expand=>'1');
	$fbot3 = $fbot->Frame()->pack(-side=>'left', -fill=>'both', -expand=>'1');
		$fbot3a = $fbot3->Frame()->pack(-side=>'top', -fill=>'both', -fill=>'both', -expand=>'1');
		$fbot3b = $fbot3->Frame()->pack(-side=>'top', -fill=>'both', -fill=>'both', -expand=>'1');
		$fbot3c = $fbot3->Frame()->pack(-side=>'top', -fill=>'both', -fill=>'both', -expand=>'1');
	$fbot4 = $fbot->Frame()->pack(-side=>'right', -anchor=>'n');
		$fbotl = $fbot4->Frame()->pack(-side=>'top', -anchor=>'n');
		$fbotn = $fbot4->Frame()->pack(-side=>'right', -anchor=>'ne', -fill=>'both');
		$fbotm = $fbot4->Frame()->pack(-side=>'right', -anchor=>'ne');
		$fboto = $fbot4->Frame()->pack(-side=>'right', -anchor=>'ne');

$canvasl = $ftopl->Canvas(-background=>'white', -width=>'800', -height=>'600')->pack();

$plot_frame_width = 5000;

$box = [0, 0, $plot_frame_width, 370]; # may need to adjust $plot_frame_width to accomodate extra digests
$canvasr  = $ftopr->Scrolled(qw/
       Canvas -bg #FFFFFF
       -xscrollincrement 1
       -confine 1
       -scrollbars s
       -width 480
       -height 370/,
       -scrollregion => $box,
       );
$canvasr->pack(-fill=>'both');
$realcanvasr = $canvasr->Subwidget('scrolled');

#$f_plot_dataset = $ftopr->Frame()->pack(-side=>'top', -fill=>'both');
	#$plot_30S = $f_plot_dataset->Radiobutton(-text=>'30S', -variable=>\$plotframe, -value=>'0', -command=>\&draw_frame)->pack(-side=>'left');
	#$plot_50S = $f_plot_dataset->Radiobutton(-text=>'50S', -variable=>\$plotframe, -value=>'1', -command=>\&draw_frame)->pack(-side=>'left');
	#$plot_70S = $f_plot_dataset->Radiobutton(-text=>'70S', -variable=>\$plotframe, -value=>'2', -command=>\&draw_frame)->pack(-side=>'left');
	#$plot_30Schlamy = $f_plot_dataset->Radiobutton(-text=>'30Schlamy', -variable=>\$plotframe, -value=>'3', -command=>\&draw_frame)->pack(-side=>'left');
	#$plot_40Shuman = $f_plot_dataset->Radiobutton(-text=>'40Shuman', -variable=>\$plotframe, -value=>'4', -command=>\&draw_frame)->pack(-side=>'left');
	#$plot_60Shuman = $f_plot_dataset->Radiobutton(-text=>'60Shuman', -variable=>\$plotframe, -value=>'5', -command=>\&draw_frame)->pack(-side=>'left');
	#$plot_cof = $f_plot_dataset->Radiobutton(-text=>'cof', -variable=>\$plotframe, -value=>'6', -command=>\&draw_frame)->pack(-side=>'left');

$f_plot_width = $ftopr->Frame()->pack(-side=>'top', -fill=>'both');
	$but_narrow = $f_plot_width->Radiobutton(-text=>'||', -variable=>\$xint_type, -value=>'0', -command=>\&adjust_plot_width)->pack(-side=>'left');
	$but_normal = $f_plot_width->Radiobutton(-text=>'| |', -variable=>\$xint_type, -value=>'1', -command=>\&adjust_plot_width)->pack(-side=>'left');
	$but_wide = $f_plot_width->Radiobutton(-text=>'|  |', -variable=>\$xint_type, -value=>'2', -command=>\&adjust_plot_width)->pack(-side=>'left');

	$but_abx = $f_plot_width->Checkbutton(-text=>'abx', -variable=>\$display_abx, -command=>\&plot_frac)->pack(-side=>'left');

	$plotframebutton = $f_plot_width->Button(-text=>'Sample Type', -command=>\&choose_sample, -width=>'10')->pack(-side=>'right');
	$button_plotconfig = $f_plot_width->Button(-text=>'Define frac_lab', -command=>\&config_plot, -width=>'11')->pack(-side=>'right');
	$button_displaybehavior = $f_plot_width->Button(-text=>'Display Behavior', -command=>\&display_behavior, -width=>'13')->pack(-side=>'right');

$box2 = [0, 0, 5000, 125]; # Is this wide enough for all contour plots?
$canvasr2  = $ftopr->Scrolled(qw/
       Canvas -bg #DDDDDD
       -xscrollincrement 1
	   -yscrollincrement 1
       -confine 1
       -scrollbars se
       -width 480
       -height 60/,
       -scrollregion => $box2,
       );
$canvasr2->pack(-side=>'bottom', -fill=>'both');
$realcanvasr2 = $canvasr2->Subwidget('scrolled');

$ftopr_vars = $ftopr->Frame()->pack(-side=>'bottom', -fill=>'both');
	$vframe1 = $ftopr_vars->Frame()->pack(-side=>'bottom', -fill=>'both');
		$comment_label = $vframe1->Label(-text=>'Comment:')->pack(-side=>'left');
		$comment_entry = $vframe1->Entry(-textvariable=>\$comment_display, -width=>'35')->pack(-side=>'left');
		$offset_entry = $vframe1->Entry(-textvariable=>\$offset_display, -width=>'15', -relief=>'groove', -state=>'readonly')->pack(-side=>'right');
		$offset_label = $vframe1->Label(-text=>'Offset:')->pack(-side=>'right');
	$vframe2 = $ftopr_vars->Frame()->pack(-side=>'bottom', -fill=>'both');
		$seq_label = $vframe2->Label(-text=>'Sequence:')->pack(-side=>'left');
		$seq_entry = $vframe2->Entry(-textvariable=>\$sequence_display, -state=>'readonly', -relief=>'groove', -width=>'35')->pack(-side=>'left');
		$chisq_entry = $vframe2->Entry(-textvariable=>\$chisq_display, -state=>'readonly', -relief=>'groove', -width=>'15')->pack(-side=>'right');
		$chisq_label = $vframe2->Label(-text=>'ChiSq:')->pack(-side=>'right');

	$vframe4 = $ftopr_vars->Frame()->pack(-side=>'bottom', -fill=>'both');
		$ail14mass_label = $vframe4->Label(-text=>'ail14 Mass:')->pack(-side=>'left');
		$ail14mass_entry = $vframe4->Entry(-textvariable=>\$ail14mass_display, -state=>'readonly', -relief=>'groove', -width=>'12')->pack(-side=>'left');
		$ail15mass_label = $vframe4->Label(-text=>'ail15 Mass:')->pack(-side=>'left');
		$ail15mass_entry = $vframe4->Entry(-textvariable=>\$ail15mass_display, -state=>'readonly', -relief=>'groove', -width=>'12')->pack(-side=>'left');
		$avg_wr_entry = $vframe4->Entry(-textvariable=>\$avg_wr_display, -state=>'readonly', -relief=>'groove', -width=>'15')->pack(-side=>'right');
		$avg_wr_label = $vframe4->Label(-text=>"avgWR:")->pack(-side=>'right');


	$vframe3 = $ftopr_vars->Frame()->pack(-side=>'bottom', -fill=>'both');
		$n14mass_label = $vframe3->Label(-text=>'N14 Mass:')->pack(-side=>'left');
		$n14mass_entry = $vframe3->Entry(-textvariable=>\$n14mass_display, -state=>'readonly', -relief=>'groove', -width=>'12')->pack(-side=>'left');
		$n15mass_label = $vframe3->Label(-text=>'N15 Mass:')->pack(-side=>'left');
		$n15mass_entry = $vframe3->Entry(-textvariable=>\$n15mass_display, -state=>'readonly', -relief=>'groove', -width=>'12')->pack(-side=>'left');
		$endres_entry = $vframe3->Entry(-textvariable=>\$endres_display, -state=>'readonly', -relief=>'groove', -width=>'4')->pack(-side=>'right');
		$endres_label = $vframe3->Label(-text=>'End:')->pack(-side=>'right');
		$startres_entry = $vframe3->Entry(-textvariable=>\$startres_display, -state=>'readonly', -relief=>'groove', -width=>'4')->pack(-side=>'right');
		$startres_label = $vframe3->Label(-text=>'Start:')->pack(-side=>'right');


	$vframe5 = $ftopr_vars->Frame()->pack(-side=>'bottom', -fill=>'both');

#$buttonl = $fbotn->Button(-text=>'Go(L)', -command=>\&go1)->pack();
#$buttonr = $fbotn->Button(-text=>'Go(R)', -command=>\&go2)->pack();

$iso3modelentry = $fbotl->Entry(-textvariable=>\$iso3_modelfile, -width=>'14')->pack(-side=>'right');
$iso3regenbutton = $fbotl->Button(-text=>'Reg(3)', -command=>\&iso3_regen, -width=>'4')->pack(-side=>'right');

$loadbutton = $fbotn->Button(-text=>'Load', -width=>'4', -command=>\&launch_load)->pack(-side=>'top');
$savebutton = $fbotn->Button(-text=>'Save', -width=>'4', -command=>\&launch_save)->pack(-side=>'top');
$spacer1 = $fbotn->Label(-text=>'_________', -width=>'6')->pack(-side=>'top');
$redunbutton = $fbotn->Button(-text=>'Redun', -width=>'4', -command=>\&launch_redun)->pack(-side=>'top');
$quitbutton = $fbotn->Button(-text=>'Quit', -width=>'4', -command=>\&quit)->pack(-side=>'bottom');

$refitbutton = $fbotm->Button(-text=>'Refit', -width=>'4', -command=>\&launch_refit)->pack();
$envelopebutton = $fbotm->Button(-text=>'Nvlp', -width=>'4', -command=>\&launch_envelope)->pack();
$spacer2 = $fbotm->Label(-text=>'_________', -width=>'6')->pack();
$autobutton = $fbotm->Button(-text=>'Auto', -width=>'4', -command=>\&launch_auto)->pack();
$rtcorrbutton = $fbotm->Button(-text=>'RTPic', -width=>'4', -command=>\&display_rtcorr)->pack(-side=>'top');

$shuntbutton = $fboto->Button(-text=>'<---', -width=>'4', -command=>\&shunt)->pack();
$ncbutton = $fboto->Button(-text=>'NewCon', -width=>'4', -command=>\&launch_NewContour)->pack();
$spacer3 = $fboto->Label(-text=>'_________', -width=>'6')->pack();
$deletebutton = $fboto->Button(-text=>'Delete', -width=>'4', -command=>\&delete_selection)->pack();
$undobutton = $fboto->Button(-text=>'Undo', -width=>'4', -command=>\&undo_deletion)->pack();
$purgebutton = $fboto->Button(-text=>'Purge', -width=>'4', -command=>\&purge_selection)->pack();



# Main, Leftmost Listbox
# Contains the main view, traditionally the list of all fits
$list1 = $fbot1a->Listbox(-listvariable=>\@view1, -state=>'normal', -selectmode=>'browse', -width=>'60', -height=>'10', -background=>$bg1)->pack(-side=>'left', -fill=>'both', -anchor=>'n', -expand=>'1');
$listscrolly1 = $fbot1a->Scrollbar(-orient=>'vertical')->pack(-side=>'right', -fill=>'both');
$list1->configure(-yscrollcommand=>['set'=>$listscrolly1]);
$listscrolly1->configure(-command=>['yview'=>$list1]);
$list1->bind('<<ListboxSelect>>'=>\&process_selection1);

# Sorting toggles for the main listbox
# By default, $sortstyle = 0, meaning default sort and no active button
$sort1 = $fbot1b->Frame()->pack(-side=>'bottom', -anchor=>'nw', -fill=>'both');
	$sort1_label = $sort1->Label(-text=>'Sort By: ')->pack(-side=>'left');
	$sort1_rt = $sort1->Radiobutton(-text=>'RT', -command=>\&sort_list1, -variable=>\$sortstyle, -value=>'1')->pack(-side=>'left');
	$sort1_mz = $sort1->Radiobutton(-text=>'M/Z', -command=>\&sort_list1, -variable=>\$sortstyle, -value=>'2')->pack(-side=>'left');
	$sort1_protein = $sort1->Radiobutton(-text=>'Protein', -command=>\&sort_list1, -variable=>\$sortstyle, -value=>'3')->pack(-side=>'left');
	$sort1_prochi = $sort1->Radiobutton(-text=>'Pro->ChiSq', -command=>\&sort_list1, -variable=>\$sortstyle, -value=>'4')->pack(-side=>'left');

	$sort1_displayflag = $sort1->Checkbutton(-text=>'Show this Point', -offvalue=>'0', -onvalue=>'1', -variable=>\$current_display_flag, -command=>\&set_display_flag)->pack(-side=>'right');

# Secondary Listbox
$list2label = $fbot2a->Label(-text=>'Ion:', -font=>'*-helvetica-medium-r-normal-10-*', -width=>'3', -foreground=>'#FF44FF')->pack(-side=>'left', -anchor=>'nw');
$list2 = $fbot2a->Listbox(-listvariable=>\@view2, -state=>'normal', -selectmode=>'browse', -width=>'50', -height=>'4', -font=>'*-helvetica-medium-r-normal-10-*', -background=>$bg2)->pack(-side=>'left', -fill=>'both', -expand=>'1');
$listscrolly2 = $fbot2a->Scrollbar(-orient=>'vertical')->pack(-side=>'left', -fill=>'both');
$list2->configure(-yscrollcommand=>['set'=>$listscrolly2]);
$listscrolly2->configure(-command=>['yview'=>$list2]);
$list2->bind('<<ListboxSelect>>'=>\&process_selection2);

# Secondary Listbox
$list3label = $fbot2b->Label(-text=>'Pep:', -font=>'*-helvetica-medium-r-normal-10-*', -width=>'3', -foreground=>'#009900')->pack(-side=>'left', -anchor=>'nw');
$list3 = $fbot2b->Listbox(-listvariable=>\@view3, -state=>'normal', -selectmode=>'browse', -width=>'50', -height=>'4', -font=>'*-helvetica-medium-r-normal-10-*', -background=>$bg3)->pack(-side=>'left', -fill=>'both', -expand=>'1');
$listscrolly3 = $fbot2b->Scrollbar(-orient=>'vertical')->pack(-side=>'left', -fill=>'both');
$list3->configure(-yscrollcommand=>['set'=>$listscrolly3]);
$listscrolly3->configure(-command=>['yview'=>$list3]);
$list3->bind('<<ListboxSelect>>'=>\&process_selection3);

# Secondary Listbox
$list4label = $fbot2c->Label(-text=>'Pks:', -font=>'*-helvetica-medium-r-normal-10-*', -width=>'3', -foreground=>'#0000FF')->pack(-side=>'left', -anchor=>'nw');
$list4 = $fbot2c->Listbox(-listvariable=>\@view4, -state=>'normal', -selectmode=>'browse', -width=>'50', -height=>'4', -font=>'*-helvetica-medium-r-normal-10-*', -background=>$bg4)->pack(-side=>'left', -fill=>'both', -expand=>'1');
$listscrolly4 = $fbot2c->Scrollbar(-orient=>'vertical')->pack(-side=>'left', -fill=>'both');
$list4->configure(-yscrollcommand=>['set'=>$listscrolly4]);
$listscrolly4->configure(-command=>['yview'=>$list4]);
$list4->bind('<<ListboxSelect>>'=>\&process_selection4);

# Tertiary Listbox
$list2b = $fbot3a->Listbox(-listvariable=>\@view2b, -state=>'normal', -selectmode=>'browse', -width=>'50', -height=>'4', -font=>'*-helvetica-medium-r-normal-10-*', -background=>$bg2)->pack(-side=>'left', -fill=>'both', -expand=>'1');
$listscrolly2b = $fbot3a->Scrollbar(-orient=>'vertical')->pack(-side=>'left', -fill=>'both');
$list2b->configure(-yscrollcommand=>['set'=>$listscrolly2b]);
$listscrolly2b->configure(-command=>['yview'=>$list2b]);
$list2b->bind('<<ListboxSelect>>'=>\&process_selection2b);

# Tertiary Listbox
$list3b = $fbot3b->Listbox(-listvariable=>\@view3b, -state=>'normal', -selectmode=>'browse', -width=>'50', -height=>'4', -font=>'*-helvetica-medium-r-normal-10-*', -background=>$bg3)->pack(-side=>'left', -fill=>'both', -expand=>'1');
$listscrolly3b = $fbot3b->Scrollbar(-orient=>'vertical')->pack(-side=>'left', -fill=>'both');
$list3b->configure(-yscrollcommand=>['set'=>$listscrolly3b]);
$listscrolly3b->configure(-command=>['yview'=>$list3b]);
$list3b->bind('<<ListboxSelect>>'=>\&process_selection3b);

# Tertiary Listbox
$list4b = $fbot3c->Listbox(-listvariable=>\@view4b, -state=>'normal', -selectmode=>'browse', -width=>'50', -height=>'4', -font=>'*-helvetica-medium-r-normal-10-*', -background=>$bg4)->pack(-side=>'left', -fill=>'both', -expand=>'1');
$listscrolly4b = $fbot3c->Scrollbar(-orient=>'vertical')->pack(-side=>'left', -fill=>'both');
$list4b->configure(-yscrollcommand=>['set'=>$listscrolly4b]);
$listscrolly4b->configure(-command=>['yview'=>$list4b]);
$list4b->bind('<<ListboxSelect>>'=>\&process_selection4b);

######################
# Keyboard Shortcuts #
######################
$mw->bind('<Up>' => \&minus1);
$mw->bind('<Down>' => \&plus1);
$mw->bind('<Left>' => \&shunt);
$mw->bind('<Control-d>' => \&delete_selection);
$mw->bind('<Control-z>' => \&undo_deletion);
$mw->bind('<Control-n>' => \&launch_NewContour);

##############################
# Stuff to Execute on Launch #
##############################
print "Launching massive...\n";

# Set up gnuplot
my $gnupid1 = open3( \*GIN1, \*GOUT1, \*GERR1, "/usr/bin/gnuplot" ) || die "Can't open gnuplot\n";
$mw->fileevent( \*GOUT1, readable => \&read_gout1 );
$mw->fileevent( \*GERR1, readable => \&read_gerr1 );
print GIN1 "set term X11 noraise nopersist\n";
#$gnuctr1 = 0;

###########################
# Setup the frac_lab plot #
###########################
@yticloc = ( 310, 280, 250, 220, 190, 160, 130, 100, 70, 40, 10 );
$xoffset = 20; # distance in from left that the graph starts;
$yoffset = 10;
$xint = 15; # x-interval between proteins
$xint_type = 1;
draw_frame();

##########################################
# If an input file is specified, load it #
##########################################
if($infile ne "")
{
	load();
}

#########################################
# interacting with the plot in $canvasr #
#########################################
$canvasr->CanvasBind('<1>' => sub
{
	my $x = $canvasr->canvasx($Tk::event->x);
	my $y = $canvasr->canvasy($Tk::event->y);

	print "($x, $y)\n";

	if($#ovals >= 0)
	{
		for(my $ctr = 0; $ctr <= $#ovals; $ctr++)
		{
			$dist = ($x - $ovals[$ctr][0])**2 + ($y - $ovals[$ctr][1])**2;
			if($ctr == 0)
			{
				$min = $ctr;
				$mindist = $dist;
			}
			else
			{
				if($dist < $mindist)
				{
					$min = $ctr;
					$mindist = $dist;
				}
			}
		}
	}

	$list1->selectionClear(0, $#view1);
	$list1->selectionSet($min);
	$list1->see($min);
	process_selection1();
});

#################
# Enter the GUI #
#################
MainLoop();

sub minus1
{
	@selarray1 = $list1->curselection;
	$selindex1 = $selarray1[0];

	$selindex1 -= 1;
	$selindex1 = $lastsel1 - 1;

	if($selindex1 < 0)
	{
		$selindex1 = 0;
	}

	$list1->selectionClear(0, $#view1);
	$list1->selectionSet($selindex1);
	#$list1->yview($selindex1);
	$list1->see($selindex1);
	process_selection1();
}

sub plus1
{
	@selarray1 = $list1->curselection;
	$selindex1 = $selarray1[0];

	$selindex1 += 1;
	$selindex1 = $lastsel1 + 1;

	if($selindex1 > $#view1)
	{
		$selindex1 = $#view1;
	}

	$list1->selectionClear(0, $#view1);
	$list1->selectionSet($selindex1);
	#$list1->yview($selindex1);
	$list1->see($selindex1);
	process_selection1();
}

sub draw_frame
{
	my $ctr;

	if($plotframe == 0)
	{
		%pro2loc = %pro2loc30S;
	}
	elsif($plotframe == 1)
	{
		%pro2loc = %pro2loc50S;
	}
	elsif($plotframe == 2)
	{
		%pro2loc = %pro2loc70S;
	}
	elsif($plotframe == 3)
	{
		%pro2loc = %pro2loc30Schlamy;
	}
	elsif($plotframe == 4)
	{
		%pro2loc = %pro2loc40Shuman;
	}
	elsif($plotframe == 5)
	{
		%pro2loc = %pro2loc60Shuman;
	}
	elsif($plotframe == 6)
	{
		%pro2loc = %pro2loccof;
	}
	elsif($plotframe == 7)
	{
		%pro2loc = %pro2loc40Syeast_old;
	}
	elsif($plotframe == 8)
	{
		%pro2loc = %pro2loc40Syeast;
	}
	elsif($plotframe == 9)
	{
		%pro2loc = %pro2loc60Syeast;
	}
	elsif($plotframe == 10)
	{
		%pro2loc = %pro2loc80Syeast;
	}
	elsif($plotframe == 11)
	{
		%pro2loc = %pro2locRNA;
	}
	elsif($plotframe == 12) # Proteomic, generate %pro2loc dynamically
	{
		%pro2loc = ();
		my %protein_hash = ();
		my @protein_array = ();

		my $ctr_d1;
		for($ctr_d1 = 0; $ctr_d1 <= $#data1; $ctr_d1++)
		{
			#my $temp_protein = $data1[$ctr_d1]{"protein"};
			#print "Protein? $temp_protein\n";
			$protein_hash{ $data1[$ctr_d1]{"protein"} } += 1;
		}

		foreach my $protein (keys %protein_hash)
		{
			push(@protein_array, $protein);
		}

		# Sort old way based on pro2loc
		if(!defined($name2col{"pks_acc_n14"}))
		{
			@protein_array = sort proteomic (@protein_array);
		}
		# Sort new way based on knowing what is an accession number
		else
		{
			@protein_array = sort proteomic2 (@protein_array);
		}

		for($ctr_d1 = 0; $ctr_d1 <= $#protein_array; $ctr_d1++)
		{
			#print "Protein! $protein_array[$ctr_d1]\n";
			$pro2loc{$protein_array[$ctr_d1]} = $ctr_d1;
		}
	}
	elsif($plotframe == 13)
	{
		%pro2loc = %pro2loc80Shuman;
	}

	%loc2pro = ();
	foreach my $protein (keys %pro2loc)
	{
		$loc2pro{$pro2loc{$protein}} = $protein;
	}

	$canvasr->delete('frame');

	for($ctr = 0; $ctr <= $#yticloc; $ctr++)
	{
		$canvasr->createLine($xoffset, $yticloc[$ctr], $plot_frame_width-$xoffset, $yticloc[$ctr], -fill=>'#CCCCCC', -tags=>['frame']);
		$textval = sprintf("%.1f", $ctr*0.1);
		$canvasr->createText(1, $yticloc[$ctr], -text=>$textval, -font=>'*-courier-medium-r-normal-10-*', -anchor=>'w', -tags=>['frame']);
	}

	$canvasr->createRectangle($xoffset, $yoffset, $plot_frame_width-$xoffset, $yticloc[0], -tags=>['frame']);

	foreach $protein (keys %pro2loc)
	{
		$canvasr->createLine($xoffset+$xint*$pro2loc{$protein}, $yticloc[0], $xoffset+$xint*$pro2loc{$protein}, $yticloc[0] + 5, -tags=>['frame']);
		$canvasr->createText($xoffset+$xint*$pro2loc{$protein}, $yticloc[0] + 5, -text=>$protein, -font=>'*-courier-medium-r-normal-10-*', -anchor=>'n', -width=>'1', -tags=>['frame']);
	}

	$canvasr->delete('pts');
	plot_frac();
	process_selection1();
}

sub launch_auto
{
	$aw = MainWindow->new();
	$aw->configure(-title=>'Automatic Filtering');
	$awframe1 = $aw->Frame(-border=>'5')->pack(-side=>'top');
	$awframe1b = $aw->Frame(-border=>'5')->pack(-side=>'top');
	$awframe1c = $aw->Frame(-border=>'5')->pack(-side=>'top');
	$awframe1d = $aw->Frame(-border=>'5')->pack(-side=>'top');
	$awframe1e = $aw->Frame(-border=>'5')->pack(-side=>'top');
	$awframe1f = $aw->Frame(-border=>'5')->pack(-side=>'top');
	$awframe1g = $aw->Frame(-border=>'5')->pack(-side=>'top');
	$awframe2 = $aw->Frame()->pack(-side=>'top');
	$awframe3 = $aw->Frame()->pack(-side=>'top');
	$awframe4 = $aw->Frame()->pack(-side=>'top', -fill=>'both');
	$awframe6 = $aw->Frame(-border=>'10')->pack(-side=>'top', -fill=>'both');
	$awframe6z = $awframe6->Frame(-border=>'10')->pack(-side=>'left', -fill=>'both');
	$awframe6a = $awframe6->Frame(-border=>'10')->pack(-side=>'left', -fill=>'both');
	$awframe6b = $awframe6->Frame(-border=>'10')->pack(-side=>'left', -fill=>'both');
	$awframe6c = $awframe6->Frame(-border=>'10')->pack(-side=>'left', -fill=>'both');
	$awframe6d = $awframe6->Frame(-border=>'10')->pack(-side=>'left', -fill=>'both');
	$awframe6e = $awframe6->Frame(-border=>'10')->pack(-side=>'left', -fill=>'both');
	$awframe7 = $aw->Frame(-border=>'5')->pack(-side=>'top', -fill=>'both');
	$awframe8 = $aw->Frame(-border=>'5')->pack(-side=>'top', -fill=>'both');
	$awframe9 = $aw->Frame(-border=>'5')->pack(-side=>'top', -fill=>'both');
	$awframe5 = $aw->Frame(-border=>'10')->pack(-side=>'top', -fill=>'both');

	$auto_count = 0;

	@colarray = ();

	foreach $column_name (keys %name2col)
	{
		#if(substr($column_name, 0, 5) eq "chisq" || substr($column_name, 0, 6) eq "avg_wr")
		if(substr($column_name, 0, 5) eq "chisq" || $column_name eq "avg_wr13c" || $column_name eq "avg_wr4d" || $column_name eq "avg_wr14d" || $column_name eq "avg_wr12c")
		{
			push(@colarray, $column_name);
		}
	}

	@colarray = sort @colarray;

	for($auto_ctr = 0; $auto_ctr <= $#colarray; $auto_ctr++)
	{
		$column_name = $colarray[$auto_ctr];

		$auto_type_array[$auto_count] = $column_name;

		if($auto_count < 6)
		{
			$aw_radio[$auto_count] = $awframe1->Radiobutton(-textvariable=>\$auto_type_array[$auto_count], -variable=>\$auto_type, -value=>$auto_count, -command=>\&define_auto_type)->pack(-side=>'left');
		}
		elsif($auto_count >= 6 && $auto_count < 12)
		{
			$aw_radio[$auto_count] = $awframe1b->Radiobutton(-textvariable=>\$auto_type_array[$auto_count], -variable=>\$auto_type, -value=>$auto_count, -command=>\&define_auto_type)->pack(-side=>'left');
		}
		elsif($auto_count >= 12 && $auto_count < 18)
		{
			$aw_radio[$auto_count] = $awframe1c->Radiobutton(-textvariable=>\$auto_type_array[$auto_count], -variable=>\$auto_type, -value=>$auto_count, -command=>\&define_auto_type)->pack(-side=>'left');
		}
		elsif($auto_count >= 18 && $auto_count < 24)
		{
			$aw_radio[$auto_count] = $awframe1d->Radiobutton(-textvariable=>\$auto_type_array[$auto_count], -variable=>\$auto_type, -value=>$auto_count, -command=>\&define_auto_type)->pack(-side=>'left');
		}
		elsif($auto_count >= 24 && $auto_count < 30)
		{
			$aw_radio[$auto_count] = $awframe1e->Radiobutton(-textvariable=>\$auto_type_array[$auto_count], -variable=>\$auto_type, -value=>$auto_count, -command=>\&define_auto_type)->pack(-side=>'left');
		}
		else
		{
			$aw_radio[$auto_count] = $awframe1f->Radiobutton(-textvariable=>\$auto_type_array[$auto_count], -variable=>\$auto_type, -value=>$auto_count, -command=>\&define_auto_type)->pack(-side=>'left');
		}

		$auto_count++;
	}
	$aw_viewbutton = $awframe1g->Button(-text=>'View', -command=>sub{\&prep_view_variable()})->pack();

	if(!$auto_type)
	{
		$auto_type = -1;
	}

#	if(defined($name2col{"chisq"}) && defined($name2col{"avg_wr"}))
#	{
#		$aw_radio1 = $awframe1->Radiobutton(-text=>'avgWR', -variable=>\$auto_type, -value=>'0', -command=>\&define_auto_type)->pack(-side=>'left');
#		$aw_radio2 = $awframe1->Radiobutton(-text=>'ChiSq', -variable=>\$auto_type, -value=>'1', -command=>\&define_auto_type)->pack(-side=>'left');
#	}
#	elsif(defined($name2col{"chisq"}) && !defined($name2col{"avg_wr"}))
#	{
#		$aw_radio2 = $awframe1->Radiobutton(-text=>'ChiSq', -variable=>\$auto_type, -value=>'1', -command=>\&define_auto_type)->pack(-side=>'left');
#	}
#	elsif(!defined($name2col{"chisq"}) && defined($name2col{"avg_wr"}))
#	{
#		$aw_radio1 = $awframe1->Radiobutton(-text=>'avgWR', -variable=>\$auto_type, -value=>'0', -command=>\&define_auto_type)->pack(-side=>'left');
#	}

#	if(defined($name2col{"chisq"}) && defined($name2col{"avg_wr"}))
#	{
#		$aw_chisqview = $awframe1->Button(-text=>'View ChiSq', -command=>sub{\&view_variable("chisq")})->pack(-side=>'left');
#		$aw_avg_wrview = $awframe1->Button(-text=>'View AvgWR', -command=>sub{\&view_variable("avg_wr")})->pack(-side=>'left');
#	}
#	elsif(defined($name2col{"chisq"}) && !defined($name2col{"avg_wr"}))
#	{
#		$aw_chisqview = $awframe1->Button(-text=>'View ChiSq', -command=>sub{\&view_variable("chisq")})->pack(-side=>'left');
#	}
#	elsif(!defined($name2col{"chisq"}) && defined($name2col{"avg_wr"}))
#	{
#		$aw_avg_wrview = $awframe1->Button(-text=>'View AvgWR', -command=>sub{\&view_variable("avg_wr")})->pack(-side=>'left');
#	}



	$awscale = $awframe2->Scale(-orient=>'horizontal', -from=>'0', -to=>'100', -variable=>\$auto_cutoff_percent, -command=>\&auto_convert, -showvalue=>'1', -length=>'480', -resolution=>'0.1')->pack(-side=>'left');

	$awleftscale = $awframe3->Label(-text=>'0')->pack(-side=>'left');
	$awrightscale = $awframe3->Label(-textvariable=>\$auto_max)->pack(-side=>'right');
	$awval3 = $awframe3->Label(-text=>')')->pack(-side=>'right');
	$awval2 = $awframe3->Label(-textvariable=>\$auto_cutoff)->pack(-side=>'right');
	$awval1 = $awframe3->Label(-text=>'(')->pack(-side=>'right');

	$awnewentry = $awframe4->Entry(-textvariable=>\$auto_new, -width=>'12')->pack(-side=>'right');
	$awnewbutton = $awframe4->Button(-text=>'Best X %', -command=>\&auto_apply)->pack(-side=>'right');

	$awthresholdentry = $awframe4->Entry(-textvariable=>\$auto_threshold, -width=>'12')->pack(-side=>'left');
	$awthresholdbutton = $awframe4->Button(-text=>'Cutoff Value', -command=>\&auto_apply_threshold)->pack(-side=>'left');

	$awdeletebutton = $awframe5->Button(-text=>'Delete Filtered Points', -command=>\&auto_delete)->pack(-side=>'left');
	$awquitbutton = $awframe5->Button(-text=>'Close', -command=>\&quit_auto)->pack(-side=>'right');

	$auto_selection = -1;

	$awproteinall = $awframe6z->Radiobutton(-text=>'All', -variable=>\$auto_selection, -value=>'-1', -command=>\&auto_protein)->pack();

	# Nobody uses the individual proteins feature 1/6/11 - removing it
	# It also messes up the screen with proteomic data
	@auto_proteins = ();
	#foreach $loc (keys %loc2pro)
	#{
	#	$auto_proteins[$loc] = $loc2pro{$loc};
	#}

	for($auto_ctr = 0; $auto_ctr <= $#auto_proteins; $auto_ctr++)
	{
		$auto_radio_values[$auto_ctr] = $auto_ctr;

		if($auto_ctr % 5 == 0)
		{
			$awprotein[$auto_ctr] = $awframe6a->Radiobutton(-textvariable=>\$loc2pro{$auto_ctr}, -variable=>\$auto_selection, -value=>$auto_radio_values[$auto_ctr], -command=>\&auto_protein)->pack(-side=>'top');
		}
		elsif($auto_ctr % 5 == 1)
		{
			$awprotein[$auto_ctr] = $awframe6b->Radiobutton(-textvariable=>\$loc2pro{$auto_ctr}, -variable=>\$auto_selection, -value=>$auto_radio_values[$auto_ctr], -command=>\&auto_protein)->pack(-side=>'top');
		}
		elsif($auto_ctr % 5 == 2)
		{
			$awprotein[$auto_ctr] = $awframe6c->Radiobutton(-textvariable=>\$loc2pro{$auto_ctr}, -variable=>\$auto_selection, -value=>$auto_radio_values[$auto_ctr], -command=>\&auto_protein)->pack(-side=>'top');
		}
		elsif($auto_ctr % 5 == 3)
		{
			$awprotein[$auto_ctr] = $awframe6d->Radiobutton(-textvariable=>\$loc2pro{$auto_ctr}, -variable=>\$auto_selection, -value=>$auto_radio_values[$auto_ctr], -command=>\&auto_protein)->pack(-side=>'top');
		}
		elsif($auto_ctr % 5 == 4)
		{
			$awprotein[$auto_ctr] = $awframe6e->Radiobutton(-textvariable=>\$loc2pro{$auto_ctr}, -variable=>\$auto_selection, -value=>$auto_radio_values[$auto_ctr], -command=>\&auto_protein)->pack(-side=>'top');
		}
	}

	$aw_b_view = $awframe7->Button(-text=>'View', -width=>'4', -command=>sub{\&view_variable("b")} )->pack(-side=>'right');
	$aw_b_max = $awframe7->Entry(-textvariable=>\$auto_b_max, -width=>'10')->pack(-side=>'right');
	$aw_b_maxlabel = $awframe7->Label(-text=>'Max:')->pack(-side=>'right');
	$aw_b_min = $awframe7->Entry(-textvariable=>\$auto_b_min, -width=>'10')->pack(-side=>'right');
	$aw_b_minlabel = $awframe7->Label(-text=>'Min:')->pack(-side=>'right');
	$aw_b_label = $awframe7->Label(-text=>'Baseline:')->pack(-side=>'right');

	$aw_gw_view = $awframe8->Button(-text=>'View', -width=>'4', -command=>sub{\&view_variable("gw")} )->pack(-side=>'right');
	$aw_gw_max = $awframe8->Entry(-textvariable=>\$auto_gw_max, -width=>'10')->pack(-side=>'right');
	$aw_gw_maxlabel = $awframe8->Label(-text=>'Max:')->pack(-side=>'right');
	$aw_gw_min = $awframe8->Entry(-textvariable=>\$auto_gw_min, -width=>'10')->pack(-side=>'right');
	$aw_gw_minlabel = $awframe8->Label(-text=>'Min:')->pack(-side=>'right');
	$aw_gw_label = $awframe8->Label(-text=>'Gaussian Width:')->pack(-side=>'right');

	$aw_off_view = $awframe9->Button(-text=>'View', -width=>'4', -command=>sub{\&view_variable("off")} )->pack(-side=>'right');
	$aw_off_max = $awframe9->Entry(-textvariable=>\$auto_off_max, -width=>'10')->pack(-side=>'right');
	$aw_off_maxlabel = $awframe9->Label(-text=>'Max:')->pack(-side=>'right');
	$aw_off_min = $awframe9->Entry(-textvariable=>\$auto_off_min, -width=>'10')->pack(-side=>'right');
	$aw_off_minlabel = $awframe9->Label(-text=>'Min:')->pack(-side=>'right');
	$aw_off_label = $awframe9->Label(-text=>'Offset:')->pack(-side=>'right');

}

sub prep_view_variable
{
	view_variable($auto_type_array[$auto_type]);
}

sub view_variable
{
	my @vars = @_;
	my $variable = shift(@vars);
	#my $variable2 = shift(@vars);

	print "Variable: |$variable|\n";
	#print "Variable2: |$variable2|\n";

	open TEMPDATA, ">$pid.tempdata" or die "Can't open $pid.tempdata\n";
	my @tempdata = ();
	for(my $local_ctr = 0; $local_ctr <= $#data1; $local_ctr++)
	{
		push(@tempdata, $data1[$local_ctr]{$variable});
	}
	@tempdata = sort numerically @tempdata;
	for($local_ctr = 0; $local_ctr <= $#tempdata; $local_ctr++)
	{
		print TEMPDATA "$tempdata[$local_ctr]\n";
	}
	close TEMPDATA;

	print GIN1 "set autoscale\n";
	print GIN1 "plot \"$pid.tempdata\"\n";

	`rm $pid.tempdata`;
}

sub auto_protein
{
	print "$auto_selection\n";
	$auto_cutoff_percent = 100;
	auto_convert();
}

sub auto_delete
{
	for($auto_ctr = 0; $auto_ctr <= $#data1; $auto_ctr++)
	{
		if($data1[$auto_ctr]{"display_flag"} == 0)
		{
			push(@deleted, splice(@data1, $auto_ctr, 1));
			$auto_ctr--;
		}
	}

	regen_view1();
	plot_frac();
	$lastsel1 = -1;
}

sub do_auto_filter
{
	#print "do_auto_filter\n";

	if($auto_type < 0)
	{
		return;
	}

	# Tie into the display flag system
	$obey_display_flag = 1;

	for($auto_ctr = 0; $auto_ctr <= $#data1; $auto_ctr++)
	{
		$auto_val = $data1[$auto_ctr]{$auto_type_array[$auto_type]};

		#if($auto_type == 0)
		#{
		#	$auto_val = $data1[$auto_ctr]{"avg_wr"};
		#}
		#elsif($auto_type == 1)
		#{
		#	$auto_val = $data1[$auto_ctr]{"chisq"};
		#}

		if($auto_val > $auto_cutoff && ($auto_selection == -1 || $pro2loc{$data1[$auto_ctr]{"protein"}} == $auto_selection))
		{
			$data1[$auto_ctr]{"display_flag"} = 0;
		}
		elsif($auto_selection == -1 || $pro2loc{$data1[$auto_ctr]{"protein"}} == $auto_selection)
		{
			$data1[$auto_ctr]{"display_flag"} = 1;
		}
	}

	plot_frac();
}

# Filter based on an actual cutoff value
sub auto_apply_threshold
{
	# Just go through @auto_array until we reach a value > $auto_threshold
	for($auto_ctr = 0; $auto_ctr <= $#auto_array; $auto_ctr++)
	{
		if($auto_array[$auto_ctr] > $auto_threshold)
		{
			$auto_cutoff_percent = 100*($auto_ctr/$#auto_array);
			last;
		}
	}

	$auto_cutoff = $auto_threshold;
	do_auto_filter();
}

# Filter based on a percentage of your peaks
sub auto_apply
{
	if($auto_new > 100)
	{
		$auto_new = 100;
	}
	if($auto_new < 0)
	{
		$auto_new = 0;
	}

	$auto_cutoff_percent = $auto_new;
	auto_convert();
}

sub auto_convert
{
	$auto_cutoff_index = ($auto_cutoff_percent/100) * $#auto_array;
	if($auto_cutoff_index > $#auto_array)
	{
		$auto_cutoff_index = $#auto_array;
	}
	elsif($auto_cutoff_index < 0)
	{
		$auto_cutoff_index = 0;
	}
	$auto_cutoff = $auto_array[$auto_cutoff_index];
	do_auto_filter();
}

sub define_auto_type
{
	print "auto_type: $auto_type -> $auto_type_array[$auto_type]\n";
	$auto_max = 0;
	@auto_array = ();

	for($auto_ctr = 0; $auto_ctr <= $#data1; $auto_ctr++)
	{
		$auto_val = $data1[$auto_ctr]{$auto_type_array[$auto_type]};

		$auto_val = trim($auto_val);

		# Check for "nan" and "inf", then set to big number
		# Go back and change the value stored in @data1 as well
		if(lc($auto_val) eq "inf" || lc($auto_val) eq "nan")
		{
			$auto_val = 1000000000;
			$data1[$auto_ctr]{$auto_type_array[$auto_type]} = $auto_val;
		}

		if($auto_val > $auto_max)
		{
			$auto_max = $auto_val;
		}

		push(@auto_array, $auto_val);
	}

	@auto_array = sort numerically @auto_array;

	$auto_cutoff_percent = 100;
	auto_convert();
}

sub quit_auto
{
	$aw->destroy();
}

sub launch_load
{
	# Try to guess at $infile based on the directory
	@infilearray = split /\//, `pwd`;
	chomp(@infilearray);
	$infile = join '', $infilearray[$#infilearray], "_iso.csv";

	$lw = MainWindow->new();
	$lw->configure(-title=>'Load New File');
	$lw->minsize( qw(600 50) );
	$lwframe = $lw->Frame(-border=>'10', -background=>'white')->pack(-fill=>both);
	$lwentry = $lwframe->Entry(-textvariable=>\$infile)->pack(-side=>'left', -expand=>'1', -fill=>x);
	$lwbutton = $lwframe->Button(-text=>'Load File', -command=>\&load)->pack(-side=>'left');

	$lw_exists = 1;
}

sub load
{
	if($lw_exists)
	{
		$lw->destroy();
		$lw_exists = 0;
	}

	clear1();

	if(!(-e $infile))
	{
		print "Error: $infile does not exist\n";
		return;
	}

	open IN, "$infile" or die "Can't open $infile\n";
	@inarray = <IN>;
	close IN;

	chomp(@inarray);
	$header = shift(@inarray);
	@headerarray = split /\,/, $header;

	for($ctr = 0; $ctr <= $#headerarray; $ctr++)
	{
		$headerarray[$ctr] = lc($headerarray[$ctr]);
		$name2col{ $headerarray[$ctr] } = $ctr;
		$col2name{ $ctr } = $headerarray[$ctr];

		# Collect the array indices of all the calculated amplitudes
		if(substr($headerarray[$ctr], 0, 3) eq "amp")
		{
			push(@amp_array, $ctr);
		}
	}

	# Check for "frac_lab" column - this is missing in isodist3
	# If missing, we will tack it onto the end
	if(defined($name2col{"frac_lab"}))
	{
		$was_frac_lab = 1;
	}
	else
	{
		$was_frac_lab = 0;
		$name2col{ "frac_lab" } = $ctr;
		$col2name{ $ctr } = "frac_lab";
	}

	if(defined($name2col{"display_flag"}))
	{
		$was_display_flag = 1;
	}
	else
	{
		$was_display_flag = 0;
	}

	# Process the file line by line
	$data_ctr = 0; # Use this to assign values to @data.  It is incremented separately as it needs to be manipulated when lines are skipped
	INPUT_LINE: for($ctr = 0; $ctr <= $#inarray; $ctr++)
	{
		@array = split /\,/, $inarray[$ctr];

		# Check the amplitudes.  If any equal "nan" or "inf", skip the entry
		for($amp_ctr = 0; $amp_ctr <= $#amp_array; $amp_ctr++)
		{
			# Remove whitespace
			@amp_test_array = split ' ', lc($array[$amp_array[$amp_ctr]]);
			$amp_test = $amp_test_array[0];

			if($amp_test eq "nan" || $amp_test eq "inf")
			{
				print "NAN/INF: Skipping $array[$name2col{\"file\"}]\n";
				next INPUT_LINE;
			}
		}

		if($data_ctr == 0)
		{
			$matchmode = $array[ $name2col{"matchmode"} ];
		}

		for($ctr2 = 0; $ctr2 <= $#array; $ctr2++)
		{
			$data1[$data_ctr]{ $col2name{$ctr2} } = $array[$ctr2];
		}

		if(defined($name2col{"isofile"}))
		{
		}
		else
		{
			$data1[$data_ctr]{"isofile"} = $array[ $name2col{"file"} ];
		}

		if($was_display_flag == 0)
		{
			$data1[$data_ctr]{"display_flag"} = 1;
		}

		# No precalculated frac_lab, so calculate based on AMP_U and AMP_L
		if($was_frac_lab == 0)
		{
			if($frac_lab_type == 0)
			{
				if( ($data1[$data_ctr]{"amp_l"} + $data1[$data_ctr]{"amp_u"}) != 0)
				{
					$data1[$data_ctr]{"frac_lab"} = $data1[$data_ctr]{"amp_l"} / ($data1[$data_ctr]{"amp_l"} + $data1[$data_ctr]{"amp_u"});
				}
				else
				{
					$data1[$data_ctr]{"frac_lab"} = -0.1;
				}
			}
			elsif($frac_lab_type == 1)
			{
				if( ($data1[$data_ctr]{"amp_u"} + $data1[$data_ctr]{"amp_l"} + $data1[$data_ctr]{"amp_f"}) != 0)
				{
					$data1[$data_ctr]{"frac_lab"} = ($data1[$data_ctr]{"amp_u"} + $data1[$data_ctr]{"amp_l"}) / ($data1[$data_ctr]{"amp_u"} + $data1[$data_ctr]{"amp_l"} + $data1[$data_ctr]{"amp_f"})
				}
				else
				{
					$data1[$data_ctr]{"frac_lab"} = -0.1;
				}
			}
			else
			{
				$data1[$data_ctr]{"frac_lab"} = -0.1;
				printf("Error: \$frac_lab_type not supported\n");
			}
		}

		# Accession numbers included, proteomic style MS/MS data
		if(defined($name2col{"pks_acc_n14"}))
		{
			$pro2acc{$array[$name2col{"protein"}]} = $array[$name2col{"pks_acc_n14"}];
		}

		# 3/7/2011
		# Each _iso.csv file includes an isoid entry, which is supposed to be a unique id
		# Unfortunately, it ceases being unique when multiple files are concatenated
		# Instead of isoid, just overwrite it with $data_ctr, which really is unique
		$data1[$data_ctr]{"isoid"} = $data_ctr;

		$data_ctr++;
	}

	regen_view1();

	plot_frac();

	# Joan has cat'd together multiple datasets
	# This means there are multiple "isofile" and "file" entries
	# So need to get all of them, and store them in these hashes
	# Must do this here as hashes may be incomplete if all entries deleted later

	%isofile_hash = ();
	%file_hash = ();

	for(my $ctr1 = 0; $ctr1 <= $#data1; $ctr1++)
	{
		@dirarray = split /\//, $data1[$ctr1]{"isofile"};
		@dirarray2 = split /\//, $data1[$ctr1]{"file"};

		$isofile_dir = $dirarray[0];
		$file_dir = $dirarray2[0];

		if(!defined($isofile_hash{$isofile_dir}))
		{
			$isofile_hash{$isofile_dir} = 1;
		}
		if(!defined($file_hash{$file_dir}))
		{
			$file_hash{$file_dir} = 1;
		}
	}

	#color_list1_bydisplayflag();
}

sub launch_save
{
	# Try to anticipate a useful name for $outfile
	$outfile = $infile;
	if(substr($outfile, -4, 4) eq ".csv")
	{
		substr($outfile, -4, 4) = "_filt.csv";
	}
	else
	{
		join '', $outfile, ".filt";
	}

	$sw = MainWindow->new();
	$sw->configure(-title=>'Save File As');
	$sw->minsize( qw(600 50) );
	$swframe = $sw->Frame(-border=>'10', -background=>'white')->pack(-fill=>both);
	$swentry = $swframe->Entry(-textvariable=>\$outfile)->pack(-side=>'left', -expand=>'1', -fill=>x);
	$swbutton = $swframe->Button(-text=>'Save File', -command=>\&save)->pack(-side=>'left');
}

sub save
{
	my $ctr, $hctr;

	# Joan has cat'd together multiple datasets
	# This means there are multiple "isofile" and "file" entries
	# So need to get all of them, and restore for each

	#%isofile_hash = ();
	#%file_hash = ();

	#for(my $ctr1 = 0; $ctr1 <= $#data1; $ctr1++)
	#{
	#	@dirarray = split /\//, $data1[$ctr1]{"isofile"};
	#	@dirarray2 = split /\//, $data1[$ctr1]{"file"};
	#
	#	$isofile_dir = $dirarray[0];
	#	$file_dir = $dirarray2[0];
	#
	#	if(!defined($isofile_hash{$isofile_dir}))
	#	{
	#		$isofile_hash{$isofile_dir} = 1;
	#	}
	#	if(!defined($file_hash{$file_dir}))
	#	{
	#		$file_hash{$file_dir} = 1;
	#	}
	#}

	foreach $key (keys %isofile_hash)
	{
		# Remove .msvbkp files (png, fit and dat) and .del files
		#@dirarray = split /\//, $data1[0]{"isofile"}; # fits dir
		@dirarray = split /\//, $key;
		`rm $dirarray[0]/*.msvbkp`; # These are refits, so only occur in fits dir
		`rm $dirarray[0]/*.del`;
	}

	foreach $key (keys %file_hash)
	{
		#@dirarray2 = split /\//, $data1[0]{"file"}; # peaks dir (typically same)
		@dirarray2 = split /\//, $key;
		`rm $dirarray2[0]/*.del`;
	}

	$sw->destroy();
	open OUT, ">$outfile" or die "Can't open $outfile\n";
	print "Writing $outfile...\n";
	print OUT "$header";
	if($was_frac_lab == 0)
	{
		print OUT ",frac_lab";
	}
	if($was_display_flag == 0)
	{
		print OUT ",display_flag";
	}
	print OUT "\n";
	for($ctr = 0; $ctr <= $#data1; $ctr++)
	{
		if($data1[$ctr]{"display_flag"} == 1 || $obey_display_flag == 0 || $only_save_displayed == 0)
		{
			for($hctr = 0; $hctr < $#headerarray; $hctr++)
			{
				print OUT "$data1[$ctr]{$headerarray[$hctr]},";
			}
			print OUT "$data1[$ctr]{$headerarray[$#headerarray]}";

			if($was_frac_lab == 0)
			{
				print OUT ",$data1[$ctr]{'frac_lab'}";
			}

			if($was_display_flag == 0)
			{
				print OUT ",$data1[$ctr]{'display_flag'}";
			}

			print OUT "\n";
		}
	}
	close OUT;

	my $typestring = "";
	# Now automatically run massiso2plotdata.pl on the saved .csv file
	if($plotframe == 0)
	{
		#%pro2loc = %pro2loc30S;
		$typestring = "30S";
	}
	elsif($plotframe == 1)
	{
		#%pro2loc = %pro2loc50S;
		$typestring = "50S";
	}
	elsif($plotframe == 2)
	{
		#%pro2loc = %pro2loc70S;
		$typestring = "70S";
	}
	elsif($plotframe == 3)
	{
		#%pro2loc = %pro2loc30Schlamy;
		$typestring = "30Schlamy";
	}
	elsif($plotframe == 4)
	{
		#%pro2loc = %pro2loc40Shuman;
		$typestring = "40Shuman";
	}
	elsif($plotframe == 5)
	{
		#%pro2loc = %pro2loc60Shuman;
		$typestring = "60Shuman";
	}
	elsif($plotframe == 6)
	{
		#%pro2loc = %pro2loccof;
		$typestring = "cof";
	}
	elsif($plotframe == 7)
	{
		#%pro2loc = %pro2loc40Syeast;
		$typestring = "40Syeast_old";
	}
	elsif($plotframe == 8)
	{
		#%pro2loc = %pro2loc40Syeast;
		$typestring = "40Syeast";
	}
	elsif($plotframe == 9)
	{
		#%pro2loc = %pro2loc40Syeast;
		$typestring = "60Syeast";
	}
	elsif($plotframe == 10)
	{
		#%pro2loc = %pro2loc40Syeast;
		$typestring = "80Syeast";
	}
	elsif($plotframe == 11)
	{
		$typestring = "RNA";
	}
	elsif($plotframe == 12)
	{
		$typestring = "proteomic";
	}
	elsif($plotframe == 13)
	{
		$typestring = "80Shuman";
	}

	`$ENV{'MASSACRE_PATH'}/bin/massiso2plotdata.pl --data=$typestring < $outfile > $outfile.stats`;
}

sub launch_NewContour
{
	if($nc_machine eq "")
	{
		$ncw = MainWindow->new();
		$ncw->configure(-title=>'-> NewContour');
		$ncw_frame = $ncw->Frame(-border=>'5')->pack(-fill=>'both');

		$ncw_userlabel = $ncw_frame->Label(-text=>'User:')->pack(-side=>'top', -anchor=>'w');
		$ncw_user = $ncw_frame->Entry(-textvariable=>\$nc_user, -width=>'8')->pack(-side=>'top', -anchor=>'w');
		$ncw_machinelabel = $ncw_frame->Label(-text=>'Computer:')->pack(-side=>'top', -anchor=>'w');
		$ncw_machine = $ncw_frame->Entry(-textvariable=>\$nc_machine, -width=>'8')->pack(-side=>'top', -anchor=>'w');
		$ncw_go = $ncw_frame->Button(-text=>'Set', -command=>\&set_NewContour)->pack(-side=>'top', -anchor=>'e');
	}
	else
	{
		go_NewContour();
	}
}

sub go_NewContour
{
	isoid_select();
	$nc_isoid = $isoid_curselection;
	$nc_idx = $iso2idx{$nc_isoid};

	if($nc_isoid == -1)
	{
		print "Nothing selected to send to NewContour\n";
		return;
	}

	$nc_file = $data1[$nc_idx]{'file'};
	if(substr($nc_file, -3, 3) eq "txt")
	{
		substr($nc_file, -3, 3) = "newcon";
	}
	else
	{
		print "Weird - file name not ending with .txt\n";
		return;
	}

	`scp $nc_file $nc_user\@$nc_machine:.`;
	@nc_filearray = split /\//, $nc_file;
	`ssh $nc_user\@$nc_machine open -a NewContour $nc_filearray[$#nc_filearray]`;
	`ssh $nc_user\@$nc_machine rm $nc_filearray[$#nc_filearray]`;
}

sub set_NewContour
{
	$ncw->destroy();
}

sub launch_envelope
{
	$nvlp_local = -1;

	$nvlpw = MainWindow->new();
	$nvlpw->configure(-title=>'-> Envelope');

	$nvlp_frame = $nvlpw->Frame(-border=>'5')->pack(-fill=>'both');

	isoid_select();
	$nvlp_isoid = $isoid_curselection;
	$nvlp_idx = $iso2idx{$nvlp_isoid};

	if($nvlp_isoid == -1)
	{
		$nvlpw->destroy();
		print "Nothing selected to send to Envelope\n";
		return;
	}

	$nvlp_start = $name2col{"tim"};
	$nvlp_end = $name2col{"symb"};

	# Detect iso3 style output
	if($nvlp_end == ($nvlp_start + 2))
	{
		$nvlp_end = $name2col{"isoid"};
	}

	$nvlp_sendframe0 = $nvlp_frame->Frame()->pack(-side=>'top', -fill=>'both');
	$nvlp_remotemachine = $nvlp_sendframe0->Entry(-textvariable=>\$nvlp_machine, -width=>'8')->pack(-side=>'right');
	$nvlp_remoteat = $nvlp_sendframe0->Label(-text=>'@')->pack(-side=>'right');
	$nvlp_remoteuser = $nvlp_sendframe0->Entry(-textvariable=>\$nvlp_user, -width=>'8')->pack(-side=>'right');

	$nvlp_sendframe1 = $nvlp_frame->Frame()->pack(-side=>'top', -fill=>'both');
	$nvlp_remotebutton = $nvlp_sendframe1->Button(-text=>'Send (Remote)', -width=>'12', -command=>\&envelope_remote)->pack(-side=>'right');

	$nvlp_sendframe2 = $nvlp_frame->Frame()->pack(-side=>'top', -fill=>'both');
	$nvlp_localbutton = $nvlp_sendframe2->Button(-text=>'Send (Local)', -width=>'12', -command=>\&envelope_local)->pack(-side=>'right');

	$nvlp_buttons_frame = $nvlp_frame->Frame()->pack(-side=>'top', -fill=>'both');

	$nvlp_closebutton = $nvlp_buttons_frame->Button(-text=>'Close', -command=>\&close_envelope)->pack(-side=>'right');
	$nvlp_savebutton = $nvlp_buttons_frame->Button(-text=>'Save', -command=>\&save_envelope)->pack(-side=>'right');

	$nvlp_spacer1 = $nvlp_frame->Frame(-background=>'black', -height=>'2', -borderwidth=>'5')->pack(-side=>'top', -fill=>'both');


	for($nvlp_ctr = $nvlp_start + 1; $nvlp_ctr < $nvlp_end; $nvlp_ctr++)
	{
		$nvlp_var = $col2name{$nvlp_ctr};

		$nvlp_ctr2 = $nvlp_ctr - $nvlp_start - 1;

		if(substr($nvlp_var, -3, 3) ne "err")
		{
			$nvlp_frames[$nvlp_ctr2] = $nvlp_frame->Frame()->pack(-side=>'top', -fill=>'both');

			$nvlp_vars[$nvlp_ctr2] = sprintf("%f", $data1[$nvlp_idx]{$nvlp_var});
			$nvlp_varnames[$nvlp_ctr2] = $nvlp_var;

			$nvlp_entries[$nvlp_ctr2] = $nvlp_frames[$nvlp_ctr2]->Entry(-textvariable=>\$nvlp_vars[$nvlp_ctr2])->pack(-side=>'right');
			$nvlp_labels[$nvlp_ctr2] = $nvlp_frames[$nvlp_ctr2]->Label(-text=>$nvlp_varnames[$nvlp_ctr2])->pack(-side=>'right');
		}
	}

	if($was_frac_lab == 0)
	{
		$nvlp_var = "frac_lab";

		$nvlp_ctr2 = $nvlp_ctr - $nvlp_start - 1;

		if(substr($nvlp_var, -3, 3) ne "err")
		{
			$nvlp_frames[$nvlp_ctr2] = $nvlp_frame->Frame()->pack(-side=>'top', -fill=>'both');

			$nvlp_vars[$nvlp_ctr2] = sprintf("%f", $data1[$nvlp_idx]{$nvlp_var});
			$nvlp_varnames[$nvlp_ctr2] = $nvlp_var;

			$nvlp_entries[$nvlp_ctr2] = $nvlp_frames[$nvlp_ctr2]->Entry(-textvariable=>\$nvlp_vars[$nvlp_ctr2])->pack(-side=>'right');
			$nvlp_labels[$nvlp_ctr2] = $nvlp_frames[$nvlp_ctr2]->Label(-text=>$nvlp_varnames[$nvlp_ctr2])->pack(-side=>'right');
		}
	}


}

sub save_envelope
{
	for($nvlp_ctr3 = 0; $nvlp_ctr3 < $nvlp_ctr2; $nvlp_ctr3++)
	{
		#$nvlp_ctr4 = $nvlp_ctr3 + $nvlp_start + 1;

		$nvlp_var = $nvlp_varnames[$nvlp_ctr3];

		#print "$nvlp_var = $nvlp_vars[$nvlp_ctr3]\n";
		$data1[$nvlp_idx]{$nvlp_var} = $nvlp_vars[$nvlp_ctr3];
	}

	if($was_frac_lab == 0)
	{
		$nvlp_var = "frac_lab";

		$data1[$nvlp_idx]{$nvlp_var} = $nvlp_vars[$nvlp_ctr3];
	}

	if($nvlp_local == 1)
	{
		$nvlp_png = join '', $data1[$nvlp_idx]{'isofile'}, ".fit.png";

		`cp $nvlp_png $nvlp_png.msvbkp`;
		`cp ~/Desktop/massive.png $nvlp_png`;
		`$ENV{'MASSACRE_PATH'}/bin/resizepng.pl --in=$nvlp_png --out=$nvlp_png --w=800 --h=600`;

		# Put up the new png
		$lastphoto->delete;
		$photo = $canvasl->Photo(-format=>'png', -file=>$nvlp_png);
		$canvasl->delete('plotimage');
		$canvasl->createImage(0,0, -image=>$photo, -anchor=>'nw', -tag=>'plotimage');
		$lastphoto = $photo;
	}
	elsif($nvlp_local == 0)
	{
		$nvlp_png = join '', $data1[$nvlp_idx]{'isofile'}, ".fit.png";
		`cp $nvlp_png $nvlp_png.msvbkp`;
		`scp $nvlp_user\@$nvlp_machine:Desktop/massive.png ./$nvlp_png`;
		`ssh $nvlp_user\@$nvlp_machine rm Desktop/massive.png`;

		`$ENV{'MASSACRE_PATH'}/bin/resizepng.pl --in=$nvlp_png --out=$nvlp_png --w=800 --h=600`;

		# Put up the new png
		$lastphoto->delete;
		$photo = $canvasl->Photo(-format=>'png', -file=>$nvlp_png);
		$canvasl->delete('plotimage');
		$canvasl->createImage(0,0, -image=>$photo, -anchor=>'nw', -tag=>'plotimage');
		$lastphoto = $photo;
	}

	$nvlpw->destroy();
}

sub close_envelope
{
	if(-e "nvlp_script")
	{
		`rm nvlp_script`;
	}
	$nvlpw->destroy();
}

sub envelope_local
{
	$nvlp_local = 1;

	open NVLP, ">nvlp_script" or print "Error: Can't open nvlp_script\n";

	if(defined($name2col{"seqmod"}))
	{
		print NVLP "sequence $data1[$nvlp_idx]{'seqmod'}\n";
	}
	else
	{
		print NVLP "sequence $data1[$nvlp_idx]{'seq'}\n";
	}
	print NVLP "charge $data1[$nvlp_idx]{'charge'}\n";
	print NVLP "type 0\n"; # 0 for protein, 1 for RNA -> protein only for now
	printf(NVLP "baseline %f\n", $data1[$nvlp_idx]{'b'});
	printf(NVLP "gwidth %f\n", $data1[$nvlp_idx]{'gw'});
	printf(NVLP "offset %f\n", $data1[$nvlp_idx]{'off'});

	close NVLP;

	`open -a /Applications/Envelope.app nvlp_script`;
	`open -a /Applications/Envelope.app $data1[$nvlp_idx]{'file'}`;
}

sub envelope_remote
{
	$nvlp_local = 0;

	open NVLP, ">nvlp_script" or print "Error: Can't open nvlp_script\n";

	if(defined($name2col{"seqmod"}))
	{
		print NVLP "sequence $data1[$nvlp_idx]{'seqmod'}\n";
	}
	else
	{
		print NVLP "sequence $data1[$nvlp_idx]{'seq'}\n";
	}
	print NVLP "charge $data1[$nvlp_idx]{'charge'}\n";
	print NVLP "type 0\n"; # 0 for protein, 1 for RNA -> protein only for now
	printf(NVLP "baseline %f\n", $data1[$nvlp_idx]{'b'});
	printf(NVLP "gwidth %f\n", $data1[$nvlp_idx]{'gw'});
	printf(NVLP "offset %f\n", $data1[$nvlp_idx]{'off'});

	close NVLP;

	`scp nvlp_script $nvlp_user\@$nvlp_machine:Desktop/`;
	# peaks_dir/peaks_file.txt -> nvlp_machine:Desktop/peaks_file.txt
	`scp $data1[$nvlp_idx]{'file'} $nvlp_user\@$nvlp_machine:Desktop/`;

	@nvlp_file_array = split /\//, $data1[$nvlp_idx]{'file'};
	$nvlp_file = $nvlp_file_array[$#nvlp_file_array];

	`ssh $nvlp_user\@$nvlp_machine open -a /Applications/Envelope.app Desktop/nvlp_script`;
	#`ssh $nvlp_user\@$nvlp_machine open -a /Applications/Envelope.app Desktop/$data1[$nvlp_idx]{'file'}`;
	`ssh $nvlp_user\@$nvlp_machine open -a /Applications/Envelope.app Desktop/$nvlp_file`;

	# Will these remove files before fully loaded?  Seems to work okay for ssh
	# Put in a sleep command to be safe
	sleep(2);
	`ssh $nvlp_user\@$nvlp_machine rm Desktop/nvlp_script`;
	`ssh $nvlp_user\@$nvlp_machine rm Desktop/$nvlp_file`;
}

sub launch_refit
{
	$rw = MainWindow->new();
	$rw->configure(-title=>'Refit Peak');
	$rw->minsize( qw(480 640) );

	$rwframe = $rw->Frame(-border=>'5')->pack(-side=>'top', -fill=>'both');
	$rwframe2 = $rw->Frame(-border=>'5')->pack(-side=>'top', -fill=>'both');

	isoid_select();
	$rw_isoid = $isoid_curselection;
	$rw_idx = $iso2idx{$rw_isoid};

	if($rw_isoid == -1)
	{
		$rw->destroy();
		print "Nothing selected to refit\n";
		return;
	}

	$rwlabelframe = $rwframe->Frame()->pack(-side=>'top', -fill=>'both');
		$rw_label = $rwlabelframe->Label(-text=>$view1[$rw_idx], -font=>'*-helvetica-medium-r-normal-16-*', -border=>'5')->pack(-side=>'top', -anchor=>'w');
		$rw_label2 = $rwlabelframe->Label(-text=>$data1[$rw_idx]{"seqmod"}, -font=>'*-courier-medium-r-normal-16-*', -border=>'5')->pack(-side=>'top', -anchor=>'w');

	$rwspaceframe1 = $rwframe->Frame(-background=>'black', -height=>'2')->pack(-side=>'top', -fill=>'both');

	$rwvarframe = $rwframe->Frame()->pack(-side=>'top', -fill=>'both');
	$rwvarframe1 = $rwvarframe->Frame()->pack(-side=>'left', -fill=>'both');
	$rwvarframe1b = $rwvarframe->Frame()->pack(-side=>'left');
	$rwvarframe2 = $rwvarframe->Frame()->pack(-side=>'left', -fill=>'both');

		$rwf1 = $rwvarframe1->Frame()->pack(-side=>'top', -fill=>'both');
			$b_entry = $rwf1->Entry(-textvariable=>\$b_in)->pack(-side=>'left');
			$b_label = $rwf1->Label(-text=>'b')->pack(-side=>'left');
		$rwf2 = $rwvarframe1->Frame()->pack(-side=>'top', -fill=>'both');
			$off_entry = $rwf2->Entry(-textvariable=>\$off_in)->pack(-side=>'left');
			$off_label = $rwf2->Label(-text=>'off')->pack(-side=>'left');
		$rwf3 = $rwvarframe1->Frame()->pack(-side=>'top', -fill=>'both');
			$gw_entry = $rwf3->Entry(-textvariable=>\$gw_in)->pack(-side=>'left');
			$gw_label = $rwf3->Label(-text=>'gw')->pack(-side=>'left');
		$rwf4 = $rwvarframe1->Frame()->pack(-side=>'top', -fill=>'both');
			$ul_amp_entry = $rwf4->Entry(-textvariable=>\$ul_amp_in)->pack(-side=>'left');
			$ul_amp_label = $rwf4->Label(-text=>'ul_amp')->pack(-side=>'left');
		$rwf5 = $rwvarframe1->Frame()->pack(-side=>'top', -fill=>'both');
			$l_amp_entry = $rwf5->Entry(-textvariable=>\$l_amp_in)->pack(-side=>'left');
			$l_amp_label = $rwf5->Label(-text=>'l_amp')->pack(-side=>'left');
		$rwf6 = $rwvarframe1->Frame()->pack(-side=>'top', -fill=>'both');
			$frac_n_entry = $rwf6->Entry(-textvariable=>\$frac_n_in)->pack(-side=>'left');
			$frac_n_label = $rwf6->Label(-text=>'frac_n')->pack(-side=>'left');
		$rwf7 = $rwvarframe1->Frame()->pack(-side=>'top', -fill=>'both');
			$frac_lab_entry = $rwf7->Entry(-textvariable=>\$frac_lab_in, -state=>'readonly', -relief=>'groove')->pack(-side=>'left');
			$frac_lab_label = $rwf7->Label(-text=>'frac_lab')->pack(-side=>'left');
		$rwf8 = $rwvarframe1->Frame()->pack(-side=>'top', -fill=>'both');
			$chisq_entry = $rwf8->Entry(-textvariable=>\$chisq_in, -state=>'readonly', -relief=>'groove')->pack(-side=>'left');
			$chisq_label = $rwf8->Label(-text=>'chisq')->pack(-side=>'left');

		$rw_arrow = $rwvarframe1b->Label(-text=>'--->')->pack(-fill=>'both');
		$rw_copybutton = $rwvarframe1b->Button(-text=>'<---', -command=>\&copy_fitvars)->pack(-fill=>'both');


		$rwf1new = $rwvarframe2->Frame()->pack(-side=>'top', -fill=>'both');
			$b_entry_new = $rwf1new->Entry(-textvariable=>\$b_out, -state=>'readonly', -relief=>'groove')->pack(-side=>'left');
		$rwf2new = $rwvarframe2->Frame()->pack(-side=>'top', -fill=>'both');
			$off_entry_new = $rwf2new->Entry(-textvariable=>\$off_out, -state=>'readonly', -relief=>'groove')->pack(-side=>'left');
		$rwf3new = $rwvarframe2->Frame()->pack(-side=>'top', -fill=>'both');
			$gw_entry_new = $rwf3new->Entry(-textvariable=>\$gw_out, -state=>'readonly', -relief=>'groove')->pack(-side=>'left');
		$rwf4new = $rwvarframe2->Frame()->pack(-side=>'top', -fill=>'both');
			$ul_amp_entry_new = $rwf4new->Entry(-textvariable=>\$ul_amp_out, -state=>'readonly', -relief=>'groove')->pack(-side=>'left');
		$rwf5new = $rwvarframe2->Frame()->pack(-side=>'top', -fill=>'both');
			$l_amp_entry_new = $rwf5new->Entry(-textvariable=>\$l_amp_out, -state=>'readonly', -relief=>'groove')->pack(-side=>'left');
		$rwf6new = $rwvarframe2->Frame()->pack(-side=>'top', -fill=>'both');
			$frac_n_entry_new = $rwf6new->Entry(-textvariable=>\$frac_n_out, -state=>'readonly', -relief=>'groove')->pack(-side=>'left');
			$stephenbutton_fracn = $rwf6new->Button(-text=>'Stephen', -command=>\&set_global_frac_n)->pack(-side=>'left');
		$rwf7new = $rwvarframe2->Frame()->pack(-side=>'top', -fill=>'both');
			$frac_lab_entry_new = $rwf7new->Entry(-textvariable=>\$frac_lab_out, -state=>'readonly', -relief=>'groove')->pack(-side=>'left');
		$rwf8new = $rwvarframe2->Frame()->pack(-side=>'top', -fill=>'both');
			$chisq_entry_new = $rwf8new->Entry(-textvariable=>\$chisq_out, -state=>'readonly', -relief=>'groove')->pack(-side=>'left');


	$b_in = sprintf("%f", $data1[$rw_idx]{"b"});
	$off_in = sprintf("%f", $data1[$rw_idx]{"off"});
	$gw_in = sprintf("%f", $data1[$rw_idx]{"gw"});
	$ul_amp_in = sprintf("%f", $data1[$rw_idx]{"ul_amp"});
	$l_amp_in = sprintf("%f", $data1[$rw_idx]{"l_amp"});

	if(defined($data1[$rw_idx]{"frac_n"}))
	{
		$frac_n_in = sprintf("%f", $data1[$rw_idx]{"frac_n"});
	}
	elsif(defined($data1[$rw_idx]{"frc_nx"}))
	{
		$frac_n_in = sprintf("%f", $data1[$rw_idx]{"frc_nx"});
	}
	else
	{
		$frac_n_in = $global_frac_n;
	}
	$frac_lab_in = sprintf("%f", $data1[$rw_idx]{"frac_lab"});
	$chisq_in = sprintf("%f", $data1[$rw_idx]{"chisq"});

	$b_out = "";
	$off_out = "";
	$gw_out = "";
	$ul_amp_out = "";
	$l_amp_out = "";
	$frac_n_out = "";
	$frac_lab_out = "";
	$chisq_out = "";

	$rwgoframe = $rwframe->Frame(-pady=>'5')->pack(-side=>'top', -anchor=>'nw');
		#$fitgobutton = $rwgoframe->Button(-text=>'Fit!', -command=>\&refit, -background=>'green', -activebackground=>'lightgreen')->pack(-side=>'left');
		$fititbutton = $rwgoframe->Button(-text=>'Fit It', -command=>\&call_fitit, -background=>'green', -activebackground=>'lightgreen')->pack(-side=>'left');
		$tryitbutton = $rwgoframe->Button(-text=>'Try It', -command=>\&call_tryit)->pack(-side=>'left');



		$keeparrow = $rwgoframe->Label(-text=>'--->')->pack(-side=>'left');
		$keepbutton = $rwgoframe->Button(-text=>'Keep New Fit', -command=>\&keep_new_fit)->pack(-side=>'left');

	$rwplotvarframe = $rwframe->Frame(-border=>'1', -background=>'black')->pack(-side=>'top', -anchor=>'w');
		$datframe = $rwplotvarframe->Frame()->pack(-side=>'top', -fill=>'both');
			$datlab = $datframe->Label(-text=>'dat:', -width=>'3')->pack(-side=>'left');
			$datradio1 = $datframe->Radiobutton(-text=>'pts', -variable=>\$datplotstyle, -value=>'0')->pack(-side=>'left');
			$datradio2 = $datframe->Radiobutton(-text=>'lines', -variable=>\$datplotstyle, -value=>'1')->pack(-side=>'left');
		$fitframe = $rwplotvarframe->Frame()->pack(-side=>'top', -fill=>'both');
			$fitlab = $fitframe->Label(-text=>'fit:', -width=>'3')->pack(-side=>'left');
			$fitradio1 = $fitframe->Radiobutton(-text=>'pts', -variable=>\$fitplotstyle, -value=>'0')->pack(-side=>'left');
			$fitradio2 = $fitframe->Radiobutton(-text=>'lines', -variable=>\$fitplotstyle, -value=>'1')->pack(-side=>'left');
		$xrangeframe = $rwplotvarframe->Frame()->pack(-side=>'top', -fill=>'both');
			$xrangelab = $xrangeframe->Label(-text=>'xrange:', -width=>'7')->pack(-side=>'left');
			$xrangeentrymin = $xrangeframe->Entry(-textvariable=>\$xrangemin, -width=>'8')->pack(-side=>'left');
			$xrangeentrymin = $xrangeframe->Entry(-textvariable=>\$xrangemax, -width=>'8')->pack(-side=>'left');
			$xrangelab2 = $xrangeframe->Label(-text=>'(0s for auto)')->pack(-side=>'left');
			$stephenbutton = $xrangeframe->Button(-text=>'Stephen', -command=>\&pop_xrange)->pack(-side=>'left');
		$yrangeframe = $rwplotvarframe->Frame()->pack(-side=>'top', -fill=>'both');
			$yrangelab = $yrangeframe->Label(-text=>'yrange:', -width=>'7')->pack(-side=>'left');
			$yrangeentrymin = $yrangeframe->Entry(-textvariable=>\$yrangemin, -width=>'8')->pack(-side=>'left');
			$yrangeentrymax = $yrangeframe->Entry(-textvariable=>\$yrangemax, -width=>'8')->pack(-side=>'left');
			$yrangelab2 = $yrangeframe->Label(-text=>'(0s for auto)')->pack(-side=>'left');

	$rwaddframe = $rwframe->Frame(-pady=>'5')->pack(-side=>'top', -fill=>'both');
		$addroundbutton = $rwaddframe->Button(-text=>'Add Round', -command=>\&add_round)->pack(-side=>'left');
		$defaultbutton = $rwaddframe->Button(-text=>'Default', -command=>\&default_schedule)->pack(-side=>'left');
		$ftframe = $rwaddframe->Frame(-background=>'black', -border=>'1')->pack(-side=>'left');
		#$fititbutton = $ftframe->Radiobutton(-text=>'Fit It', -variable=>\$tryitfitit, -value=>'0')->pack(-side=>'left');
		#$tryitbutton = $ftframe->Radiobutton(-text=>'Try It', -variable=>\$tryitfitit, -value=>'1')->pack(-side=>'left');
		#$fititbutton = $ftframe->Button(-text=>'Fit It', -command=>\&call_fitit)->pack(-side=>'left');
		#$tryitbutton = $ftframe->Button(-text=>'Try It', -command=>\&call_tryit)->pack(-side=>'left');

	$rwspaceframe2 = $rwframe2->Frame(-background=>'black', -height=>'2')->pack(-side=>'top', -fill=>'both');
	$rwcloseframe = $rwframe2->Frame()->pack(-side=>'top', -fill=>'both');
		$rw_closebutton = $rwcloseframe->Button(-text=>'Close', -command=>\&close_refit)->pack(-side=>'left');

	build_schedule();
}

sub set_global_frac_n
{
	$global_frac_n = $frac_n_in;
}

sub call_tryit
{
	$tryitfitit = 1;
	refit();
}

sub call_fitit
{
	$tryitfitit = 0;
	refit();
}

sub pop_xrange
{
	my $n14value = $data1[$rw_idx]{"n14mass"};
	#$n14value = ($n14value - (1.5 / $data1[$rw_idx]{"charge"}))*$data1[$rw_idx]{"charge"};
	$n14value = $n14value * $data1[$rw_idx]{"charge"};
	my $n15value = $data1[$rw_idx]{"n15mass"};
	#$n15value = ($n15value*1.005 + (2.5/$data1[$rw_idx]{"charge"}))*$data1[$rw_idx]{"charge"};
	$n15value = $n15value * $data1[$rw_idx]{"charge"};
	my $midpoint = ($n14value + $n15value) / 2;
	my $fraction = ($n15value - $n14value)/4;
	$xrangemin = $midpoint - $fraction;
	$xrangemax = $midpoint + $fraction + 1;
}

sub copy_fitvars
{
	$b_in = $b_out;
	$gw_in = $gw_out;
	$off_in = $off_out;
	$ul_amp_in = $ul_amp_out;
	$l_amp_in = $l_amp_out;
	$frac_n_in = $frac_n_out;
	$frac_lab_in = $frac_lab_out;
	$chisq_in = $chisq_out;
}

#####
#fitit (specify routine to run)
#          15N_pulse_short.batch (batchfile)
#          3 (fit_model 1=Natural Abundance, 2=Fractional Label, 3=Unlabeled + Labeled)
#          4 (niter -> number of iterations?)
#          100.0 (sig_global)
#          1.0 (baseline)
#          1.2 (unlabeled_amp)
#          0.5 (labeled_amp)
#          0.01 (offset)
#          0.0003 (gaussian_width)
#          0.82 (fraction_labeled)
#5 (nround -> number of fitting rounds, below)
#6 (nfitpar -> number of parameters to fit, below)
#B OFF GW UL_AMP LAB_AMP FRAC (fit schedule - next lines indicate whether to fit each variable on that pass
#1 0 0 1 1 0
#1 0 0 0 0 1
#1 0 1 1 1 1
#1 1 1 0 0 1
#1 1 1 1 1 1
#####
sub refit
{
	regen_schedule();

	# Generate the input parameter file
	open IN, ">isodist_batch.in" or die "Can't open isodist_batch.in\n";
	if($tryitfitit == 0)
	{
		print IN "fitit\n";
	}
	elsif($tryitfitit == 1)
	{
		print IN "tryit\n";
	}
	print IN "          massive.batch\n";
	print IN "          3\n";
	print IN "          4\n";
	print IN "          100.0\n";
	printf(IN "          %f\n", $b_in);
	printf(IN "          %f\n", $ul_amp_in);
	printf(IN "          %f\n", $l_amp_in);
	printf(IN "          %f\n", $off_in);
	printf(IN "          %f\n", $gw_in);
	printf(IN "          %f\n", $frac_n_in);
	print IN "$nround\n";
	print IN "$nvar\n";
	print IN "B OFF GW UL_AMP LAB_AMP FRAC\n";
	my @array = split ' ', $fitschedule;
	for(my $ctr1 = 0; $ctr1 < $nround; $ctr1++)
	{
		for(my $ctr2 = 0; $ctr2 < $nvar; $ctr2++)
		{
			print IN "$array[$ctr1*$nvar + $ctr2] ";
		}
		print IN "\n";
	}

	close IN;

	# Generate the batch file
	open BATCH, ">massive.batch" or die "Can't open massive.batch\n";
	print BATCH "$data1[$rw_idx]{'seqmod'} $data1[$rw_idx]{'charge'} refit.txt\n";
	close BATCH;

	# Copy the appropriate file into refit.txt
	$refit_txt = $data1[$rw_idx]{"file"}; # Comes from peaks directory (msc04)
	$refit_base = $data1[$rw_idx]{"isofile"}; # Possible alternate fits directory (msc05)
	if(substr($refit_txt, -3, 3) eq "txt")
	{
	}
	else # Should not occur as all version of massacre output have {"file"}
	{
		$refit_txt = join '', $refit_txt, ".txt";
	}

	# Alternate case for early massacre output which lacked {"isofile"}
	if(substr($refit_base, -3, 3) eq "txt")
	{
		$refit_dat = $refit_base;
		substr($refit_dat, -3, 3) = "dat";
		$refit_fit = $refit_base;
		substr($refit_fit, -3, 3) = "fit";
		$refit_png = $refit_base;
		substr($refit_png, -3, 3) = "fit.png";
	}
	# Proper {"isofile"} output should trigger this
	else
	{
		$refit_dat = join '', $refit_base, ".dat";
		$refit_fit = join '', $refit_base, ".fit";
		$refit_png = join '', $refit_base, ".fit.png";
	}

	`cp $refit_txt refit.txt`;
	`$ENV{'MASSACRE_PATH'}/bin/isodist_batch`;

	print GIN1 "set autoscale\n";

	if($xrangemin != 0 || $xrangemax != 0)
	{
		print GIN1 "set xrange[$xrangemin to $xrangemax]\n";
	}

	if($yrangemin != 0 || $yrangemax != 0)
	{
		print GIN1 "set yrange[$yrangemin to $yrangemax]\n";
	}

	if($datplotstyle == 0)
	{
		print GIN1 "plot \"refit.dat\", ";
	}
	elsif($datplotstyle == 1)
	{
		print GIN1 "plot \"refit.dat\" with lines ,";
	}

	if($fitplotstyle == 0)
	{
		print GIN1 "\"refit.fit\"\n";
	}
	elsif($fitplotstyle == 1)
	{
		print GIN1 "\"refit.fit\" with lines\n";
	}

	#print GIN1 "plot \"refit.dat\", \"refit.fit\" with lines\n";

	# Read in new parameters
	open CSV, "massive.csv" or die "Can't open massive.csv\n";
	$csvheader = <CSV>;
	chomp($csvheader);
	@csvarray = split /\,/, $csvheader;
	for(my $ctr = 0; $ctr <= $#csvarray; $ctr++)
	{
		$csv_name2col{$csvarray[$ctr]} = $ctr;
	}
	$newparamline = <CSV>;
	chomp($newparamline);
	@newparamarray = split /\,/, $newparamline;

	close CSV;

	$b_out = sprintf("%f", $newparamarray[ $csv_name2col{"b"} ]);
	$gw_out = sprintf("%f", $newparamarray[ $csv_name2col{"gw"} ]);
	$off_out = sprintf("%f", $newparamarray[ $csv_name2col{"off"} ]);
	$ul_amp_out = sprintf("%f", $newparamarray[ $csv_name2col{"ul_amp"} ]);
	$l_amp_out = sprintf("%f", $newparamarray[ $csv_name2col{"l_amp"} ]);
	$frac_n_out = sprintf("%f", $newparamarray[ $csv_name2col{"frac_n"} ]);
	$frac_lab_out = sprintf("%f", $newparamarray[ $csv_name2col{"frac_lab"} ]);
	$chisq_out = sprintf("%f", $newparamarray[ $csv_name2col{"chisq"} ]);
}

sub keep_new_fit
{
	if($b_out eq "")
	{
		print "No new fit\n";
		return;
	}

	# Move variables from *_out -> *_in
	$b_in = $b_out;
	$gw_in = $gw_out;
	$off_in = $off_out;
	$ul_amp_in = $ul_amp_out;
	$l_amp_in = $l_amp_out;
	$frac_n_in = $frac_n_out;

	# Copy relevant variables into @data1
	$data1[$rw_idx]{"b"} = $b_out;
	$data1[$rw_idx]{"b_err"} = $newparamarray[ $csv_name2col{"b_err"} ];
	$data1[$rw_idx]{"gw"} = $gw_out;
	$data1[$rw_idx]{"gw_err"} = $newparamarray[ $csv_name2col{"gw_err"} ];
	$data1[$rw_idx]{"off"} = $off_out;
	$data1[$rw_idx]{"off_err"} = $newparamarray[ $csv_name2col{"off_err"} ];
	$data1[$rw_idx]{"ul_amp"} = $ul_amp_out;
	$data1[$rw_idx]{"ul_amp_err"} = $newparamarray[ $csv_name2col{"ul_amp_err"} ];
	$data1[$rw_idx]{"l_amp"} = $l_amp_out;
	$data1[$rw_idx]{"l_amp_err"} = $newparamarray[ $csv_name2col{"l_amp_err"} ];
	$data1[$rw_idx]{"frac_n"} = $frac_n_out;
	$data1[$rw_idx]{"frac_n_err"} = $newparamarray[ $csv_name2col{"frac_n_err"} ];
	$data1[$rw_idx]{"frac_lab"} = $newparamarray[ $csv_name2col{"frac_lab"} ];
	$data1[$rw_idx]{"frac_lab_err"} = $newparamarray[ $csv_name2col{"frac_lab_err"} ];

	$data1[$rw_idx]{"tim"} = $newparamarray[ $csv_name2col{"tim"} ];
	$data1[$rw_idx]{"chisq"} = $newparamarray[ $csv_name2col{"chisq"} ];

	# Copy variables into secondary/tertiary arrays (overkill, but just in case)
	for(my $ctr = 0; $ctr <= $#data2; $ctr++)
	{
		if($data2[$ctr]{"isoid"} == $rw_isoid)
		{
			$data2[$ctr] = $data1[$rw_idx];
		}
	}
	for(my $ctr = 0; $ctr <= $#data3; $ctr++)
	{
		if($data3[$ctr]{"isoid"} == $rw_isoid)
		{
			$data3[$ctr] = $data1[$rw_idx];
		}
	}
	for(my $ctr = 0; $ctr <= $#data4; $ctr++)
	{
		if($data4[$ctr]{"isoid"} == $rw_isoid)
		{
			$data4[$ctr] = $data1[$rw_idx];
		}
	}
	for(my $ctr = 0; $ctr <= $#data2b; $ctr++)
	{
		if($data2b[$ctr]{"isoid"} == $rw_isoid)
		{
			$data2b[$ctr] = $data1[$rw_idx];
		}
	}
	for(my $ctr = 0; $ctr <= $#data3b; $ctr++)
	{
		if($data3b[$ctr]{"isoid"} == $rw_isoid)
		{
			$data3b[$ctr] = $data1[$rw_idx];
		}
	}
	for(my $ctr = 0; $ctr <= $#data4b; $ctr++)
	{
		if($data4b[$ctr]{"isoid"} == $rw_isoid)
		{
			$data4b[$ctr] = $data1[$rw_idx];
		}
	}

	# Generate new png of fit
	open GNU_IN, ">massive.gnuplot" or die "Can't open massive.gnuplot\n";
	$hostname = `hostname`;
	chomp($hostname);
	if($hostname eq "xtopher")
	{
		print GNU_IN "set term png picsize 800 600\n";
		print GNU_IN "set output \"massive.png\"\n";
		print GNU_IN "plot \"refit.dat\" with points pointtype 4 pointsize 1.5 title \"$refit_dat\", ";
		print GNU_IN "\"refit.fit\" with lines title \"$refit_fit\"\n";
	}
	else
	{
		print GNU_IN "set term png medium size 800,600\n";
		print GNU_IN "set output \"massive.png\"\n";
		print GNU_IN "plot \"refit.dat\" with points pointtype 2 pointsize 1.0 title \"$refit_dat\", ";
		print GNU_IN "\"refit.fit\" with lines title \"$refit_fit\"\n";
	}
	#print GNU_IN "set output \"massive.png\"\n";
	#print GNU_IN "plot \"refit.dat\" with points pointtype 4 pointsize 1.5, ";
	#print GNU_IN "\"refit.fit\" with lines\n";
	close GNUPLOT;

	$gnuline = `/usr/bin/gnuplot < massive.gnuplot`;

	# Backup old png dat fit in case the changes are not saved
	`cp $refit_png $refit_png.msvbkp`;
	# For space reasons we no longer keep .fit and .dat files
	#`cp $refit_fit $refit_fit.msvbkp`;
	#`cp $refit_dat $refit_dat.msvbkp`;

	# Copy png etc
	`cp massive.png $refit_png`;
	# No longer keep .fit and .dat files
	#`cp refit.fit $refit_fit`;
	#`cp refit.dat $refit_dat`;
	# screw the contour.png

	# Put up the new png
	$lastphoto->delete;
	$photo = $canvasl->Photo(-format=>'png', -file=>$refit_png);
	$canvasl->delete('plotimage');
	$canvasl->createImage(0,0, -image=>$photo, -anchor=>'nw', -tag=>'plotimage');

	$lastphoto = $photo;

	# Update the frac_lab plot
	plot_frac();

	# Update the highlighted frac_lab points as locations may have changed
	$previous_isoid = $data1[$lastsel1]{"isoid"};
	$previous_tag = join '', "pt", $previous_isoid;
	$canvasr->delete($selectrect);
	@selectcoor = $canvasr->coords($previous_tag);
	$selectcoor[0] -= 1;
	$selectcoor[1] -= 1;
	$selectcoor[2] += 2;
	$selectcoor[3] += 2;
	$selectrect = $canvasr->createRectangle(@selectcoor, -outline=>'#FF0000', -width=>'2');

	color_pts234();
}

sub build_schedule
{
	#print "$fitschedule\n";
	@fitarray = split ' ', $fitschedule;
	$nround = ($#fitarray + 1) / $nvar;
	for(my $ctr = 0; $ctr < $nround; $ctr++)
	{
		$fitframe[$ctr] = $rwframe->Frame()->pack(-side=>'top', -fill=>'both');
		$roundnum[$ctr] = join '', "Round ", $ctr + 1, ":";
		$roundlabel[$ctr] = $fitframe[$ctr]->Label(-text=>$roundnum[$ctr])->pack(-side=>'left');

		for(my $ctr2 = 0; $ctr2 < $nvar; $ctr2++)
		{
			$fitarray_2d[$ctr][$ctr2] = $fitarray[$nvar*$ctr + $ctr2];

			$fitchecks[$ctr][$ctr2] = $fitframe[$ctr]->Checkbutton(-text=>$fitvars[$ctr2], -offvalue=>'0', -onvalue=>'1', -variable=>\$fitarray_2d[$ctr][$ctr2])->pack(-side=>'left');
		}

		$fitbuttons[$ctr] = $fitframe[$ctr]->Button(-text=>'Delete', -command=>[\&delete_round, $ctr])->pack(-side=>'left');
	}

	regen_schedule();
	#print "$fitschedule\n";
}

sub default_schedule
{
	for(my $ctr = 0; $ctr < $nround; $ctr++)
	{
		$fitframe[$ctr]->destroy();
	}

	$fitschedule = $default_fitschedule;

	build_schedule();
}

sub delete_round
{
	my $ctr1;
	my $ctr2;

	my @params = @_;
	my $which_round = $params[0];

	@fitarray = ();
	for($ctr1 = 0; $ctr1 < $nround; $ctr1++)
	{
		if($ctr1 == $which_round)
		{
			next;
		}

		for($ctr2 = 0; $ctr2 < $nvar; $ctr2++)
		{
			push(@fitarray, $fitarray_2d[$ctr1][$ctr2]);
		}
	}

	$fitschedule = join ' ', @fitarray;

	for($ctr1 = 0; $ctr1 < $nround; $ctr1++)
	{
		$fitframe[$ctr1]->destroy();
	}

	build_schedule();
}

sub add_round
{
	my $ctr;

	for($ctr = 0; $ctr < $nvar; $ctr++)
	{
		$fitschedule = join ' ', $fitschedule, "0";
	}

	for($ctr = 0; $ctr < $nround; $ctr++)
	{
		$fitframe[$ctr]->destroy();
	}

	build_schedule();
}

sub regen_schedule
{
	my $ctr1;
	my $ctr2;

	$fitschedule = "";
	for($ctr1 = 0; $ctr1 < $nround; $ctr1++)
	{
		for($ctr2 = 0; $ctr2 < $nvar; $ctr2++)
		{
			$fitschedule = join ' ', $fitschedule, $fitarray_2d[$ctr1][$ctr2];
		}
	}
}

sub close_refit
{
	regen_schedule();

	$rw->destroy();
}

sub show_all_plot
{
	my $ctr;

	for($ctr = 0; $ctr <= $#data1; $ctr++)
	{
		$data1[$ctr]{"display_flag"} = 1;
	}

	plot_frac();
	#color_list1_bydisplayflag();
	$current_display_flag = 1;
}

sub hide_all_plot
{
	my $ctr;

	for($ctr = 0; $ctr <= $#data1; $ctr++)
	{
		$data1[$ctr]{"display_flag"} = 0;
	}

	plot_frac();
	#color_list1_bydisplayflag();
	$current_display_flag = 0;
}

sub plot_frac
{
	my $ctr;
	my $rad = 3;

	$canvasr->delete('pts');
	@ovals = ();

	for($ctr = 0; $ctr <= $#data1; $ctr++)
	{
		if($data1[$ctr]{"display_flag"} == 0 && $obey_display_flag == 1)
		{
			next;
		}

		$protein = $data1[$ctr]{"protein"};
		$frac_lab = $data1[$ctr]{"frac_lab"};
		$isoid = $data1[$ctr]{"isoid"};
		$tag = join '', "pt", $isoid;

		#print "-$#data1-$ctr-$protein-$isoid-$frac_lab-\n";

		if(!defined($pro2loc{$protein}))
		{
			print "ERROR: $protein location not defined\n";
		}

		$x = $xoffset + $xint*$pro2loc{$protein};
		$y = $yticloc[0] + $frac_lab*($yticloc[$#yticloc] - $yticloc[0]);

		$x1 = $x-$rad;
		$x2 = $x+$rad;
		$y1 = $y-$rad;
		$y2 = $y+$rad;

		@{ $ovals[$ctr] } = ($x, $y);
		if($display_abx)
		{
			$prot_abx = substr($protein, -1, 1);

			if($prot_abx eq "A")
			{
				$canvasr->createRectangle($x1, $y1, $x2, $y2, -outline=>'#000000', -tags=>[$tag, 'pts']);
			}
			elsif($prot_abx eq "B")
			{

				$canvasr->createPolygon($x1, $y, $x, $y1, $x2, $y, $x, $y2, -outline=>'#000000', -tags=>[$tag, 'pts']);
			}
			else
			{
				$canvasr->createOval($x1, $y1, $x2, $y2, -outline=>'#000000', -tags=>[$tag, 'pts']);
			}
		}
		else
		{
			$canvasr->createOval($x1, $y1, $x2, $y2, -outline=>'#000000', -tags=>[$tag, 'pts']);
		}
	}

	print "----------\n";
}

sub adjust_plot_width
{
	if($xint_type == 0)
	{
		$xint = 10;
	}
	elsif($xint_type == 1)
	{
		$xint = 15;
	}
	elsif($xint_type == 2)
	{
		$xint = 20;
	}

	draw_frame();
	plot_frac();
	process_selection1();
}

# Launch a new window with potential plot configuration options
# Initially, just allow the user to toggle what is displayed as frac_lab
# So scan through the header, identify all the "amp" columns, and offer up as options
sub config_plot
{
	@amps = ();

	foreach $key (keys %name2col)
	{
		if($key eq "ul_amp" || $key eq "l_amp" || substr($key, 0, 4) eq "amp_")
		{
			push(@amps, $key);
		}
	}

	@amps = sort @amps;

	$pw = MainWindow->new();
	$pw->configure(-title=>'Plot Configuration', -background=>'black', -border=>'2');

	$pw_top = $pw->Frame()->pack(-side=>'top');
	$pw_bot = $pw->Frame()->pack(-side=>'bottom', -fill=>'both');

	$pw_l = $pw_top->Frame(-border=>'2')->pack(-side=>'left');
	$pw_r = $pw_top->Frame(-border=>'2')->pack(-side=>'right');

	$pw_l_label = $pw_l->Label(-text=>'These:')->pack(-side=>'top');
	$pw_r_label = $pw_r->Label(-text=>'/ These:')->pack(-side=>'top');

	for($ampctr = 0; $ampctr <= $#amps; $ampctr++)
	{
		$frac_lab_num[$ampctr] = 0;
		$frac_lab_den[$ampctr] = 0;

		$pw_l_check[$ampctr] = $pw_l->Checkbutton(-text=>$amps[$ampctr], -offvalue=>'0', -onvalue=>'1', -variable=>\$frac_lab_num[$ampctr])->pack(-side=>'top');
		$pw_r_check[$ampctr] = $pw_r->Checkbutton(-text=>$amps[$ampctr], -offvalue=>'0', -onvalue=>'1', -variable=>\$frac_lab_den[$ampctr])->pack(-side=>'top');
		print "$amps[$ampctr]\n";
	}

	$pw_applybutton = $pw_bot->Button(-text=>'Apply', -command=>\&apply_config_plot)->pack(-side=>'right');
	$pw_cancelbutton = $pw_bot->Button(-text=>'Cancel', -command=>\&cancel_config_plot)->pack(-side=>'right');
}

sub apply_config_plot
{
	for($datactr = 0; $datactr <= $#data1; $datactr++)
	{
		$num = 0;
		$den = 0;

		for($ampctr = 0; $ampctr <= $#amps; $ampctr++)
		{
			if($frac_lab_num[$ampctr] == 1)
			{
				$num += $data1[$datactr]{ $amps[$ampctr] };
			}
			if($frac_lab_den[$ampctr] == 1)
			{
				$den += $data1[$datactr]{ $amps[$ampctr] };
			}
		}

		if($den != 0)
		{
			$data1[$datactr]{"frac_lab"} = $num / $den;
		}
		else
		{
			$data1[$datactr]{"frac_lab"} = -0.1
		}
	}

	$pw->destroy();

	plot_frac();
}

sub cancel_config_plot
{
	$pw->destroy();
}

sub choose_sample
{
	$sample_w = MainWindow->new();
	$sample_w->configure(-title=>'Sample Type', -background=>'black', -border=>'2');

	$sample_frame = $sample_w->Frame()->pack(-fill=>'both');

	$plot_30S = $sample_frame->Radiobutton(-text=>'30S ecoli', -variable=>\$plotframe, -value=>'0', -command=>\&draw_frame)->pack(-side=>'top', -anchor=>'w');
	$plot_50S = $sample_frame->Radiobutton(-text=>'50S ecoli', -variable=>\$plotframe, -value=>'1', -command=>\&draw_frame)->pack(-side=>'top', -anchor=>'w');
	$plot_70S = $sample_frame->Radiobutton(-text=>'70S ecoli', -variable=>\$plotframe, -value=>'2', -command=>\&draw_frame)->pack(-side=>'top', -anchor=>'w');
	$plot_30Schlamy = $sample_frame->Radiobutton(-text=>'30S chlamy', -variable=>\$plotframe, -value=>'3', -command=>\&draw_frame)->pack(-side=>'top', -anchor=>'w');
	$plot_40Shuman = $sample_frame->Radiobutton(-text=>'40S human', -variable=>\$plotframe, -value=>'4', -command=>\&draw_frame)->pack(-side=>'top', -anchor=>'w');
	$plot_60Shuman = $sample_frame->Radiobutton(-text=>'60S human', -variable=>\$plotframe, -value=>'5', -command=>\&draw_frame)->pack(-side=>'top', -anchor=>'w');
	$plot_80Shuman = $sample_frame->Radiobutton(-text=>'80S human', -variable=>\$plotframe, -value=>'13', -command=>\&draw_frame)->pack(-side=>'top', -anchor=>'w');
	$plot_cof = $sample_frame->Radiobutton(-text=>'cofactors', -variable=>\$plotframe, -value=>'6', -command=>\&draw_frame)->pack(-side=>'top', -anchor=>'w');
	$plot_40Syeast_old = $sample_frame->Radiobutton(-text=>'40S yeast (old)', -variable=>\$plotframe, -value=>'7', -command=>\&draw_frame)->pack(-side=>'top', -anchor=>'w');
	$plot_40Syeast = $sample_frame->Radiobutton(-text=>'40S yeast', -variable=>\$plotframe, -value=>'8', -command=>\&draw_frame)->pack(-side=>'top', -anchor=>'w');
	$plot_60Syeast = $sample_frame->Radiobutton(-text=>'60S yeast', -variable=>\$plotframe, -value=>'9', -command=>\&draw_frame)->pack(-side=>'top', -anchor=>'w');
	$plot_80Syeast = $sample_frame->Radiobutton(-text=>'80S yeast', -variable=>\$plotframe, -value=>'10', -command=>\&draw_frame)->pack(-side=>'top', -anchor=>'w');
	$plot_RNA = $sample_frame->Radiobutton(-text=>'RNA', -variable=>\$plotframe, -value=>'11', -command=>\&draw_frame)->pack(-side=>'top', -anchor=>'w');
	$plot_proteomic = $sample_frame->Radiobutton(-text=>'Proteomic', -variable=>\$plotframe, -value=>'12', -command=>\&draw_frame)->pack(-side=>'top', -anchor=>'w');

	$sample_donebutton = $sample_frame->Button(-text=>'Done', -command=>\&done_choose_sample)->pack(-side=>'top', -anchor=>'e');
}

sub done_choose_sample
{
	$sample_w->destroy();
}

sub display_behavior
{
	$dbw = MainWindow->new();
	$dbw->configure(-title=>'Plot Display Behavior', -background=>'black', -border=>'2');
	#$dbw->minsize( qw(400 400) );


	$dbw_frame = $dbw->Frame()->pack(-fill=>'both');
	$dbw_frame1 = $dbw_frame->Frame()->pack(-side=>'top', -fill=>'both');
		$dbw_obeydisplayflag = $dbw_frame1->Checkbutton(-text=>'Only Show Checked Entries', -variable=>\$obey_display_flag, -command=>\&plot_frac)->pack(-side=>'top');
		$dbw_savedisplayed = $dbw_frame1->Checkbutton(-text=>'Only Save Displayed Entries', -variable=>\$only_save_displayed)->pack(-side=>'top');

	$dbw_frame2 = $dbw_frame->Frame()->pack(-side=>'top', -fill=>'both');
		$dbw_showall = $dbw_frame2->Button(-text=>'Reset All Entries (Check)', -command=>\&show_all_plot)->pack(-side=>'right');
	$dbw_frame3 = $dbw_frame->Frame()->pack(-side=>'top', -fill=>'both');
		$dbw_hideall = $dbw_frame3->Button(-text=>'Reset All Entries (Uncheck)', -command=>\&hide_all_plot)->pack(-side=>'right');

	$dbw_frame4 = $dbw_frame->Frame(-height=>'15')->pack(-side=>'top');

	$dbw_frameX = $dbw_frame->Frame()->pack(-side=>'top', -fill=>'both');
		$dbw_close = $dbw_frameX->Button(-text=>'Close', -command=>\&close_dbw)->pack(-side=>'right');


}

sub close_dbw
{
	$dbw->destroy();
}

sub color_list1_bydisplayflag
{
	my $ctr;

	for($ctr = 0; $ctr <= $#data1; $ctr++)
	{
		if($data1[$ctr]{"display_flag"} == 0)
		{
			$list1->itemconfigure($ctr, -foreground=>'#888888');
		}
		elsif($data1[$ctr]{"display_flag"} == 1)
		{
			$list1->itemconfigure($ctr, -foreground=>'#000000');
		}
	}
}

sub process_selection1
{
	#print "1 $#data1\n";
	save_comment();
	#print "2 $#data1\n";
	foreground();
	#print "3 $#data1\n";

	# Restore the appropriate colors to $list1
	if(defined($lastsel1) && $lastsel1 >= 0)
	{
		$list1->itemconfigure($lastsel1, -foreground=>'#000000');
		$list1->itemconfigure($lastsel1, -selectforeground=>'#000000');
		$list1->itemconfigure($lastsel1, -background=>$bg1);
		$list1->itemconfigure($lastsel1, -selectbackground=>$bg1);
	}
	# Restore default plot colors and fill
	$canvasr->itemconfigure('pts', -outline=>'#000000', -width=>'1');

	#print "4 $#data1\n";
	clear234();
	clear234b();
	#print "5 $#data1\n";

	@selarray1 = $list1->curselection;
	$selindex1 = $selarray1[0];

	if(!defined($selindex1) || $selindex1 > $#view1)
	{
		return;
	}

	populate_list234();
	#print "6 $#data1\n";

	#plot_fit($selindex1, 1); # Seg Fault?!
	put_photo($selindex1, 1);
	#print "7 $#data1\n";

	# Color the plot point appropriately
	$isoid = $data1[$selindex1]{"isoid"};
	$tag = join '', "pt", $isoid;
	$canvasr->itemconfigure($tag, -outline=>'#FF0000', -width=>'3');
	$realcanvasr->raise($tag, 'all');

	$canvasr->delete($selectrect);
	@selectcoor = $canvasr->coords($tag);
	$selectcoor[0] -= 1;
	$selectcoor[1] -= 1;
	$selectcoor[2] += 2;
	$selectcoor[3] += 2;
	$selectrect = $canvasr->createRectangle(@selectcoor, -outline=>'#FF0000', -width=>'2');
	#print "8 $#data1\n";
	#####
	#####

	$lastsel1 = $selindex1;
	$last_isoid = $data1[$selindex1]{"isoid"};

	#color_list1_bydisplayflag();

	$list1->itemconfigure($selindex1, -foreground=>'#FF0000');
	$list1->itemconfigure($selindex1, -selectforeground=>'#FF0000');
	$list1->itemconfigure($selindex1, -background=>$bg);
	$list1->itemconfigure($selindex1, -selectbackground=>$bg);

	# Variable Display Stuff
	$n14mass_display = $data1[$selindex1]{"n14mass"}*$data1[$selindex1]{"charge"};
	$n15mass_display = $data1[$selindex1]{"n15mass"}*$data1[$selindex1]{"charge"};
	$ail14mass_display = $data1[$selindex1]{"mz_n14"}*$data1[$selindex1]{"charge"};
	$ail15mass_display = $data1[$selindex1]{"mz_n15"}*$data1[$selindex1]{"charge"};
	$sequence_display = $data1[$selindex1]{"seqmod"};
	$chisq_display = $data1[$selindex1]{"chisq"};
	$avg_wr_display = $data1[$selindex1]{"avg_wr"};
	$comment_display = $data1[$selindex1]{"comment"};
	$current_display_flag = $data1[$selindex1]{"display_flag"};
	$startres_display = $data1[$selindex1]{"startres"};
	$endres_display = $data1[$selindex1]{"endres"};
	$offset_display = sprintf("%.6f", $data1[$selindex1]{"off"});
	#print "9 $#data1\n";
}

sub process_selection2
{
	save_comment();
	foreground();

	clear234b();

	@selarray2 = $list2->curselection;
	$selindex2 = $selarray2[0];

	if(!defined($selindex2) || $selindex2 > $#view2)
	{
		return;
	}
	if($view2[$selindex2] eq "---")
	{
		return;
	}

	put_photo($selindex2, 2);

	$isoid = $data2[$selindex2]{"isoid"};
	$tag = join '', "pt", $isoid;
	$canvasr->delete($selectrect);
	@selectcoor = $canvasr->coords($tag);
	$selectcoor[0] -= 1;
	$selectcoor[1] -= 1;
	$selectcoor[2] += 2;
	$selectcoor[3] += 2;
	$selectrect = $canvasr->createRectangle(@selectcoor, -outline=>'#FF00FF', -width=>'2');

	populate_list234b($isoid);

	#####
	#####

	$lastsel2 = $selindex2;
	$last_isoid = $data2[$selindex2]{"isoid"};

	$list2->itemconfigure($selindex2, -foreground=>'#FF0000');
	$list2->itemconfigure($selindex2, -selectforeground=>'#FF0000');

	# Variable Display Stuff
	$n14mass_display = $data2[$selindex2]{"n14mass"}*$data2[$selindex2]{"charge"};
	$n15mass_display = $data2[$selindex2]{"n15mass"}*$data2[$selindex2]{"charge"};
	$ail14mass_display = $data2[$selindex2]{"mz_n14"}*$data2[$selindex2]{"charge"};
	$ail15mass_display = $data2[$selindex2]{"mz_n15"}*$data2[$selindex2]{"charge"};
	$sequence_display = $data2[$selindex2]{"seqmod"};
	$chisq_display = $data2[$selindex2]{"chisq"};
	$avg_wr_display = $data2[$selindex2]{"avg_wr"};
	$comment_display = $data2[$selindex2]{"comment"};
	$current_display_flag = $data2[$selindex2]{"display_flag"};
	$startres_display = $data2[$selindex2]{"startres"};
	$endres_display = $data2[$selindex2]{"endres"};
	$offset_display = sprintf("%.6f", $data2[$selindex2]{"off"});
}

sub process_selection3
{
	save_comment();
	foreground();

	clear234b();

	@selarray3 = $list3->curselection;
	$selindex3 = $selarray3[0];

	#print "selindex3: $selindex3\n";

	if(!defined($selindex3) || $selindex3 > $#view3)
	{
		return;
	}
	if($view3[$selindex3] eq "---")
	{
		return;
	}

	put_photo($selindex3, 3);

	$isoid = $data3[$selindex3]{"isoid"};
	$tag = join '', "pt", $isoid;
	$canvasr->delete($selectrect);
	@selectcoor = $canvasr->coords($tag);
	$selectcoor[0] -= 1;
	$selectcoor[1] -= 1;
	$selectcoor[2] += 2;
	$selectcoor[3] += 2;
	$selectrect = $canvasr->createRectangle(@selectcoor, -outline=>'#00FF00', -width=>'2');

	populate_list234b($isoid);

	#####
	#####

	$lastsel3 = $selindex3;
	$last_isoid = $data3[$selindex3]{"isoid"};

	$list3->itemconfigure($selindex3, -foreground=>'#FF0000');
	$list3->itemconfigure($selindex3, -selectforeground=>'#FF0000');

	# Variable Display Stuff
	$n14mass_display = $data3[$selindex3]{"n14mass"}*$data3[$selindex3]{"charge"};
	$n15mass_display = $data3[$selindex3]{"n15mass"}*$data3[$selindex3]{"charge"};
	$ail14mass_display = $data3[$selindex3]{"mz_n14"}*$data3[$selindex3]{"charge"};
	$ail15mass_display = $data3[$selindex3]{"mz_n15"}*$data3[$selindex3]{"charge"};
	$sequence_display = $data3[$selindex3]{"seqmod"};
	$chisq_display = $data3[$selindex3]{"chisq"};
	$avg_wr_display = $data3[$selindex3]{"avg_wr"};
	$comment_display = $data3[$selindex3]{"comment"};
	$current_display_flag = $data3[$selindex3]{"display_flag"};
	$startres_display = $data3[$selindex3]{"startres"};
	$endres_display = $data3[$selindex3]{"endres"};
	$offset_display = sprintf("%.6f", $data3[$selindex3]{"off"});
}

sub process_selection4
{
	save_comment();
	foreground();

	clear234b();

	@selarray4 = $list4->curselection;
	$selindex4 = $selarray4[0];

	if(!defined($selindex4) || $selindex4 > $#view4)
	{
		return;
	}
	if($view4[$selindex4] eq "---")
	{
		return;
	}

	put_photo($selindex4, 4);

	$isoid = $data4[$selindex4]{"isoid"};
	$tag = join '', "pt", $isoid;
	$canvasr->delete($selectrect);
	@selectcoor = $canvasr->coords($tag);
	$selectcoor[0] -= 1;
	$selectcoor[1] -= 1;
	$selectcoor[2] += 2;
	$selectcoor[3] += 2;
	$selectrect = $canvasr->createRectangle(@selectcoor, -outline=>'#0000FF', -width=>'2');

	populate_list234b($isoid);

	#####
	#####

	$lastsel4 = $selindex4;
	$last_isoid = $data4[$selindex4]{"isoid"};

	$list4->itemconfigure($selindex4, -foreground=>'#FF0000');
	$list4->itemconfigure($selindex4, -selectforeground=>'#FF0000');

	# Variable Display Stuff
	$n14mass_display = $data4[$selindex4]{"n14mass"}*$data4[$selindex4]{"charge"};
	$n15mass_display = $data4[$selindex4]{"n15mass"}*$data4[$selindex4]{"charge"};
	$ail14mass_display = $data4[$selindex4]{"mz_n14"}*$data4[$selindex4]{"charge"};
	$ail15mass_display = $data4[$selindex4]{"mz_n15"}*$data4[$selindex4]{"charge"};
	$sequence_display = $data4[$selindex4]{"seqmod"};
	$chisq_display = $data4[$selindex4]{"chisq"};
	$avg_wr_display = $data4[$selindex4]{"avg_wr"};
	$comment_display = $data4[$selindex4]{"comment"};
	$current_display_flag = $data4[$selindex4]{"display_flag"};
	$startres_display = $data4[$selindex4]{"startres"};
	$endres_display = $data4[$selindex4]{"endres"};
	$offset_display = sprintf("%.6f", $data4[$selindex4]{"off"});
}

sub process_selection2b
{
	save_comment();
	foreground();

	@selarray2b = $list2b->curselection;
	$selindex2b = $selarray2b[0];

	if(!defined($selindex2b) || $selindex2b > $#view2b)
	{
		return;
	}
	if($view2b[$selindex2b] eq "---")
	{
		return;
	}

	put_photo($selindex2b, 5);

	#####
	#####

	$lastsel2b = $selindex2b;
	$last_isoid = $data2b[$selindex2b]{"isoid"};

	$list2b->itemconfigure($selindex2b, -foreground=>'#FF0000');
	$list2b->itemconfigure($selindex2b, -selectforeground=>'#FF0000');

	# Variable Display Stuff
	$n14mass_display = $data2b[$selindex2b]{"n14mass"}*$data2b[$selindex2b]{"charge"};
	$n15mass_display = $data2b[$selindex2b]{"n15mass"}*$data2b[$selindex2b]{"charge"};
	$ail14mass_display = $data2b[$selindex2b]{"mz_n14"}*$data2b[$selindex2b]{"charge"};
	$ail15mass_display = $data2b[$selindex2b]{"mz_n15"}*$data2b[$selindex2b]{"charge"};
	$sequence_display = $data2b[$selindex2b]{"seqmod"};
	$chisq_display = $data2b[$selindex2b]{"chisq"};
	$avg_wr_display = $data2b[$selindex2b]{"avg_wr"};
	$comment_display = $data2b[$selindex2b]{"comment"};
	$current_display_flag = $data2b[$selindex2b]{"display_flag"};
	$startres_display = $data2b[$selindex2b]{"startres"};
	$endres_display = $data2b[$selindex2b]{"endres"};
	$offset_display = sprintf("%.6f", $data2b[$selindex2b]{"off"});
}

sub process_selection3b
{
	save_comment();
	foreground();

	@selarray3b = $list3b->curselection;
	$selindex3b = $selarray3b[0];

	if(!defined($selindex3b) || $selindex3b > $#view3b)
	{
		return;
	}
	if($view3b[$selindex3b] eq "---")
	{
		return;
	}

	put_photo($selindex3b, 6);

	#####
	#####

	$lastsel3b = $selindex3b;
	$last_isoid = $data3b[$selindex3b]{"isoid"};

	$list3b->itemconfigure($selindex3b, -foreground=>'#FF0000');
	$list3b->itemconfigure($selindex3b, -selectforeground=>'#FF0000');

	# Variable Display Stuff
	$n14mass_display = $data3b[$selindex3b]{"n14mass"}*$data3b[$selindex3b]{"charge"};
	$n15mass_display = $data3b[$selindex3b]{"n15mass"}*$data3b[$selindex3b]{"charge"};
	$ail14mass_display = $data3b[$selindex3b]{"mz_n14"}*$data3b[$selindex3b]{"charge"};
	$ail15mass_display = $data3b[$selindex3b]{"mz_n15"}*$data3b[$selindex3b]{"charge"};
	$sequence_display = $data3b[$selindex3b]{"seqmod"};
	$chisq_display = $data3b[$selindex3b]{"chisq"};
	$avg_wr_display = $data3b[$selindex3b]{"avg_wr"};
	$comment_display = $data3b[$selindex3b]{"comment"};
	$current_display_flag = $data3b[$selindex3b]{"display_flag"};
	$startres_display = $data3b[$selindex3b]{"startres"};
	$endres_display = $data3b[$selindex3b]{"endres"};
	$offset_display = sprintf("%.6f", $data3b[$selindex3b]{"off"});
}

sub process_selection4b
{
	save_comment();
	foreground();

	@selarray4b = $list4b->curselection;
	$selindex4b = $selarray4b[0];

	if(!defined($selindex4b) || $selindex4b > $#view4b)
	{
		return;
	}
	if($view4b[$selindex4b] eq "---")
	{
		return;
	}

	put_photo($selindex4b, 7);

	#####
	#####

	$lastsel4b = $selindex4b;
	$last_isoid = $data4b[$selindex4b]{"isoid"};

	$list4b->itemconfigure($selindex4b, -foreground=>'#FF0000');
	$list4b->itemconfigure($selindex4b, -selectforeground=>'#FF0000');

	# Variable Display Stuff
	$n14mass_display = $data4b[$selindex4b]{"n14mass"}*$data4b[$selindex4b]{"charge"};
	$n15mass_display = $data4b[$selindex4b]{"n15mass"}*$data4b[$selindex4b]{"charge"};
	$ail14mass_display = $data4b[$selindex4b]{"mz_n14"}*$data4b[$selindex4b]{"charge"};
	$ail15mass_display = $data4b[$selindex4b]{"mz_n15"}*$data4b[$selindex4b]{"charge"};
	$sequence_display = $data4b[$selindex4b]{"seqmod"};
	$chisq_display = $data4b[$selindex4b]{"chisq"};
	$avg_wr_display = $data4b[$selindex4b]{"avg_wr"};
	$comment_display = $data4b[$selindex4b]{"comment"};
	$current_display_flag = $data4b[$selindex4b]{"display_flag"};
	$startres_display = $data4b[$selindex4b]{"startres"};
	$endres_display = $data4b[$selindex4b]{"endres"};
	$offset_display = sprintf("%.6f", $data4b[$selindex4b]{"off"});
}

sub save_comment
{
	$comment_display =~ s/\,//g;
	my $comment_id = $last_isoid;

	if($comment_id == -1)
	{
		return;
	}

	if(!defined($iso2idx{$comment_id}))
	{
		#print "Deletion\n";
		return;
	}

	$comment_idx = $iso2idx{$comment_id};

	# Apply comment to the main array
	$data1[$comment_idx]{"comment"} = $comment_display;

	# Apply comment to the secondary/tertiary arrays as needed
	for(my $ctr = 0; $ctr <= $#data2; $ctr++)
	{
		if($data2[$ctr]{"isoid"} == $comment_id)
		{
			$data2[$ctr]{"comment"} = $comment_display;
		}
	}
	for(my $ctr = 0; $ctr <= $#data3; $ctr++)
	{
		if($data3[$ctr]{"isoid"} == $comment_id)
		{
			$data3[$ctr]{"comment"} = $comment_display;
		}
	}
	for(my $ctr = 0; $ctr <= $#data4; $ctr++)
	{
		if($data4[$ctr]{"isoid"} == $comment_id)
		{
			$data4[$ctr]{"comment"} = $comment_display;
		}
	}
	for(my $ctr = 0; $ctr <= $#data2b; $ctr++)
	{
		if($data2b[$ctr]{"isoid"} == $comment_id)
		{
			$data2b[$ctr]{"comment"} = $comment_display;
		}
	}
	for(my $ctr = 0; $ctr <= $#data3b; $ctr++)
	{
		if($data3b[$ctr]{"isoid"} == $comment_id)
		{
			$data3b[$ctr]{"comment"} = $comment_display;
		}
	}
	for(my $ctr = 0; $ctr <= $#data4b; $ctr++)
	{
		if($data4b[$ctr]{"isoid"} == $comment_id)
		{
			$data4b[$ctr]{"comment"} = $comment_display;
		}
	}
}

sub set_display_flag
{
	my $display_flag_id = $last_isoid;

	if($display_flag_id == -1)
	{
		return;
	}

	if(!defined($iso2idx{$display_flag_id}))
	{
		#print "Deletion\n";
		return;
	}

	$display_flag_idx = $iso2idx{$display_flag_id};

	# Apply comment to the main array
	$data1[$display_flag_idx]{"display_flag"} = $current_display_flag;

	# Apply comment to the secondary/tertiary arrays as needed
	for(my $ctr = 0; $ctr <= $#data2; $ctr++)
	{
		if($data2[$ctr]{"isoid"} == $display_flag_id)
		{
			$data2[$ctr]{"display_flag"} = $current_display_flag;
		}
	}
	for(my $ctr = 0; $ctr <= $#data3; $ctr++)
	{
		if($data3[$ctr]{"isoid"} == $display_flag_id)
		{
			$data3[$ctr]{"display_flag"} = $current_display_flag;
		}
	}
	for(my $ctr = 0; $ctr <= $#data4; $ctr++)
	{
		if($data4[$ctr]{"isoid"} == $display_flag_id)
		{
			$data4[$ctr]{"display_flag"} = $current_display_flag;
		}
	}
	for(my $ctr = 0; $ctr <= $#data2b; $ctr++)
	{
		if($data2b[$ctr]{"isoid"} == $display_flag_id)
		{
			$data2b[$ctr]{"display_flag"} = $current_display_flag;
		}
	}
	for(my $ctr = 0; $ctr <= $#data3b; $ctr++)
	{
		if($data3b[$ctr]{"isoid"} == $display_flag_id)
		{
			$data3b[$ctr]{"display_flag"} = $current_display_flag;
		}
	}
	for(my $ctr = 0; $ctr <= $#data4b; $ctr++)
	{
		if($data4b[$ctr]{"isoid"} == $display_flag_id)
		{
			$data4b[$ctr]{"display_flag"} = $current_display_flag;
		}
	}

	plot_frac();

	#color_list1_bydisplayflag();
}

sub delete_selection
{
	isoid_select();
	if($isoid_curselection == -1)
	{
		print "Nothing to delete...\n";
		return;
	}

	my $delindex = $iso2idx{$isoid_curselection};

	if($delindex > $#data1)
	{
		print "ERROR: Deletion outside of range\n";
	}

	#print "$delindex: $#data1 -> ";
	# Delete the entry from @data1
	print "$#data1 $#deleted\n";
	push(@deleted, splice(@data1, $delindex, 1));
	#splice(@data1, $delindex, 1);
	print "$#data1 $#deleted---\n";
	#print "$#data1\n";

	# Regenerate @view1 and %iso2idx
	regen_view1();
	#print "A $#data1\n";

	# Delete the item highlighted in list1
	# Selection index is unchanged, as the next item increments to the current selection index
	if($delindex == $lastsel1)
	{
		$newselindex = $lastsel1;
	}
	# Delete an item which occurs earlier in list1
	# Must decrement the selection so that we still select the same item in list1
	elsif($delindex < $lastsel1)
	{
		$newselindex = $lastsel1 - 1;
	}
	# Delete an item which occurs later in list1
	# Selection index is unchanged
	elsif($delindex > $lastsel1)
	{
		$newselindex = $lastsel1;
	}

	# Check bounds of selection
	if($newselindex < 0)
	{
		$newselindex = 0;
	}
	elsif($newselindex > $#view1)
	{
		$newselindex = $#view1;
	}

	$lastsel1 = $newselindex; # Set this or foreground() might complain

	# Regenerate the plot
	plot_frac();
	#print "B $#data1\n";

	# Reprocess the selection
	$list1->selectionClear(0, $#view1);
	$list1->selectionSet($newselindex);
	$list1->see($newselindex);
	#print "B2 $#data1\n";
	process_selection1();

	#print "C $#data1\n";
}

# This is just like delete_selection() (above)...
# ...but we tack ".del" onto the end of the relevant files (.txt, .png and .contour.png)
sub purge_selection
{
	isoid_select();
	if($isoid_curselection == -1)
	{
		print "Nothing to delete...\n";
		return;
	}

	my $delindex = $iso2idx{$isoid_curselection};

	if($delindex > $#data1)
	{
		print "ERROR: Deletion outside of range\n";
	}

	#print "$delindex: $#data1 -> ";
	# Delete the entry from @data1
	print "$#data1 $#deleted\n";
	push(@deleted, splice(@data1, $delindex, 1));
	#splice(@data1, $delindex, 1);
	print "$#data1 $#deleted---\n";
	#print "$#data1\n";

	#####
	# Additional Purge Section
	#####
	$purge_txt = $deleted[$#deleted]{"file"}; # Comes from peaks directory (msc04)
	$purge_cnt = $purge_txt;
	substr($purge_cnt, -3, 3) = "contour.png";
	$purge_base = $deleted[$#deleted]{"isofile"}; # Possible alt fits dir (msc05)
	$purge_png = join '', $purge_base, ".fit.png";

	$del_txt = join '', $purge_txt, ".del";
	$del_cnt = join '', $purge_cnt, ".del";
	$del_png = join '', $purge_png, ".del";

	`mv $purge_txt $del_txt`;
	`mv $purge_cnt $del_cnt`;
	`mv $purge_png $del_png`;
	#####
	# /Additonal Purge Section
	#####

	# Regenerate @view1 and %iso2idx
	regen_view1();
	#print "A $#data1\n";

	# Delete the item highlighted in list1
	# Selection index is unchanged, as the next item increments to the current selection index
	if($delindex == $lastsel1)
	{
		$newselindex = $lastsel1;
	}
	# Delete an item which occurs earlier in list1
	# Must decrement the selection so that we still select the same item in list1
	elsif($delindex < $lastsel1)
	{
		$newselindex = $lastsel1 - 1;
	}
	# Delete an item which occurs later in list1
	# Selection index is unchanged
	elsif($delindex > $lastsel1)
	{
		$newselindex = $lastsel1;
	}

	# Check bounds of selection
	if($newselindex < 0)
	{
		$newselindex = 0;
	}
	elsif($newselindex > $#view1)
	{
		$newselindex = $#view1;
	}

	$lastsel1 = $newselindex; # Set this or foreground() might complain

	# Regenerate the plot
	plot_frac();
	#print "B $#data1\n";

	# Reprocess the selection
	$list1->selectionClear(0, $#view1);
	$list1->selectionSet($newselindex);
	$list1->see($newselindex);
	#print "B2 $#data1\n";
	process_selection1();

	#print "C $#data1\n";
}

sub undo_deletion
{
	if($#deleted < 0)
	{
		print "Error: No deletions to undo\n";
		return;
	}

	print "$#data1 $#deleted\n";
	#push(@data1, shift(@deleted)); # Oops, off wrong end!
	push(@data1, pop(@deleted));
	print "$#data1 $#deleted\n---\n";

	# Now we need to check for .del files, and restore if necessary
	$purge_txt = $data1[$#data1]{"file"}; # Comes from peaks directory (msc04)
	$purge_cnt = $purge_txt;
	substr($purge_cnt, -3, 3) = "contour.png";
	$purge_base = $data1[$#data1]{"isofile"}; # Possible alt fits dir (msc05)
	$purge_png = join '', $purge_base, ".fit.png";

	$del_txt = join '', $purge_txt, ".del";
	$del_cnt = join '', $purge_cnt, ".del";
	$del_png = join '', $purge_png, ".del";

	if(-e $del_txt)
	{
		`mv $del_txt $purge_txt`;
	}
	if(-e $del_cnt)
	{
		`mv $del_cnt $purge_cnt`;
	}
	if(-e $del_png)
	{
		`mv $del_png $purge_png`;
	}

	sort_list1();
	#plot_frac(); # Not needed - included when re-sorting
	process_selection1();
}

sub foreground
{
	if(defined($lastsel1) && $lastsel1 >= 0)
	{
		$list1->itemconfigure($lastsel1, -foreground=>'#000000');
		$list1->itemconfigure($lastsel1, -selectforeground=>'#000000');
	}
	for(my $ctr = 0; $ctr <= $#view2; $ctr++)
	{
		$list2->itemconfigure($ctr, -foreground=>'#000000');
		$list2->itemconfigure($ctr, -selectforeground=>'#000000');
	}
	for(my $ctr = 0; $ctr <= $#view3; $ctr++)
	{
		$list3->itemconfigure($ctr, -foreground=>'#000000');
		$list3->itemconfigure($ctr, -selectforeground=>'#000000');
	}
	for(my $ctr = 0; $ctr <= $#view4; $ctr++)
	{
		$list4->itemconfigure($ctr, -foreground=>'#000000');
		$list4->itemconfigure($ctr, -selectforeground=>'#000000');
	}
	for(my $ctr = 0; $ctr <= $#view2b; $ctr++)
	{
		$list2b->itemconfigure($ctr, -foreground=>'#000000');
		$list2b->itemconfigure($ctr, -selectforeground=>'#000000');
	}
	for(my $ctr = 0; $ctr <= $#view3b; $ctr++)
	{
		$list3b->itemconfigure($ctr, -foreground=>'#000000');
		$list3b->itemconfigure($ctr, -selectforeground=>'#000000');
	}
	for(my $ctr = 0; $ctr <= $#view4b; $ctr++)
	{
		$list4b->itemconfigure($ctr, -foreground=>'#000000');
		$list4b->itemconfigure($ctr, -selectforeground=>'#000000');
	}
}

sub clear1
{
	for(my $ctr = 0; $ctr <= $#view1; $ctr++)
	{
		$list1->itemconfigure($ctr, -background=>$bg1);
		$list1->itemconfigure($ctr, -selectbackground=>$bg1);
		$list1->itemconfigure($ctr, -foreground=>'#000000');
		$list1->itemconfigure($ctr, -selectforeground=>'#000000');
	}

	@view1 = ("---");
	@data1 = ("---");

	$last_isoid = -1;
}

sub clear234
{
	for(my $ctr = 0; $ctr <= $#view2; $ctr++)
	{
		$list2->itemconfigure($ctr, -background=>$bg2);
		$list2->itemconfigure($ctr, -selectbackground=>$bg2);
		$list2->itemconfigure($ctr, -foreground=>'#000000');
		$list2->itemconfigure($ctr, -selectforeground=>'#000000');
	}
	for(my $ctr = 0; $ctr <= $#view3; $ctr++)
	{
		$list3->itemconfigure($ctr, -background=>$bg3);
		$list3->itemconfigure($ctr, -selectbackground=>$bg3);
		$list3->itemconfigure($ctr, -foreground=>'#000000');
		$list3->itemconfigure($ctr, -selectforeground=>'#000000');
	}
	for(my $ctr = 0; $ctr <= $#view4; $ctr++)
	{
		$list4->itemconfigure($ctr, -background=>$bg4);
		$list4->itemconfigure($ctr, -selectbackground=>$bg4);
		$list4->itemconfigure($ctr, -foreground=>'#000000');
		$list4->itemconfigure($ctr, -selectforeground=>'#000000');
	}

	@view2 = ("---");
	@data2 = ("---");
	@view3 = ("---");
	@data3 = ("---");
	@view4 = ("---");
	@data4 = ("---");
}

sub clear234b
{
	for(my $ctr = 0; $ctr <= $#view2b; $ctr++)
	{
		$list2b->itemconfigure($ctr, -background=>$bg2);
		$list2b->itemconfigure($ctr, -selectbackground=>$bg2);
		$list2b->itemconfigure($ctr, -foreground=>'#000000');
		$list2b->itemconfigure($ctr, -selectforeground=>'#000000');
	}
	for(my $ctr = 0; $ctr <= $#view3b; $ctr++)
	{
		$list3b->itemconfigure($ctr, -background=>$bg3);
		$list3b->itemconfigure($ctr, -selectbackground=>$bg3);
		$list3b->itemconfigure($ctr, -foreground=>'#000000');
		$list3b->itemconfigure($ctr, -selectforeground=>'#000000');
	}
	for(my $ctr = 0; $ctr <= $#view4b; $ctr++)
	{
		$list4b->itemconfigure($ctr, -background=>$bg4);
		$list4b->itemconfigure($ctr, -selectbackground=>$bg4);
		$list4b->itemconfigure($ctr, -foreground=>'#000000');
		$list4b->itemconfigure($ctr, -selectforeground=>'#000000');
	}

	@view2b = ("---");
	@data2b = ("---");
	@view3b = ("---");
	@data3b = ("---");
	@view4b = ("---");
	@data4b = ("---");
}

sub regen_view1
{
	@view1 = ();
	%iso2idx = ();
	undef(%iso2idx);
	for(my $ctr = 0; $ctr <= $#data1; $ctr++)
	{
		my $temp = $data1[$ctr]{"isofile"};
		my @tarray = split /\//, $temp;
		my $isoid = $data1[$ctr]{"isoid"};
		$view1[$ctr] = $tarray[$#tarray];
		$iso2idx{$isoid} = $ctr;
	}
	$list1->update;
}

sub regen_view234
{
	my @params = @_;
	$highlight2 = $params[0];
	$highlight3 = $params[1];
	$highlight4 = $params[2];

	@view2 = ();
	for(my $ctr = 0; $ctr <= $#data2; $ctr++)
	{
		my $temp = $data2[$ctr]{"isofile"};
		my @tarray = split /\//, $temp;
		$view2[$ctr] = $tarray[$#tarray];
	}
	$list2->itemconfigure($highlight2, -background=>$bg);
	$list2->itemconfigure($highlight2, -selectbackground=>$bg);

	@view3 = ();
	for(my $ctr = 0; $ctr <= $#data3; $ctr++)
	{
		my $temp = $data3[$ctr]{"isofile"};
		my @tarray = split /\//, $temp;
		$view3[$ctr] = $tarray[$#tarray];
	}
	$list3->itemconfigure($highlight3, -background=>$bg);
	$list3->itemconfigure($highlight3, -selectbackground=>$bg);

	@view4 = ();
	for(my $ctr = 0; $ctr <= $#data4; $ctr++)
	{
		my $temp = $data4[$ctr]{"isofile"};
		my @tarray = split /\//, $temp;
		$view4[$ctr] = $tarray[$#tarray];
	}
	$list4->itemconfigure($highlight4, -background=>$bg);
	$list4->itemconfigure($highlight4, -selectbackground=>$bg);
}

sub regen_view234b
{
	my @params = @_;
	$highlight2b = $params[0];
	$highlight3b = $params[1];
	$highlight4b = $params[2];

	@view2b = ();
	for(my $ctr = 0; $ctr <= $#data2b; $ctr++)
	{
		my $temp = $data2b[$ctr]{"isofile"};
		my @tarray = split /\//, $temp;
		$view2b[$ctr] = $tarray[$#tarray];
	}
	$list2b->itemconfigure($highlight2b, -background=>$bg);
	$list2b->itemconfigure($highlight2b, -selectbackground=>$bg);

	@view3b = ();
	for(my $ctr = 0; $ctr <= $#data3b; $ctr++)
	{
		my $temp = $data3b[$ctr]{"isofile"};
		my @tarray = split /\//, $temp;
		$view3b[$ctr] = $tarray[$#tarray];
	}
	$list3b->itemconfigure($highlight3b, -background=>$bg);
	$list3b->itemconfigure($highlight3b, -selectbackground=>$bg);

	@view4b = ();
	for(my $ctr = 0; $ctr <= $#data4b; $ctr++)
	{
		my $temp = $data4b[$ctr]{"isofile"};
		my @tarray = split /\//, $temp;
		$view4b[$ctr] = $tarray[$#tarray];
	}
	$list4b->itemconfigure($highlight4b, -background=>$bg);
	$list4b->itemconfigure($highlight4b, -selectbackground=>$bg);
}

sub sort_list1
{
	#print "$#data1 -> ";
	if($sortstyle == 1)
	{
		print "Sorting by RT\n";
		@data1 = sort byrt @data1;
	}
	elsif($sortstyle == 2)
	{
		print "Sorting by M/Z\n";
		@data1 = sort bymz @data1;
	}
	elsif($sortstyle == 3)
	{
		print "Sorting by Protein\n";
		@data1 = sort byprotein @data1;
	}
	elsif($sortstyle == 4)
	{
		print "Sort by Protein then Chi Squared\n";
		@data1 = sort byprochi @data1;
	}
	#print "$#data1 -> ";
	plot_frac(); # Needed to regenerate the coordinates associated with points
	regen_view1();
	#print "$#data1\n";
}

sub byrt
{
	my $rta = $$a{"rt_n14"} + $$a{"rt_n15"};
	my $rtb = $$b{"rt_n14"} + $$b{"rt_n15"};
	$rta <=> $rtb;
}

sub bymz
{
	$$a{"n14mass"} <=> $$b{"n14mass"};
}

sub byprotein
{
	$pro2loc{ $$a{"protein"} } <=> $pro2loc{ $$b{"protein"} }
		or
	$$a{"startres"} <=> $$b{"startres"}
		or
	$$a{"endres"} <=> $$b{"endres"}
		or
	$$a{"n14mass"} <=> $$b{"n14mass"};
}

sub byprochi
{
	$pro2loc{ $$a{"protein"} } <=> $pro2loc{ $$b{"protein"} }
		or
	$$a{"chisq"} <=> $$b{"chisq"};
}


sub quit
{
	# Joan has cat'd together multiple datasets
	# This means there are multiple "isofile" and "file" entries
	# So need to get all of them, and restore for each

	#%isofile_hash = ();
	#%file_hash = ();

	#for(my $ctr1 = 0; $ctr1 <= $#data1; $ctr1++)
	#{
	#	@dirarray = split /\//, $data1[$ctr1]{"isofile"};
	#	@dirarray2 = split /\//, $data1[$ctr1]{"file"};
	#
	#	$isofile_dir = $dirarray[0];
	#	$file_dir = $dirarray2[0];
	#
	#	if(!defined($isofile_hash{$isofile_dir}))
	#	{
	#		$isofile_hash{$isofile_dir} = 1;
	#	}
	#	if(!defined($file_hash{$file_dir}))
	#	{
	#		$file_hash{$file_dir} = 1;
	#	}
	#}

	foreach $key (keys %isofile_hash)
	{
		# If .msvbkp files exist, these need to be replaced as originals, since this indicates changes were not saved
		#@dirarray = split /\//, $data1[0]{"isofile"}; # fits directory
		@dirarray = split /\//, $key;

		@msvbkp_ls = `ls $dirarray[0]/*.msvbkp`;

		if($#msvbkp_ls >= 0)
		{
			chomp(@msvbkp_ls);
			for(my $ctr = 0; $ctr <= $#msvbkp_ls; $ctr++)
			{
				$origfile = $msvbkp_ls[$ctr];
				substr($origfile, -7, 7) = "";
				print "Restoring $origfile\n";
				`cp $msvbkp_ls[$ctr] $origfile`;
			}
			# After replacement, delete .msvbkp files
			print "Removing $dirarray[0]/*.msvbkp\n";
			`rm $dirarray[0]/*.msvbkp`;
		}

		# If .del files exist, these also need to be replaced - originally marked for deletion, but changes were not saved
		@del_ls = `ls $dirarray[0]/*.del`;

		if($#del_ls >= 0)
		{
			chomp(@del_ls);
			for(my $ctr = 0; $ctr <= $#del_ls; $ctr++)
			{
				$origfile = $del_ls[$ctr];
				substr($origfile, -4, 4) = "";
				print "Restoring $origfile\n";
				`cp $del_ls[$ctr] $origfile`;
			}
			print "Removing $dirarray[0]/*.del\n";
			`rm $dirarray[0]/*.del`;
		}
	}

	foreach $key (keys %file_hash)
	{
		#@dirarray2 = split /\//, $data1[0]{"file"}; # txt directory (often same)
		@dirarray2 = split /\//, $key;

		# Repeat for second directory
		@del_ls = `ls $dirarray2[0]/*.del`;

		if($#del_ls >= 0)
		{
			chomp(@del_ls);
			for(my $ctr = 0; $ctr <= $#del_ls; $ctr++)
			{
				$origfile = $del_ls[$ctr];
				substr($origfile, -4, 4) = "";
				print "Restoring $origfile\n";
				`cp $del_ls[$ctr] $origfile`;
			}
			print "Removing $dirarray2[0]/*.del\n";
			`rm $dirarray2[0]/*.del`;
		}
	}

	print "Goodbye\n";
	exit(0);
}

sub populate_list234
{
	$idflag2a = "digestionid";
	$idflag2b = "digestionid";

	$idflag3a = "digestpepid";
	$idflag3b = "digestpepid";

	$idflag4a = "ailid_n14";
	$idflag4b = "ailid_n15";

	if($matchmode == 0)
	{
		# Do nothing, default of N14/N15 pairs
	}
	elsif($matchmode == 1)
	{
		# N14 only
		$idflag4b = "ailid_n14";
	}
	elsif($matchmode == 2)
	{
		# N15 only
		$idflag4a = "ailid_n15";
	}

	$id = $data1[$selindex1]{"isoid"};

	$id2a = $data1[$selindex1]{$idflag2a};
	$id2b = $data1[$selindex1]{$idflag2b};
	$id3a = $data1[$selindex1]{$idflag3a};
	$id3b = $data1[$selindex1]{$idflag3b};
	$id4a = $data1[$selindex1]{$idflag4a};
	$id4b = $data1[$selindex1]{$idflag4b};

	# Detect MS/MS IDs (digestionid & digestpepid = -1)
	if($id2a == -1 && $id2b == -1)
	{
		$data2[0] = $data1[$selindex1];
		$data3[0] = $data1[$selindex1];
		#$data4[0] = $data1[$selindex1];
		$bold2 = 0;
		$bold3 = 0;
		$bold4 = 0;
	}

	for(my $ctr = 0; $ctr <= $#data1; $ctr++)
	{
		$idflag = 0;
		$test2a = $data1[$ctr]{$idflag2a}; # digestionid
		$test2b = $data1[$ctr]{$idflag2b};
		$test3a = $data1[$ctr]{$idflag3a}; # digestpepid
		$test3b = $data1[$ctr]{$idflag3b};
		$test4a = $data1[$ctr]{$idflag4a}; # ailid_n14
		$test4b = $data1[$ctr]{$idflag4b}; # ailid_n15

		$testid = $data1[$ctr]{"isoid"};
		if($testid == $id)
		{
			$idflag = 1;
		}

		if($test2a >= 0 && $test2b >= 0) # MS/MS IDs generate digestionid/digestpepid = -1
		{
			if($test2a == $id2a || $test2b == $id2b)
			{
				if($data2[0] eq "---")
				{
					$data2[0] = $data1[$ctr];
				}
				else
				{
					push(@data2, $data1[$ctr]);
				}
				if($idflag == 1)
				{
					$bold2 = $#data2;
				}
			}
		}
		# Here we want to include other charges of the same peptide
		# But we want to ignore other versions of same charge (which go in @data2)
		if($test3a >= 0 && $test3b >= 0)
		{
			if($test3a == $id3a || $test3b == $id3b)
			{
				# If it's not the identical isoid, but it is the same ion, skip
				if($idflag != 1 && ( $test2a == $id2a || $test2b == $id2b ) )
				{
					# Do Nothing
				}
				else
				{
					if($data3[0] eq "---")
					{
						$data3[0] = $data1[$ctr];
					}
					else
					{
						push(@data3, $data1[$ctr]);
					}
					if($idflag == 1)
					{
						$bold3 = $#data3;
					}
				}
			}
		}
		if($test4a == $id4a || $test4b == $id4b || $test4a == $id4b || $test4b == $id4a)
		{
			if($data4[0] eq "---")
			{
				$data4[0] = $data1[$ctr];
			}
			else
			{
				push(@data4, $data1[$ctr]);
			}
			if($idflag == 1)
			{
				$bold4 = $#data4;
			}
		}
	}

	color_pts234();
	regen_view234($bold2, $bold3, $bold4);
}

sub populate_list234b
{
	my @params = @_;
	my $id = $params[0];

	$idflag2a = "digestionid";
	$idflag2b = "digestionid";

	$idflag3a = "digestpepid";
	$idflag3b = "digestpepid";

	$idflag4a = "ailid_n14";
	$idflag4b = "ailid_n15";

	if($matchmode == 0)
	{
		# Do nothing, default of N14/N15 pairs
	}
	elsif($matchmode == 1)
	{
		# N14 only
		$idflag4b = "ailid_n14";
	}
	elsif($matchmode == 2)
	{
		# N15 only
		$idflag4a = "ailid_n15";
	}

	for(my $ctr = 0; $ctr <= $#data1; $ctr++)
	{
		$testid = $data1[$ctr]{"isoid"};
		if($testid == $id)
		{
			$id2a = $data1[$ctr]{$idflag2a};
			$id2b = $data1[$ctr]{$idflag2b};
			$id3a = $data1[$ctr]{$idflag3a};
			$id3b = $data1[$ctr]{$idflag3b};
			$id4a = $data1[$ctr]{$idflag4a};
			$id4b = $data1[$ctr]{$idflag4b};

			# Detect MS/MS IDs (digestionid & digestpepid = -1)
			if($id2a == -1 && $id2b == -1)
			{
				$data2b[0] = $data1[$ctr];
				$data3b[0] = $data1[$ctr];
				#$data4b[0] = $data1[$selindex1];
				$bold2 = 0;
				$bold3 = 0;
				$bold4 = 0;
			}

		}
	}



	for(my $ctr = 0; $ctr <= $#data1; $ctr++)
	{
		$idflag = 0;
		$test2a = $data1[$ctr]{$idflag2a};
		$test2b = $data1[$ctr]{$idflag2b};
		$test3a = $data1[$ctr]{$idflag3a};
		$test3b = $data1[$ctr]{$idflag3b};
		$test4a = $data1[$ctr]{$idflag4a};
		$test4b = $data1[$ctr]{$idflag4b};

		$testid = $data1[$ctr]{"isoid"};
		if($testid == $id)
		{
			$idflag = 1;
		}

		if($test2a >= 0 && $test2b >= 0)
		{
			if($test2a == $id2a || $test2b == $id2b)
			{
				if($data2b[0] eq "---")
				{
					$data2b[0] = $data1[$ctr];
				}
				else
				{
					push(@data2b, $data1[$ctr]);
				}
				if($idflag == 1)
				{
					$bold2 = $#data2b;
				}
			}
		}
		if($test3a >= 0 && $test3b >= 0)
		{
			if($test3a == $id3a || $test3b == $id3b)
			{
				if($idflag != 1 && ( $test2a == $id2a || $test2b == $id2b ) )
				{
					# Do Nothing
				}
				else
				{
					if($data3b[0] eq "---")
					{
						$data3b[0] = $data1[$ctr];
					}
					else
					{
						push(@data3b, $data1[$ctr]);
					}
					if($idflag == 1)
					{
						$bold3 = $#data3b;
					}
				}
			}
		}
		if($test4a == $id4a || $test4b == $id4b || $test4a == $id4b || $test4b == $id4a)
		{
			if($data4b[0] eq "---")
			{
				$data4b[0] = $data1[$ctr];
			}
			else
			{
				push(@data4b, $data1[$ctr]);
			}
			if($idflag == 1)
			{
				$bold4 = $#data4b;
			}
		}
	}

	regen_view234b($bold2, $bold3, $bold4);
}

sub color_pts234
{
	my $ctr2, $ctr3, $ctr4;

	# Color in a specific order
	# 1) Using the same peaks (purple)
	# 2) Same peptide (green)
	# 3) Same ion (blue)
	# 4) The actual point itself is colored elsewhere
	for($ctr2 = 0; $ctr2 <= $#data2; $ctr2++)
	{
		$isoid = $data2[$ctr2]{"isoid"};
		$tag = join '', "pt", $isoid;
		$canvasr->itemconfigure($tag, -outline=>'#FF00FF', -width=>'3');
		$realcanvasr->raise($tag, 'all');
	}
	for($ctr3 = 0; $ctr3 <= $#data3; $ctr3++)
	{
		$isoid = $data3[$ctr3]{"isoid"};
		$tag = join '', "pt", $isoid;
		$canvasr->itemconfigure($tag, -outline=>'#00FF00', -width=>'3');
		$realcanvasr->raise($tag, 'all');
	}
	for($ctr4 = 0; $ctr4 <= $#data4; $ctr4++)
	{
		$isoid = $data4[$ctr4]{"isoid"};
		$tag = join '', "pt", $isoid;
		$canvasr->itemconfigure($tag, -outline=>'#0000FF', -width=>'3');
		$realcanvasr->raise($tag, 'all');
	}
}

#
# sets the variable $isoid_curselection
#
sub isoid_select
{
	@sa1 = $list1->curselection;
	$s1 = $sa1[0];
	@sa2 = $list2->curselection;
	$s2 = $sa2[0];
	@sa3 = $list3->curselection;
	$s3 = $sa3[0];
	@sa4 = $list4->curselection;
	$s4 = $sa4[0];
	@sa2b = $list2b->curselection;
	$s2b = $sa2b[0];
	@sa3b = $list3b->curselection;
	$s3b = $sa3b[0];
	@sa4b = $list4b->curselection;
	$s4b = $sa4b[0];

	if(defined($s1))
	{
		$isoid_curselection = $data1[$s1]{"isoid"};
		$array_curselection = 1;
	}
	elsif(defined($s2))
	{
		$isoid_curselection = $data2[$s2]{"isoid"};
		$array_curselection = 2;
	}
	elsif(defined($s3))
	{
		$isoid_curselection = $data3[$s3]{"isoid"};
		$array_curselection = 3;
	}
	elsif(defined($s4))
	{
		$isoid_curselection = $data4[$s4]{"isoid"};
		$array_curselection = 4;
	}
	elsif(defined($s2b))
	{
		$isoid_curselection = $data2b[$s2b]{"isoid"};
		$array_curselection = 5;
	}
	elsif(defined($s3b))
	{
		$isoid_curselection = $data3b[$s3b]{"isoid"};
		$array_curselection = 6;
	}
	elsif(defined($s4b))
	{
		$isoid_curselection = $data4b[$s4b]{"isoid"};
		$array_curselection = 7;
	}
	else
	{
		$isoid_curselection = -1;
		print "No current selection\n";
	}

	#print "|$s1-$s2-$s3-$s4-$s2b-$s3b-$s4b|\n";
}

sub shunt
{
	isoid_select();
	$shuntid = $isoid_curselection;

	if($shuntid == -1)
	{
		print "Nothing to move\n";
		return;
	}

	$list1->selectionClear(0, $#view1);
	$list1->selectionSet($iso2idx{$isoid});
	$list1->see($iso2idx{$isoid});
	process_selection1();
}

sub put_photo
{
	my @params = @_;
	my $fitindex = $params[0];
	my $arrayid = $params[1];
	my $isofile = "";
	my $pngfile = "";
	my $contourfile = "";

	if($arrayid == 1)
	{
		#print "1\n";
		$isofile = $data1[$fitindex]{"isofile"};
		$contourfile = $data1[$fitindex]{"file"};
		#print "isofile: $isofile\n";
	}
	elsif($arrayid == 2)
	{
		#print "2\n";
		$isofile = $data2[$fitindex]{"isofile"};
		$contourfile = $data2[$fitindex]{"file"};
	}
	elsif($arrayid == 3)
	{
		#print "3\n";
		$isofile = $data3[$fitindex]{"isofile"};
		$contourfile = $data3[$fitindex]{"file"};
	}
	elsif($arrayid == 4)
	{
		#print "4\n";
		$isofile = $data4[$fitindex]{"isofile"};
		$contourfile = $data4[$fitindex]{"file"};
	}
	elsif($arrayid == 5)
	{
		$isofile = $data2b[$fitindex]{"isofile"};
		$contourfile = $data2b[$fitindex]{"file"};
	}
	elsif($arrayid == 6)
	{
		$isofile = $data3b[$fitindex]{"isofile"};
		$contourfile = $data3b[$fitindex]{"file"};
	}
	elsif($arrayid == 7)
	{
		$isofile = $data4b[$fitindex]{"isofile"};
		$contourfile = $data4b[$fitindex]{"file"};
	}

	if(substr($isofile, -3, 3) eq "txt")
	{
		$pngfile = $isofile;
		substr($pngfile, -3, 3) = "fit.png";
	}
	else
	{
		$pngfile = join '', $isofile, ".fit.png";
	}

	substr($contourfile, -3, 3) = "contour.png";

	if($lastphoto)
	{
		$lastphoto->delete;
		if(defined($lastphoto2))
		{
			$lastphoto2->delete;
		}
	}

	if(-e $pngfile)
	{
		$photo = $canvasl->Photo(-format=>'png', -file=>$pngfile);
		$canvasl->delete('plotimage');
		$canvasl->createImage(0,0, -image=>$photo, -anchor=>'nw', -tag=>'plotimage');
	}

	if(-e $contourfile)
	{
		$photo2 = $canvasr2->Photo(-format=>'png', -file=>$contourfile);
		$canvasr2->delete('contourimage');
		$canvasr2->createImage(0,0, -image=>$photo2, -anchor=>'nw', -tag=>'contourimage');
	}

	$lastphoto = $photo;
	$lastphoto2 = $photo2;
}

##############################
# unused gnuplot stuff below #
##############################
sub plot_fit
{
	my @params = @_;
	my $fitindex = $params[0];
	my $arrayid = $params[1];
	my $isofile = "";
	my $fitfile = "";
	my $datfile = "";

	print "Plotting... ($fitindex, $arrayid)\n";

	if($arrayid == 1)
	{
		print "1\n";
		$isofile = $data1[$fitindex]{"isofile"};
		print "isofile: $isofile\n";
	}
	elsif($arrayid == 2)
	{
		print "2\n";
		$isofile = $data2[$fitindex]{"isofile"};
	}
	elsif($arrayid == 3)
	{
		print "3\n";
		$isofile = $data3[$fitindex]{"isofile"};
	}
	elsif($arrayid == 4)
	{
		print "4\n";
		$isofile = $data4[$fitindex]{"isofile"};
	}

	print "||$isofile||\n";

	if(substr($isofile, -3, 3) eq "txt")
	{
		$fitfile = $isofile;
		substr($fitfile, -3, 3) = "fit";
		$datfile = $isofile;
		substr($datfile, -3, 3) = "dat";
	}
	else
	{
		$fitfile = join '', $isofile, ".fit";
		$datfile = join '', $isofile, ".dat";
	}

	print "plot \"$fitfile\", \"$datfile\"\n";
	print GIN1 "plot \"$fitfile\", \"$datfile\"\n";
}

sub go1
{
	print GIN1 "a=$gnuctr1\n";
	print GIN1 "plot a\n";
	$gnuctr1++;
}

sub go2
{
	print GIN2 "a=$gnuctr2\n";
	print GIN2 "plot a\n";
	$gnuctr2++;
}

sub read_gerr1
{
	my $num = sysread(GERR1, my $buffer, 1024);
	print "$num\n";
	print "$buffer\n";
}

sub read_gerr2
{
	my $num = sysread(GERR2, my $buffer, 1024);
}

sub read_gout1
{
	my $buffer = <GOUT1>;
	my $can = $canvasl;
	eval($buffer);
	$mw->update;
}

sub read_gout2
{
	my $buffer = <GOUT2>;
	my $can = $canvasr;
	eval($buffer);
	$mw->update;
}

sub iso3_regen
{
	isoid_select();
	$regen_isoid = $isoid_curselection;
	$regen_idx = $iso2idx{$regen_isoid};

	if($regen_isoid == -1)
	{
		print "Nothing selected to regenerate\n";
		return;
	}

	my @regarray = split /\./, $iso3_modelfile;

	# If the user has specified a massparam file, extract the model file from it
	# Also extract the atom file, allowing for overriding of hidden default
	if($regarray[$#regarray] eq "massparam")
	{
		open MASSPARAM, "$iso3_modelfile" or die "Can't open $iso3_modelfile\n";
		MP: while(defined($mpinput=<MASSPARAM>))
		{
			chomp($mpinput);
			my @mparray = split ' ', $mpinput;
			if($mparray[0] eq "model")
			{
				$iso3_modelfile = $mparray[1];
			}
			if($mparray[0] eq "atoms")
			{
				$iso3_atomsfile = $mparray[1];
			}
		}
		close MASSPARAM;
	}

	# Original atoms file: $iso3_atomsfile
	# Original model file: $iso3_modelfile

	# Keep original atoms file
	# Generate a new model file with appropriate parameters
	# Generate a .in file
	# Generate a .batch file
	# Capture the .fit and .dat
	# Capture the .csv

	# Retrive B/GW/OFF
	my $rg_b = sprintf("%f", $data1[$regen_idx]{"b"});
	my $rg_gw = sprintf("%f", $data1[$regen_idx]{"gw"});
	my $rg_off = sprintf("%f", $data1[$regen_idx]{"off"});
	# Search for all AMP and FRC entries

	@rg_amps = ();
	@rg_frcs = ();

	foreach $rg_key (keys %{$data1[$regen_idx]})
	{
		if(substr($rg_key, 0, 4) eq "amp_")
		{
			push(@rg_amps, $data1[$regen_idx]{$rg_key});
		}
		if(substr($rg_key, 0, 4) eq "frc_")
		{
			push(@rg_frcs, $data1[$regen_idx]{$rg_key});
		}
	}

	open MODEL, "$iso3_modelfile" or die "Can't open $iso3_modelfile\n";
	open RMODEL, ">regen.model" or die "Can't open regen.model\n";
	my $rline = <MODEL>;
	print RMODEL "$rline";
	$rline = <MODEL>;
	chomp($rline);
	@rarray = split ' ', $rline;
	if(($#rg_amps+1) != $rarray[0])
	{
		print "Error, num_amplitudes != num_species\n";
		return;
	}
	print RMODEL "$rline\n";

	# Deal with species and amplitudes
	for($rctr = 0; $rctr <= $#rg_amps; $rctr++)
	{
		$rline = <MODEL>;
		chomp($rline);
		@rarray = split ' ', $rline;
		$rarray[1] = sprintf("%f", $rg_amps[$rctr]);
		$rline = join ' ', @rarray;
		print RMODEL "$rline\n";
	}

	# Deal with atom types
	$r_num_var = 0; # Number of atoms and/or residues set to "variable"
	$rline = <MODEL>;
	print RMODEL "$rline";
	chomp($rline);
	@rarray = split ' ', $rline;
	$r_num_atom = $rarray[0];
	for($rctr = 0; $rctr < $r_num_atom; $rctr++)
	{
		$rline = <MODEL>;
		chomp($rline);
		@rarray = split ' ', $rline;
		if($rarray[1] eq "variable")
		{
			$rarray[2] = sprintf("%f", $rg_frcs[$r_num_var]);
			$r_num_var += 1;
			$rline = join ' ', @rarray;
		}
		print RMODEL "$rline\n";
	}

	# Deal with residues
	while(defined($rline=<MODEL>))
	{
		chomp($rline);
		@rarray = split ' ', $rline;
		if($rarray[1] eq "variable")
		{
			$rarray[2] = $rg_frcs[$r_num_var];
			$r_num_var += 1;
			$rline = join ' ', @rarray;
		}
		print RMODEL "$rline\n";
		for($rctr = 0; $rctr <= $#rg_amps; $rctr++)
		{
			$rline = <MODEL>;
			print RMODEL "$rline";
		}
	}

	close MODEL;
	close RMODEL;

	open RIN, ">regen.in" or die "Can't open regen.in\n";
	print RIN "fitit\n";
	print RIN "regen.batch\n";
	print RIN "$iso3_atomsfile\n";
	print RIN "regen.model\n";
	print RIN "5\n";
	print RIN "100.0\n";
	print RIN "$rg_b fixed\n";
	print RIN "$rg_off fixed\n";
	print RIN "$rg_gw fixed\n";
	close RIN;

	@regen_array = split /\//, $data1[$regen_idx]{'file'};
	$regen_base = $regen_array[$#regen_array];
	substr($regen_base, -4, 4) = "";

	open RBATCH, ">regen.batch" or die "Can't open regen.batch\n";
	#print RBATCH "$data1[$regen_idx]{'seq'} $data1[$regen_idx]{'charge'} $data1[$regen_idx]{'file'}\n";
	print RBATCH "$data1[$regen_idx]{'seq'} $data1[$regen_idx]{'charge'} regen_peaks/$regen_base.txt\n";
	close RBATCH;

	`mkdir regen_peaks`;
	`cp $data1[$regen_idx]{'file'} regen_peaks/`;

	$hostname = `hostname`;
	chomp($hostname);
	if($hostname eq "xtopher")
	{
		`$ENV{'MASSACRE_PATH'}/bin/isodist3_ppc regen.in`;
	}
	else
	{
		`$ENV{'MASSACRE_PATH'}/bin/isodist3_x86 regen.in`;
	}

	#$rfitfile = $data1[$regen_idx]{'file'};
	#$rdatfile = $data1[$regen_idx]{'file'};
	$rfitfile = join '', "regen_peaks/", $regen_base, ".fit";
	$rdatfile = join '', "regen_peaks/", $regen_base, ".dat";

	#substr($rfitfile, -3, 3) = "fit";
	#substr($rdatfile, -3, 3) = "dat";

	`mv $rfitfile ./regen.fit`;
	`mv $rdatfile ./regen.dat`;

	#`rm regen.batch`;
	#`rm regen.model`;
	#`rm regen.in`;

	#open RGNU, ">regen.gnuplot" or die "Can't open regen.gnuplot\n";
	#print RGNU "set term x11\n";
	print GIN1 "plot \"regen.dat\" with points, \"regen.fit\" with lines\n";
	#close RGNU;

	`gnuplot < regen.gnuplot`;

	#`rm regen.fit`;
	#`rm regen.dat`;
	`rm regen_peaks/$regen_base.txt`;
	`rmdir regen_peaks`;
}

sub launch_redun
{
	# Automatically perform the redundancy checking on launch
	flag_redundant();

	$redunw = MainWindow->new();
	$redunw->configure(-title=>'Locate Redundancies');
	$redunframe = $redunw->Frame(-border=>'5')->pack(-side=>'top');
	$redunframe1 = $redunframe->Frame()->pack(-side=>'top');
	$redunframe2 = $redunframe->Frame()->pack(-side=>'top');
	$redunframe3 = $redunframe->Frame()->pack(-side=>'top');
	$redunframe4 = $redunframe->Frame()->pack(-side=>'top', -fill=>'both', -expand=>'1');

	$obey_display_flag = 1;

	$hidepeaksbutton = $redunframe1->Checkbutton(-text=>'Show Peptides With Shared Peaks', -variable=>\$show_redunpeaks, -onvalue=>'1', -offvalue=>'0', -command=>\&redun_togglepeaks)->pack(-side=>'left');

	$redundelbutton = $redunframe4->Button(-text=>'Delete', -command=>\&redun_delete)->pack(-side=>'left');
	$redunquitbutton = $redunframe4->Button(-text=>'Quit', -command=>\&redun_quit)->pack(-side=>'right');
}

sub redun_delete
{
	for(my $ctr = 0; $ctr <= $#data1; $ctr++)
	{
		if($data1[$ctr]{"peak_redun"} == 1)
		{
			push(@deleted, splice(@data1, $ctr, 1));
			$ctr--;
		}
	}

	regen_view1();
	plot_frac();
	$lastsel1 = -1;
}

sub redun_togglepeaks
{
	for(my $ctr = 0; $ctr <= $#data1; $ctr++)
	{
		if($show_redunpeaks == 1)
		{
			$data1[$ctr]{"display_flag"} = 1;
		}
		elsif($data1[$ctr]{"peak_redun"} == 1)
		{
			$data1[$ctr]{"display_flag"} = 0;
		}
	}

	plot_frac();
}

sub redun_quit
{
	$redunw->destroy();
}

sub flag_redundant
{
	for(my $ctr = 0; $ctr <= $#data1; $ctr++)
	{
		$data1[$ctr]{"peak_redun"} = 0; # Set to 0 for shared peaks
		$data1[$ctr]{"seq_redun"} = 0; # Set to 0 for multiple IDs of same sequence
	}

	RONE: for(my $ctr1 = 0; $ctr1 <= $#data1; $ctr1++)
	{
		my $n14mass1 = $data1[$ctr1]{"n14mass"}; # Exact N14 Mass
		my $n15mass1 = $data1[$ctr1]{"n15mass"}; # Exact N15 Mass
		my $ailn141 = $data1[$ctr1]{"ailid_n14"}; # AIL N14 Peak ID
		my $ailn151 = $data1[$ctr1]{"ailid_n15"}; # AIL N15 Peak ID
		my $pep1 = $data1[$ctr1]{"isopep"}; # Peptide Sequence
		my $pep1sort = join('', sort(split '', $pep1)); # Sorted Peptide Sequence
		my $charge1 = $data1[$ctr1]{"charge"};

		RTWO: for(my $ctr2 = ($ctr1+1); $ctr2 <= $#data1; $ctr2++)
		{
			my $n14mass2 = $data1[$ctr2]{"n14mass"}; # Exact N14 Mass
			my $n15mass2 = $data1[$ctr2]{"n15mass"}; # Exact N15 Mass
			my $ailn142 = $data1[$ctr2]{"ailid_n14"}; # AIL N14 Peak ID
			my $ailn152 = $data1[$ctr2]{"ailid_n15"}; # AIL N15 Peak ID
			my $pep2 = $data1[$ctr2]{"isopep"}; # Peptide Sequence
			my $pep2sort = join('', sort(split '', $pep2)); # Sorted Peptide Sequence
			my $charge2 = $data1[$ctr2]{"charge"};

			# N14 Only Search
			if($ailn151 == 0 && $ailn152 == 0)
			{
				if($ailn141 == $ailn142)
				{
					$data1[$ctr1]{"peak_redun"} = 1;
					$data1[$ctr2]{"peak_redun"} = 1;
				}
			}
			# N15 Only Search
			elsif($ailn141 == 0 && $ailn142 == 0)
			{
				if($ailn151 == $ailn152)
				{
					$data1[$ctr1]{"peak_redun"} = 1;
					$data1[$ctr2]{"peak_redun"} = 1;
				}
			}
			# This happens when pairs are searched for
			else
			{
				if($ailn141 == $ailn142 ||
					$ailn151 == $ailn152 ||
					$ailn141 == $ailn152 ||
					$ailn151 == $ailn142 )
				{
					$data1[$ctr1]{"peak_redun"} = 1;
					$data1[$ctr2]{"peak_redun"} = 1;
				}
			}

			# NOTE: In practice, all peptide redundancy cases are subsets of shared peaks
			# Things with the same mass will always map to the same peaks

			# Same identical peptide
			# This is a lesser sequence condition as it is often resolved in RT
			# You compare to RT of other ion states
			if($pep1 eq $pep2 && $charge1 == $charge2)
			{
				# Add this check so as not to downgrade a 2 to a 1
				if($data1[$ctr1]{"seq_redun"} == 0)
				{
					$data1[$ctr1]{"seq_redun"} = 1;
				}
				if($data1[$ctr2]{"seq_redun"} == 0)
				{
					$data1[$ctr2]{"seq_redun"} = 1;
				}
			}
			# Peptides are different, but formula is the same
			# This is bad, as impossible to distinguish
			elsif($pep1sort eq $pep2sort && $charge1 == $charge2)
			{
				$data1[$ctr1]{"seq_redun"} = 2;
				$data1[$ctr2]{"seq_redun"} = 2;
			}
		}
	}
}

sub numerically
{
	return $a <=> $b;
}

sub trim
{
	my $string = shift;
	$string =~ s/^\s+//; #remove leading spaces
	$string =~ s/\s+$//; #remove trailing spaces
	return $string;
}

sub proteomic
{
	my $loca1;
	my $locb1;
	my $loca2;
	my $locb2;
	my $loca3;
	my $locb3;

	if(!defined($pro2loc70S{$a}))
	{
		$loca1 = 1000000000;
	}
	else
	{
		$loca1 = $pro2loc70S{$a};
	}

	if(!defined($pro2loc70S{$b}))
	{
		$locb1 = 1000000000;
	}
	else
	{
		$locb1 = $pro2loc70S{$b};
	}

	if(!defined($pro2loc80Syeast{$a}))
	{
		$loca2 = 1000000000;
	}
	else
	{
		$loca2 = $pro2loc80Syeast{$a};
	}

	if(!defined($pro2loc80Syeast{$b}))
	{
		$locb2 = 1000000000;
	}
	else
	{
		$locb2 = $pro2loc80Syeast{$b};
	}

	if(!defined($pro2loc80Shuman{$a}))
	{
		$loca3 = 1000000000;
	}
	else
	{
		$loca3 = $pro2loc80Shuman{$a};
	}

	if(!defined($pro2loc80Shuman{$b}))
	{
		$locb3 = 1000000000;
	}
	else
	{
		$locb3 = $pro2loc80Shuman{$b};
	}

	$loca1 <=> $locb1
		or
	$loca2 <=> $locb2
		or
	$loca3 <=> $locb3
		or
	$a cmp $b;


#	$pro2loc70S{$a} <=> $pro2loc70S{$b}
#		or
#	$pro2loc80Syeast{$a} <=> $pro2loc80Syeast{$b}
#		or
#	$pro2loc80Shuman{$a} <=> $pro2loc80Shuman{$b}
#		or
#	$a cmp $b
#	;
}

sub proteomic2
{
	my $a_is_acc = 0;
	my $b_is_acc = 0;
	my $a_type = 0;
	my $b_type = 0;
	my $a_num = 0;
	my $b_num = 0;
	my $a_mod = 0;
	my $b_mod = 0;

	my $astring = $a;
	my $bstring = $b;

	my $type;

	# If these match, $a/$b are accession numbers
	if($a eq $pro2acc{$a})
	{
		#print "a eq\n";
		$a_is_acc = 1;
	}
	if($b eq $pro2acc{$b})
	{
		#print "b eq\n";
		$b_is_acc = 1;
	}

	# Check for S--, L--, RS--, RL--, RT--, RM--

	if($a =~ /^S/)
	{
		$a_type = 1;
		substr($astring, 0, 1) = "";
	}
	elsif($a =~ /^L/)
	{
		$a_type = 2;
		substr($astring, 0, 1) = "";
	}
	elsif($a =~ /^RS/)
	{
		$a_type = 3;
		substr($astring, 0, 2) = "";
	}
	elsif($a =~ /^RL/)
	{
		$a_type = 4;
		substr($astring, 0, 2) = "";
	}
	elsif($a =~ /^RT/)
	{
		$a_type = 5;
		substr($astring, 0, 2) = "";
	}
	elsif($a =~ /^RM/)
	{
		$a_type = 6;
		substr($astring, 0, 2) = "";
	}

	if($astring =~ /(\d+)(\D+)/)
	{
		$a_num = $1;
		$a_mod = $2;
	}
	elsif($astring =~ /(\D+)(\d+)/)
	{
		$a_num = $1;
		$a_mod = $2;
	}
	else
	{
		$a_num = $astring;
	}

	if($b =~ /^S/)
	{
		$b_type = 1;
		substr($bstring, 0, 1) = "";
	}
	elsif($b =~ /^L/)
	{
		$b_type = 2;
		substr($bstring, 0, 1) = "";
	}
	elsif($b =~ /^RS/)
	{
		$b_type = 3;
		substr($bstring, 0, 2) = "";
	}
	elsif($b =~ /^RL/)
	{
		$b_type = 4;
		substr($bstring, 0, 2) = "";
	}
	elsif($b =~ /^RT/)
	{
		$b_type = 5;
		substr($bstring, 0, 2) = "";
	}
	elsif($b =~ /^RM/)
	{
		$b_type = 6;
		substr($bstring, 0, 2) = "";
	}

	if($bstring =~ /(\d+)(\D+)/)
	{
		$b_num = $1;
		$b_mod = $2;
	}
	elsif($bstring =~ /(\D+)(\d+)/)
	{
		$b_num = $1;
		$b_mod = $2;
	}
	else
	{
		$b_num = $bstring;
	}

	$a_is_acc <=> $b_is_acc
		or
	$a_type <=> $b_type
		or
	$a_num <=> $b_num
		or
	$a_mod <=> $b_mod
		or
	$pro2acc{$a} cmp $pro2acc{$b};
}

sub display_rtcorr
{
	#print "display_rtcorr\n";
	isoid_select();
	$rt_isoid = $isoid_curselection;
	$rt_idx = $iso2idx{$rt_isoid};

	if($rt_isoid == -1)
	{
		print "Nothing selected to display RT correction plots for\n";
		return;
	}

	my @rtfiles = split /\#/, $data1[$rt_idx]{'pks_rtcorrfile_n14'};

	for(my $ctr = 0; $ctr <= $#rtfiles; $ctr++)
	{
		print "$rtfiles[$ctr]\n";
		my $rtwindow = MainWindow->new();
		my $rtcanvas = $rtwindow->Canvas(-background=>'white', -width=>'640', -height=>'480')->pack();
		my $rtphoto = $rtcanvas->Photo(-format=>'png', -file=>$rtfiles[$ctr]);
		$rtcanvas->createImage(0,0, -image=>$rtphoto, -anchor=>'nw');
	}
}
