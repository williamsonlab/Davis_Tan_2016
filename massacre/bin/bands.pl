#!/usr/bin/perl

use strict;
use warnings;
use Tk;
require Tk::Image;
use Tk::Image;
use Tk::PNG;
use Tk::JPEG;

my $invert = 0;

my @rectangles;
my @rectcoor;
my @rectscale;
my @rectlabels;
my @rectdisplay = ("---");
#my @rectdisplay;
my $numrect = 0; # Just a display/label thing - always increments up, delete rectangles have no effect on this value
my $nextrect = 0; # Index of a new rectangle - replacing the old $numrect usage, equivalent to ($#rectdisplay + 1)
my $infilename = "test.png";
my $outfilename = "bands.txt";
my $statefilename = "mystate.bands";
my $photo;
my $photo2;
my $image;
my $image_orig;
my $radioselect;

my $bg;
my @bgcoor;
my $bgscale;
my $bglabel;

my @listselarray; 
my $listselect;
my $prevselect;

my $area;
my $intensity;
my $intensity_nobgcorr;
my $intensity_display;
my $intensity_display_nobgcorr;
my $background;
my $background_display;

my $xprofilebase;
my $yprofilebase;
my $xprofile;
my $yprofile;
my @xprofilecoor;
my @yprofilecoor;

my $scale = 1;

my $lanewidth;
my $lw;
my @lwcoor;
my $lwscale;
my $lwlabel;

##########
##########

my $mw = new MainWindow;
$mw->configure(-title=>'bands v1.1');
$mw->minsize( qw(1024 768));

my $frame = $mw->Frame(-width=>'1024', -height=>'768')->pack(-expand=>'1', -fill=>'both');
my $frame_left = $frame->Frame(-width=>'768', -height=>'768')->pack(-side =>'left', -fill=>'both', -expand=>'1');
my $frame_right = $frame->Frame(-width=>'256', -height=>'768')->pack(-side=>'top', -fill=>'y');

my $frame_right_A = $frame_right->Frame(-width=>'256', -height=>'100')->pack(-side=>'top', -fill=>'x', -expand=>'1');
my $frame_right_A2 = $frame_right->Frame(-height=>'10')->pack(-side=>'top');
my $frame_right_B = $frame_right->Frame(-width=>'256', -height=>'100')->pack(-side=>'top', -fill=>'x', -expand=>'1');
my $frame_right_B2 = $frame_right->Frame(-height=>'10')->pack(-side=>'top');
my $frame_right_Ca = $frame_right->Frame(-width=>'256', -height=>'100')->pack(-side=>'top', -fill=>'x', -expand=>'1');
my $frame_right_Cb = $frame_right->Frame(-width=>'256', -height=>'100')->pack(-side=>'top', -fill=>'x', -expand=>'1');
my $frame_right_Cc = $frame_right->Frame(-width=>'256', -height=>'100')->pack(-side=>'top', -fill=>'x', -expand=>'1');
my $frame_right_Cd = $frame_right->Frame(-width=>'256', -height=>'100')->pack(-side=>'top', -fill=>'x', -expand=>'1');
my $frame_right_Ce = $frame_right->Frame(-width=>'256', -height=>'100')->pack(-side=>'top', -fill=>'x', -expand=>'1');

my $frame_right_C2 = $frame_right->Frame(-height=>'10')->pack(-side=>'top');
my $frame_right_D = $frame_right->Frame(-width=>'256')->pack(-side=>'top', -fill=>'both', -expand=>'1');
	my $frame_right_DL = $frame_right_D->Frame()->pack(-side=>'left', -fill=>'both');
	my $frame_right_DR = $frame_right_D->Frame()->pack(-side=>'right', -fill=>'both');
my $frame_right_D2 = $frame_right->Frame(-height=>'10')->pack(-side=>'top');
my $frame_right_E = $frame_right->Frame(-width=>'256', -height=>'100')->pack(-side=>'top', -fill=>'x', -expand=>'1');
my $frame_right_E1 = $frame_right->Frame(-width=>'256', -height=>'100')->pack(-side=>'top', -fill=>'x', -expand=>'1');
my $frame_right_E2 = $frame_right->Frame(-height=>'10')->pack(-side=>'top');

my $frame_right_F0 = $frame_right->Frame(-height=>'10')->pack(-side=>'top');
my $frame_right_F1 = $frame_right->Frame(-width=>'256', -height=>'100')->pack(-side=>'top', -fill=>'x', -expand=>'1');
#my $frame_right_F2 = $frame_right->Frame(-height=>'10')->pack(-side=>'top');

my $frame_right_G = $frame_right->Frame(-width=>'256', -height=>'100')->pack(-side=>'top', -fill=>'x', -expand=>'1');
my $frame_right_G2 = $frame_right->Frame(-height=>'10')->pack(-side=>'top');

my $frame_right_Z = $frame_right->Frame(-width=>'256', -height=>'100')->pack(-side=>'bottom', -fill=>'x', -expand=>'1');

#####

my $infilebox = $frame_right_A->Entry(-textvariable=>\$infilename)->pack(-side=>'left');
my $infilebutton = $frame_right_A->Button(-text=>'Load', -command=>\&load_gel)->pack(-side=>'right');

my $invert_blackpos = $frame_right_A2->Radiobutton(-text=>'Black Positive', -variable=>\$invert, -value=>'0')->pack(-side=>'top');
my $invert_whitepos = $frame_right_A2->Radiobutton(-text=>'White Positive', -variable=>\$invert, -value=>'1')->pack(-side=>'top');

my $zoominbutton = $frame_right_B->Button(-text=>'Zoom In', -command=>\&zoom_in)->pack(-side=>'left');
my $zoomoutbutton = $frame_right_B->Button(-text=>'Zoom Out', -command=>\&zoom_out)->pack(-side=>'left');

my $selector0 = $frame_right_Ca->Radiobutton(-text=>'Set Background', -variable=>\$radioselect, -value=>'0')->pack(-side=>'left');
$selector0->select;
my $selector1 = $frame_right_Cb->Radiobutton(-text=>'Define Band', -variable=>\$radioselect, -value=>'1')->pack(-side=>'left');
my $selector2 = $frame_right_Cc->Radiobutton(-text=>'Repeat Band', -variable=>\$radioselect, -value=>'2')->pack(-side=>'left');
my $selector3 = $frame_right_Cd->Radiobutton(-text=>'Move Band', -variable=>\$radioselect, -value=>'3')->pack(-side=>'left');
my $selector4 = $frame_right_Ce->Radiobutton(-text=>'Gel Profile', -variable=>\$radioselect, -value=>'4')->pack(-side=>'left');

my $bgval = $frame_right_Ca->Entry(-textvariable=>\$background_display, -state=>'readonly', -width=>'5')->pack(-side=>'left');

my $rectlist = $frame_right_DL->Listbox(-listvariable=>\@rectdisplay, -state=>'normal', -selectmode=>'browse', -height=>'20', -font=>'Courier', -width=>'4')->pack(-side=>'left', -fill=>'none');
my $listscrolly = $frame_right_DL->Scrollbar(-orient=>'vertical')->pack(-side=>'left', -fill=>'y');
$rectlist->configure(-yscrollcommand=>['set'=>$listscrolly]);
$listscrolly->configure(-command=>['yview'=>$rectlist]);
$rectlist->bind('<<ListboxSelect>>'=>\&process_selection);

my $arealabel = $frame_right_DR->Label(-text=>'Area:')->pack(-side=>'top');
my $areaval = $frame_right_DR->Entry(-textvariable=>\$area, -state=>'readonly')->pack(-side=>'top');
my $intensitylabel0 = $frame_right_DR->Label(-text=>'Intensity (Uncorr):')->pack(-side=>'top');
my $intensityval0 = $frame_right_DR->Entry(-textvariable=>\$intensity_display_nobgcorr, -state=>'readonly')->pack(-side=>'top');
my $intensitylabel = $frame_right_DR->Label(-text=>'Intensity (Corr):')->pack(-side=>'top');
my $intensityval = $frame_right_DR->Entry(-textvariable=>\$intensity_display, -state=>'readonly')->pack(-side=>'top');
my $deletebutton = $frame_right_DR->Button(-text=>'Delete', -command=>\&delete_band)->pack(-side=>'top');

my $exportbox = $frame_right_E->Entry(-textvariable=>\$outfilename, -width=>'18')->pack(-side=>'left');
my $exportbutton = $frame_right_E->Button(-text=>'Export', -command=>\&export_bands)->pack(-side=>'right');

my $savestatebox = $frame_right_E1->Entry(-textvariable=>\$statefilename, -width=>'14')->pack(-side=>'left');
my $savestatebutton = $frame_right_E1->Button(-text=>'Save', -command=>\&save_state)->pack(-side=>'right');
my $loadstatebutton = $frame_right_E1->Button(-text=>'Load', -command=>\&load_state)->pack(-side=>'right');

my $selector5 = $frame_right_F1->Radiobutton(-text=>'Define Lane Width', -variable=>\$radioselect, -value=>'5')->pack(-side=>'left');
my $lwval = $frame_right_F1->Entry(-textvariable=>\$lanewidth, -state=>'readonly', -width=>'3')->pack(-side=>'left');

my $bgcoorbutton = $frame_right_G->Button(-text=>'Correct Background', -command=>\&bg_correct)->pack(-side=>'left');

my $quitbutton = $frame_right_Z->Button(-text=>'Quit', -command=>\&quit)->pack(-side=>'right');

##########
##########

my $box = [0, 0, 500, 500];

my $c  = $frame_left->Scrolled(qw/
       Canvas -bg #EEFFEE
       -xscrollincrement 1
       -yscrollincrement 1
       -confine 1
       -scrollbars se
       -width 700
       -height 700/,
       -scrollregion => $box,
       );
$c->pack(-fill=>'both', -expand=>'1');
my $realc = $c->Subwidget('scrolled'); 

$c->CanvasBind('<1>' => sub
{
	my $x = $c->canvasx($Tk::event->x);
	my $y = $c->canvasy($Tk::event->y);
	
	my $xctr;
	my $yctr;
				
	print "Before: $x $y\n";
	if($x % $scale != 0)
	{
		$x += $x % $scale;
	}
	if($y % $scale != 0)
	{
		$y += $y % $scale;
	}
	print "After: $x $y\n";
	
	if($radioselect == 0)
	{
		if(defined($bg))
		{
			$c->delete($bg);
			$c->delete($bglabel);
			#$bglabel->delete;
		}
		@bgcoor = ($x,$y,$x,$y);
		$bg = $c->createRectangle(@bgcoor, -outline=>'blue', -tags=>['ontop']);	
		$bgscale = $scale;
	}
	elsif($radioselect == 1)
	{
		if($rectdisplay[0] eq "---")
		{
			$nextrect = 0;
		}
		else
		{
			$nextrect = $#rectdisplay + 1;
		}
		print "NEXT: $nextrect\n";
		@{$rectcoor[$nextrect]} = ($x,$y,$x,$y);
		$rectangles[$nextrect] = $c->createRectangle(@{$rectcoor[$nextrect]}, -outline=>'red', -tags=>['ontop']);
		$rectscale[$nextrect] = $scale;
		$rectlist->selectionClear(0,$#rectcoor);
		$area = "";
		$intensity_display_nobgcorr = "";
		$intensity_display = "";
	}
	elsif($radioselect == 2)
	{
		if($#rectdisplay < 0 || $rectdisplay[0] eq "---")
		{
			print STDERR "Error: No bands to replicate\n";
		}
		else
		{
			$nextrect = $#rectdisplay + 1;
			print "NEXT: $nextrect\n";
			
			@listselarray = $rectlist->curselection;
			$listselect = $listselarray[0];		
		
			if(!defined($listselect))
			{
				$rectlist->selectionSet($#rectcoor);
				$listselect = $#rectcoor;
			}
			
			my $xshift = $x - $rectcoor[$listselect][0];
			my $yshift = $y - $rectcoor[$listselect][1];
			
			my $repscale = $scale / $rectscale[$listselect];
			
			@{$rectcoor[$nextrect]} = ($rectcoor[$listselect][0] + $xshift, $rectcoor[$listselect][1] + $yshift, $rectcoor[$listselect][2] + $xshift, $rectcoor[$listselect][3] + $yshift);
			$rectcoor[$nextrect][0] *= $repscale;
			$rectcoor[$nextrect][1] *= $repscale;
			$rectcoor[$nextrect][2] *= $repscale;
			$rectcoor[$nextrect][3] *= $repscale;
			
			$rectangles[$nextrect] = $c->createRectangle(@{$rectcoor[$nextrect]}, -outline=>'red', -tags=>['ontop']);

			$rectscale[$nextrect] = $scale;
			$rectlist->selectionClear(0,$#rectcoor);
			$area = "";
			$intensity = "";
		}
	}
	elsif($radioselect == 3)
	{
		if($#rectdisplay < 0 || $rectdisplay[0] eq "---")
		{
			print STDERR "Error: No bands to move\n";
		}
		else
		{		
			@listselarray = $rectlist->curselection;
			$listselect = $listselarray[0];		
		
			if(!defined($listselect))
			{
				$rectlist->selectionSet($#rectcoor);
				$listselect = $#rectcoor;
			}
			
			my $xshift = $x - $rectcoor[$listselect][0];
			my $yshift = $y - $rectcoor[$listselect][1];
			
			my $movescale = $scale / $rectscale[$listselect];

			$rectcoor[$listselect][0] += $xshift;
			$rectcoor[$listselect][1] += $yshift;
			$rectcoor[$listselect][2] += $xshift;
			$rectcoor[$listselect][3] += $yshift;
			$rectcoor[$listselect][0] *= $movescale;
			$rectcoor[$listselect][1] *= $movescale;
			$rectcoor[$listselect][2] *= $movescale;
			$rectcoor[$listselect][3] *= $movescale;
			
			$rectscale[$listselect] = $scale;

			$c->coords($rectangles[$listselect] => ($rectcoor[$listselect][0], $rectcoor[$listselect][1], $rectcoor[$listselect][2], $rectcoor[$listselect][3]) );
			$c->coords($rectlabels[$listselect] => ( ($rectcoor[$listselect][0] + $rectcoor[$listselect][2])/2, ($rectcoor[$listselect][1] + $rectcoor[$listselect][3])/2 ));
			
			$intensity_display_nobgcorr = "";
			$intensity_display = "";
		}
	}
	elsif($radioselect == 4)
	{
		my @colorarray;
		my $norm_intensity;
	
		if(defined($xprofilebase))
		{
			$c->delete($xprofilebase);
		}
		if(defined($yprofilebase))
		{
			$c->delete($yprofilebase);
		}
		if(defined($xprofile))
		{
			$c->delete($xprofile);
			@xprofilecoor = ();
		}
		if(defined($yprofile))
		{
			$c->delete($yprofile);
			@yprofilecoor = ();
		}
		
		$xprofilebase = $c->createLine(0, $y, $photo->width, $y, -fill=>'purple', -tags=>['ontop']);
		$yprofilebase = $c->createLine($x, 0, $x, $photo->height, -fill=>'orange', -tags=>['ontop']);

		for($xctr = 0; $xctr < $photo->width; $xctr++)
		{
			push(@xprofilecoor, $xctr);
			@colorarray = $photo->get($xctr, $y);
			if($invert == 0)
			{
				$norm_intensity = (( (255-$colorarray[0]) + (255-$colorarray[1]) + (255-$colorarray[2]) ) / 3) / 1;
			}
			elsif($invert == 0)
			{
				$norm_intensity = (($colorarray[0] + $colorarray[1] + $colorarray[2]) / 3) / 1;
			}
			$norm_intensity *= $scale;
			push(@xprofilecoor, $y - $norm_intensity);
		}
		$xprofile = $c->createLine( @xprofilecoor, -fill=>'purple', -tags=>['ontop']);

		for($yctr = 0; $yctr < $photo->height; $yctr++)
		{
			@colorarray = $photo->get($x, $yctr);
			if($invert == 0)
			{
				$norm_intensity = (( (255-$colorarray[0]) + (255-$colorarray[1]) + (255-$colorarray[2]) ) / 3) / 1;
			}
			elsif($invert == 0)
			{
				$norm_intensity = (($colorarray[0] + $colorarray[1] + $colorarray[2]) / 3) / 1;
			}
			$norm_intensity *= $scale;
			push(@yprofilecoor, $x - $norm_intensity);
			push(@yprofilecoor, $yctr);
		}
		$yprofile = $c->createLine( @yprofilecoor, -fill=>'orange', -tags=>['ontop']);
	}
	elsif($radioselect == 5)
	{
		if(defined($lw))
		{
			$c->delete($lw);
			$c->delete($lwlabel);
		}
		@lwcoor = ($x,$y,$x,$y);
		$lw = $c->createLine(@lwcoor, -fill=>'orange', -tags=>['ontop']);	
		$lwscale = $scale;	
		$lanewidth = abs($lwcoor[2] - $lwcoor[0]) + $lwscale;
		$lanewidth = $lanewidth / $lwscale;
	}
		
	print "$x $y\n";
	#print "Box: @$box[0] @$box[1] @$box[2] @$box[3]\n";
});

$c->CanvasBind('<B1-Motion>' => sub
{
	my $x = $c->canvasx($Tk::event->x);
	my $y = $c->canvasy($Tk::event->y);

	if($x % $scale != 0)
	{
		$x += $x % $scale;
	}
	if($y % $scale != 0)
	{
		$y += $y % $scale;
	}
	
	if($radioselect == 0)
	{
		$bgcoor[2] = $x;
		$bgcoor[3] = $y;
		$c->coords($bg => @bgcoor);
	}
	elsif($radioselect == 1)
	{
		$rectcoor[$nextrect][2] = $x;
		$rectcoor[$nextrect][3] = $y;
		$c->coords($rectangles[$nextrect] => @{$rectcoor[$nextrect]} );
				
		@rectcoor = @rectcoor;
	}
	elsif($radioselect == 2)
	{
		if($#rectdisplay < 0 || $rectdisplay[0] eq "---")
		{
			print STDERR "Error: No bands to replicate\n";
		}
		else
		{		
			my $xshift = $x - $rectcoor[$nextrect][0];
			my $yshift = $y - $rectcoor[$nextrect][1];
			$c->coords($rectangles[$nextrect] => ($rectcoor[$nextrect][0] + $xshift, $rectcoor[$nextrect][1] + $yshift, $rectcoor[$nextrect][2] + $xshift, $rectcoor[$nextrect][3] + $yshift) );

			$rectcoor[$nextrect][0] += $xshift;
			$rectcoor[$nextrect][1] += $yshift;
			$rectcoor[$nextrect][2] += $xshift;
			$rectcoor[$nextrect][3] += $yshift;
			
			@rectcoor = @rectcoor;
		}
	}
	elsif($radioselect == 3)
	{
		if($#rectdisplay < 0 || $rectdisplay[0] eq "---")
		{
			print STDERR "Error: No bands to move\n";
		}
		else
		{	
			my $xshift = $x - $rectcoor[$listselect][0];
			my $yshift = $y - $rectcoor[$listselect][1];

			$rectcoor[$listselect][0] += $xshift;
			$rectcoor[$listselect][1] += $yshift;
			$rectcoor[$listselect][2] += $xshift;
			$rectcoor[$listselect][3] += $yshift;

			$c->coords($rectangles[$listselect] => ($rectcoor[$listselect][0], $rectcoor[$listselect][1], $rectcoor[$listselect][2], $rectcoor[$listselect][3]) );
			$c->coords($rectlabels[$listselect] => ( ($rectcoor[$listselect][0]+$rectcoor[$listselect][2])/2, ($rectcoor[$listselect][1]+$rectcoor[$listselect][3])/2 ));
		}
	}
	elsif($radioselect == 4)
	{
		my @colorarray;
		my $norm_intensity;
		my $xctr;
		my $yctr;
	
		$c->delete($xprofilebase);
		$c->delete($yprofilebase);
		$c->delete($xprofile);
		@xprofilecoor = ();
		$c->delete($yprofile);
		@yprofilecoor = ();
		
		$xprofilebase = $c->createLine(0, $y, $photo->width, $y, -fill=>'purple', -tags=>['ontop']);
		$yprofilebase = $c->createLine($x, 0, $x, $photo->height, -fill=>'orange', -tags=>['ontop']);

		for($xctr = 0; $xctr < $photo->width; $xctr++)
		{
			push(@xprofilecoor, $xctr);
			@colorarray = $photo->get($xctr, $y);
			if($invert == 0)
			{
				$norm_intensity = (( (255-$colorarray[0]) + (255-$colorarray[1]) + (255-$colorarray[2]) ) / 3) / 1;
			}
			elsif($invert == 0)
			{
				$norm_intensity = (($colorarray[0] + $colorarray[1] + $colorarray[2]) / 3) / 1;
			}
			$norm_intensity *= $scale;
			push(@xprofilecoor, $y - $norm_intensity);
		}
		$xprofile = $c->createLine( @xprofilecoor, -fill=>'purple', -tags=>['ontop']);

		for($yctr = 0; $yctr < $photo->height; $yctr++)
		{
			@colorarray = $photo->get($x, $yctr);
			if($invert == 0)
			{
				$norm_intensity = (( (255-$colorarray[0]) + (255-$colorarray[1]) + (255-$colorarray[2]) ) / 3) / 1;
			}
			elsif($invert == 0)
			{
				$norm_intensity = (($colorarray[0] + $colorarray[1] + $colorarray[2]) / 3) / 1;
			}
			$norm_intensity *= $scale;
			push(@yprofilecoor, $x - $norm_intensity);
			push(@yprofilecoor, $yctr);
		}
		$yprofile = $c->createLine( @yprofilecoor, -fill=>'orange', -tags=>['ontop']);
	}
	elsif($radioselect == 5)
	{
		$lwcoor[2] = $x;
		$lwcoor[3] = $y;
		$c->coords($lw => @lwcoor);
		
		$lanewidth = abs($lwcoor[2] - $lwcoor[0]) + $lwscale;
		$lanewidth = $lanewidth / $lwscale;
	}

	print "--$x $y\n";
});

$c->CanvasBind('<B1-ButtonRelease>' => sub
{
	my $x = $c->canvasx($Tk::event->x);
	my $y = $c->canvasy($Tk::event->y);

	if($x % $scale != 0)
	{
		$x += $x % $scale;
	}
	if($y % $scale != 0)
	{
		$y += $y % $scale;
	}
	
	if($radioselect == 0)
	{
		my $fontx = ($bgcoor[0] + $bgcoor[2]) / 2;
		my $fonty = ($bgcoor[1] + $bgcoor[3]) / 2;
		
		$bglabel = $c->createText($fontx, $fonty, -text=>'bg', -fill=>'blue', -font=>'Courier', -tags=>['ontop']);
	
		calc_background();
		
		if($#rectdisplay < 0 || $rectdisplay[0] eq "---")
		{
			#
		}
		else
		{
			process_selection();
		}
	}
	elsif($radioselect == 1 || $radioselect == 2)
	{
		if( ($#rectdisplay < 0 || $rectdisplay[0] eq "---") && $radioselect == 2)
		{
			print STDERR "Error: No bands to replicate\n";
		}
		else
		{				
			my $fontx = ($rectcoor[$nextrect][0] + $rectcoor[$nextrect][2]) / 2;
			my $fonty = ($rectcoor[$nextrect][1] + $rectcoor[$nextrect][3]) / 2;
				
			$rectlabels[$nextrect] = $c->createText($fontx,$fonty, -text=>$numrect, -fill=>'red', -font=>'Courier', -tags=>['ontop']);
			$rectdisplay[$nextrect] = "$numrect";
			
			$rectlist->selectionSet($#rectcoor);
						
			$numrect++;
			
			print "----$x $y\n";
			print "$numrect rectangles\n";

			process_selection();
		}
	}
	elsif($radioselect == 3)
	{
		if($#rectdisplay < 0 || $rectdisplay[0] eq "---")
		{
			print STDERR "Error: No bands to move\n";
		}
		else
		{	
			process_selection();
		}
	}
	elsif($radioselect == 5)
	{
		my $fontx = ($lwcoor[0] + $lwcoor[2]) / 2;
		my $fonty = ($lwcoor[1] + $lwcoor[3]) / 2;
		
		$lwlabel = $c->createText($fontx, $fonty, -text=>'lw', -fill=>'orange', -font=>'Courier', -tags=>['ontop']);	
	}
});

sub bg_correct
{
	# Only perform the operation if the lane width has been defined
	if(!defined($lanewidth))
	{
		printf("Error: Lane Width not Defined\n");
	}
	else
	{
		my $pid = $$;
		my $bgcoor_infile = $infilename;
		my $bgcoor_outfile = join '', "tmp.bands.bgcoor.", $pid, ".png";

		# Call the external program
		`./bgcoor.exe $bgcoor_infile $bgcoor_outfile $lanewidth`;
		
		if(defined($image))
		{
			$c->delete($image);
		}

		$photo = $c->Photo(-file=>$bgcoor_outfile);

		$image = $c->createImage(0,0, -image=>$photo, -anchor=>'nw', -tags=>['gelimage']);
		
		my $iwidth = $photo->width;
		my $iheight = $photo->height;
		print "($iheight, $iwidth)\n";
		
		$c->configure(-scrollregion=>[0,0,$iwidth,$iheight]);

		$realc->raise('ontop', 'gelimage');
		
		# Check to make sure background has already been defined
		# Similar checking for selection is done within that subroutine
		if(defined($bg))
		{
			calc_background();
		}
		if(defined($rectlist->curselection))
		{
			process_selection();
		}
		
		# Remove the corrected image from the disk
		`rm $bgcoor_outfile`;
	}
}

sub zoom_in
{
	my $sf = 2;
	my $w_i;
	my $h_i;
	my $w_f;
	my $h_f;

	$scale = $scale * $sf;

	$w_i = $photo->width;
	$h_i = $photo->height;
	print "Before: $w_i, $h_i\n";

	$w_f = $sf*$w_i;
	$h_f = $sf*$h_i;

	$c->scale('all' => 0, 0, $sf, $sf);

	$photo2 = $c->Photo(-width=>$w_f, -height=>$h_f);
	$photo2->copy($photo, -zoom=>$sf,$sf, -from=>0,0,$w_i,$h_i, -to=>0,0,$w_f,$h_f);
	$photo->configure(-width=>$w_f, -height=>$h_f);
	$photo->copy($photo2);
	$image = $c->createImage(0,0, -image=>$photo, -anchor=>'nw', -tags=>['gelimage']);

	$_ *= $sf for @$box;
	$c->configure(-scrollregion=>[0,0,$w_f,$h_f]);
		
	$w_i = $photo->width;
	$h_i = $photo->height;
	print "After: $w_i, $h_i\n";

	$realc->raise('ontop', 'gelimage');
	
	@rectcoor = @rectcoor;
}

sub zoom_out
{
	my $sf = 0.5;
	my $w_i;
	my $h_i;
	my $w_f;
	my $h_f;

	$scale = $scale * $sf;

	if($scale < 1)
	{
		$scale = $scale * 2;
	}
	else
	{
		$w_i = $photo->width;
		$h_i = $photo->height;
		print "Before: $w_i, $h_i\n";

		$w_f = $sf*$w_i;
		$h_f = $sf*$h_i;
		
		$c->scale('all' => 0, 0, $sf, $sf);
		
		$photo2 = $c->Photo(-width=>$w_f, -height=>$h_f);
		$photo2->copy($photo, -shrink, -from=>0,0,$w_i,$h_i, -to=>0,0,$w_f,$h_f, -subsample=>2,2);
		$photo->configure(-width=>$w_f, -height=>$h_f);
		$photo->copy($photo2);
		$image = $c->createImage(0,0, -image=>$photo, -anchor=>'nw', -tags=>['gelimage']);	
		
		$_ *= $sf for @$box;
		$c->configure(-scrollregion=>[0,0,$w_f,$h_f]);

		$w_i = $photo->width;
		$h_i = $photo->height;
		print "After: $w_i, $h_i\n";

		$realc->raise('ontop', 'gelimage');
	}
}

sub load_gel
{
	if(defined($image))
	{
		$c->delete($image);
	}

	$photo = $c->Photo(-file=>$infilename);

	$image = $c->createImage(0,0, -image=>$photo, -anchor=>'nw', -tags=>['gelimage']);
	my $iwidth = $photo->width;
	my $iheight = $photo->height;
	print "($iheight, $iwidth)\n";
	$c->configure(-scrollregion=>[0,0,$iwidth,$iheight]);

	$realc->raise('ontop', 'gelimage');
}

sub process_selection
{
	my $sf;
	my $ctrx;
	my $ctry;
	my @colorarray;
	my $iheight;
	my $iwidth;
	my $xi;
	my $xf;
	my $yi;
	my $yf;

	print "Processing Selection\n";

	if(!defined($rectlist->curselection))
	{
		print "Error: No selection\n";
		return;
	}

	if(defined($prevselect) && $prevselect >= 0)
	{
		$rectlist->itemconfigure($prevselect, -foreground=>'#000000');
		$rectlist->itemconfigure($prevselect, -selectforeground=>'#000000');
		$c->itemconfigure($rectangles[$prevselect], -outline=>'red');
	}

	@colorarray = $photo->get(99,99);
	print "Test: @colorarray\n";
	@colorarray = $photo->get(100,99);
	print "Test: @colorarray\n";

	$iwidth = $photo->width;
	$iheight = $photo->height;
	print "Calc: $iwidth $iheight\n";

	@listselarray = $rectlist->curselection;
	$listselect = $listselarray[0];	

	print "Select: $listselect\n";

	$area = (abs($rectcoor[$listselect][2]/$rectscale[$listselect] - $rectcoor[$listselect][0]/$rectscale[$listselect]) + 1) * (abs($rectcoor[$listselect][3]/$rectscale[$listselect] - $rectcoor[$listselect][1]/$rectscale[$listselect]) + 1);

	if($rectcoor[$listselect][0] <= $rectcoor[$listselect][2])
	{
		$xi = $rectcoor[$listselect][0];
		$xf = $rectcoor[$listselect][2];
	}
	else
	{
		$xf = $rectcoor[$listselect][0];
		$xi = $rectcoor[$listselect][2];	
	}
	if($rectcoor[$listselect][1] <= $rectcoor[$listselect][3])
	{
		$yi = $rectcoor[$listselect][1];
		$yf = $rectcoor[$listselect][3];
	}
	else
	{
		$yf = $rectcoor[$listselect][1];
		$yi = $rectcoor[$listselect][3];	
	}
	
	$intensity = 0;
	for($ctrx = $xi*($scale/$rectscale[$listselect]); $ctrx <= $xf*($scale/$rectscale[$listselect]); $ctrx += $scale)
	{
		for($ctry = $yi*($scale/$rectscale[$listselect]); $ctry <= $yf*($scale/$rectscale[$listselect]); $ctry += $scale)
		{
			@colorarray = $photo->get($ctrx, $ctry);
			#print "@colorarray\n";
			if($invert == 0)
			{
				$intensity += ( (255-$colorarray[0])+(255-$colorarray[1])+(255-$colorarray[2]) ) / 3;
			}
			elsif($invert == 1)
			{
				$intensity += ($colorarray[0] + $colorarray[1] + $colorarray[2]) / 3;
			}
		}
	}
	
	$intensity_nobgcorr = $intensity;
	$intensity_display_nobgcorr = sprintf("%.1f", $intensity_nobgcorr);
	if(defined($background))
	{
		$intensity = $intensity - $area*$background;
	}
	$intensity_display = sprintf("%.1f", $intensity);

	print "Select: $listselect\n";

	$rectlist->yview($listselect);
	
	$prevselect = $listselect;
	
	$rectlist->itemconfigure($listselect, -foreground=>'forestgreen');
	$rectlist->itemconfigure($listselect, -selectforeground=>'forestgreen');
	$c->itemconfigure($rectangles[$listselect], -outline=>'forestgreen');
}

sub calc_background
{
	my $ctrx;
	my $ctry;
	my @colorarray;
	my $bgarea;
	my $xi;
	my $xf;
	my $yi;
	my $yf;
	
	$bgarea = (abs($bgcoor[2]/$bgscale - $bgcoor[0]/$bgscale) + 1) * (abs($bgcoor[3]/$bgscale - $bgcoor[1]/$bgscale) + 1);
	
	if($bgcoor[0] <= $bgcoor[2])
	{
		$xi = $bgcoor[0];
		$xf = $bgcoor[2];
	}
	else
	{
		$xf = $bgcoor[0];
		$xi = $bgcoor[2];
	}
	if($bgcoor[1] <= $bgcoor[3])
	{
		$yi = $bgcoor[1];
		$yf = $bgcoor[3];
	}
	else
	{
		$yf = $bgcoor[1];
		$yi = $bgcoor[3];
	}	
	
	$background = 0;
	for($ctrx = $xi*($scale/$bgscale); $ctrx <= $xf*($scale/$bgscale); $ctrx += $scale)
	{
		for($ctry = $yi*($scale/$bgscale); $ctry <= $yf*($scale/$bgscale); $ctry += $scale)
		{	
			@colorarray = $photo->get($ctrx, $ctry);
			#print "@colorarray\n";
			if($invert == 0)
			{
				$background += ( (255-$colorarray[0])+(255-$colorarray[1])+(255-$colorarray[2]) ) / 3;
			}
			elsif($invert == 1)
			{
				$background += ($colorarray[0] + $colorarray[1] + $colorarray[2]) / 3;
			}
		}
	}
	$background = $background / $bgarea;
	$background_display = sprintf("%.1f", $background);

	$selector1->select;
}

sub export_bands
{
	my $ctr;
	my $outline;
	
	open OUTFILE, ">$outfilename" or die "Can't open $outfilename\n";
	
	for($ctr = 0; $ctr <= $#rectdisplay; $ctr++)
	{
		$rectlist->selectionClear(0,$#rectdisplay);
		$rectlist->selectionSet($ctr);
		process_selection();
		$outline = sprintf("%5d %15d %20.2f %20.2f", $rectdisplay[$ctr], $area, $intensity_nobgcorr, $intensity);
		print OUTFILE "$outline\n";
	}
	
	close OUTFILE;
}

sub save_state
{
	my $ctr;
	
	open OUTFILE, ">$statefilename" or die "Can't open $statefilename\n";
	
	print OUTFILE "$scale\n";
	print OUTFILE "@bgcoor\n";
	
	for($ctr = 0; $ctr <= $#rectcoor; $ctr++)
	{
		print OUTFILE "@{$rectcoor[$ctr]}\n";
	}
	
	close OUTFILE;
}

sub load_state
{
	open INFILE, "$statefilename" or die "Can't open $statefilename\n";
	my $inline;
	my $ctr = 0;
	my $fontx;
	my $fonty;
	
	# Zeroth line is the scale
	$inline = <INFILE>;
	chomp($inline);
	$scale = $inline;
	
	# First line is the background
	$inline = <INFILE>;
	chomp($inline);
	@bgcoor = split ' ', $inline;
	$bgscale = $scale;
	$bg = $c->createRectangle(@bgcoor, -outline=>'blue', -tags=>['ontop']);	
	$fontx = ($bgcoor[0] + $bgcoor[2]) / 2;
	$fonty = ($bgcoor[1] + $bgcoor[3]) / 2;
	$bglabel = $c->createText($fontx, $fonty, -text=>'bg', -fill=>'blue', -font=>'Courier', -tags=>['ontop']);
	calc_background();

	# Subsequent lines are bands
	while(defined($inline=<INFILE>))
	{
		chomp($inline);
		@{$rectcoor[$ctr]} = split ' ', $inline;
		$rectangles[$ctr] = $c->createRectangle(@{$rectcoor[$ctr]}, -outline=>'red', -tags=>['ontop']);
		$rectscale[$ctr] = $scale;
		$fontx = ($rectcoor[$ctr][0] + $rectcoor[$ctr][2]) / 2;
		$fonty = ($rectcoor[$ctr][1] + $rectcoor[$ctr][3]) / 2;
		$rectlabels[$ctr] = $c->createText($fontx,$fonty, -text=>$ctr, -fill=>'red', -font=>'Courier', -tags=>['ontop']);
		$rectdisplay[$ctr] = "$ctr";
		$rectlist->selectionClear(0,$#rectdisplay);
		$rectlist->selectionSet($ctr);
		process_selection();
		$numrect++;
		
		$ctr++;
	}
}

sub delete_band
{
	@listselarray = $rectlist->curselection;
	$listselect = $listselarray[0];	

	print "DELETE: $listselect\n";

	if(!defined($listselect))
	{
		print STDERR "Error: No band selected\n";
	}
	else
	{	
		print "Clearing Selection\n";
		$rectlist->selectionClear($listselect);
	
		print "Splicing Arrays\n";
		$c->delete($rectangles[$listselect]);
		$c->delete($rectlabels[$listselect]);
		splice(@rectangles, $listselect, 1);
		splice(@rectcoor, $listselect, 1);
		splice(@rectscale, $listselect, 1);
		splice(@rectlabels, $listselect, 1);
		splice(@rectdisplay, $listselect, 1);

		print "list1 $listselect\n";
		print "num $#rectdisplay\n";

		if($listselect > $#rectdisplay)
		{
			$listselect = $#rectdisplay;
			$prevselect = $#rectdisplay;
		}

		print "list2 $listselect\n";

		if($#rectdisplay < 0)
		{
			print "NO ITEMS\n";
			$rectdisplay[0] = "---";
			$area = "";
			$intensity_display_nobgcorr = "";
			$intensity_display = "";
			$rectlist->itemconfigure('0', -foreground=>'#000000');
			$rectlist->itemconfigure('0', -selectforeground=>'#000000');
		}
		else
		{
			print "YES ITEMS\n";
			@rectdisplay = @rectdisplay;
			$rectlist->selectionSet($listselect);
			process_selection();
		}
	}
	
	print "Finished Deletion\n";
}

sub quit
{
	print "Goodbye\n";
	exit(0);
}

MainLoop();