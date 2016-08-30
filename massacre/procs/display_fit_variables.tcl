proc display_fit_variables { } {

##### Fit Variables Window

global fitlist
global var_list
global val_list
global lab_list
global fit_mod

set w4 ".t4"
if {[winfo exists $w4] !=1} { 
	set w4 [toplevel .t4]
	wm title $w4 {Initial Fit Variables}
	wm geometry $w4 -20+600
	set frv [frame $w4.frv] 
	set nb 0
	foreach f $fitlist {
		set var $f
		global $var
		set val [lindex $val_list [lsearch $var_list $var]]
		set lab [lindex $lab_list [lsearch $var_list $var]]
		set vbl [frame $frv.$nb]
		label $vbl.lab -text $lab -width 15 -justify left
		entry $vbl.ent -textvariable $var -width 12
		pack $vbl.lab  -side left -anchor e -expand 10
		pack $vbl.ent -side left -anchor w -padx 10 -fill x
		pack $vbl -side top
#		set $var $val
		incr nb
	}
	pack $frv
	set frvb [frame $w4.frvb]
	button $frvb.b1 -text "Close" -default active -command "destroy $w4"
	pack $frvb.b1
	pack $frvb
}
}

