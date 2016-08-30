proc display_fit_schedule { } {

global nround
global nvar
global slist
global fitlist
global cb

set bl " "
set col ":"

set w5 ".t5"

if {[winfo exists $w5] !=1} { 
	set w5 [toplevel .t5]
	wm title $w5 {Fit Schedule}
	wm geometry $w5 -20+600
	set frs [frame $w5.frs -relief groove -bd 2] 
	set nb 0
	set sbfr [frame $frs.sbfr -relief groove -bd 2]
#	label $sbfr.lab -text "Number of Rounds:"
#	entry $sbfr.ent  -textvariable nround -width 4
	button $sbfr.but -text Add_Round -default active -command {add_round}
	button $sbfr.but2 -text Default_Schedule -default active -command {default_fit_schedule;display_fit_schedule}
	button $sbfr.but3 -text Close -default active -command "destroy $w5"
	pack $sbfr.but $sbfr.but2 $sbfr.but3 -side left
	pack $sbfr -side top
	set i 1
	while {$i <= $nround} {
		set sfr [frame $frs.$i -relief groove -bd 2]
		set slab "Round "
		label $sfr.lab -text $slab$i$col -width 10
		pack $sfr.lab -side left
		set j 0
		while {$j < $nvar} {
			checkbutton $sfr.$j -text [lindex $fitlist $j] -variable cb($i$j)
			pack $sfr.$j -side left
			incr j
		}
		button $sfr.but -text "Delete" -default active -command "delete_round $i"
		pack $sfr.but -side left
		pack $sfr -side top
		incr i
	}
	pack $frs
}
}
