proc browse_batch { } {

global batchfile
global blist
global w7
global lb7

set bop [open $batchfile r]

set w7 ".t7"
if {[winfo exists $w7] == 1} { destroy $w7 }

set w7 [toplevel .t7]
wm title $w7 {Browse Batch File}
wm geometry $w7 +144+144
set f7 [frame $w7.fr]
set lb7 [listbox $f7.lb -width 80 -height 40 -yscrollcommand "$f7.sb set"]
scrollbar $f7.sb -command "$f7.lb yview"
button $f7.b1 -text Select -default active -command {set_input [lindex $blist [$lb7 curselection]]}
button $f7.b2 -text Close -default active -command {destroy $w7}
pack $f7.lb -side left
pack $f7.sb -side left -fill y
pack $f7.b1 -side top
pack $f7.b2 -side top
pack $f7 -side left


set i 0
set blist {}
while {[eof $bop] == 0} {
	set line [gets $bop]
	lappend blist $line
	$f7.lb insert $i $line
	incr i
}
}

