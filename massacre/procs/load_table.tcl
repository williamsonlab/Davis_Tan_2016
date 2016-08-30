proc load_table { } {

global batchfile 
global w7
global lb7
global blist
global mfile
global tab_list
global table


if {[winfo exists $w7] == 0} {browse_batch}

set mop [open $mfile r]

set line [gets $mop]
set tab_list [split $line ","]

set i 0

while {[eof $mop] == 0} {
	set line [gets $mop]
	set val_list [split $line ","]
	set j 0
	foreach var $tab_list {
		set table($i,$j) [lindex $val_list $j]
		incr j
	}
	incr i
}
}

		
	
