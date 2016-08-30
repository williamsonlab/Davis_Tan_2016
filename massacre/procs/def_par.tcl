proc def_par { } {
global fitlist
global fit_mod


if {$fit_mod == 1 } {set def {0.0 1.0 0.0003 0.0}}
if {$fit_mod == 2 } {set def {0.0 1.0 0.0 0.0003 0.5}}
if {$fit_mod == 3 } {set def {0.0 1.0 1.0 0.0 0.0003 0.5}}

set i 0
foreach var $fitlist {
	global $var
	set $var [lindex $def $i]
	incr i
}
}	
