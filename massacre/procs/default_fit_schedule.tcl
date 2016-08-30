proc default_fit_schedule { } {
global nround
global nvar
global fit_mod
global fs_lab_list
global cb

if {$fit_mod == 1 } {
	set nround 4
	set nvar 4
	set fs_lab_list {B UL_AMP OFF GW}
	set slist(1) {0 1 0 0}
	set slist(2) {0 1 0 1}
	set slist(3) {0 1 1 1}
	set slist(4) {1 1 1 1}
}

if {$fit_mod == 2 } {
	set nround 5
	set nvar 5
	set fs_lab_list {B LAB_AMP OFF GW FRAC}
	set slist(1) {0 1 0 0 0}
	set slist(2) {0 0 0 0 1}
	set slist(3) {0 1 0 1 1}
	set slist(4) {0 0 1 1 1}
	set slist(5) {1 1 1 1 1}
}


if {$fit_mod == 3 } {
	set nround 5
	set nvar 6
	set fs_lab_list {B UL_AMP LAB_AMP OFF GW FRAC}
	set slist(1) {0 1 1 0 0 0}
	set slist(2) {0 0 0 0 0 1}
	set slist(3) {0 1 1 0 1 1}
	set slist(4) {0 0 0 1 1 1}
	set slist(5) {1 1 1 1 1 1}
}

for {set i 1} {$i <= $nround} {incr i} {
	for {set j 0} {$j < $nvar} {incr j} {
		set cb($i$j) [lindex $slist($i) $j]
	}
}
set w5 ".t5"
if {[winfo exists $w5] ==1} { destroy $w5 }

#display_fit_schedule

}
