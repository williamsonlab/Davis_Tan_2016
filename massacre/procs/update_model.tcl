proc update_model { } {

global fitlist
global fit_mod
global var_list
global nvar

if {$fit_mod == 1} {
	set fitlist {b ul_amp off gw}
}
if {$fit_mod == 2} {
	set fitlist {b lab_amp off gw frac}
}
if {$fit_mod == 3} {
	set fitlist {b ul_amp lab_amp off gw frac}
}
set nvar [llength $fitlist]
default_fit_schedule
if {[winfo exists ".t4"] == 1} {destroy ".t4";display_fit_variables}
if {[winfo exists ".t5"] == 1} {destroy ".t5";display_fit_schedule}

}
