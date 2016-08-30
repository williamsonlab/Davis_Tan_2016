proc add_round { } {
global nround
global nvar
global cb

incr nround
for {set j 0} {$j < $nvar} {incr j} {set cb($nround$j) 0}
set w5 ".t5"
destroy $w5
display_fit_schedule
}
