proc delete_round { i } {
global nround
global nvar
global cb

set bl " "
for {set j [expr $i + 1]} {$j <= $nround} {incr j} {
	set k [expr $j - 1]
	puts $k$j
	for {set m 0} {$m < $nvar} {incr m} {set cb($k$m) $cb($j$m)}	
}
for {set j 0} {$j < $nvar} {incr j} {set cb($nround$j) 0}
set nround [expr $nround - 1]
set w5 ".t5"
destroy $w5
display_fit_schedule

}
