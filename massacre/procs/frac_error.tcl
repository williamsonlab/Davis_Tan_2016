proc frac_error {a14 a15 s14 s15} {

set f [expr $a15/($a14+$a15)]

set n [expr $s14*$s14+$s15*$s15]
set d [expr ($a14+$a15)*($a14+$a15)]
set t [expr $s15*$s15/($a15*$a15)]

set e [expr sqrt($f*$f*($t+$n/$d))]

return $e

}

