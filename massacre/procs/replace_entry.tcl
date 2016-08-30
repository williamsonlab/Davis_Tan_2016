proc replace_entry { } {

global w7 lb7 tab_list blist table

set j [$lb7 curselection]

set  bj [lindex $blist $j]

set val1 $table($j,0)
set val2 $table($j,1)
set val3 $table($j,2)

set olist [list $val1 $val2 $val3]

puts ""
puts $bj
puts $olist

}
