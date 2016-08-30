proc cur_par { } {

global fitlist
global fitvars

set i 0
foreach var $fitlist {
	set cur [lindex $fitvars $i]
	global $cur
	global $var
	set $var [expr $$cur]
	incr i
}
}
