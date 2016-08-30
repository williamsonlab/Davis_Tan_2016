#!/bin/csh
#
~/bin/diffmap/diffmap.exe << eof
M
./00_source_maps/map1.mrc
./00_source_maps/map2.mrc
1.45
./00_source_maps/map1-map2.mrc
./00_source_maps/map2_scaled.mrc
eof
#
