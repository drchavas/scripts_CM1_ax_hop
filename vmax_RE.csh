#!/bin/csh

cd vmax

#pgf90 -o vmax.exe vmax_RE.F -Mfree

./vmax.exe << EOF
$1
$2
$3
$4

EOF

cd ..


