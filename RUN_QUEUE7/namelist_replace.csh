#!/bin/csh

## INITIAL CONDITIONS ##########
set nx          =  $1
set dx          =  $2
set dy          =  $3
set nz          =  $4
set dz          =  $5
set dtl         =  $6
set l_h         =  $7
set xhd         =  $8

set l_v         =  $9
set timax       =  ${10}
set tapfrq      =  ${11}
set rstfrq      =  ${12}
set irst        =  ${13}
set rstnum      =  ${14}
set rterm       =  ${15}
set irandp      =  ${16}
set imove       =  ${17}
set umove       =  ${18}
set vmove       =  ${19}
set radopt      =  ${20}
set tsk0        =  ${21}
set tmn0        =  ${22}
set cecd        =  ${23}
set cnstce      =  ${24}
set cnstcd      =  ${25}
set fcor        =  ${26}
set zd          =  ${27}
################################################################################

#adjust namelist.input parameters in namelist.input_drc

sed 's| nx .*=.*\,| nx            =     '"${nx}"'\,|g' namelist.input > temp_file
mv temp_file namelist.input
sed 's| dx .*=.*\,| dx      =    '"${dx}"'\,|g' namelist.input > temp_file
mv temp_file namelist.input
sed 's| dy .*=.*\,| dy      =    '"${dy}"'\,|g' namelist.input > temp_file
mv temp_file namelist.input
sed 's| nz .*=.*\,| nz            =     '"${nz}"'\,|g' namelist.input > temp_file
mv temp_file namelist.input
sed 's| dz .*=.*\,| dz      =    '"${dz}"'\,|g' namelist.input > temp_file
mv temp_file namelist.input
sed 's| dtl .*=.*\,| dtl     =      '"${dtl}"'\,|g' namelist.input > temp_file
mv temp_file namelist.input
sed 's| l_h .*=.*\,| l_h      =    '"${l_h}"'\,|g' namelist.input > temp_file
mv temp_file namelist.input
sed 's| l_v .*=.*\,| l_v      =    '"${l_v}"'\,|g' namelist.input > temp_file
mv temp_file namelist.input
sed 's| xhd .*=.*\,| xhd      =    '"${xhd}"'\,|g' namelist.input > temp_file
mv temp_file namelist.input
sed 's| timax .*=.*\,| timax   = '"${timax}"'\,|g' namelist.input > temp_file
mv temp_file namelist.input
sed 's| tapfrq .*=.*\,| tapfrq  =   '"${tapfrq}"'\,|g' namelist.input > temp_file
mv temp_file namelist.input
sed 's| rstfrq .*=.*\,| rstfrq  = '"${rstfrq}"'\,|g' namelist.input > temp_file
mv temp_file namelist.input
sed 's| irst .*=.*\,| irst       =  '"${irst}"'\,|g' namelist.input > temp_file
mv temp_file namelist.input
sed 's| rstnum .*=.*\,| rstnum     =  '"${rstnum}"'\,|g' namelist.input > temp_file
mv temp_file namelist.input
sed 's| rterm .*=.*\,| rterm      =  '"${rterm}"'\,|g' namelist.input > temp_file
mv temp_file namelist.input
sed 's| irandp .*=.*\,| irandp     =  '"${irandp}"'\,|g' namelist.input > temp_file
mv temp_file namelist.input
sed 's| imove .*=.*\,| imove      =  '"${imove}"'\,|g' namelist.input > temp_file
mv temp_file namelist.input
sed 's| umove .*=.*\,| umove    =  '"${umove}"'\,|g' namelist.input > temp_file
mv temp_file namelist.input
sed 's| vmove .*=.*\,| vmove    =  '"${vmove}"'\,|g' namelist.input > temp_file
mv temp_file namelist.input
sed 's| radopt .*=.*\,| radopt   =        '"${radopt}"'\,|g' namelist.input > temp_file
mv temp_file namelist.input
sed 's| tsk0 .*=.*\,| tsk0        = '"${tsk0}"'\,|g' namelist.input > temp_file
mv temp_file namelist.input
sed 's| tmn0 .*=.*\,| tmn0        = '"${tmn0}"'\,|g' namelist.input > temp_file
mv temp_file namelist.input
sed 's| cecd .*=.*\,| cecd        = '"${cecd}"'\,|g' namelist.input > temp_file
mv temp_file namelist.input
sed 's| cnstce .*=.*\,| cnstce        = '"${cnstce}"'\,|g' namelist.input > temp_file
mv temp_file namelist.input
sed 's| cnstcd .*=.*\,| cnstcd        = '"${cnstcd}"'\,|g' namelist.input > temp_file
mv temp_file namelist.input
sed 's| fcor .*=.*\,| fcor        = '"${fcor}"'\,|g' namelist.input > temp_file
mv temp_file namelist.input
sed 's| zd .*=.*\,| zd        = '"${zd}"'\,|g' namelist.input > temp_file
mv temp_file namelist.input

#copy updated namelist.input_drc file to namelist.input
#cp namelist.input_drc namelist.input


