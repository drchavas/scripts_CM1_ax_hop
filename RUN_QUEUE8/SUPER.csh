#!/bin/csh

#########################################################################################
## USER INPUT: INITIAL CONDITIONS #######################################################
set runnum = 8  #IMPORTANT!  This is the number of the directory where everything will be copied and saved in output

set r0s          = ('400000.0') #400000.0 outer radius of vortex, where V=0 [m]
set r0drmaxs     = ('5.0')      #5.0 r0/rmax [-]
set vmaxs        = ('0.0')      #12.5 max wind speed in vortex [m/s]
set zcs          = ('4375.0')   #4375.0 height of center of vortex above ground [m]; NOTE:dz=1.25km, z0=.625km
set dz0s         = ('2875.0')   #2875.0 vertical scale of vortex edge from center [m]
set rc_qvs       = ('0.0')      #0.0 radius of center of anomaly [m]

set r0_qvs       = ('200000.0') #200000.0 outer radius of anomaly [m]
set rhpert_maxs  = ('1.9')      #0.3 max relative humidity perturbation (%); 0.01455 kg/kg --> RH=80% at lowest level
set zc_qvs       = ('4375.0')   #4375.0 height of center of bubble above ground [m]
set dz0_qvs      = ('5000.0')   #2875.0 vertical scale of bubble edge from center [m]
set shape_qvs    = ('2')        #2 1=gaussian, 2=constant

set z_bltop      = '1500.0'     #1500.0 CONSTANT height of boundary layer top [m] (initial perturbations = 0 below this height)
set fcor         = '.0001'     #.00005 CONSTANT coriolis parameter [s-1]
set sfcphys      = '2'          #2 1=original; 2=gust+u'; 21=gust+u'+ring modification; 3=gust only; 4=azimuthal mean magnitude (AX ONLY)
 set U_gust       = '3'         #4.15 gustiness added to surface flux wind speed [m/s]; only used if sfcphys = 2 or 3
 set r_inner      = '900'       #radius [km] of inner edge of ring
 set r_outer      = '1000'      #radius [km] of outer edge of ring
set radtype      = '3'          #3 1=RE87 relaxation; 2=interactive; 3=constant cooling (DRC 05-24-11)
 set dTdt_halfday = '-0.5'      #-0.5 radiative HEATING rate [K/half-day]; only used if radtype = 3
 set T_tpp = '150'              #200; [K] temperature below which newtonian relaxation kicks in
##########################################################################################
## namelist.input PARAMETERS ##########
set l_h         =  12000.0       #1500.0

set dx          =  4000.0       #4000.0
set dy          =  4000.0       #4000.0
set dtl         =  8.0         #16.0 -- DOUBLE/HALVE WITH dx/dy!
set xhd         =  400000.0     #100000.0 -- DOUBLE/HALVE WITH dx/dy!

set nx          =  3072          #384 -- DO NOT DOUBLE/HALVE!
set nz          =  48           #40
set zd          =  25000.00     #20000.00
set dz          =  625.0        #625.0
set l_v         =  100.0        #100.0
set timax       =  12960000.0    #8640000.0
set tapfrq      =  21600.0      #21600.0
set rstfrq      =  86400.0      #86400.0
set irst        =  1            #0
set rstnum      =  72            #1
set irandp      =  1            #1
set imove       =  0            #0
set umove       =  5.0          #5.0
set vmove       =  0.0          #0.0
set tsk0        =  300.00       #299.28
set tmn0        =  297.28       #297.28
set cecd        =  1            #1
set cnstce      =  .0015        #.0015
set cnstcd      =  .0015        #.0015
################################################################################

##########################################################################################

#Control parameter values (for purpose of naming output directories)
set CTRL_r0          = '400000.0'  # outer radius of vortex, where V=0 [km]
set CTRL_r0drmax     = '5.0'  # r0/rmax [-]
set CTRL_vmax        = '12.5'  # max wind speed in vortex [m/s]
set CTRL_zc          = '4375.0'  # height of center of vortex above ground [m]
set CTRL_dz0         = '2875.0'  # vertical scale of vortex edge from center [m]
set CTRL_rc_qv       = '0.0'  # radius of center of anomaly [km]
set CTRL_r0_qv       = '200000.0'  # outer radius of anomaly [km]
set CTRL_rhpert_max  = '0.3'  # max relative humidity perturbation (%); 0.01455 kg/kg --> RH=80% at lowest level
set CTRL_zc_qv       = '4375.0'  # height of center of bubble above ground [m]
set CTRL_dz0_qv      = '2875.0'  # vertical scale of bubble edge from center [m]
set CTRL_shape_qv    = '2'    # 1=gaussian, 2=constant

##########################################################################################
##########################################################################################

#Determine settings for radiation scheme and newtonian relaxation term
if ($radtype == 2) then
  set radopt   = 1
  set rterm    = 0
else
  set radopt   = 0
  set rterm    = 1
endif

set SCRATCH = "/data/drchavas/"

cd ~/scripts_CM1_ax_hop/

#define some variables
set model = CM1
echo $model
set version = v15
echo $version


#loop over everything
foreach r0 (${r0s})
foreach r0drmax (${r0drmaxs})
foreach vmax (${vmaxs})
foreach zc (${zcs})
foreach dz0 (${dz0s})
foreach rc_qv (${rc_qvs})
foreach r0_qv (${r0_qvs})
foreach rhpert_max (${rhpert_maxs})
foreach zc_qv (${zc_qvs})
foreach dz0_qv (${dz0_qvs})
foreach shape_qv (${shape_qvs})


## Keep line of sub-directories with model input files

cd ~/scripts_CM1_ax_hop/


set files_dir = '~/scripts_CM1_ax_hop/RUN_QUEUE'
set files_dir = "${files_dir}${runnum}"
if (! -d ${files_dir}) then
  mkdir ${files_dir}	#make the sub-directory into which CM1 input files will be saved
else
  echo 'ALREADY A RUN WITH THAT QUEUE NUMBER!'
  set files_dir = 'FAIL'
endif

#name output simulation directory
set outdir = 'CTRL'

if ($r0 != $CTRL_r0) then
  set temp = `echo "$r0"| cut -d . -f 1`	#take the integer only
  set outdir = "${outdir}ro${temp}"
else
endif
if ($r0drmax != $CTRL_r0drmax) then
  set temp = `echo "$r0drmax"| cut -d . -f 1`	#take the integer only
  set outdir = "${outdir}rodrm${temp}"
else
endif
if ($vmax != $CTRL_vmax) then
  set temp = `echo "$vmax"| cut -d . -f 1`	#take the integer only
  set outdir = "${outdir}v${temp}"
else
endif
if ($zc != $CTRL_zc) then
  set temp = `echo "$zc"| cut -d . -f 1`	#take the integer only
  set outdir = "${outdir}zc${temp}"
else
endif
if ($dz0 != $CTRL_dz0) then
  set temp = `echo "$dz0"| cut -d . -f 1`	#take the integer only
  set outdir = "${outdir}dz${temp}"
else
endif
if ($rc_qv != $CTRL_rc_qv) then
  set temp = `echo "$rc_qv"| cut -d . -f 1`	#take the integer only
  set outdir = "${outdir}qrc${temp}"
else
endif
if ($r0_qv != $CTRL_r0_qv) then
  set temp = `echo "$r0_qv"| cut -d . -f 1`	#take the integer only
  set outdir = "${outdir}qro${temp}"
else
endif
if ($rhpert_max != $CTRL_rhpert_max) then
  set temp = `echo "$rhpert_max"| cut -d . -f 2`	#take the decimal only
  set outdir = "${outdir}qrh${temp}"
else
endif
if ($zc_qv != $CTRL_zc_qv) then
  set temp = `echo "$zc_qv"| cut -d . -f 1`	#take the integer only
  set outdir = "${outdir}qzc${temp}"
else
endif
if ($dz0_qv != $CTRL_dz0_qv) then
  set temp = `echo "$dz0_qv"| cut -d . -f 1`	#take the integer only
  set outdir = "${outdir}qdz${temp}"
else
endif
if ($shape_qv != $CTRL_shape_qv) then
  set temp = `echo "$shape_qv"| cut -d . -f 1`	#take the integer only
  set outdir = "${outdir}qsh${temp}"
else
endif

echo $outdir > run.name

###########################################################
## SET UP SIMULATION ######################################

cd ~/scripts_CM1_ax_hop/

##Copy input/model files to be modified into sub-directory
cp namelist.input ${files_dir}/.
cp input_sounding ${files_dir}/.
cp vmax_RE.csh ${files_dir}/.
cp -r vmax ${files_dir}/.
cp namelist_replace.csh ${files_dir}/.
cp run_CM1 ${files_dir}/run_CM1_${runnum}

##Copy these as well for reference
cp SUPER.csh ${files_dir}/.
mv run.name ${files_dir}/.

#UPDATE run_CM1_${runnum} to have the same runnum within
cd ${files_dir}

sed 's|set runnum = 1 |set runnum = '"${runnum}"' |g' run_CM1_${runnum} > temp_file
mv temp_file run_CM1_${runnum}

#UPDATE namelist.input
echo "Updating namelist.input"
./namelist_replace.csh ${nx} ${dx} ${dy} ${nz} ${dz} ${dtl} ${l_h} ${xhd} ${l_v} ${timax} ${tapfrq} ${rstfrq} ${irst} ${rstnum} ${rterm} ${irandp} ${imove} ${umove} ${vmove} ${radopt} ${tsk0} ${tmn0} ${cecd} ${cnstce} ${cnstcd} ${fcor} ${zd}


## Calculate inputs for RE87 initial vortex that will give desired vmax and rmax
cd ${files_dir}
./vmax_RE.csh ${vmax} ${r0} ${r0drmax} ${fcor} > vm_RE.out
set rmax_RE = `head -1 vm_RE.out`
set vmax_RE = `tail -1 vm_RE.out`
rm vm_RE.out vmax_RE.csh
rm -r vmax

## Adjust initial conditions in copied files
echo "Adjusting initial conditions in init3d.F, solve.F, sfcphys.F" 

cd ~/scripts_CM1_ax_hop/CM1files_drc

cp init3d_drc_backup.F init3d_drc.F
cp init3d_drc.F ${files_dir}/.
cd ${files_dir}

sed 's|r0 .*=.* !|r0 = '"${r0}"' !|g' init3d_drc.F > temp_file
mv temp_file init3d_drc.F
sed 's|rmax .*=.* !|rmax = '"${rmax_RE}"' !|g' init3d_drc.F > temp_file
mv temp_file init3d_drc.F
sed 's|vmax .*=.* !|vmax = '"${vmax_RE}"' !|g' init3d_drc.F > temp_file
mv temp_file init3d_drc.F
sed 's|zc .*=.* !|zc = '"${zc}"' !|g' init3d_drc.F > temp_file
mv temp_file init3d_drc.F
sed 's|dz0 .*=.* !|dz0 = '"${dz0}"' !|g' init3d_drc.F > temp_file
mv temp_file init3d_drc.F
sed 's|rc_qv .*=.* !|rc_qv = '"${rc_qv}"' !|g' init3d_drc.F > temp_file
mv temp_file init3d_drc.F
sed 's|r0_qv .*=.* !|r0_qv = '"${r0_qv}"' !|g' init3d_drc.F > temp_file
mv temp_file init3d_drc.F
sed 's|rhpert_max .*=.* !|rhpert_max = '"${rhpert_max}"' !|g' init3d_drc.F > temp_file
mv temp_file init3d_drc.F
sed 's|zc_qv .*=.* !|zc_qv = '"${zc_qv}"' !|g' init3d_drc.F > temp_file
mv temp_file init3d_drc.F
sed 's|dz0_qv .*=.* !|dz0_qv = '"${dz0_qv}"' !|g' init3d_drc.F > temp_file
mv temp_file init3d_drc.F
sed 's|shape_qv .*=.* !|shape_qv = '"${shape_qv}"' !|g' init3d_drc.F > temp_file
mv temp_file init3d_drc.F
sed 's|z_bltop .*=.* !|z_bltop = '"${z_bltop}"' !|g' init3d_drc.F > temp_file
mv temp_file init3d_drc.F

## Copy updated init3d file
mv init3d_drc.F init3d.F


## Copy updated sfcphys file into files directory
cd ~/scripts_CM1_ax_hop/CM1files_drc

if ($sfcphys == 1) then #original sfc flux formulation
   cp sfcphys_orig.F ${files_dir}/sfcphys.F
else if ($sfcphys == 2) then #sfc gustiness + (umove,vmove)
  cp sfcphys_DRCgust_backup.F sfcphys_DRCgust.F
  cp sfcphys_DRCgust.F ${files_dir}/.

  cd ${files_dir}

  sed 's|w1(i,j)=.* +|w1(i,j)= '"${U_gust}"' +|g' sfcphys_DRCgust.F > temp_file
  mv temp_file sfcphys_DRCgust.F

  mv sfcphys_DRCgust.F sfcphys.F
else if ($sfcphys == 222) then #sfc gustiness + (umove,vmove) in BOTH enthalpy AND drag formulae
  cp sfcphys_DRCgust_drag_backup.F sfcphys_DRCgust_drag.F
  cp sfcphys_DRCgust_drag.F ${files_dir}/.

  cd ${files_dir}

  sed 's|w1(i,j)=.* +|w1(i,j)= '"${U_gust}"' +|g' sfcphys_DRCgust_drag.F > temp_file
  mv temp_file sfcphys_DRCgust_drag.F
  sed 's|wspd =.* +|wspd = '"${U_gust}"' +|g' sfcphys_DRCgust_drag.F > temp_file
  mv temp_file sfcphys_DRCgust_drag.F

  mv sfcphys_DRCgust_drag.F sfcphys.F
else if ($sfcphys == 21) then #sfc gustiness + (umove,vmove) + reduction in ring
  cp sfcphys_DRCgust_ringmod_backup.F sfcphys_DRCgust_ringmod.F
  cp sfcphys_DRCgust_ringmod.F ${files_dir}/.

  cd ${files_dir}

  sed 's|w1(i,j)=.* +|w1(i,j)= '"${U_gust}"' +|g' sfcphys_DRCgust_ringmod.F > temp_file
  mv temp_file sfcphys_DRCgust_ringmod.F
  sed 's|r_inner =.*|r_inner = '"${r_inner}"' |g' sfcphys_DRCgust_ringmod.F > temp_file
  mv temp_file sfcphys_DRCgust_ringmod.F
  sed 's|r_outer =.*|r_outer = '"${r_outer}"' |g' sfcphys_DRCgust_ringmod.F > temp_file
  mv temp_file sfcphys_DRCgust_ringmod.F

  mv sfcphys_DRCgust_ringmod.F sfcphys.F
else if ($sfcphys == 3) then #sfc gustiness only (i.e. no WISHE)
  cp sfcphys_DRCgustonly_backup.F sfcphys_DRCgustonly.F
  cp sfcphys_DRCgustonly.F ${files_dir}/.

  cd ${files_dir}

  sed 's|w1(i,j)=.* +|w1(i,j)= '"${U_gust}"' +|g' sfcphys_DRCgustonly.F > temp_file
  mv temp_file sfcphys_DRCgustonly.F

  mv sfcphys_DRCgustonly.F sfcphys.F
else	#sfc flux wind magnitudes increased by constant value (check file!)
   cp sfcphys_DRCax.F ${files_dir}/sfcphys.F
endif

## Copy updated solve.F, which contains radiation scheme, into files directory
#if using constant cooling, update solve.F with input value for dTdt_halfday
cd ~/scripts_CM1_ax_hop/CM1files_drc

if ($radtype == 3) then 
   cp solve_DRCradconst_backup.F solve_DRCradconst.F
   cp solve_DRCradconst.F ${files_dir}/.

   cd ${files_dir}

   sed 's|thrad_drc =.* |thrad_drc = '"${dTdt_halfday}"' / |g' solve_DRCradconst.F > temp_file
   mv temp_file solve_DRCradconst.F
   sed 's|T_tpp =.*|T_tpp = '"${T_tpp}"'|g' solve_DRCradconst.F > temp_file
   mv temp_file solve_DRCradconst.F


   echo "Selecting input_sounding file"
   rm input_sounding
   set radstr = `echo $dTdt_halfday | cut -c2-10`
   
   if(${radstr} == 0.5) then
      if (-e ~/scripts_CM1_ax_hop/sounding_files/input_sounding_3dRCE_nx48_SST${tsk0}K_Tthresh${T_tpp}K_usfc${U_gust}) then
        cp ~/scripts_CM1_ax_hop/sounding_files/input_sounding_3dRCE_nx48_SST${tsk0}K_Tthresh${T_tpp}K_usfc${U_gust} input_sounding
        echo input_sounding_3dRCE_nx48_SST${tsk0}K_Tthresh${T_tpp}K_usfc${U_gust} 
      else
        cp ~/scripts_CM1_ax_hop/sounding_files/input_sounding_RE87_T200K_orig input_sounding
        echo 'DID NOT FIND SOUNDING FILE ASSOCIATED WITH GIVEN PARAMETERS!  USING ORIGINAL FILE'
      endif
   else
      if (-e ~/scripts_CM1_ax_hop/sounding_files/input_sounding_3dRCE_nx48_SST${tsk0}K_Tthresh${T_tpp}K_usfc${U_gust}_rad${radstr}K) then
        cp ~/scripts_CM1_ax_hop/sounding_files/input_sounding_3dRCE_nx48_SST${tsk0}K_Tthresh${T_tpp}K_usfc${U_gust}_rad${radstr}K input_sounding
        echo input_sounding_3dRCE_nx48_SST${tsk0}K_Tthresh${T_tpp}K_usfc${U_gust}_rad${radstr}K
      else
        cp ~/scripts_CM1_ax_hop/sounding_files/input_sounding_RE87_T200K_orig input_sounding
        echo 'DID NOT FIND SOUNDING FILE ASSOCIATED WITH GIVEN PARAMETERS!  USING ORIGINAL FILE'
      endif
   endif

else if ($radtype == 2) then
  
   cd ${files_dir}

   if (-e ~/scripts_CM1_ax_hop/sounding_files/input_sounding_3dRCE_nx48_SST${tsk0}K_radfull_usfc${U_gust}) then
     cp ~/scripts_CM1_ax_hop/sounding_files/input_sounding_3dRCE_nx48_SST${tsk0}K_radfull_usfc${U_gust} input_sounding
   else
     cp ~/scripts_CM1_ax_hop/sounding_files/input_sounding_RE87_T200K_orig input_sounding
     echo 'DID NOT FIND SOUNDING FILE ASSOCIATED WITH GIVEN PARAMETERS!  USING ORIGINAL FILE'
   endif

endif

if ($radtype <= 2) then #original solve.F (either newtonian relax or interactive rad)
   cp ~/scripts_CM1_ax_hop/CM1files_drc/solve_orig.F solve.F
else  #constant radiational cooling
   mv solve_DRCradconst.F solve.F
endif

cd ${files_dir}

echo "Initialization complete, submitting model to queue!" 

## Run the PGI program #########
#cd ${files_dir}
cd ~/scripts_CM1_ax_hop/
qsub -T hyperthread ${files_dir}/run_CM1_${runnum}


end
end
end
end
end
end
end
end
end
end
end
