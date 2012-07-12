#!/bin/csh

set JOB_NAME = most_recent_run

#define some variables
set model = CM1
echo $model
set version = v15
echo $version

#Go to output CM1 directory
cd ~/${model}/${version}/${model}_pgi/run
echo -n "I am currently in directory "; pwd



############################################################################
#Post-Processing

#extract output directory name from file outdir.info
cd ~/scripts_CM1/
set outdir = `less outdir.info`
rm outdir.info

cd ~/${model}/${version}/

#move output to archive directory
mkdir ./${model}_output/$JOB_NAME

mv ./${model}_pgi/run/cm1out* ./${model}_pgi/run/proc* ./${model}_output/$JOB_NAME/.

cp ./${model}_pgi/run/namelist.input ./${model}_output/$JOB_NAME/.

cp ./${model}_pgi/run/input_sounding ./${model}_output/$JOB_NAME/.

less ./${model}_pgi/src/init3d.F | grep '\!Initial vortex' -A 7 > ./${model}_output/$JOB_NAME/ic.info
echo '' >> ./${model}_output/$JOB_NAME/ic.info
less ./${model}_pgi/src/init3d.F | grep '\!Initial relative humidity anomaly' -A 7 >> ./${model}_output/$JOB_NAME/ic.info

#Go to output directory
cd ./${model}_output/$JOB_NAME
echo -n "I am currently in directory "; pwd

#output text file listing all differences between namelist.input and namelist.input_CTRL
diff namelist.input ../namelist.input_CTRL > params_diff.info


## Code from POST_PROCESS.csh ###############################

#DETERMINE NUMBER OF TILES
set num_tiles = `ls cm1out_????_0001.* | wc -l`

#DETERMINE NUMBER OF TIMES
set num_times = `ls cm1out_0001_????.* | wc -l`

#copy in the fortran-generated executable that combines netcdf tile files into single file for each time (these must be in same dir as files to be run)
cp ../combine_r12_template_DRC.exe combine_r12.exe

#add '.nc' suffix to end of each output netcdf file for each time
@ t = 1
while( $t <= $num_times)	#loop over times
echo $t

	@ i = 0
	if ( $t < 10 ) then
		while ( $i < $num_tiles )	#loop over tiles
			if ( $i < 10 ) then
				mv cm1out_000${i}_000${t}.cdf cm1out_000${i}_000${t}.nc
			else
				mv cm1out_00${i}_000${t}.cdf cm1out_00${i}_000${t}.nc
			endif
			@ i = $i + 1
		end
		set OUTPUT_FILE = 'cm1out_t000'${t}'.nc'
	else if ( $t < 100 ) then
		while ( $i < $num_tiles )
			if ( $i < 10 ) then
				mv cm1out_000${i}_00${t}.cdf cm1out_000${i}_00${t}.nc
			else
				mv cm1out_00${i}_00${t}.cdf cm1out_00${i}_00${t}.nc
			endif
			@ i = $i + 1
		end
		set OUTPUT_FILE = 'cm1out_t00'${t}'.nc'
	else if ( $t < 1000 ) then
		while ( $i < $num_tiles )
			if ( $i < 10 ) then
				mv cm1out_000${i}_0${t}.cdf cm1out_000${i}_0${t}.nc
			else
				mv cm1out_00${i}_0${t}.cdf cm1out_00${i}_0${t}.nc
			endif
			@ i = $i + 1
		end
		set OUTPUT_FILE = 'cm1out_t0'${t}'.nc'
	else
		while ( $i < $num_tiles )
			if ( $i < 10 ) then
				mv cm1out_000${i}_${t}.cdf cm1out_000${i}_${t}.nc
			else
				mv cm1out_00${i}_${t}.cdf cm1out_00${i}_${t}.nc
			endif
			@ i = $i + 1
		end
		set OUTPUT_FILE = 'cm1out_t'${t}'.nc'
	endif

#run combine_r12.exe to combine the files: can do one time each execution
#echo ${num_tiles}
./combine_r12.exe << EOF
${num_tiles}
${t}
${OUTPUT_FILE}

EOF

	@ t = $t + 1

end	#time

rm combine_r12.exe

#output text file listing variables in files
module add pgi
module add netcdf
ncdump -h cm1out_t0001.nc > filecontent.info

#move files for transfer to new folder "TRANSFER"
mkdir ./TRANSFER
mv ./cm1out_t*.nc ./cm1out_stats.nc ./params_diff.info ./namelist.input ./input_sounding ./filecontent.info ./ic.info ./TRANSFER/.
rm *.*

#rename directory
cd ..
mv $JOB_NAME /data/drchavas/CM1_output/$outdir

