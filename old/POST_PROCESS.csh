#!/bin/csh

#Created: 2 Nov 2010, Dan Chavas

#Purpose: This script does post-processing.

#Syntax: './POST_PROCESS.csh [new subdirectory name]'

### USER INPUT #########################
set JOB_NAME = most_recent_run  #is name of input subdirectory in CM1_output/
########################################

#switch to output directory
cd '/home/drchavas/CM1/v15/CM1_output/'$JOB_NAME

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
cp ./cm1out_t*.nc ./cm1out_stats.nc ./namelist.input ./input_sounding ./filecontent.info ./ic.info ./TRANSFER/.

#rename directory
cd ..
mv $JOB_NAME $1

#transfer files to zapata
#cd ..
#scp -r ./${JOB_NAME}/cm1out_t*.nc ./${JOB_NAME}/namelist.input ./${JOB_NAME}/filecontent.info drchavas@zapata.mit.edu:"/media/Iomega\ HDD/Research/CM1/pgi_output/${JOB_NAME}"
