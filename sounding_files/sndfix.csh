#!/bin/csh

set ssts = (295 300 305)
set tthreshs = (150 200 250)
set usfcs = (1 3 5)
#set rads = (.25 .5 2 4)

foreach sst (${ssts})
  foreach tthresh (${tthreshs})
    foreach usfc (${usfcs})

#input_sounding_3dRCE_nx48_SST${sst}.00K_Tthresh${tthresh}K_usfc${usfc}_rad${rad}K 
#input_sounding_3dRCE_nx48_SST${sst}.00K_Tthresh${tthresh}K_usfc${usfc}

echo '   40000.00000    1064.6500       0.00000      0.00000     0.00000' >> input_sounding_3dRCE_nx48_SST${sst}.00K_Tthresh${tthresh}K_usfc${usfc}

    end
  end
end
