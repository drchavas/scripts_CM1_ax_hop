#!/bin/csh

cp sfcphys_DRCgust.F sfcphys_test.F
set U_gust = '6.78'
sed 's|w1(i,j)=.* +|w1(i,j)= '"${U_gust}"' +|g' sfcphys_test.F > temp_file
mv temp_file sfcphys_test.F


