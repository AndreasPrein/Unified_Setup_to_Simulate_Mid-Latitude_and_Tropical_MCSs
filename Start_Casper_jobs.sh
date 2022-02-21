#!/bin/bash
# ======================================
#
# Start_Casper_jobs.sh
#
# This scripts submits individual ASR MCS simulations
# to perform feature tracking on the cloud and precipitation
# field
#
#
# Dependencies:
# none
#
# ======================================
# ======================================
# USER Input

DataDir='/glade/campaign/mmm/c3we/prein/Projects/2019_ASR-MCS/data/Physics-Impacts_CAPE-CIN-Clouds/'

# declare -a Physics=('Thomson_YSU' 'Thomson_MYJ' 'Thomson_MYNN2.5' 'Morrison_YSU' 'Morrison_MYJ' 'Morrison_MYNN2.5' 'P351_YSU' 'P351_MYJ' 'P351_MYNN2.5')

declare -a Simulations=('sgp_20120531_04:00:00_L' 'sgp_20120615_07:00:00_L' 'sgp_20130509_07:00:00_L' 'sgp_20130605_09:00:00_L' 'sgp_20130617_07:00:00_L' 'sgp_20140602_04:00:00_L' 'sgp_20140605_12:00:00_L' 'sgp_20140612_06:00:00_L' 'sgp_20140628_16:00:00_L' 'sgp_20140710_10:00:00_L' 'sgp_20160308_15:00:00_L' 'sgp_20160618_10:00:00_L' 'sgp_20160729_09:00:00_L' 'mao_20140816_14:00:00_L' 'mao_20140917_17:00:00_L' 'mao_20150621_14:00:00_L' 'mao_20140401_15:00:00_L' 'mao_20141210_14:00:00_L' 'mao_20150328_15:00:00_L' 'mao_20150412_12:00:00_L' 'mao_20141004_13:00:00_L' 'mao_20141018_14:00:00_L' 'mao_20141117_18:00:00_L' 'mao_20151106_12:00:00_L')

# declare -a Simulations=('sgp_20120531_04:00:00_L' 'sgp_20120615_07:00:00_L' 'sgp_20130509_07:00:00_L' 'sgp_20130605_09:00:00_L' 'sgp_20130617_07:00:00_L' 'sgp_20140602_04:00:00_L' 'sgp_20140605_12:00:00_L' 'sgp_20140612_06:00:00_L' 'sgp_20140628_16:00:00_L' 'sgp_20140710_10:00:00_L' 'sgp_20160308_15:00:00_L' 'sgp_20160618_10:00:00_L' 'sgp_20160729_09:00:00_L')



# declare -a Simulations=('mao_20140319_09:00:00_L' 'mao_20140320_11:00:00_L' 'mao_20140323_18:00:00_L' 'mao_20140331_06:00:00_L' 'mao_20140401_15:00:00_L' 'mao_20140414_07:00:00_L' 'mao_20140421_07:00:00_L' 'mao_20140509_12:00:00_L' 'mao_20140616_18:00:00_L' 'mao_20140619_14:00:00_L' 'mao_20140703_10:00:00_L' 'mao_20140709_15:00:00_L' 'mao_20140712_17:00:00_L' 'mao_20140816_14:00:00_L' 'mao_20140917_17:00:00_L' 'mao_20141004_13:00:00_L' 'mao_20141015_08:00:00_L' 'mao_20141016_13:00:00_L' 'mao_20141018_14:00:00_L' 'mao_20141109_14:00:00_L' 'mao_20141117_18:00:00_L' 'mao_20141120_10:00:00_L' 'mao_20141127_13:00:00_L' 'mao_20141210_14:00:00_L' 'mao_20141211_11:00:00_L' 'mao_20141214_15:00:00_L' 'mao_20150109_11:00:00_L' 'mao_20150110_13:00:00_L' 'mao_20150111_14:00:00_L' 'mao_20150119_19:00:00_L' 'mao_20150222_18:00:00_L' 'mao_20150303_20:00:00_L' 'mao_20150311_11:00:00_L' 'mao_20150328_15:00:00_L' 'mao_20150329_16:00:00_L' 'mao_20150331_13:00:00_L' 'mao_20150412_12:00:00_L' 'mao_20150425_15:00:00_L' 'mao_20150607_16:00:00_L' 'mao_20150621_14:00:00_L' 'mao_20151106_12:00:00_L' 'mao_20151126_07:00:00_L' 'sgp_20120531_04:00:00_L' 'sgp_20120615_07:00:00_L' 'sgp_20130509_07:00:00_L' 'sgp_20130605_09:00:00_L' 'sgp_20130617_07:00:00_L' 'sgp_20140602_04:00:00_L' 'sgp_20140605_12:00:00_L' 'sgp_20140612_06:00:00_L' 'sgp_20140628_16:00:00_L' 'sgp_20140710_10:00:00_L' 'sgp_20160308_15:00:00_L' 'sgp_20160618_10:00:00_L' 'sgp_20160729_09:00:00_L')

# # MAO selected
# declare -a Simulations=('mao_20140816_14:00:00_L' 'mao_20140917_17:00:00_L' 'mao_20150621_14:00:00_L'
# 'mao_20140401_15:00:00_L' 'mao_20141210_14:00:00_L' 'mao_20150328_15:00:00_L' 'mao_20150412_12:00:00_L' 'mao_20141004_13:00:00_L' 'mao_20141018_14:00:00_L' 'mao_20141117_18:00:00_L' 'mao_20151106_12:00:00_L')

# 'sgp_20120531_04:00:00_L' 'sgp_20120615_07:00:00_L' 'sgp_20130509_07:00:00_L' 'sgp_20130605_09:00:00_L' 'sgp_20130617_07:00:00_L' 'sgp_20140602_04:00:00_L' 'sgp_20140605_12:00:00_L' 'sgp_20140612_06:00:00_L' 'sgp_20140628_16:00:00_L' 'sgp_20140710_10:00:00_L' 'sgp_20160308_15:00:00_L' 'sgp_20160618_10:00:00_L' 'sgp_20160729_09:00:00_L')


# loop over the simulations and submitt progressing scripts to Casper

for i in "${Simulations[@]}"
do
   Logfile=$DataDir'P351_MYNN2.5_'$i'_4_CAPE-CIN-LCL-Hydromet.npz'
   if [ ! -f $Logfile ]; then 
       echo $Logfile
       # start the casper submit script for the selected month
       sed "s/RUN/$i/g" CasperSubmit.sh > CasperSubmit_fin.sh
       # sed -i "s/RUN/$i/g" CasperSubmit_fin.sh
       # sed -i "s/TempDicpl/$TempDicpl/g" CasperSubmit_fin.sh
       # sed -i "s/SpatialDispl/$SpatialDispl/g" CasperSubmit_fin.sh
       # sed -i "s/DeltaDeg/$DeltaDeg/g" CasperSubmit_fin.sh
       # sed -i "s/Zlevel/$Zlevel/g" CasperSubmit_fin.sh
       qsub CasperSubmit_fin.sh
       sleep 1s
   fi
done

exit



