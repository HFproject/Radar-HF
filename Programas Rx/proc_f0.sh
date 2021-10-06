#!/bin/bash
echo "PROC_F0"

while true;do python hf_proc_gdf2hdf5.py  --station-name=HF --freq=f0 -N 90 -x 10 --rawdata-from-doy=/media/igp-114/RAWDATA/ --procdata=/media/igp-114/PROCDATA/  --graphics=file/media/igp-114/PROCDATA/ --send-graphs=0 --auto-delete-rawdata=1 -I $2 --mode-campaign=$1 --stationtx-codes=$3;done
