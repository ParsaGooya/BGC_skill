
'''

bash script to extract the newly created assimilation runs from disc to netcdf files

'''


#!/bin/bash
export CCCMA_REF=/home/scrd102/cccma_libs/cccma/latest/   
source $CCCMA_REF/src/CCCma_tools/generic/u3_site_profile
load_cccma_env
# export PATH=$CCCMA_REF/scr/CCCma_tools/tools:$PATH  
export DATAPATH_DB=/space/hall7/sitestore/eccc/crd/cccma/users/rpg002/datapath_local.db 
umask 022
set -a



#=============================== user specifies the configs below ===================================
vars=("mlotst" "so" "thetao")
realms=("Omon" "Omon" "Omon")
start_year=2017
#====================================================================================================



if [ "${#vars[@]}" -ne "${#realms[@]}" ]; then
    echo "Error: vars and realms must have the same length."
    exit 1
fi

for i in "${!vars[@]}"; do
var="${vars[$i]}"
realm="${realms[$i]}"
echo "============================================================================ 
  extracting assimilation runs for ${var} : 
============================================================================"
workdir="/space/hall5/sitestore/eccc/crd/ccrn/users/rpg002/data/${var}/assimilation/CanESM5/extentions"
mkdir -p ${workdir}
cd ${workdir}

echo "Extracting ${var} from ${realm}"
 
for year in $(seq "${start_year}" 2021); do


  hpcarcuser="cpd102"
 
  for e in `seq 1 10` ; do
    e2ch=`echo ${e} | awk '{print($1 +100)}' | cut -c2-3`
    hpcarchive_pfx="nc_d2a-asm-e${e2ch}-re_${year}01_${year}12"
    hpcarchive=`hpcarchive -p crd_cccma -u ${hpcarcuser} -L -x -s -c "^$hpcarchive_pfx" | grep $hpcarchive_pfx | head -1`
    hpcarcname=${hpcarchive##* }
    echo $hpcarcname
    # hpcarchive -p crd_cccma -u ${hpcarcuser} -r -b -x -c ${hpcarcname} -f ${var}_${realm}_CanESM5_dcppA-assim_r${e}i1p2f1_gn_${year}01-${year}12.nc
    hpcarchive -p crd_cccma -u ${hpcarcuser} -r -x -c ${hpcarcname} -f ${var}_${realm}_CanESM5_dcppA-assim_r${e}i1p2f1_gn_${year}01-${year}12.nc

  done
 
done
 
for year in 2022 2023; do
  hpcarcuser="scrd107"
 
  for e in `seq 1 10` ; do
    e2ch=`echo ${e} | awk '{print($1 +100)}' | cut -c2-3`
    hpcarchive_pfx="nc_d2a-asm-e${e2ch}-u2_${year}01_${year}12"
    hpcarchive=`hpcarchive -p crd_cccma -u ${hpcarcuser} -L -x -s -c "^$hpcarchive_pfx" | grep $hpcarchive_pfx | head -1`
    hpcarcname=${hpcarchive##* }
    echo $hpcarcname
    # hpcarchive -p crd_cccma -u ${hpcarcuser} -r -b -x -c ${hpcarcname} -f ${var}_${realm}_CanESM5_dcppA-assim_r${e}i1p2f1_gn_${year}01-${year}12.nc
    hpcarchive -p crd_cccma -u ${hpcarcuser} -r -x -c ${hpcarcname} -f ${var}_${realm}_CanESM5_dcppA-assim_r${e}i1p2f1_gn_${year}01-${year}12.nc

  done
 
done
done




# # #!/usr/bin/env bash
# var='mlotst'
# realm='Omon'
# year=2023
# for ens_id in {1..10}; do
#     e2ch=`echo ${ens_id} | awk '{print($1 +100)}' | cut -c2-3`
#     mv \
#     /space/hall5/sitestore/eccc/crd/ccrn/users/rpg002/data/${var}/assimilation/CanESM5/extentions/nc_d2a-asm-e${e2ch}-u2_${year}01_${year}12_*/CMIP6/DCPP/CCCma/CanESM5/dcppA-assim/r${ens_id}i1p2f1/${realm}/${var}/gn/v20190429/${var}_${realm}_CanESM5_dcppA-assim_r${ens_id}i1p2f1_gn_${year}01-${year}12.nc \
#     /space/hall5/sitestore/eccc/crd/ccrn/users/rpg002/data/${var}/assimilation/CanESM5/extentions
# done



# var='mlotst'
# realm='Omon'
# year=2021
# for ens_id in {1..10}; do
#     e2ch=`echo ${ens_id} | awk '{print($1 +100)}' | cut -c2-3`
#     mv \
#     /space/hall5/sitestore/eccc/crd/ccrn/users/rpg002/data/${var}/assimilation/CanESM5/extentions/output/CMIP6/DCPP/CCCma/CanESM5/dcppA-assim/r${ens_id}i1p2f1/${realm}/${var}/gn/v20190429/${var}_${realm}_CanESM5_dcppA-assim_r${ens_id}i1p2f1_gn_${year}01-${year}12.nc \
#     /space/hall5/sitestore/eccc/crd/ccrn/users/rpg002/data/${var}/assimilation/CanESM5/extentions
# done
