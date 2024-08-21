#!/bin/bash


export CCCMA_REF=/home/scrd102/cccma_libs/cccma/latest/                  # defines the lib version to use
export PATH=$CCCMA_REF/CanESM_source_link/CCCma_tools/tools:$PATH        # Access to setup/s scripts
source $CCCMA_REF/CanESM_source_link/CCCma_tools/generic/u2_site_profile # basic ordenv setup & ssm loads of maestro etc
 
alias load_cccma_env='source $CCCMA_REF/env_setup_file'                  # command to activate a full env with binaries like ggstat
umask 022                                                                 # Default read permissions on new files for group



set -a
 
workdir="/space/hall5/sitestore/eccc/crd/ccrn/users/rpg002/data/chlos/assimilation/extentions"
mkdir -p ${workdir}
cd ${workdir}
 
for year in 2021 ; do
  hpcarcuser="cpd102"
 
  for e in `seq 1 40` ; do
    e2ch=`echo ${e} | awk '{print($1 +100)}' | cut -c2-3`
    hpcarchive_pfx="nc_d2a-asm-e${e2ch}-re_${year}01_${year}12"
    hpcarchive=`hpcarchive -p crd_cccma -u ${hpcarcuser} -L -x -s -c "^$hpcarchive_pfx" | grep $hpcarchive_pfx | head -1`
    hpcarcname=${hpcarchive##* }
    echo $hpcarcname
    hpcarchive -p crd_cccma -u ${hpcarcuser} -r -b -x -c ${hpcarcname} -f chlos_Omon_CanESM5_dcppA-assim_r${e}i1p2f1_gn_${year}01-${year}12.nc
  done
 
done
 
for year in 2022 2023 ; do
  hpcarcuser="scrd107"
 
  for e in `seq 1 40` ; do
    e2ch=`echo ${e} | awk '{print($1 +100)}' | cut -c2-3`
    hpcarchive_pfx="nc_d2a-asm-e${e2ch}-u2_${year}01_${year}12"
    hpcarchive=`hpcarchive -p crd_cccma -u ${hpcarcuser} -L -x -s -c "^$hpcarchive_pfx" | grep $hpcarchive_pfx | head -1`
    hpcarcname=${hpcarchive##* }
    echo $hpcarcname
    hpcarchive -p crd_cccma -u ${hpcarcuser} -r -b -x -c ${hpcarcname} -f chlos_Omon_CanESM5_dcppA-assim_r${e}i1p2f1_gn_${year}01-${year}12.nc
  done
 
done