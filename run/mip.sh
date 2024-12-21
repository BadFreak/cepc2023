#! /bin/bash
isTrackFit='1' 
source ./global_config

#------ I/O directories for each processor -------
# Decode: Raw to ROOT
decode_i=${raw_path}
decode_o=${decode_path}

# Pedestal
pedestal_i=${decode_path}
pedestal_o_1=${ped_hist}
pedestal_o_2=${ped_fit}

# MIP Fit
mip_i=${decode_path}
mip_o_1=${mip_hist}
mip_o_2=${mip_fit}

#------ Create Directores -----

mkdir -p ${mip_o_1}
mkdir -p ${mip_o_2}

mkdir -p ${log_dir}

#------ Decode -----
prefix=${file_name}
#prefix=${file_name%.dat}
#prefix=${prefix%.root}

echo "Fitting MIP: ${file_name}"
echo " ... ... "
${bin_dir}/mip ${mip_i}/${prefix}.root ${pede_factor} ${mip_factor} "${mip_o_1}/" "${mip_o_2}/" ${isTrackFit} # > ${log_dir}/${prefix}.mip.log
echo "    Output histogram ===> ${mip_o_1}/${prefix}.root"
echo "    Output fit result ===> ${mip_o_2}/${prefix}.root"

echo "=============================="

# if [[ -e ${decode_o}/${date_tag}.root ]]
# then
# 	rm ${decode_o}/${date_tag}.root
# 	hadd ${decode_o}/${date_tag}.root ${decode_o}/*
# 	rm -i ${decode_o}/*_*.root
# else
# 	hadd ${decode_o}/${date_tag}.root ${decode_o}/*
# 	rm -i ${decode_o}/*_*.root
# fi
