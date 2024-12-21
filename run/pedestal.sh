#! /bin/bash
source ./global_config

#------ I/O directories for each processor -------
# Decode: Raw to ROOT
decode_i=${raw_path}
decode_o=${decode_path}

# Pedestal
#pedestal_i=${decode_path}
pedestal_i=${hl_decode_path}
pedestal_o_1=${ped_hist}
pedestal_o_2=${ped_fit}

# MIP Fit
mip_i=${decode_path}
mip_o_1=${mip_hist}
mip_o_2=${mip_fit}

#------ Create Directores -----

mkdir -p ${pedestal_o_1}
mkdir -p ${pedestal_o_2}

mkdir -p ${log_dir}

#------ Decode -----
prefix=${file_name}
#prefix=${prefix%.root}

#------ Diagnose -----
echo "Fitting Pedestal: ${file_name}"
echo " ... ... "
${bin_dir}/pedestal ${pedestal_i}/${prefix}.root "${pedestal_o_1}/" "${pedestal_o_2}/" # > ${log_dir}/${prefix}.ped.log
echo "    Output histogram ===> ${pedestal_o_1}/${prefix}.root"
echo "    Output fit result ===> ${pedestal_o_2}/${prefix}.root"

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
#！/bin/sh

