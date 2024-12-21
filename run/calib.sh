#! /bin/bash
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

# Calibration
calib_o=${calib_path}

#------ Create Directores -----
mkdir -p ${calib_o}

prefix=${file_name%.dat}
prefix=${prefix%.root}

echo "Calibration ${prefix}"
echo " ... ... "
${bin_dir}/calib ${pede_factor} ${inter_factor} ${mip_factor} ${decode_o}/${prefix}.root ${calib_o}/${prefix}.root #> ${log_dir}/${prefix}.log
echo "    Output dir ===> ${calib_o}/${prefix}.root "

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
