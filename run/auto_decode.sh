#! /bin/bash
source ./global_config

#------ I/O directories for each processor -------
# Decode: Raw to ROOT
decode_i=${raw_path}
decode_o=${decode_path}

# Auto Decode: Raw to ROOT
auto_decode_i=${raw_path}
auto_decode_o=${auto_decode_path}

# Pedestal
pedestal_i=${decode_path}
pedestal_o_1=${ped_hist}
pedestal_o_2=${ped_fit}

# MIP Fit
mip_i=${decode_path}
mip_o_1=${mip_hist}
mip_o_2=${mip_fit}

#------ Create Directores -----
mkdir -p ${auto_decode_o}

#------ Decode -----
prefix=${file_name}

echo "Decoding ${prefix}"
echo " ... ... "
${bin_dir}/auto_decode ${auto_decode_i}/${file_name}.dat ${auto_decode_o}/${prefix}.root #> ${log_dir}/${prefix}.log
echo "    Output dir ===> ${auto_decode_o}/${prefix}.root "

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
