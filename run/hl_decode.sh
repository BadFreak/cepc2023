#! /bin/bash
source ./global_config

#------ I/O directories for each processor -------
# Decode: Raw to ROOT
decode_i=${raw_path}
decode_o=${decode_path}

# HL Decode: Raw to ROOT
hl_decode_i=${raw_path}
hl_decode_o=${hl_decode_path}

#------ Create Directores -----
mkdir -p ${hl_decode_o}

#------ Decode -----
prefix=${file_name}

echo "Decoding ${prefix}"
echo " ... ... "
${bin_dir}/hl_decode ${hl_decode_i}/${file_name}.dat ${hl_decode_o}/${prefix}.root #> ${log_dir}/${prefix}.log
echo "    Output dir ===> ${hl_decode_o}/${prefix}.root "

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
