#！/bin/sh
source ./global_config

#------I/O dircetories ----
# HL Ratio
hl_ratio_i=${decode_path}
hl_ratio_o_1=${hl_ratio_fit}
hl_ratio_o_2=${hl_ratio_hist}

#------Create dircetories----
mkdir -p ${hl_ratio_o_1}
mkdir -p ${hl_ratio_o_2}

#------HL ratio ------
prefix=${file_name}

echo "Running HL ratio ${prefix}"
echo "... ..."
${bin_dir}/hl_ratio ${hl_ratio_i}/${prefix}.root ${hl_ratio_o_2}/${prefix}.root ${hl_ratio_o_1}/${prefix}.root  ${pede_factor}
echo "		Output fit dir ===> ${hl_ratio_o_1}/${prefix}.root"
echo "		Output hist dir ===> ${hl_ratio_o_2}/${prefix}.root"
echo "=================="

# path="$PWD"
# datafile=20201216_20210308
# inputfile="/mnt2/ScECAL_CR/results/trackFit/"
# caliFile="/mnt2/ScECAL_CR/results/calibration/"
# DIR_IN=${inputfile}${datafile}".root"
# DIR_OUT1=${caliFile}
# 
# ./main ${DIR_IN} ${DIR_OUT1}
# echo "all over"
