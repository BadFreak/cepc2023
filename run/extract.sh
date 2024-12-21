#! /bin/bash
source ./global_config

#------ I/O directories for each processor -------

echo " Extract Parameter && Dead Channels"
echo " ... ... "
${bin_dir}/extract ${pede_factor} ${inter_factor} ${mip_factor}

echo "    Output file ===> DeadList.txt "
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
