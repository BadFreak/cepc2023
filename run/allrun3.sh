#!/bin/sh
##################

DATADIR="/home/wangjx/HEP/CEPC2023/Result_MC/decode/e-"

##################
#cd ${DATADIR}
echo "execute all run file..."
function allrun3(){
for ENERGY in `ls $1`; do
	if [[ $ENERGY = *"GeV"* ]]; then
		for FILENAME in `ls $1/${ENERGY}`; do 	
			# echo " ====== NOW file is ${FILENAME} ======"
			if [[ $FILENAME = *"GeV"* ]]; then
				# echo " ====== NOW file is ${FILENAME} ======"
				echo " PROCESSING : modifying the global_config ..."
				sed -i "9c particle_tag=e-" global_config
				sed -i "10c energy_tag=$ENERGY" global_config
				sed -i "11c file_name=$FILENAME" global_config
				source calib.sh
				echo " ====== NOW ${FILENAME} finished ====="
			else
				echo "SKIPING ${FILENAME} ... ..."
			fi
		done
	fi
done
}

allrun3 ${DATADIR}

