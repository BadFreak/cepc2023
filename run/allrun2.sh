#!/bin/sh
##################

DATADIR="/home/wangjx/HEP/CEPC2023/Result_SPS/decode/mu-"

##################
#cd ${DATADIR}
echo "execute all run file..."
function allrun2(){
for ENERGY in `ls $1`
do
	for FILENAME in `ls $1/${ENERGY}`
	do 	
		fileno=${FILENAME#*Run}
		fileno=${fileno%%_*}
		let filenum=fileno
		#if [ $filenum -le 107 -a $filenum -ge 96 ] || [ $filenum -le 138 -a $filenum -ge 134 ]
		if [ $filenum -ge 139 -a $filenum -le 166 ]
		#if [ $filenum -le 133 -a $filenum -ge 133 ] 
		then
			echo " ====== NOW file is ${FILENAME%.root} ======"
			echo " modifying the global_config ..."
			sed -i "9c particle_tag=mu-" global_config
			sed -i "10c energy_tag=$ENERGY" global_config
			sed -i "11c file_name=${FILENAME%.root}" global_config
			source calib.sh
			echo " ====== NOW ${FILENAME} finished ====="
		fi
	done
done
}

allrun2 ${DATADIR}

