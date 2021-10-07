#!/bin/bash
echo "Borrando Archivos Incompletos..."
cd /home/igp-114/
screen -S "delete" -d -m ./delete_gdf_size0.sh

echo "Dont Forget sudo sysctl -w net.core.rmem_max=100000000"

echo "Iniciando PRE-Procesamientos..."
cd /home/igp-114/proc_GDF2HDF5
echo "Iniciando F0"

###############################################################################################################################
# NOTA 
#  * VALOR 1 :    EL PRIMER VALOR  PUEDE SER
#	          O (MODO NORMAL DE OPERACION 100 PERFILES)  permite escoger 10 o 60 segundos con el numero de integraciones incoherentes.
#       	  1 (MODO CAMPAÑA 600 PERFILES)  10 segundos seteado las integraciones incoherentes siempre a 0.
# 
#  * VALOR 2:
#  		  0-6   EL SEGUNDO VALOR PUEDE SER UN NUMERO DE  0 A 6 INDICA  EL NUMERO DE INTEGRACIONES INCOHERENTES 
#           	  PARA GENERAR LOS SPECTROS.
#               
#                 screen -S "PROC_F0" -d -m ./proc_f0.sh VALOR1 VALOR2
#
#  * VALOR 3:     EL TERCER VALOR PUEDE SER  
#		  "0", "0,1" ,"0,2" , "0,1,2"  ESTE PARAMETRO REPRESENTA  EL NUMERO DE TX PRESENTES Y EL NUMERO DE DECODIFICACIONES QUE EL PROGRAMA
#		  REALIZA PARA GENERAR LOS ARCHIVOS hdf5, EN EL CASO DE '0' SE REALIZA UNA DECODIFICACION CON CODIGO-0 QUE PROVIENEN DEL TX UBICADO EN ANCON,
#		  PARA EL CASO "1" LE CORRESPONTE EL CODIGO-1 Y EL TX UBICADO EN HUANCAYO Y "2" LE CORRESPONDE EL CODIGO-2 Y EL TX UBICADO EN ICA
#
#
#
#
################################################################################################################################

#screen -S "PROC_F0" -d -m ./proc_f0.sh 1 0 "0,1,2" #CAMPAIGN MODE
screen -S "PROC_F0" -d -m ./proc_f0.sh 0 6 "0,1,2" # NORMAL MODE
for i in {1..4}
do
   echo -n "."
   sleep 1s
done
echo ""
echo "PROC_F0 Ejecutado"

echo "Iniciando F1"
#screen -S "PROC_F1" -d -m ./proc_f1.sh 1 0 "0,1,2"   # CAMPAIGN MODE
screen -S "PROC_F1" -d -m ./proc_f1.sh 0 6 "0,1,2"   # NORMAL MODE
for i in {1..4}
do
   echo -n "."
   sleep 1s
done
echo ""
echo "PROC_F1 Ejecutado"