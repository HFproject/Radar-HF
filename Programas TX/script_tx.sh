#!/bin/bash
echo "Iniciando Secuencia de Transmision"
echo "Ubicando el Directorio de Datos"

cd /home/igp-114/HF/gr-filesink/trunk/apps/hfradar/
echo ""
export DISPLAY=:0
for i in {1..3}
do 
       echo -n "."
       sleep 1s
done
echo ""
echo "Iniciando Tx F0"
screen -S "TX_F0" -d -m ./hftx_f2.72.sh
for i in {1..10}
do
       echo -n "."
       sleep 1s
done
echo ""
echo "Iniciando Tx F1"
screen -S "TX_F1" -d -m ./hftx_f3.64.sh
for i in {1..10}
do 
       echo -n "."
       sleep 1s
done
screen -S "PING" -d -m ping 8.8.8.8
echo ""
