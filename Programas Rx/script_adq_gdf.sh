#!/bin/bash
PIDFILE=/tmp/SoftRX.pid

case $1 in
        start)
                cd /home/igp-114/proc_GDF2HDF5/
                python /home/igp-114/proc_GDF2HDF5/hfrx_gdf.py 2>/dev/null &
                # Get its PID and store it
                echo $! > ${PIDFILE}
        ;;
        stop)
                kill 'cat ${PIDFILE}'
                # Now that its killed, don't  forget to remote the PID file
                rm ${PIDFILE}
        ;;
        *)
                echo "usage: Software RX {start|stop}";;
esac
exit 0
