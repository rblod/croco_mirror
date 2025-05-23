#!/bin/bash
#===================================

set -u
###set -x
###DO NOT USE SET -E !

source CONFIGURE_GLOBAL
source configure_file

#cd $SUBMIT_DIR
#echo "$SUBMIT_DIR"

#===================================

par1='MPI'
# Compile
msg1="- Compilation failure for ${TEST_NAME} : ${par1}..."
msg2="${FMT_REDBLD}${msg1}${FMT_ORD}"
./jobcomp_rvtk.bash Compile_$par1 > jobcomp_${par1}_${TEST_NAME}.log  2>&1 || { echo -e "   $msg2" | tee -a mylog.txt ; echo -e $msg1 ; exit 1 ; }
/bin/mv croco croco_${par1}.exe

# Run
msg1="- Execution failure for ${TEST_NAME} : ${par1}..."
msg2="${FMT_REDBLD}${msg1}${FMT_ORD}"

msg3="- Parallel repro. failure for ${TEST_NAME} : ${par1}..."
msg4="${FMT_REDBLD}${msg3}${FMT_ORD}"

$MPIRUN --allow-run-as-root -np $NBPROCS ./croco_${par1}.exe $CROCOIN > mpi_${NBPROCS}_${TEST_NAME}.log 2>&1  || { echo -e "   $msg2" | tee -a mylog.txt ; echo -e $msg1 ; exit 2 ; }
exec_status=$?
echo "execution_status is "$exec_status
#  =0, OK or  clean stop before the end (bugbin or blow up)
# !=0, KO and bad stop before the end (input problem)

grep 'BUGBIN' mpi_${NBPROCS}_${TEST_NAME}.log > /dev/null 2>&1
bugbin_detec=$?
echo "bugbin detection flag is " $bugbin_detec

if [ $exec_status -eq 0 ] && [ $bugbin_detec -eq 0 ] ; then
    # sortie car prb  exec CAR erreur de repro (bugbin detection)
    echo -e "   $msg4" | tee -a mylog.txt
    echo -e $msg3
    exit 3 
fi
#echo 'Output message status is' $?

#===
# # FAKE CHANGE #2
# # Additional check in case of clean stop before the end, mainly Parallel repro failure case
# #     Clean stop in case of blow up or bugbin => exit status =0
# #     Bad stop in case of bad input or else  => exit status !=0

# SUCCESS_TMP=1
# grep 'MAIN: DONE'  mpi_${NBPROCS}_${TEST_NAME}.log || SUCCESS_TMP=0
# if [  "$SUCCESS_TMP" -eq 0 ]; then
#   echo -e "   $msg4" | tee -a mylog.txt
#   echo -e $msg3 
#   exit 2 
# fi	
#===
