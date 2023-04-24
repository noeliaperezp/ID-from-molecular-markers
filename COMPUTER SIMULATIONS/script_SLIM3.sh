#!/bin/bash
#$ -cwd

rm script_SLIM3_PRUEBA.sh.*

if [ $# -ne 1 ]
then
	echo "Usage: $0 <n>"
	exit 1
fi

#Set arguments
n=$1

WDIR=$PWD

#Output directory
if [ -d "$WDIR/OUTPUT_SLIM_$n" ]
then
rm -r $WDIR/OUTPUT_SLIM_$n
fi

mkdir -p $WDIR/OUTPUT_SLIM_$n

mkdir -p /state/partition1/armandoSLIM$n/$SLURM_JOBID/

###################### TRANSFER TO state/partition1 #########################

cp SNP_BP_SLIM3_2 /state/partition1/armandoSLIM$n/$SLURM_JOBID/
cp slim3INPUT /state/partition1/armandoSLIM$n/$SLURM_JOBID/

touch $WDIR/$SLURM_JOBID.`hostname`.`date +%HH%MM`
cd /state/partition1/armandoSLIM$n/$SLURM_JOBID

###################### run SLIM3 #########################

START=$(date +%s)
module load SLiM/3.3.2
slim slim3INPUT > slimout

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "SLIM3 took 	$DIFF seconds" >> timefile

cp -r /state/partition1/armandoSLIM$n/$SLURM_JOBID/slimout $WDIR/OUTPUT_SLIM_$n

###################### SNP_BP_SLIM3_2 #########################

num=$RANDOM
echo "$num" > seedfile

START=$(date +%s)

./SNP_BP_SLIM3_2<<@
-99
1000	N
100000000 chrlen
@

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "SNP_BP took 	$DIFF seconds" >> timefile

######################## TRANSFER OF FILES TO DIRECTORY #########################

cp -r /state/partition1/armandoSLIM$n/$SLURM_JOBID/dataBP.ped $WDIR/OUTPUT_SLIM_$n
cp -r /state/partition1/armandoSLIM$n/$SLURM_JOBID/dataBP.map $WDIR/OUTPUT_SLIM_$n
cp -r /state/partition1/armandoSLIM$n/$SLURM_JOBID/list* $WDIR/OUTPUT_SLIM_$n
cp -r /state/partition1/armandoSLIM$n/$SLURM_JOBID/timefile $WDIR/OUTPUT_SLIM_$n
cp -r /state/partition1/armandoSLIM$n/$SLURM_JOBID/slimout $WDIR/OUTPUT_SLIM_$n
cp -r /state/partition1/armandoSLIM$n/$SLURM_JOBID/collect* $WDIR/OUTPUT_SLIM_$n

######################## state/partition1 CLEANING #########################

rm -r /state/partition1/armandoSLIM$n/$SLURM_JOBID/
rm $WDIR/$SLURM_JOBID.*
