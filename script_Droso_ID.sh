#!/bin/bash
#$ -cwd

rm script_Droso_ID.sh.*
rm data.F
rm dataF.het
rm data.hom.indiv
rm data.hom.summary
rm data.ibc
rm outfileID
rm timefile

#Check number of arguments
if [ $# -ne 4 ]  
then
	echo "Usage: $0 <NIND> <n> <dataname> <phen>" 
	exit 1
fi

### Set arguments

nind=$1
n=$2
dataname=$3
phen=$4

WDIR=$PWD 
mkdir -p /state/partition1/noelia_IDroso$n/$SLURM_JOBID/

################## Output directory
mkdir -p $WDIR/$dataname/

###################### TRANSFER TO state/partition1 #########################

cp shell_F_values /state/partition1/noelia_IDroso$n/$SLURM_JOBID/
cp plink1.9 /state/partition1/noelia_IDroso$n/$SLURM_JOBID/
cp gcta64 /state/partition1/noelia_IDroso$n/$SLURM_JOBID/
cp $dataname.map /state/partition1/noelia_IDroso$n/$SLURM_JOBID/data.map
cp $dataname.ped /state/partition1/noelia_IDroso$n/$SLURM_JOBID/data.ped
cp $phen.phe /state/partition1/noelia_IDroso$n/$SLURM_JOBID/qt.phe
cp SNP_ID_Droso /state/partition1/noelia_IDroso$n/$SLURM_JOBID/

touch $WDIR/$SLURM_JOBID.`hostname`.`date +%HH%MM`
cd /state/partition1/noelia_IDroso$n/$SLURM_JOBID

START=$(date +%s)

./SNP_ID_Droso<<@
$nind	nind
@
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "SNP_ID took 	$DIFF seconds" >> timefile

######################## TRANSFER OF FILES TO DIRECTORY #########################

cp -r /state/partition1/noelia_IDroso$n/$SLURM_JOBID/data.F $WDIR/$dataname/data.F
cp -r /state/partition1/noelia_IDroso$n/$SLURM_JOBID/outfileID $WDIR/$dataname/outfileID
cp -r /state/partition1/noelia_IDroso$n/$SLURM_JOBID/timefile $WDIR/$dataname/

######################## state/partition1 CLEANING #########################

rm -r /state/partition1/noelia_IDroso$n/$SLURM_JOBID/
rm $WDIR/$SLURM_JOBID.*
