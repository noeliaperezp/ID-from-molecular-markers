#!/bin/bash
#$ -cwd

rm script_Droso_ID_randomSNPs.sh.*
rm data.F
rm dataF.het
rm data.hom.indiv
rm data.hom.summary
rm data.ibc
rm outfileID
rm timefile

#Check number of arguments
if [ $# -ne 5 ]  
then
	echo "Usage: $0 <NIND> <dataname> <phen> <n> <SNPS> " 
	exit 1
fi

### Set arguments

nind=$1
dataname=$2
phen=$3
n=$4
SNPS=$5

WDIR=$PWD 
mkdir -p /state/partition1/noelia_IDroso$n/$SLURM_JOBID/

################## Output directory
dirname="${dataname}_${SNPS}"
mkdir -p $WDIR/$dirname/

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

cut -f 2 data.map > snps.map
shuf -n $SNPS snps.map > snps.subset.map
./plink1.9 --file data --aec --extract snps.subset.map --recode --out data

START=$(date +%s)

./SNP_ID_Droso<<@
$nind	nind
@
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "SNP_ID took 	$DIFF seconds" >> timefile

######################## TRANSFER OF FILES TO DIRECTORY #########################

cp -r /state/partition1/noelia_IDroso$n/$SLURM_JOBID/data.F $WDIR/$dirname/s$n.data.F
cp -r /state/partition1/noelia_IDroso$n/$SLURM_JOBID/outfileID $WDIR/$dirname/s$n.outfileID
cp -r /state/partition1/noelia_IDroso$n/$SLURM_JOBID/timefile $WDIR/$dirname/
#cp -r /state/partition1/noelia_IDroso$n/$SLURM_JOBID/data.map $WDIR/$dirname/s$n.data.map
#
cp -r /state/partition1/noelia_IDroso$n/$SLURM_JOBID/ROH100_data.hom $WDIR/$dirname/
cp -r /state/partition1/noelia_IDroso$n/$SLURM_JOBID/ROH100_data.hom.indiv $WDIR/$dirname/
cp -r /state/partition1/noelia_IDroso$n/$SLURM_JOBID/ROH1000_data.hom $WDIR/$dirname/
cp -r /state/partition1/noelia_IDroso$n/$SLURM_JOBID/ROH1000_data.hom.indiv $WDIR/$dirname/
cp -r /state/partition1/noelia_IDroso$n/$SLURM_JOBID/ROH2000_data.hom $WDIR/$dirname/
cp -r /state/partition1/noelia_IDroso$n/$SLURM_JOBID/ROH2000_data.hom.indiv $WDIR/$dirname/

done

######################## state/partition1 CLEANING #########################

rm -r /state/partition1/noelia_IDroso$n/$SLURM_JOBID/
rm $WDIR/$SLURM_JOBID.*
