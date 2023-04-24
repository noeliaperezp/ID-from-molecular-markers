#!/bin/bash
#$ -cwd

rm script_SimIDDrosoPneu.sh.*

if [ $# -ne 8 ]  
then
	echo "Usage: $0 <TNIND> <NIND> <FS> <G> <genBP> <REPS> <type> <n>" 
	exit 1
fi

#Set arguments
TNIND=$1
NIND=$2
FS=$3
G=$4
genBP=$5
REPS=$6
type=$7
n=$8

WDIR=$PWD 
if [ -d "RESULTS.T$TNIND.N$NIND.neu.$n" ]
then
rm -r RESULTS.T$TNIND.N$NIND.neu.$n
fi
mkdir -p $WDIR/RESULTS.T$TNIND.N$NIND.neu.$n

mkdir -p /state/partition1/EstMolNe$n/$SLURM_JOBID/
cp SimIDDrosoPneu /state/partition1/EstMolNe$n/$SLURM_JOBID/
cp dataBP.map /state/partition1/EstMolNe$n/$SLURM_JOBID/dataBP.map
cp dataBP.ped /state/partition1/EstMolNe$n/$SLURM_JOBID/dataBP.ped
cp list_allsnps /state/partition1/EstMolNe$n/$SLURM_JOBID/list_allsnps
touch $WDIR/$SLURM_JOBID.`hostname`.`date +%HH%MM`
cd /state/partition1/EstMolNe$n/$SLURM_JOBID

num=$RANDOM
echo "$num" > seedfile

time ./SimIDDrosoPneu>>out<<@
1
-99
1000	NP
$TNIND	TNIND
$NIND	NIND
0.1	var_effects
1.0	ha (positive effect)
2.0	VE
20.0	L (99: free recombination)
$G	generations
$genBP		genBP
$FS	RM (0), FS (1)
1	neutral (1) for fitness
0	INITIALFREQ05
0.0	MAF
$REPS	replicates
$type	type
@

echo "gen ind W    Fped     Fibd     Fvr1      Fvr2      Fyang1   Fyang2   FLH1     FLH2     Fhom" > kk1
mv list_individual_Fs kk2
cat kk1 kk2 > list_individual_Fs
mv list_individual_FsBP kk2
cat kk1 kk2 > list_individual_FsBP
mv list_individual_Fs_pairs kk3
cat kk1 kk3 > list_individual_Fs_pairs

cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/genfile.dat $WDIR/RESULTS.T$TNIND.N$NIND.neu.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/outfileF.dat $WDIR/RESULTS.T$TNIND.N$NIND.neu.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/outfileID.dat $WDIR/RESULTS.T$TNIND.N$NIND.neu.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/outfileNe.dat $WDIR/RESULTS.T$TNIND.N$NIND.neu.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/outfileR.dat $WDIR/RESULTS.T$TNIND.N$NIND.neu.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/outfileVF.dat $WDIR/RESULTS.T$TNIND.N$NIND.neu.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/outfileFBP.dat $WDIR/RESULTS.T$TNIND.N$NIND.neu.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/outfileIDBP.dat $WDIR/RESULTS.T$TNIND.N$NIND.neu.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/outfileNeBP.dat $WDIR/RESULTS.T$TNIND.N$NIND.neu.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/outfileRBP.dat $WDIR/RESULTS.T$TNIND.N$NIND.neu.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/outfileVFBP.dat $WDIR/RESULTS.T$TNIND.N$NIND.neu.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/outfileFLASTGEN.dat $WDIR/RESULTS.T$TNIND.N$NIND.neu.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/outfileVFLASTGEN.dat $WDIR/RESULTS.T$TNIND.N$NIND.neu.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/OUTFILEDROSO $WDIR/RESULTS.T$TNIND.N$NIND.neu.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/frequencies.dat $WDIR/RESULTS.T$TNIND.N$NIND.neu.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/list* $WDIR/RESULTS.T$TNIND.N$NIND.neu.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/GRAPHS $WDIR/RESULTS.T$TNIND.N$NIND.neu.$n
#cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/repfile.dat $WDIR/RESULTS.T$TNIND.N$NIND.neu.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/dfilename* $WDIR/RESULTS.T$TNIND.N$NIND.neu.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/data.map $WDIR/RESULTS.T$TNIND.N$NIND.neu.$n/data$n.map
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/data.ped $WDIR/RESULTS.T$TNIND.N$NIND.neu.$n/data$n.ped

cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/out $WDIR/RESULTS.T$TNIND.N$NIND.neu.$n

rm -r /state/partition1/EstMolNe$n/$SLURM_JOBID/
rm $WDIR/$SLURM_JOBID.*