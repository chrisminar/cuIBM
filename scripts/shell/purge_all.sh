#!/bin/sh
#deletes all file output made during any run for all cases
#WARNING USE THIS WITH CAUTION, MAKE A BACKUP BEFORE YOU RUN THIS
CUIBM_DIR=/scratch/src/cuIBM

CAVITY_DIR=$CUIBM_DIR/validation/lidDrivenCavity
CYLINDER_DIR=$CUIBM_DIR/validation/cylinder

for filename in 100 1000 10000
do
	echo "Cleaning lid driven cavity for Reynolds number " $filename
	cd $CAVITY_DIR/Re$filename
	pwd
	find . ! -name 'bodies.yaml' ! -name 'domain.yaml' ! -name 'flow.yaml' ! -name 'simParams.yaml' -type f -exec rm -f {} +
	find . ! -name '.' -type d -exec rm -r {} +
	mkdir output
	cd ..
done

for filename in 40 75 100 150 200 550 3000
do
	echo "Cleaning cylinder for Reynolds number " $filename
	cd $CYLINDER_DIR/Re$filename
	pwd
	find . ! -name 'bodies.yaml' ! -name 'domain.yaml' ! -name 'flow.yaml' ! -name 'simParams.yaml' -type f -exec rm -f {} +
	find . ! -name '.' -type d -exec rm -r {} +
	mkdir output
	cd ..
done
