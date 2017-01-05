#!/bin/sh
#deletes run files made for all cases, keeps folders, pngs and .yamls
CUIBM_DIR=/scratch/src/cuIBM

CAVITY_DIR=$CUIBM_DIR/validation/lidDrivenCavity
CYLINDER_DIR=$CUIBM_DIR/validation/cylinder

for filename in 100 1000 10000
do
	echo "Cleaning lid driven cavity files for Reynolds number " $filename
	cd $CAVITY_DIR/Re$filename
	find . ! -name 'bodies.yaml' ! -name 'domain.yaml' ! -name 'flow.yaml' ! -name 'simParams.yaml' ! -name '*.png' -type f -exec rm -f {} +
	cd ..
done

for filename in 40 75 100 150 200 550 3000
do
	echo "Cleaning cylinder files for Reynolds number " $filename
	cd $CYLINDER_DIR/Re$filename
	find . ! -name 'bodies.yaml' ! -name 'domain.yaml' ! -name 'flow.yaml' ! -name 'simParams.yaml' ! -name '*.png' -type f -exec rm -f {} +
	cd ..
done
