#!/bin/bash

DATASET_NAME=$1
SERVER=$2
REF_NAME=$3
ANALYZER=$4

REF_PATH="/nfs/data8/chuyu/data/my_register/refbrain/${REF_NAME}_crop_16bit.nrrd"
DATA_PATH="/nfs/data${SERVER}/${ANALYZER}/data/${DATASET_NAME}"

hs=`hostname`
hs=${hs:5}
if [ "$hs" == "$SERVER" ]
then
    DATA_PATH="/data/${ANALYZER/data}/${DATASET_NAME}"
fi

if [ "$hs" == 8 ]
then
    REF_PATH="/data/chuyu/data/my_register/refbrain/${REF_NAME}_crop_16bit.nrrd"
fi


IMG_PATH="${DATA_PATH}/${DATASET_NAME}_ds_stack_mean.nrrd"
RESULT_PATH="${DATA_PATH}/${DATASET_NAME}_${REF_NAME}"
OUTPUT_NAME="${DATASET_NAME}_${REF_NAME}_registered.nrrd"


OUTPUT_NAME_INIT="${DATASET_NAME}_${REF_NAME}_registered_initial.nrrd"
OUTPUT_NAME_RIGID="${DATASET_NAME}_${REF_NAME}_registered_rigid.nrrd"
OUTPUT_NAME_WARP="${DATASET_NAME}_${REF_NAME}_registered_warp.nrrd"

OUTPUT_NAME_INIT_INV="${DATASET_NAME}_${REF_NAME}_registered_initial_inv.nrrd"
OUTPUT_NAME_RIGID_INV="${DATASET_NAME}_${REF_NAME}_registered_rigid_inv.nrrd"
OUTPUT_NAME_WARP_INV="${DATASET_NAME}_${REF_NAME}_registered_warp_inv.nrrd"

echo "make_initial_affine"
cmtk make_initial_affine --centers-of-mass $IMG_PATH $REF_PATH $RESULT_PATH/initial.list

echo "registration"
cmtk registration --initial $RESULT_PATH/initial.list --nmi --threads 64 --dofs 6 --dofs 12 --nmi --exploration 8 --accuracy 0.8 -o $RESULT_PATH/affine.list $IMG_PATH $REF_PATH

echo "warp"
cmtk warp --nmi --threads 64 --jacobian-weight 0 --fast -e 18 --grid-spacing 100 --energy-weight 1e-1 --refine 4 --coarsest 10 --ic-weight 0 --output-intermediate --accuracy 0.5 -o $RESULT_PATH/warp.list $RESULT_PATH/affine.list

echo "reformatx"
cmtk reformatx --pad-out 0 -o $RESULT_PATH/$OUTPUT_NAME_INIT --floating $REF_PATH $IMG_PATH $RESULT_PATH/initial.list
cmtk reformatx --pad-out 0 -o $RESULT_PATH/$OUTPUT_NAME_RIGID --floating $REF_PATH $IMG_PATH $RESULT_PATH/affine.list
cmtk reformatx --pad-out 0 -o $RESULT_PATH/$OUTPUT_NAME_WARP --floating $REF_PATH $IMG_PATH $RESULT_PATH/warp.list

# cmtk reformatx --pad-out 0 -o $RESULT_PATH/$OUTPUT_NAME_INIT_INV --floating $IMG_PATH $REF_PATH --inverse $RESULT_PATH/initial.list
# cmtk reformatx --pad-out 0 -o $RESULT_PATH/$OUTPUT_NAME_RIGID_INV --floating $IMG_PATH $REF_PATH --inverse $RESULT_PATH/affine.list
# cmtk reformatx --pad-out 0 -o $RESULT_PATH/$OUTPUT_NAME_WARP_INV --floating $IMG_PATH $REF_PATH --inverse $RESULT_PATH/warp.list




echo "streamxform"
cmtk streamxform $RESULT_PATH/warp.list < $DATA_PATH/${DATASET_NAME}_nmf_loc.txt > $RESULT_PATH/${DATASET_NAME}_nmf_loc_${REF_NAME}_registered.txt

echo "done!"