#!/bin/bash

IMG_PATH=$1
REF_PATH=$2
RESULT_PATH=$3



OUTPUT_NAME_INIT="image_registered_initial.nrrd"
OUTPUT_NAME_RIGID="image_registered_rigid.nrrd"
OUTPUT_NAME_WARP="image_registered_warp.nrrd"


echo "make_initial_affine"
cmtk make_initial_affine --identity $REF_PATH $IMG_PATH $RESULT_PATH/initial.list

# echo "registration"
# cmtk registration --initial $RESULT_PATH/initial.list  --threads 64 --nmi --dofs 6,9 --exploration 1 --accuracy 0.8 -o $RESULT_PATH/affine.list $REF_PATH $IMG_PATH

# echo "warp"
# cmtk warp --nmi --threads 64 --jacobian-weight 0 --fast -e 18 --grid-spacing 100 --energy-weight 1e-1 --refine 4 --coarsest 10 --ic-weight 0 --output-intermediate --accuracy 0.5 -o $RESULT_PATH/warp.list $RESULT_PATH/affine.list

echo "reformatx"
cmtk reformatx --pad-out 0 -o $RESULT_PATH/$OUTPUT_NAME_INIT --floating $IMG_PATH $REF_PATH $RESULT_PATH/initial.list
# cmtk reformatx --pad-out 0 -o $RESULT_PATH/$OUTPUT_NAME_RIGID --floating $IMG_PATH $REF_PATH $RESULT_PATH/affine.list
# cmtk reformatx --pad-out 0 -o $RESULT_PATH/$OUTPUT_NAME_WARP --floating $IMG_PATH $REF_PATH $RESULT_PATH/warp.list

echo "done!"