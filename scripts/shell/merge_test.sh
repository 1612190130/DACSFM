

#num_images_ub=$2
#log_folder=$3
#completeness_ratio=$4
#VOC_TREE_PATH=$5
# image_overlap=$3
# max_num_cluster_pairs=$4

DATASET_PATH=$1

/home/yzc/Projects/DACSFM_bc/colmap/build/src/exe/colmap hierarchical_mapper \
--database_path=$DATASET_PATH/database.db \
--image_path=$DATASET_PATH/images \
--output_path=$DATASET_PATH/cut \
--branching=2 \
--merge_useba=true \
--Mapper.ba_refine_focal_length=false \
--Mapper.ba_refine_principal_point=false \
--leaf_max_num_images=550 \
| tee $DATASET_PATH/cut/log.txt \

<<!
#--image_overlap=0 \
#--leaf_max_num_images=150 \

#--Mapper.ba_refine_extra_params=false \
!


