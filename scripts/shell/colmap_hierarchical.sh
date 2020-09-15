DATASET_PATH=$1

#num_images_ub=$2
#log_folder=$3
#completeness_ratio=$4
#VOC_TREE_PATH=$5
# image_overlap=$3
# max_num_cluster_pairs=$4


colmap hierarchical_mapper \
--database_path=$DATASET_PATH/database.db \
--image_path=$DATASET_PATH/images \
--output_path=$DATASET_PATH/cut_colmap_for_compare \
#--leaf_max_num_images=50 \
#--image_overlap=50 \

