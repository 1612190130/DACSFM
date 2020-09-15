DATASET_PATH=$1

#num_images_ub=$2
#log_folder=$3
#completeness_ratio=$4
#VOC_TREE_PATH=$5
# image_overlap=$3
# max_num_cluster_pairs=$4


/home/intern/colmap/build/src/exe/colmap feature_extractor \
--database_path=$DATASET_PATH/database.db \
--image_path=$DATASET_PATH/images \
--SiftExtraction.num_threads=8 \
--SiftExtraction.use_gpu=1 \
--SiftExtraction.gpu_index=-1 \
--ImageReader.single_camera=true \
--ImageReader.camera_params=718.856,607.1928,185.2157,0

<<!
/home/intern/colmap/build/src/exe/colmap exhaustive_matcher \
--database_path=$DATASET_PATH/database.db \
--SiftMatching.num_threads=8 \
--SiftMatching.use_gpu=1 \
--SiftMatching.gpu_index=-1
!

## Or use vocabulary tree matcher
/home/intern/colmap/build/src/exe/colmap vocab_tree_matcher \
--database_path=$DATASET_PATH/database.db \
--SiftMatching.num_threads=8 \
--SiftMatching.use_gpu=1 \
--SiftMatching.gpu_index=-1 \
--VocabTreeMatching.num_images=100 \
--VocabTreeMatching.num_nearest_neighbors=5 \
--VocabTreeMatching.vocab_tree_path=/home/intern/SFM/vocab_tree_flickr100K_words32K.bin

/home/intern/colmap/build/src/exe/colmap hierarchical_mapper \
--database_path=$DATASET_PATH/database.db \
--image_path=$DATASET_PATH/images \
--output_path=$DATASET_PATH/cut \
#--Mapper.ba_refine_focal_length=false \
#--Mapper.ba_refine_principal_point=false \
#--Mapper.ba_refine_extra_params=false \
#--leaf_max_num_images=500 \
#--image_overlap=50 \

