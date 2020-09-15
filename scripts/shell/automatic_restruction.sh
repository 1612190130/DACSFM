DATASET_PATH=$1

/home/intern/colmap/build/src/exe/colmap mapper \
--input_path=$DATASET_PATH \
--database_path=$DATASET_PATH/database.db \
--image_path=$DATASET_PATH \
--image_list_path=$DATASET_PATH \
--output_path=$DATASET_PATH/sparse \