#ifndef SRC_CONTROLLERS_DIVIDEANDCONQUER_MAPPER_H_
#define SRC_CONTROLLERS_DIVIDEANDCONQUER_MAPPER_H_

#include "util/threading.h"
#include "base/reconstruction_manager.h"
#include "cut/scenecut.h"
#include "controllers/incremental_mapper_controller.h"
#include "cut/scenecut.h"


using namespace colmap;

namespace DACSFM{

class DivideAndConquerMapperController:public Thread
{
enum class DataType { INDIVIDUAL, VIDEO, INTERNET };
  enum class Quality { LOW, MEDIUM, HIGH, EXTREME };
  enum class Mesher { POISSON, DELAUNAY };

  struct Options {
    // The path to the workspace folder in which all results are stored.
    std::string workspace_path;

    // The path to the image folder which are used as input.
    std::string image_path;

    // The path to the mask folder which are used as input.
    std::string mask_path;

    // The path to the vocabulary tree for feature matching.
    std::string vocab_tree_path;

    // The type of input data used to choose optimal mapper settings.
    DataType data_type = DataType::INDIVIDUAL;

    // Whether to perform low- or high-quality reconstruction.
    Quality quality = Quality::HIGH;

    // Whether to use shared intrinsics or not.
    bool single_camera = false;

    // Which camera model to use for images.
    std::string camera_model = "SIMPLE_RADIAL";

    // Whether to perform sparse mapping.
    bool sparse = true;

// Whether to perform dense mapping.
#ifdef CUDA_ENABLED
    bool dense = true;
#else
    bool dense = false;
#endif

    // The meshing algorithm to be used.
    Mesher mesher = Mesher::POISSON;

    // The number of threads to use in all stages.
    int num_threads = -1;

    // Whether to use the GPU in feature extraction and matching.
    bool use_gpu = true;

    // Index of the GPU used for GPU stages. For multi-GPU computation,
    // you should separate multiple GPU indices by comma, e.g., "0,1,2,3".
    // By default, all GPUs will be used in all stages.
    std::string gpu_index = "-1";
    /************************************************************/
    int num_workers = -1;

    int init_num_trials = 10;

    std::String ouput_path;

  };
    DivideAndConquerMapperController(
      const Options& options, ReconstructionManager* reconstruction_manager);
private:
    
    void Run() override;
    void RunFeatureExtraction();
    void RunFeatureMatching();

    const Options options_;
    OptionManager option_manager_;
    Thread* active_thread_;
    SceneCut::Options cut_options_;

    std::unique_ptr<Thread> feature_extractor_;
    std::unique_ptr<Thread> exhaustive_matcher_;
    std::unique_ptr<Thread> sequential_matcher_;
    std::unique_ptr<Thread> vocab_tree_matcher_;

    /*************************************/
    /*************************************/
    /*************************************/
    void LoadDatase(std::vector<image_t>& image_ids,std::vector<std::pair<image_t, image_t>>& image_pairs,
                                                      std::vector<int>& num_inliers);
    std::vector<ImageCluster> CutScenes(const std::vector<image_t>& image_ids,const std::vector<std::pair<image_t, image_t>>& image_pairs,
                                                      const std::vector<int>& num_inliers)

};

};


#endif