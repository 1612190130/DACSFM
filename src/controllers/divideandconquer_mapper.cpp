#include "controllers/divideandconquer_mapper.h"
#include "feature/extraction.h"
#include "feature/matching.h"

#include "util/option_manager.h"
#include "util/misc.h"
#include "base/database.h"

#include <iostream>
#include <utility>
#include <vector>
#include <fstream>
#include <ceres/rotation.h>
#include <boost/filesystem.hpp>


using namespace colmap;

namespace DACSFM
{
    DivideAndConquerMapperController::DivideAndConquerMapperController(const Options& options
    , ReconstructionManager* reconstruction_manager,SceneCut::Options cut_options):
        options_(options),reconstruction_manager_(reconstruction_manager),active_thread_(nullptr),
        cut_options_(cut_options)
        
        
        
    {
        CHECK(ExistsDir(options_.workspace_path));
        CHECK(ExistsDir(options_.image_path));
        CHECK_NOTNULL(reconstruction_manager_);

        option_manager_.AddAllOptions();

        *option_manager_.image_path = options_.image_path;
        *option_manager_.database_path =
            JoinPaths(options_.workspace_path, "database.db");

        if (options_.data_type == DataType::VIDEO) {
            option_manager_.ModifyForVideoData();
        } else if (options_.data_type == DataType::INDIVIDUAL) {
            option_manager_.ModifyForIndividualData();
        } else if (options_.data_type == DataType::INTERNET) {
            option_manager_.ModifyForInternetData();
        } else {
            LOG(FATAL) << "Data type not supported";
        }

        CHECK(ExistsCameraModelWithName(options_.camera_model));

        if (options_.quality == Quality::LOW) {
            option_manager_.ModifyForLowQuality();
        } else if (options_.quality == Quality::MEDIUM) {
            option_manager_.ModifyForMediumQuality();
        } else if (options_.quality == Quality::HIGH) {
            option_manager_.ModifyForHighQuality();
        } else if (options_.quality == Quality::EXTREME) {
            option_manager_.ModifyForExtremeQuality();
        }

        option_manager_.sift_extraction->num_threads = options_.num_threads;
        option_manager_.sift_matching->num_threads = options_.num_threads;
        option_manager_.mapper->num_threads = options_.num_threads;
        option_manager_.poisson_meshing->num_threads = options_.num_threads;

        ImageReaderOptions reader_options = *option_manager_.image_reader;
        reader_options.database_path = *option_manager_.database_path;
        reader_options.image_path = *option_manager_.image_path;
        if (!options_.mask_path.empty()) {
            reader_options.mask_path = options_.mask_path;
            option_manager_.image_reader->mask_path = options_.mask_path;
        }
        reader_options.single_camera = options_.single_camera;
        reader_options.camera_model = options_.camera_model;

        option_manager_.sift_extraction->use_gpu = options_.use_gpu;
        option_manager_.sift_matching->use_gpu = options_.use_gpu;

        option_manager_.sift_extraction->gpu_index = options_.gpu_index;
        option_manager_.sift_matching->gpu_index = options_.gpu_index;
        option_manager_.patch_match_stereo->gpu_index = options_.gpu_index;

        feature_extractor_.reset(new SiftFeatureExtractor(
            reader_options, *option_manager_.sift_extraction));

        exhaustive_matcher_.reset(new ExhaustiveFeatureMatcher(
            *option_manager_.exhaustive_matching, *option_manager_.sift_matching,
            *option_manager_.database_path));

        if (!options_.vocab_tree_path.empty()) {
            option_manager_.sequential_matching->loop_detection = true;
            option_manager_.sequential_matching->vocab_tree_path =
                options_.vocab_tree_path;
        }

        sequential_matcher_.reset(new SequentialFeatureMatcher(
            *option_manager_.sequential_matching, *option_manager_.sift_matching,
            *option_manager_.database_path));

        if (!options_.vocab_tree_path.empty()) {
            option_manager_.vocab_tree_matching->vocab_tree_path =
                options_.vocab_tree_path;
            vocab_tree_matcher_.reset(new VocabTreeFeatureMatcher(
                *option_manager_.vocab_tree_matching, *option_manager_.sift_matching,
                *option_manager_.database_path));
        }
  }

    void DivideAndConquerMapperController::Run() {
        RunFeatureExtraction();

        RunFeatureMatching();


        std::vector<image_t>& image_ids;
        std::vector<std::pair<image_t, image_t>>& image_pairs;
        std::vector<int>& num_inliers;
        std::unordered_map<image_t,std::string>& image_id_to_name;
        LoadDatabase(image_ids,image_pairs,num_inliers,image_id_to_name);



        /*************************************/
        /********  cut images ****************/
        /*************************************/

        std::vector<ImageCluster> inter_clusters = CutScenes(image_ids,image_pairs,num_inliers);

        /*******************************************/
        /********reconstruct partitions*************/
        /*******************************************/

        std::unordered_map<const ImageCluster*,ReconstructionManager> reconstruction_managers;
        std::vector<Reconstruction*> reconstructions;
        ReconstructPartitions(image_id_to_name,inter_clusters,reconstruction_managers,reconstructions);

        

    }

    void DivideAndConquerMapperController::RunFeatureExtraction() {
        CHECK(feature_extractor_);
        active_thread_ = feature_extractor_.get();
        feature_extractor_->Start();
        feature_extractor_->Wait();
        feature_extractor_.reset();
        active_thread_ = nullptr;
    }

    void DivideAndConquerMapperController::RunFeatureMatching() {
        Thread* matcher = nullptr;
        if (options_.data_type == DataType::VIDEO) {
            matcher = sequential_matcher_.get();
        } else if (options_.data_type == DataType::INDIVIDUAL ||
                    options_.data_type == DataType::INTERNET) {
            Database database(*option_manager_.database_path);
            const size_t num_images = database.NumImages();
            if (options_.vocab_tree_path.empty() || num_images < 200) {
            matcher = exhaustive_matcher_.get();
            } else {
            matcher = vocab_tree_matcher_.get();
            }
        }

        CHECK(matcher);
        active_thread_ = matcher;
        matcher->Start();
        matcher->Wait();
        exhaustive_matcher_.reset();
        sequential_matcher_.reset();
        vocab_tree_matcher_.reset();
        active_thread_ = nullptr;
    }



    void DivideAndConquerMapperController::LoadDatase(std::vector<image_t>& image_ids,std::vector<std::pair<image_t, image_t>>& image_pairs,
                                                      std::vector<int>& num_inliers,std::unordered_map<image_t,std::string>& image_id_to_name )
    {
        PrintHeading1("Loading database");

        Database database(*option_manager_.database_path);

        const auto images = database.ReadAllImages();
        for(auto& image：images){
            image_id_to_name.emplace(image.ImageId(),image.Name());
            image_ids.push_back(image.ImageId());
        }

        database.ReadTwoViewGeometryNumInliers(image_pairs,num_inliers);
    }

    std::vector<ImageCluster> DivideAndConquerMapperController::CutScenes(const std::vector<image_t>& image_ids,
                                                      const std::vector<std::pair<image_t, image_t>>& image_pairs,
                                                      const std::vector<int>& num_inliers)
    {
            ImageCluster image_cluster;
            image_cluster.image_ids = image_ids;

            for(int i=0;i<image_pairs.size();i++)
            {
                image_cluster.edges[image_pairs[i]]=num_inliers[i];
            }

            scene_cut_ = std::unique_ptr<SceneCut>(new SceneCut(cut_options_,image_cluster));

            scene_cut_->Cut();
            scene_cut_->Expand();

            std::vector<ImageCluster> inter_clusters = scene_cut_->GetInterClusters();

            for(auto cluster:inter_clusters)
            {
                cluster.showInfo();
            }

            return inter_clusters;
            
    }

    void DivideAndConquerMapperController::ReconstructPartitions(const std::unordered_map<image_t,std::string>& image_id_to_name,
                                                                 std::vector<ImageCluster>& inter_clusters,
                                                                 std::unordered_map<const ImageCluster*,ReconstructManager>& reconstruction_managers,
                                                                 std::vector<Reconstruction*>& reconstructions)
    {
        const int kMaxNumThreads = -1;
        const int num_eff_threads = GetEffectiveNumThreads(kMaxNumThreads);
        const int kDefaultNumWorkers = 8;
        const int num_eff_workers = options_.num_workers < 1? std::min(static_cast<int>(inter_clusters.size()),
                                                                       std::min(kDefaultNumWorkers,num_eff_threads))
                                                                       :options_.num_workers;

        

        const int num_threads_per_worker = std::max(1,num_eff_threads / num_eff_workers);

        auto ReconstructCluster = [&,this](const ImageCluster& cluster,
                                           ReconstructionManager* reconstruction_manager)
        {
            IncrementalMapperOptions custom_options = mapper_options_;
            custom_options.max_model_overlap = 3;

            custom_options.init_num_trials = options_.init_num_trials;
            custom_options.num_threads = num_threads_per_worker ; 

            for(const auto image_id:cluster.image_ids)
            {
                custom_options.image_names.insert(image_id_to_name.at(image_id));
            }

            IncrementalMapperController mapper(&custom_options,options_.image_path,
                                               options_.database_path,reconstruction_manager);
            mapper.Start();
            mapper.Wait();

        }

        const auto cmp = [](const ImageCluster& cluster1,const ImageCluster& cluster2){
            return cluater1.image_ids.size()>cluster2.image_ids.size();
        }

        std::sort(inter_clusters.begin(),inter_cluster.end(),cmp);

        reconstruction_managers.reserve(inter_cluster.size());

        bool is_recons_exist = IsPartialReconsExist(reconstructions);
        if(is_recons_exist)
        {
                LOG(INFO) << "Loaded from previous reconstruction partitions.";
        }
        else
        {
                ThreadPool thread_pool(num_eff_workers);
                for(const auto& cluster:inter_clusters)
                {
                    thread_pool.AddTask(ReconstructCluster,cluster,&reconstruction_managers[&cluster]);
                }
                thread_pool.Wait();

                for(const auto & cluster: inter_clusters)
                {
                    auto& recon_manager = reconstruction_managers.at(cluster);
                    for(size_t i=0;i<recon_manager.Size();i++)
                    {
                            reconstructions.push_back(&recon_manager.Get(i));
                    }
                }

                for(int i=0;i<reconstructions.size();i++)
                {
                    const std::string reconstruction_path = JoinPaths(options_.output_path,"partition_"+std::to_string(i));
                    CreateDirIfNotExists(reconstruction_path);
                    reconstructions[i]->Writer(reconstruction_path);
                }
        }

    }


    bool DivideAndConquerMapperController::IsPartialReconsExist(std::vector<Reconstruction*>& recons)
    {
        // const auto sparse_path = JoinPaths(options_.workspace_path, "sparse");
        // if (ExistsDir(sparse_path)) {
        //     auto dir_list = GetDirList(sparse_path);
        //     std::sort(dir_list.begin(), dir_list.end());
        //     if (dir_list.size() > 0) {
        //     std::cout << std::endl
        //                 << "WARNING: Skipping sparse reconstruction because it is "
        //                 "already computed"
        //                 << std::endl;
        //     for (const auto& dir : dir_list) {
        //         reconstruction_manager_->Read(dir);
        //     }
        //     return;
        //     }
        // }
        return false;


    }


}