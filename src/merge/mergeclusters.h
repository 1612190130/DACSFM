#ifndef DACSFM_SRC_MERGE_MERGECLUSTERS_H_
#define DACSFM_SRC_MERGE_MERGECLUSTERS_H_


#include "controllers/hierarchical_mapper.h"

#include "base/scene_clustering.h"

#include "base/reconstruction.h"
#include "util/option_manager.h"
#include "merge/similaritytransform.h"

using namespace colmap;
namespace DACSFM
{

    struct MergeControllerOptions
    {
        double kMinInlierObservations = 0.3;
        double kMaxReprojError = 8.0;
        bool useba = true;
        std::string output_path;
        MergeControllerOptions(){}
        MergeControllerOptions(std::string outputpath)
        {
            output_path = outputpath; 
        }
    };
    class MergeController
    {
public:
    MergeController();

    MergeController(MergeControllerOptions mergeoptions);

    void MergeClustersOnTreeRun(
    const SceneClustering::Cluster& cluster,
    std::unordered_map<const SceneClustering::Cluster*, ReconstructionManager>*
        reconstruction_managers);

    void MergeClustersOnTree(
    const SceneClustering::Cluster& cluster,
    std::unordered_map<const SceneClustering::Cluster*, ReconstructionManager>*
        reconstruction_managers,const OptionManager& options);
    
    void MergeClustersIntoLargestReconstructionRun(
    std::vector<Reconstruction*>& reconstructions,
    ReconstructionManager* reconstruction_manager,
    std::vector<std::string>& cluster_name);

    void MergeClustersIntoLargestReconstruction(std::vector<Reconstruction*>& reconstructions,
                                std::vector<std::string>& cluster_name);


    void MergeClusterwithBARun(std::vector<Reconstruction*>& reconstructions,
    ReconstructionManager* reconstruction_manager,
    std::vector<std::string>& cluster_name);

    void MergeClustersintoMaximum(
    std::vector<Reconstruction*>& reconstructions,
    ReconstructionManager* reconstruction_manager,
    std::unordered_map<Reconstruction*,std::string>& cluster_name);

    void MergeClustersintoGlobal(std::vector<Reconstruction*>& reconstructions,ReconstructionManager* reconstruction_manager,
    std::unordered_map<Reconstruction*,std::string>& cluster_name);

private:
    void CallBundleAdjustment(std::vector<Reconstruction*>& reconstructions,const double kMaxReprojError,
                              std::unordered_map<Reconstruction*,const SimilarityTransform3>& similaritytransform3_to_global);
    bool EstimateSRT(Reconstruction& main_reconstruction,const Reconstruction& sub_reconstruction,
                           const double max_reproj_error,SimilarityTransform& srt);


    bool MergeReconstructionsAmbiguously(Reconstruction& main_reconstruction,const Reconstruction& sub_reconstruction,
                           const double max_reproj_error,const SimilarityTransform3& tform);

    size_t FilterPoints3DWithLargeReprojectionError(Reconstruction& reconstruction,
    const double max_reproj_error,const std::unordered_set<point3D_t>& point3D_ids);


    bool MergeImagesToObtainSRT(Reconstruction& main_reconstruction,const Reconstruction& sub_reconstruction,
                           const double max_reproj_error,SimilarityTransform* srt);

    bool MergeReconstructionsUnambiguously(Reconstruction& main_reconstruction,const Reconstruction& sub_reconstruction,
                           const double max_reproj_error,const SimilarityTransform& tform);

    bool MergeAllForSRT(Reconstruction& main_reconstruction,const Reconstruction& sub_reconstruction,
                           const double max_reproj_error,SimilarityTransform& srt);

    double MeanSquaredErrorByRecons(const Reconstruction& reconstruction);
    void WriteWithoutBA(const Reconstruction& reconstruction);
    void WriteKittiPose(const std::string filename, const Reconstruction& reconstruction);
    void WriteEuRoC(const std::string filename,const Reconstruction& reconstruction);

    
    MergeControllerOptions mergeoptions_;
    };

    //void CallBundleAdjustment(const OptionManager& options, Reconstruction* reconstruction);
}

#endif