#ifndef DACSFM_SRC_MERGE_BUNDLEADJUSTMENT_H_
#define DACSFM_SRC_MERGE_BUNDLEADJUSTMENT_H_

#include <memory>
#include <unordered_set>
#include <Eigen/Core>
#include <ceres/ceres.h>

#include "optim/bundle_adjustment.h"
#include "merge/similaritytransform.h"
#include "base/camera_rig.h"
#include "base/reconstruction.h"

using namespace colmap;
namespace  DACSFM
{
    
    struct BundleAdjustmentOptions
    {
        //enum class LossFunctionType{TRIVIAL,SOFT_L1,CAUCHY};
        colmap::BundleAdjustmentOptions::LossFunctionType loss_function_type = colmap::BundleAdjustmentOptions::LossFunctionType::TRIVIAL;

        double loss_function_scale = 1.0;

        bool refine_focal_length = true;

        bool refine_principal_point = false;

        bool refine_extra_params = true;

        bool refine_extrinsics = true;

        bool print_summary = true;

        int min_num_residuals_for_multi_threading = 50000;

        ceres::Solver::Options solver_options;

        BundleAdjustmentOptions(){
            solver_options.function_tolerance = 0.0;
            solver_options.gradient_tolerance = 0.0;
            solver_options.parameter_tolerance = 0.0;
            solver_options.minimizer_progress_to_stdout = false;
            solver_options.max_num_iterations = 100;
            solver_options.max_linear_solver_iterations = 200;
            solver_options.max_num_consecutive_invalid_steps = 10;
            solver_options.num_threads = -1;
    #if CERES_VERSION_MAJOR < 2
                solver_options.num_linear_solver_threads = -1;
    #endif
        }
        BundleAdjustmentOptions(const colmap::BundleAdjustmentOptions& ba_options)
        {
            loss_function_type = ba_options.loss_function_type;
            loss_function_scale = ba_options.loss_function_scale;
            refine_focal_length = ba_options.refine_focal_length;
            refine_principal_point = ba_options.refine_principal_point;
            refine_extra_params = ba_options.refine_extra_params;
            refine_extrinsics = ba_options.refine_extrinsics;
            print_summary = ba_options.print_summary;
            min_num_residuals_for_multi_threading = ba_options.min_num_residuals_for_multi_threading;
            solver_options = ba_options.solver_options;
        }
        ceres::LossFunction* CreateLossFunction() const;

        bool Check() const;

    };

    class BundleAdjustmentConfig{
        public:
        BundleAdjustmentConfig();


        size_t NumImages() const;
        size_t NumPoints() const;
        size_t NumConstantCameras() const;
        size_t NumConstantPoses() const;
        size_t NumConstantTvecs() const;
        size_t NumVariablePoints() const;
        size_t NumConstantPoints() const;

        // Determine the number of residuals for the given reconstruction. The number
        // of residuals equals the number of observations times two.
        size_t NumResiduals(const Reconstruction& reconstruction) const;

        // Add / remove images from the configuration.
        void AddImage(const image_t image_id);
        bool HasImage(const image_t image_id) const;
        void RemoveImage(const image_t image_id);

        // Set cameras of added images as constant or variable. By default all
        // cameras of added images are variable. Note that the corresponding images
        // have to be added prior to calling these methods.
        void SetConstantCamera(const camera_t camera_id);
        void SetVariableCamera(const camera_t camera_id);
        bool IsConstantCamera(const camera_t camera_id) const;

        // Set the pose of added images as constant. The pose is defined as the
        // rotational and translational part of the projection matrix.
        void SetConstantPose(const image_t image_id);
        void SetVariablePose(const image_t image_id);
        bool HasConstantPose(const image_t image_id) const;

        // Set the translational part of the pose, hence the constant pose
        // indices may be in [0, 1, 2] and must be unique. Note that the
        // corresponding images have to be added prior to calling these methods.
        void SetConstantTvec(const image_t image_id, const std::vector<int>& idxs);
        void RemoveConstantTvec(const image_t image_id);
        bool HasConstantTvec(const image_t image_id) const;

        // Add / remove points from the configuration. Note that points can either
        // be variable or constant but not both at the same time.
        void AddVariablePoint(const point3D_t point3D_id);
        void AddConstantPoint(const point3D_t point3D_id);
        bool HasPoint(const point3D_t point3D_id) const;
        bool HasVariablePoint(const point3D_t point3D_id) const;
        bool HasConstantPoint(const point3D_t point3D_id) const;
        void RemoveVariablePoint(const point3D_t point3D_id);
        void RemoveConstantPoint(const point3D_t point3D_id);

        // Access configuration data.
        const std::unordered_set<image_t>& Images() const;
        const std::unordered_set<point3D_t>& VariablePoints() const;
        const std::unordered_set<point3D_t>& ConstantPoints() const;
        const std::vector<int>& ConstantTvec(const image_t image_id) const;

        private:
        std::unordered_set<camera_t> constant_camera_ids_;
        std::unordered_set<image_t> image_ids_;
        std::unordered_set<point3D_t> variable_point3D_ids_;
        std::unordered_set<point3D_t> constant_point3D_ids_;
        std::unordered_set<image_t> constant_poses_;
        std::unordered_map<image_t, std::vector<int>> constant_tvecs_;
    };

    class BundleAdjustment
    {
        public :
        BundleAdjustment(const BundleAdjustmentOptions& options,
                         const std::unordered_map<Reconstruction*,BundleAdjustmentConfig*>& configs);

        bool Solve(std::vector<Reconstruction*>& reconstructions,
                   std::unordered_map<Reconstruction*,SimilarityTransform>& srtform_to_global);

        bool SolveOGBA(Reconstruction* reconstruction);

        private:
        // Largest reconstruction should be placed on the first(subscript 0)
        void SetUpInter(std::vector<Reconstruction*>& reconstructions,
                        std::unordered_map<Reconstruction*,SimilarityTransform>& srtform_to_global,
                        ceres::LossFunction* loss_function);
        void SetUpIntra(Reconstruction* reconstruction,
                   ceres::LossFunction* loss_function);
        
        void AddImageToProblem(const image_t image_id,Reconstruction* reconstruction,
                               ceres::LossFunction* loss_function,BundleAdjustmentConfig* config_);

        void AddPointToProblem(const point3D_t point3D_id,
                               Reconstruction* reconstruction,
                               ceres::LossFunction* loss_function,BundleAdjustmentConfig* config_);
        void AddOverlappingImageSingleSRT(Reconstruction* src_reconstruction,Reconstruction* ref_reconstruction,
                             SimilarityTransform* similaritytransform,ceres::LossFunction* loss_function);
        void AddOverlappingImageDoubleSRT(Reconstruction* src_reconstruction,Reconstruction* ref_reconstruction,
                                           SimilarityTransform* src_similaritytransform,SimilarityTransform* ref_similaritytransform,
                                           ceres::LossFunction* loss_function);
        
        protected:
        void ParameterizeCameras(Reconstruction* reconstruction,BundleAdjustmentConfig* config_);
        void ParameterizePoints(Reconstruction* reconstruction,BundleAdjustmentConfig* config_);

        const BundleAdjustmentOptions options_;
        std::unordered_map<Reconstruction*,BundleAdjustmentConfig*> configs_;

        std::unique_ptr<ceres::Problem> problem_;
        ceres::Solver::Summary summary_;
        
        std::unordered_set<camera_t> camera_ids_;
        std::unordered_map<point3D_t,size_t> point3D_num_observations_;
        std::unordered_set<image_t> image_ids_;

    };

}

#endif