
#include "merge/bundleadjustment_controller.h"
#include <ceres/ceres.h>
#include "util/misc.h"
#include "merge/bundle_adjustment.h"

namespace DACSFM
{

    
    class BundleAdjustmentIterationCallback: public ceres::IterationCallback{
        public:
        explicit BundleAdjustmentIterationCallback(Thread* thread)
        :thread_(thread){}

        virtual ceres::CallbackReturnType operator()(const ceres::IterationSummary& summary){
            CHECK_NOTNULL(thread_);
            thread_->BlockIfPaused();
            if(thread_->IsStopped()){
                return ceres::SOLVER_TERMINATE_SUCCESSFULLY;
            }else{
                return ceres::SOLVER_CONTINUE;
            }
        }


        private:
        Thread* thread_;
    };


    MyBundleAdjustmentController::MyBundleAdjustmentController(const OptionManager& options,
                             Reconstruction* reconstruction):
    options_(options),reconstruction_(reconstruction){}


    void MyBundleAdjustmentController::Run()
    {

        CHECK_NOTNULL(reconstruction_);

        PrintHeading1("Global bundle adjustment");
        std::unordered_map<Reconstruction*,BundleAdjustmentConfig*> ba_configs;

        const std::vector<image_t>& reg_image_ids = reconstruction_->RegImageIds();

        if (reg_image_ids.size() < 2) {
            std::cout << "ERROR: Need at least two views." << std::endl;
            return;
        }

        // Avoid degeneracies in bundle adjustment.
        //reconstruction_->FilterObservationsWithNegativeDepth();

        BundleAdjustmentOptions ba_options = *options_.bundle_adjustment;
        ba_options.solver_options.minimizer_progress_to_stdout = true;

        BundleAdjustmentIterationCallback iteration_callback(this);
        ba_options.solver_options.callbacks.push_back(&iteration_callback);

        // Configure bundle adjustment.
        BundleAdjustmentConfig ba_config;
        for (const image_t image_id : reg_image_ids) {
            ba_config.AddImage(image_id);
        }
        ba_config.SetConstantPose(reg_image_ids[0]);
        ba_config.SetConstantTvec(reg_image_ids[1], {0});
        ba_configs.insert(std::make_pair(reconstruction_ ,&ba_config));
        // Run bundle adjustment.
        BundleAdjustment bundle_adjuster(ba_options, ba_configs);
        bundle_adjuster.SolveOGBA(reconstruction_);

        GetTimer().PrintMinutes();

    }


    void MyBundleAdjustmentController::Solve(const OptionManager& options,
                             std::vector<Reconstruction*>& reconstructions,
                             std::unordered_map<Reconstruction*,SimilarityTransform>& srtforms)
    {
        auto start = std::chrono::system_clock::now();
        PrintHeading1("Global bundle adjustment");

        std::unordered_map<Reconstruction*,BundleAdjustmentConfig*> ba_configs;

        BundleAdjustmentOptions ba_options(*options.bundle_adjustment);
        ba_options.solver_options.minimizer_progress_to_stdout = true;

        //BundleAdjustmentIterationCallback  ba_iterationcallback(this);
        //ba_options.solver_options.callbacks.push_back(&ba_iterationcallback);
        

        for(auto reconstruction:reconstructions)
        {
            CHECK_NOTNULL(reconstruction);

            const std::vector<image_t> reg_image_ids = reconstruction->RegImageIds();

            if(reg_image_ids.size() < 2){
                std::cout<<"ERROR: Need at least two views."<<std::endl;
                return;
            }

            reconstruction->FilterObservationsWithNegativeDepth();

            BundleAdjustmentConfig* ba_config = new BundleAdjustmentConfig();
            for(const image_t image_id:reg_image_ids){
                ba_config->AddImage(image_id);
                ba_config->SetConstantCamera(reconstruction->Image(image_id).CameraId()); //?
            }
            
            ba_config->SetConstantPose(reg_image_ids[0]);
            ba_config->SetConstantTvec(reg_image_ids[1],{0});
            ba_configs.insert(std::make_pair(reconstruction ,ba_config));
        }

        BundleAdjustment ba_solver(ba_options,ba_configs);
        ba_solver.Solve(reconstructions,srtforms);
        PrintHeading2("bundle adjustment complete");
        for(auto config:ba_configs)
        {
            delete config.second;
        }

        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
         std::cout<< "elapsed time: " << elapsed_seconds.count()/60.0 << "m\n";
        // PrintHeading1("Normalize");

        // for(auto reconstruction:reconstructions){
        //     reconstruction_->Normalize();          //???
        // }
        

        //GetTimer().PrintMinutes();

    }
}