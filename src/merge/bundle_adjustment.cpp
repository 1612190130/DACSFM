#include "merge/bundle_adjustment.h"

#include <iomanip>

#ifdef OPENMP_ENABLED
#include <omp.h>
#endif

#include "base/camera_models.h"
#include "merge/cost_function.h"
#include "base/projection.h"
#include "util/misc.h"
#include "util/threading.h"
#include "util/timer.h"

namespace DACSFM{

    /****************************************************/
    //  BundleAdjustmentOptions 
    /****************************************************/
    ceres::LossFunction* BundleAdjustmentOptions::CreateLossFunction() const{
        ceres::LossFunction* loss_function = nullptr;

        switch(loss_function_type){
            case colmap::BundleAdjustmentOptions::LossFunctionType::TRIVIAL:
                loss_function = new ceres::TrivialLoss();
                break;
            case colmap::BundleAdjustmentOptions::LossFunctionType::SOFT_L1:
                loss_function = new ceres::SoftLOneLoss(loss_function_scale);
                break;
            case colmap::BundleAdjustmentOptions::LossFunctionType::CAUCHY:
                loss_function = new ceres::CauchyLoss(loss_function_scale);
                break;
        }
        CHECK_NOTNULL(loss_function);
        return loss_function;
    }

    bool BundleAdjustmentOptions::Check() const {
        CHECK_OPTION_GE(loss_function_scale, 0);
        return true;
    }

    ////////////////////////////////////////////////////////////////////////////////
    // BundleAdjustmentConfig
    ////////////////////////////////////////////////////////////////////////////////

    BundleAdjustmentConfig::BundleAdjustmentConfig() {}

    size_t BundleAdjustmentConfig::NumImages() const { return image_ids_.size(); }

    size_t BundleAdjustmentConfig::NumPoints() const {
    return variable_point3D_ids_.size() + constant_point3D_ids_.size();
    }

    size_t BundleAdjustmentConfig::NumConstantCameras() const {
    return constant_camera_ids_.size();
    }

    size_t BundleAdjustmentConfig::NumConstantPoses() const {
    return constant_poses_.size();
    }

    size_t BundleAdjustmentConfig::NumConstantTvecs() const {
    return constant_tvecs_.size();
    }

    size_t BundleAdjustmentConfig::NumVariablePoints() const {
    return variable_point3D_ids_.size();
    }

    size_t BundleAdjustmentConfig::NumConstantPoints() const {
    return constant_point3D_ids_.size();
    }

    size_t BundleAdjustmentConfig::NumResiduals(
        const Reconstruction& reconstruction) const {
    // Count the number of observations for all added images.
    size_t num_observations = 0;
    for (const image_t image_id : image_ids_) {
        num_observations += reconstruction.Image(image_id).NumPoints3D();
    }

    // Count the number of observations for all added 3D points that are not
    // already added as part of the images above.

    auto NumObservationsForPoint = [this,
                                    &reconstruction](const point3D_t point3D_id) {
        size_t num_observations_for_point = 0;
        const auto& point3D = reconstruction.Point3D(point3D_id);
        for (const auto& track_el : point3D.Track().Elements()) {
        if (image_ids_.count(track_el.image_id) == 0) {
            num_observations_for_point += 1;
        }
        }
        return num_observations_for_point;
    };

    for (const auto point3D_id : variable_point3D_ids_) {
        num_observations += NumObservationsForPoint(point3D_id);
    }
    for (const auto point3D_id : constant_point3D_ids_) {
        num_observations += NumObservationsForPoint(point3D_id);
    }

    return 2 * num_observations;
    }

    void BundleAdjustmentConfig::AddImage(const image_t image_id) {
    image_ids_.insert(image_id);
    }

    bool BundleAdjustmentConfig::HasImage(const image_t image_id) const {
    return image_ids_.find(image_id) != image_ids_.end();
    }

    void BundleAdjustmentConfig::RemoveImage(const image_t image_id) {
    image_ids_.erase(image_id);
    }

    void BundleAdjustmentConfig::SetConstantCamera(const camera_t camera_id) {
    constant_camera_ids_.insert(camera_id);
    }

    void BundleAdjustmentConfig::SetVariableCamera(const camera_t camera_id) {
    constant_camera_ids_.erase(camera_id);
    }

    bool BundleAdjustmentConfig::IsConstantCamera(const camera_t camera_id) const {
    return constant_camera_ids_.find(camera_id) != constant_camera_ids_.end();
    }

    void BundleAdjustmentConfig::SetConstantPose(const image_t image_id) {
    CHECK(HasImage(image_id));
    CHECK(!HasConstantTvec(image_id));
    constant_poses_.insert(image_id);
    }

    void BundleAdjustmentConfig::SetVariablePose(const image_t image_id) {
    constant_poses_.erase(image_id);
    }

    bool BundleAdjustmentConfig::HasConstantPose(const image_t image_id) const {
    return constant_poses_.find(image_id) != constant_poses_.end();
    }

    void BundleAdjustmentConfig::SetConstantTvec(const image_t image_id,
                                                const std::vector<int>& idxs) {
    CHECK_GT(idxs.size(), 0);
    CHECK_LE(idxs.size(), 3);
    CHECK(HasImage(image_id));
    CHECK(!HasConstantPose(image_id));
    CHECK(!VectorContainsDuplicateValues(idxs))
        << "Tvec indices must not contain duplicates";
    constant_tvecs_.emplace(image_id, idxs);
    }

    void BundleAdjustmentConfig::RemoveConstantTvec(const image_t image_id) {
    constant_tvecs_.erase(image_id);
    }

    bool BundleAdjustmentConfig::HasConstantTvec(const image_t image_id) const {
    return constant_tvecs_.find(image_id) != constant_tvecs_.end();
    }

    const std::unordered_set<image_t>& BundleAdjustmentConfig::Images() const {
    return image_ids_;
    }

    const std::unordered_set<point3D_t>& BundleAdjustmentConfig::VariablePoints()
        const {
    return variable_point3D_ids_;
    }

    const std::unordered_set<point3D_t>& BundleAdjustmentConfig::ConstantPoints()
        const {
    return constant_point3D_ids_;
    }

    const std::vector<int>& BundleAdjustmentConfig::ConstantTvec(
        const image_t image_id) const {
    return constant_tvecs_.at(image_id);
    }

    void BundleAdjustmentConfig::AddVariablePoint(const point3D_t point3D_id) {
    CHECK(!HasConstantPoint(point3D_id));
    variable_point3D_ids_.insert(point3D_id);
    }

    void BundleAdjustmentConfig::AddConstantPoint(const point3D_t point3D_id) {
    CHECK(!HasVariablePoint(point3D_id));
    constant_point3D_ids_.insert(point3D_id);
    }

    bool BundleAdjustmentConfig::HasPoint(const point3D_t point3D_id) const {
    return HasVariablePoint(point3D_id) || HasConstantPoint(point3D_id);
    }

    bool BundleAdjustmentConfig::HasVariablePoint(
        const point3D_t point3D_id) const {
    return variable_point3D_ids_.find(point3D_id) != variable_point3D_ids_.end();
    }

    bool BundleAdjustmentConfig::HasConstantPoint(
        const point3D_t point3D_id) const {
    return constant_point3D_ids_.find(point3D_id) != constant_point3D_ids_.end();
    }

    void BundleAdjustmentConfig::RemoveVariablePoint(const point3D_t point3D_id) {
    variable_point3D_ids_.erase(point3D_id);
    }

    void BundleAdjustmentConfig::RemoveConstantPoint(const point3D_t point3D_id) {
    constant_point3D_ids_.erase(point3D_id);
    }

    ////////////////////////////////////////////////////////////////////////////////
    // BundleAdjustment
    ////////////////////////////////////////////////////////////////////////////////
    BundleAdjustment::BundleAdjustment(const BundleAdjustmentOptions& options,
                         const std::unordered_map<Reconstruction*,BundleAdjustmentConfig*>& configs)
    :options_(options),configs_(configs){
        CHECK(options_.Check());
    }


    bool BundleAdjustment::Solve(std::vector<Reconstruction*>& reconstructions,
                                 std::unordered_map<Reconstruction*,SimilarityTransform>& srtform_to_global)
    {
            
            CHECK(!problem_)<<"Cannot use the same BundleAdjuster multiple times";

            problem_.reset(new ceres::Problem());

            //ceres::LossFunction* loss_function = options_.CreateLossFunction();
            ceres::LossFunction* loss_function = new ceres::HuberLoss(2.0);
            SetUpInter(reconstructions,srtform_to_global,loss_function);
            for(const auto reconstruction:reconstructions){
                CHECK_NOTNULL(reconstruction);
                SetUpIntra(reconstruction,loss_function);
            }
            
            if(problem_->NumResiduals() == 0){
                return false;
            }

            //ceres::Solver::Options solver_options = options_.solver_options;
            ceres::Solver::Options solver_options ;
            solver_options.minimizer_progress_to_stdout = true;
            //solver_options.max_linear_solver_iterations = 200;
            //solver_options.max_num_consecutive_invalid_steps = 10;
            solver_options.num_threads = -1;

            const size_t kMaxNumImageDirectDenseSolver = 50;
            const size_t kMaxNUmImageDirectSparseSolver = 1000;
            const size_t num_images_all = image_ids_.size();
            size_t num_images = 0;
            for(auto config:configs_)
            {
                num_images += config.second->NumImages();
            }
            std::cout<<"num_images(overlap):"<<num_images<<"\n";
            std::cout<<"num_images_all:"<<num_images_all<<"\n";
            if(num_images <= kMaxNumImageDirectDenseSolver){
                solver_options.linear_solver_type = ceres::DENSE_SCHUR;
                std::cout<<"ceres::DENSE_SCHUR \n";
            }
            else if(num_images <= kMaxNUmImageDirectSparseSolver){
                solver_options.linear_solver_type = ceres::SPARSE_SCHUR;
                std::cout<<"ceres::SPARSE_SCHUR\n";
            }
            else{
                solver_options.linear_solver_type = ceres::ITERATIVE_SCHUR;
                solver_options.preconditioner_type = ceres::SCHUR_JACOBI;
                std::cout<<" ceres::ITERATIVE_SCHUR+SCHUR_JACOBI\n";
            }
            solver_options.linear_solver_type = ceres::SPARSE_SCHUR;
            solver_options.trust_region_strategy_type = ceres::DOGLEG;
            solver_options.preconditioner_type = ceres::SCHUR_JACOBI;
            std::cout<<"ceres::SPARSE_SCHUR + ceres::DOGLEG + ceres::SCHUR_JACOBI\n";
            solver_options.max_num_iterations = 10000;
            
            std::cout<<"solver_options.max_num_iterations:"<<solver_options.max_num_iterations<<"\n";
            std::cout<<"problem_->NumResiduals():"<<problem_->NumResiduals()<<"\n";
            if(problem_->NumResiduals()<
                        options_.min_num_residuals_for_multi_threading){
                        solver_options.num_threads = 1;
            #if CERES_VERSION_MAJOR <2
                solver_options.num_linear_solver_threads = 1;
            #endif
            }
            else{
                
                solver_options.num_threads = GetEffectiveNumThreads(solver_options.num_threads);
                std::cout<<solver_options.num_threads<<" threads!!\n";
            #if CERES_VERSION_MAJOR < 2
                solver_options.num_linear_solver_threads = 
                    GetEffectiveNumThreads(solver_options.num_linear_solver_threads);
            #endif
            }

            std::string solver_error;
            CHECK(solver_options.IsValid(&solver_error)) <<solver_error;

            ceres::Solve(solver_options,problem_.get(),&summary_);

            if(solver_options.minimizer_progress_to_stdout){
                std::cout<<std::endl;
            }

            if(options_.print_summary){
                PrintHeading2("Bundle adjustment report");
                PrintSolverSummary(summary_);
            }

            //TearDown(reconstruction);

            return true;

    }

    void BundleAdjustment::SetUpInter(std::vector<Reconstruction*>& reconstructions,
                                      std::unordered_map<Reconstruction*,SimilarityTransform>& srtform_to_global,
                                      ceres::LossFunction* loss_function)
    {
        
        CHECK_GE(reconstructions.size(),1);
        // Reconstruction* largest_reconstruction = reconstructions[0];
        // // add edge from local reconstructions to global reconstruction.
        // for(size_t i = 1 ;i < reconstructions.size(); i++)
        // {
        //     if(srtform_to_global.count(reconstructions[i])== 0) continue;
        //     auto similaritytransform = srtform_to_global.at(reconstructions[i]);
        //     AddOverlappingImageSingleSRT(reconstructions[i],largest_reconstruction,similaritytransform,loss_function);
        // }
        // add edge among local reconstructions.(one-way)

        /*for(size_t i = 0 ; i < reconstructions.size();i++)
        {
            std::cout<<"Reconstruction "<<i<<":\n"
                            <<"Qvec addresss:"<<srtform_to_global.at(reconstructions[i]).Qvec.data()<<"\n"
                            <<"Qvec:"<<srtform_to_global.at(reconstructions[i]).Qvec<<"\n"
                            <<"tvec:"<<srtform_to_global.at(reconstructions[i]).tvec<<"\n"
                            <<"scale:"<<srtform_to_global.at(reconstructions[i]).scale<<"\n";
        }*/
        
        
        for(size_t i = 0;i < reconstructions.size();i++)
        {
            if(srtform_to_global.count(reconstructions[i])== 0) continue;
            SimilarityTransform* similaritytransform1 = &srtform_to_global.at(reconstructions[i]);
            for(size_t j = 0; j< i ; j++)
            {
                // projection from reconstruction i to reconstruction j 
                if(srtform_to_global.count(reconstructions[j])== 0) continue;
                SimilarityTransform* similaritytransform2 = &srtform_to_global.at(reconstructions[j]);
                AddOverlappingImageDoubleSRT(reconstructions[i],reconstructions[j],
                                             similaritytransform1,similaritytransform2,loss_function);
            }
        }

        //Parameterize the srt
        problem_->SetParameterBlockConstant(&srtform_to_global.at(reconstructions[0]).scale);
        problem_->SetParameterBlockConstant(srtform_to_global.at(reconstructions[0]).Qvec.data());
        problem_->SetParameterBlockConstant(srtform_to_global.at(reconstructions[0]).tvec.data());

        
        
        for(size_t i = 1;i < reconstructions.size();i++)
        {
            if(srtform_to_global.count(reconstructions[i])== 0) continue;
            SimilarityTransform* similaritytransform = &srtform_to_global.at(reconstructions[i]);
            ceres::LocalParameterization* quaternion_parameterization = new ceres::QuaternionParameterization;   //??
            double* srt1_Qvec_data = similaritytransform->Qvec.data();
            
            /*std::cout<<"Reconstruction "<<reconstructions[i]<<"("<<similaritytransform->Qvec.data()
                                                                    <<":"<<srt1_Qvec_data[0]
                                                                    <<","<<srt1_Qvec_data[1]
                                                                    <<","<<srt1_Qvec_data[2]
                                                                    <<","<<srt1_Qvec_data[3]
                                                                    <<") is parameterizing.\n";*/
            problem_->SetParameterization(srt1_Qvec_data,quaternion_parameterization);
            //std::cout<<"Reconstruction "<<reconstructions[i]<<"("<<srt1_Qvec_data<<") is parameterized.\n";
        }
        
    }
    // add reprojection error from src_reconstruction's 3D point to ref_reconstrcution's 2D point by overlapping image
    void BundleAdjustment::AddOverlappingImageDoubleSRT(Reconstruction* src_reconstruction,Reconstruction* ref_reconstruction,
                                           SimilarityTransform* src_similaritytransform,SimilarityTransform* ref_similaritytransform,
                                           ceres::LossFunction* loss_function)
    {
        const auto& common_image_ids =
                        src_reconstruction->FindCommonRegImageIds(*ref_reconstruction);
            if(common_image_ids.size() == 0) {
                return;
            }

            
            CHECK_NOTNULL(src_similaritytransform);
            CHECK_NOTNULL(ref_similaritytransform);
            CHECK_NOTNULL(src_reconstruction);
            CHECK_NOTNULL(ref_reconstruction);
            
            // src_similaritytransform transforms src_reconstruction to global_reconstruction
            
            double* srt1_scale_data = &src_similaritytransform->scale;
            double* srt1_Qvec_data = src_similaritytransform->Qvec.data();
            double* srt1_tvec_data = src_similaritytransform->tvec.data();
            // ref_similaritytransform transforms ref_reconstruction to global_reconstruction
            double* srt2_scale_data = &ref_similaritytransform->scale;
            double* srt2_Qvec_data = ref_similaritytransform->Qvec.data();
            double* srt2_tvec_data = ref_similaritytransform->tvec.data();
            size_t num_observations = 0;
            for(const image_t image_id : common_image_ids)
            {
                Image& src_image = src_reconstruction->Image(image_id);
                Image& ref_image = ref_reconstruction->Image(image_id);
                Camera& ref_camera = ref_reconstruction->Camera(ref_image.CameraId());
                // Qvec_data,tvec_data transforms world coordinate to image plane in ref_reconstruction
                double* Qvec_data = ref_image.Qvec().data();
                double* tvec_data = ref_image.Tvec().data();
                double* camera_params = ref_camera.ParamsData();
                for(const Point2D& point2D:src_image.Points2D())
                {
                    if(point2D.HasPoint3D())
                    {
                        num_observations += 1;
                        Point3D& point3D = src_reconstruction->Point3D(point2D.Point3DId());
                        assert(point3D.Track().Length()>1);
                        ceres::CostFunction* cost_function = nullptr;

                        switch(ref_camera.ModelId())
                        {
                        #define  CAMERA_MODEL_CASE(CameraModel)   \
                            case CameraModel::kModelId:              \
                                cost_function =                    \
                                DoubleSrtBundleAdjustmentCostFunction<CameraModel>::Create( \
                                    point2D.XY()); \
                                break;
                        CAMERA_MODEL_SWITCH_CASES

                        #undef  CAMERA_MODEL_CASE
                        }
                        
                        problem_->AddResidualBlock(cost_function,loss_function,srt1_scale_data,srt1_Qvec_data,srt1_tvec_data,  
                                                   srt2_scale_data,srt2_Qvec_data,srt2_tvec_data,     
                                                   Qvec_data,tvec_data,point3D.XYZ().data(),camera_params);
                    }   
                }
                
                // if(num_observations > 0){
                //     if(!image_ids_.count(image_id)){
                //     //std::cout<<"Qvec_data:"<<Qvec_data[0]<<","<<Qvec_data[1]<<","<<Qvec_data[2]<<","<<Qvec_data[3]<<"\n";
                //     ceres::LocalParameterization* quaternion_parameterization = new ceres::QuaternionParameterization;
                //     std::cout<<"image_id("<<image_id<<") is parameterizing.\n";
                //     problem_->SetParameterization(Qvec_data,quaternion_parameterization);
                //     std::cout<<"image_id("<<image_id<<") is parameterized.\n\n";
                //     }
                //     else
                //     {
                //         std::cout<<"image_id("<<image_id<<") has been parameterized\n\n";
                //     }
                    
                // }
            }
            if(num_observations > 0) std::cout<<"Reconstruction "<<src_reconstruction
                                                                    <<"("<<src_similaritytransform->Qvec.data()
                                                                    <<":"<<srt1_Qvec_data[0]
                                                                    <<","<<srt1_Qvec_data[1]
                                                                    <<","<<srt1_Qvec_data[2]
                                                                    <<","<<srt1_Qvec_data[3]
                                                                    <<") is added.\n";
            
            //tvec??camera_params??parameterization??
    }

    // add reprojection error from src_reconstruction's 3D point to ref_reconstrcution's 2D point by overlapping image
    void BundleAdjustment::AddOverlappingImageSingleSRT(Reconstruction* src_reconstruction,Reconstruction* ref_reconstruction,
                                           SimilarityTransform* similaritytransform,ceres::LossFunction* loss_function)
    {
            const auto& common_image_ids =
                        src_reconstruction->FindCommonRegImageIds(*ref_reconstruction);
            if(common_image_ids.size() == 0) {
                return;
            }
            CHECK_NOTNULL(similaritytransform);
            CHECK_NOTNULL(src_reconstruction);
            CHECK_NOTNULL(ref_reconstruction);
            
            // similaritytransform transforms src_reconstruction to ref_reconstruction
            double* srt_scale_data = &similaritytransform->scale;
            double* srt_Qvec_data = similaritytransform->Qvec.data();
            double* srt_tvec_data = similaritytransform->tvec.data();

            for(const image_t image_id : common_image_ids)
            {
                Image& src_image = src_reconstruction->Image(image_id);
                Image& ref_image = ref_reconstruction->Image(image_id);
                Camera& ref_camera = ref_reconstruction->Camera(ref_image.CameraId());
                // Qvec_data,tvec_data transforms world coordinate to image plane in ref_reconstruction
                double* Qvec_data = ref_image.Qvec().data();
                double* tvec_data = ref_image.Tvec().data();
                double* camera_params = ref_camera.ParamsData();
                for(const Point2D& point2D:src_image.Points2D())
                {
                    if(point2D.HasPoint3D())
                    {
                        Point3D& point3D = src_reconstruction->Point3D(point2D.Point3DId());
                        assert(point3D.Track().Length()>1);
                        ceres::CostFunction* cost_function = nullptr;

                        switch(ref_camera.ModelId())
                        {
                        #define  CAMERA_MODEL_CASE(CameraModel)   \
                            case CameraModel::kModelId:              \
                                cost_function =                    \
                                SRTBundleAdjustmentCostFunction<CameraModel>::Create( \
                                    point2D.XY()); \
                                break;
                        CAMERA_MODEL_SWITCH_CASES

                        #undef  CAMERA_MODEL_CASE
                        }
                        
                        problem_->AddResidualBlock(cost_function,loss_function,srt_scale_data,srt_Qvec_data,srt_tvec_data,Qvec_data,tvec_data,point3D.XYZ().data(),camera_params);
                    }   
                }

            }
            
            //tvec??camera_params??parameterization??

    }


     void BundleAdjustment::SetUpIntra(Reconstruction* reconstruction,
                   ceres::LossFunction* loss_function)
    {
            camera_ids_.clear();
            point3D_num_observations_.clear();
            //PrintHeading2("begin setupintra ");
            BundleAdjustmentConfig* config_=configs_.at(reconstruction);
            for(const auto image_id:config_->Images()){
                AddImageToProblem(image_id,reconstruction,loss_function,config_);
            }
            //PrintHeading2("AddImageToProblem complete");
            for(const auto point3d_id:config_->VariablePoints()){
                AddPointToProblem(point3d_id,reconstruction,loss_function,config_);
            }
            //PrintHeading2("AddVariablePointToProblem complete");
            for(const auto point3d_id:config_->ConstantPoints()){
                AddPointToProblem(point3d_id,reconstruction,loss_function,config_);
            }
            //PrintHeading2("AddConstantPointToProblem complete");
            ParameterizeCameras(reconstruction,config_);
            //PrintHeading2("ParameterizeCameras complete");
            ParameterizePoints(reconstruction,config_);
            //PrintHeading2("ParameterizePoints complete");
    }

    void BundleAdjustment::AddImageToProblem(const image_t image_id,Reconstruction* reconstruction,
                               ceres::LossFunction* loss_function,BundleAdjustmentConfig* config_)
    {
        Image& image = reconstruction->Image(image_id);
        Camera& camera = reconstruction->Camera(image.CameraId());
        

        image.NormalizeQvec();

        double* Qvec_data = image.Qvec().data();
        double* tvec_data = image.Tvec().data();
        double* camera_params = camera.ParamsData();
        

        size_t num_observation = 0;

        bool constant_pose = 
           !options_.refine_extrinsics || config_->HasConstantPose(image_id);
        
        
        for(const auto& point2D:image.Points2D())
        {
            if(!point2D.HasPoint3D()){
                continue;
            }  
            
            num_observation += 1;
            point3D_num_observations_[point2D.Point3DId()] += 1;

            Point3D& point3D = reconstruction->Point3D(point2D.Point3DId());

            assert(point3D.Track().Length()>1);
            ceres::CostFunction* cost_function = nullptr;
            if(constant_pose){
                switch(camera.ModelId())
                {
                #define  CAMERA_MODEL_CASE(CameraModel)   \
                    case CameraModel::kModelId:             \
                        cost_function =                    \
                        BundleAdjustmentConstantPoseCostFunction<CameraModel>::Create( \
                            image.Qvec(),image.Tvec(),point2D.XY()); \
                        break;
                CAMERA_MODEL_SWITCH_CASES

                #undef  CAMERA_MODEL_CASE

                }
                problem_->AddResidualBlock(cost_function,loss_function,point3D.XYZ().data(),camera_params);
            }
            else
            {
                switch(camera.ModelId())
                {
                #define  CAMERA_MODEL_CASE(CameraModel)   \
                    case CameraModel::kModelId:              \
                        cost_function =                    \
                        BundleAdjustmentCostFunction<CameraModel>::Create( \
                            point2D.XY()); \
                        break;
                CAMERA_MODEL_SWITCH_CASES

                #undef  CAMERA_MODEL_CASE
                }
                problem_->AddResidualBlock(cost_function,loss_function,Qvec_data,
                                           tvec_data,point3D.XYZ().data(),camera_params);
            }

            
        }

        if(num_observation)
        {
            camera_ids_.insert(image.CameraId());
            image_ids_.insert(image_id);
            
            if(!constant_pose)
            {
                ceres::LocalParameterization* qvec_quaternionparameterization =
                    new ceres::QuaternionParameterization;
                problem_->SetParameterization(Qvec_data,qvec_quaternionparameterization);

                if(config_->HasConstantTvec(image_id))
                {
                    const auto& constant_tvec_idxs = config_->ConstantTvec(image_id);
                    ceres::SubsetParameterization* tvec_subsetparameterization =
                    new ceres::SubsetParameterization(3,constant_tvec_idxs);
                    problem_->SetParameterization(tvec_data,tvec_subsetparameterization);
                }
            }
        }
    }


    void BundleAdjustment::AddPointToProblem(const point3D_t point3D_id,
                               Reconstruction* reconstruction,
                               ceres::LossFunction* loss_function,BundleAdjustmentConfig* config_)
    {
        Point3D& point3D = reconstruction->Point3D(point3D_id);
        

        if(point3D_num_observations_[point3D_id] == point3D.Track().Length())
        {
            return;
        }
        for(const auto& track_item:point3D.Track().Elements())
        {
            image_t image_id = track_item.image_id;
            point2D_t point2D_idx = track_item.point2D_idx;
            if(config_->HasImage(image_id))
            {
                continue;
            }

            point3D_num_observations_[point3D_id] += 1;
            Image& image = reconstruction->Image(image_id);
            Camera& camera = reconstruction->Camera(image.CameraId());
            Point2D& point2D = image.Point2D(point2D_idx);
            double* camera_params = camera.ParamsData();

            if(camera_ids_.count(image.CameraId()) == 0)
            {
                camera_ids_.insert(image.CameraId());
                config_->SetConstantCamera(image.CameraId());
            }
            std::cout<<"AddPointToProblem():image_id("<<image_id<<")"<<",point2D_idx("<<point2D_idx<<"):point3D_id("<<point3D_id
                                                                    <<")"<<" Reconstruction("<<reconstruction<<") add constantly\n";
            ceres::CostFunction* cost_function = nullptr;
            switch(camera.ModelId())
                {
                #define  CAMERA_MODEL_CASE(CameraModel)   \
                    case CameraModel::kModelId:             \
                        cost_function =                    \
                        BundleAdjustmentConstantPoseCostFunction<CameraModel>::Create( \
                            image.Qvec(),image.Tvec(),point2D.XY()); \
                        break;
                CAMERA_MODEL_SWITCH_CASES

                #undef  CAMERA_MODEL_CASE

                }
            problem_->AddResidualBlock(cost_function,loss_function,point3D.XYZ().data(),camera_params);

            
        }
    }


    void BundleAdjustment::ParameterizeCameras(Reconstruction* reconstruction,BundleAdjustmentConfig* config_)
    {
        const bool constant_camera = !options_.refine_focal_length &&
                                     !options_.refine_principal_point &&
                                     !options_.refine_extra_params;
        
        for(auto camera_id:camera_ids_)
        {
            Camera& camera = reconstruction->Camera(camera_id);
            if(constant_camera ||config_->IsConstantCamera(camera_id))
            {
                std::cout<<"ParameterizeCameras():camera_id:"<<camera_id<<" is constant\n";
                problem_->SetParameterBlockConstant(camera.ParamsData());
            }
            else
            {
                std::cout<<"ParameterizeCameras():camera_id:"<<camera_id<<")"<<" Reconstruction("<<reconstruction<<" is variable\n";
                std::vector<int> constant_camera_params;

                if(!options_.refine_focal_length)
                {
                    const std::vector<size_t>& params_idxs = camera.FocalLengthIdxs();
                    constant_camera_params.insert(constant_camera_params.end(),params_idxs.begin(),
                                                  params_idxs.end());
                }

                if(!options_.refine_principal_point)
                {
                    const std::vector<size_t>& params_idxs = camera.PrincipalPointIdxs();
                    constant_camera_params.insert(constant_camera_params.end(),params_idxs.begin(),
                                                  params_idxs.end());
                }

                if(!options_.refine_extra_params)
                {
                    const std::vector<size_t>& params_idxs = camera.ExtraParamsIdxs();
                    constant_camera_params.insert(constant_camera_params.end(),params_idxs.begin(),
                                                  params_idxs.end());
                }

                if(constant_camera_params.size() > 0){
                    ceres::SubsetParameterization* camera_params_subsetparameterization =
                    new ceres::SubsetParameterization(static_cast<int>(camera.NumParams()),constant_camera_params);
                    problem_->SetParameterization(camera.ParamsData(),camera_params_subsetparameterization);
                }
            }
        }

    }


    void BundleAdjustment::ParameterizePoints(Reconstruction* reconstruction,BundleAdjustmentConfig* config_)
    {
        
        for(const auto elem:point3D_num_observations_)
        {
            
            Point3D& point3D = reconstruction->Point3D(elem.first);
            if(elem.second < point3D.Track().Length())
            {
                std::cout<<"ParameterizePoints():point3d("<<elem.first<<")"<<" Reconstruction("<<reconstruction<<")is constant\n";
                problem_->SetParameterBlockConstant(point3D.XYZ().data());
            }
            
        }

        for(const point3D_t point3D_idx:config_->ConstantPoints())
        {
            Point3D& point3D = reconstruction->Point3D(point3D_idx);
            
            problem_->SetParameterBlockConstant(point3D.XYZ().data());
            
        }

        
    }

    bool BundleAdjustment::SolveOGBA(Reconstruction* reconstruction)
    {
            std::cout<<"/////////DACSFM::OGBA/////////\n";
            CHECK(!problem_)<<"Cannot use the same BundleAdjuster multiple times";

            problem_.reset(new ceres::Problem());

            //ceres::LossFunction* loss_function = options_.CreateLossFunction();
            ceres::LossFunction* loss_function = new ceres::HuberLoss(2.0);
           
            
            CHECK_NOTNULL(reconstruction);
            SetUpIntra(reconstruction,loss_function);
            
            
            if(problem_->NumResiduals() == 0){
                return false;
            }

            //ceres::Solver::Options solver_options = options_.solver_options;
            ceres::Solver::Options solver_options ;
            solver_options.minimizer_progress_to_stdout = true;
            //solver_options.max_linear_solver_iterations = 200;
            //solver_options.max_num_consecutive_invalid_steps = 10;
            solver_options.num_threads = -1;
            
            

            const size_t kMaxNumImageDirectDenseSolver = 50;
            const size_t kMaxNUmImageDirectSparseSolver = 1000;
            const size_t num_images_all = image_ids_.size();
            size_t num_images = 0;
            for(auto config:configs_)
            {
                num_images += config.second->NumImages();
            }
            std::cout<<"num_images(overlap):"<<num_images<<"\n";
            std::cout<<"num_images_all:"<<num_images_all<<"\n";
            if(num_images <= kMaxNumImageDirectDenseSolver){
                solver_options.linear_solver_type = ceres::DENSE_SCHUR;
                std::cout<<"ceres::DENSE_SCHUR \n";
            }
            else if(num_images <= kMaxNUmImageDirectSparseSolver){
                solver_options.linear_solver_type = ceres::SPARSE_SCHUR;
                std::cout<<"ceres::SPARSE_SCHUR\n";
            }
            else{
                solver_options.linear_solver_type = ceres::SPARSE_SCHUR;
                std::cout<<"ceres::SPARSE_SCHUR\n";
                // solver_options.linear_solver_type = ceres::ITERATIVE_SCHUR;
                // solver_options.preconditioner_type = ceres::SCHUR_JACOBI;
                // std::cout<<" ceres::ITERATIVE_SCHUR+SCHUR_JACOBI\n";
            }
             solver_options.linear_solver_type = ceres::SPARSE_SCHUR;
             solver_options.trust_region_strategy_type = ceres::DOGLEG;
             solver_options.preconditioner_type = ceres::SCHUR_JACOBI;
             std::cout<<"ceres::SPARSE_SCHUR + ceres::DOGLEG + ceres::SCHUR_JACOBI\n";
            solver_options.max_num_iterations = 1000;
            
            std::cout<<"solver_options.max_num_iterations:"<<solver_options.max_num_iterations<<"\n";
            std::cout<<"problem_->NumResiduals():"<<problem_->NumResiduals()<<"\n";
            if(problem_->NumResiduals()<
                        options_.min_num_residuals_for_multi_threading){
                        solver_options.num_threads = 1;
            #if CERES_VERSION_MAJOR <2
                solver_options.num_linear_solver_threads = 1;
            #endif
            }
            else{
                
                solver_options.num_threads = GetEffectiveNumThreads(solver_options.num_threads);
                std::cout<<solver_options.num_threads<<" threads!!\n";
            #if CERES_VERSION_MAJOR < 2
                solver_options.num_linear_solver_threads = 
                    GetEffectiveNumThreads(solver_options.num_linear_solver_threads);
            #endif
            }

            std::string solver_error;
            CHECK(solver_options.IsValid(&solver_error)) <<solver_error;

            ceres::Solve(solver_options,problem_.get(),&summary_);

            if(solver_options.minimizer_progress_to_stdout){
                std::cout<<std::endl;
            }

            if(options_.print_summary){
                PrintHeading2("Bundle adjustment report");
                PrintSolverSummary(summary_);
            }

            //TearDown(reconstruction);

            return true;

    }








}
