#include "merge/mergeclusters.h"

#include "util/option_manager.h"
#include "merge/bundleadjustment_controller.h"
#include "base/similarity_transform.h"
#include "base/projection.h"
#include <ceres/ceres.h>
#include "controllers/bundle_adjustment.h"
#include "util/misc.h"
#include <iomanip>

namespace DACSFM {

    // class BundleAdjustmentIterationCallback::public ceres::IterationCallback{
    //     public:
    //     explicit BundleAdjustmentIterationCallback(Thread* thread)
    //     :thread_(thread){}

    //     virtual ceres::CallbackReturnType operator()(const IterationSummary& summary){
    //         CHECK_NOTNULL(thread_);
    //         thread->BlockIfPaused();
    //         if(thread_->IsStopped()){
    //             return ceres::SOLVER_TERMINATE_SUCCESSFULLY;
    //         }else{
    //             return ceres::SOLVER_CONTINUE;
    //         }
    //     }


    //     private:
    //     Thread* thread_;
    // };

    MergeController::MergeController(){}

    MergeController::MergeController(MergeControllerOptions mergeoptions):
        mergeoptions_(mergeoptions){}

    void MergeController::MergeClustersOnTreeRun(
    const SceneClustering::Cluster& cluster,
    std::unordered_map<const SceneClustering::Cluster*, ReconstructionManager>*
        reconstruction_managers)
    {
         OptionManager options;
         //options.AddRequiredOption("input_path", &input_path);
         //options.AddRequiredOption("output_path", &output_path);
         options.AddBundleAdjustmentOptions();
         
         //MergeClustersOnTree(cluster,reconstruction_managers,options);
    }


    void MergeController::MergeClustersOnTree(
    const SceneClustering::Cluster& cluster,
    std::unordered_map<const SceneClustering::Cluster*, ReconstructionManager>*
        reconstruction_managers,const OptionManager& options) {
    // Extract all reconstructions from all child clusters.
    std::vector<Reconstruction*> reconstructions;
    for (const auto& child_cluster : cluster.child_clusters) {
        if (!child_cluster.child_clusters.empty()) {
        MergeClustersOnTree(child_cluster, reconstruction_managers,options);
        }

        auto& reconstruction_manager = reconstruction_managers->at(&child_cluster);
        for (size_t i = 0; i < reconstruction_manager.Size(); ++i) {
        reconstructions.push_back(&reconstruction_manager.Get(i));
        }
    }

    // Try to merge all child cluster reconstruction.
    while (reconstructions.size() > 1) {
        bool merge_success = false;
        for (size_t i = 0; i < reconstructions.size(); ++i) {
        for (size_t j = 0; j < i; ++j) {
            const double kMaxReprojError = 8.0;
            if (reconstructions[i]->Merge(*reconstructions[j], kMaxReprojError)) {
            reconstructions.erase(reconstructions.begin() + j);
            merge_success = true;
            /*****************************/
            // DACSFM::MyBundleAdjustmentController ba_controller(options, reconstructions[i]);
            // ba_controller.Start();
            // ba_controller.Wait();
            /****************************/ 

            break;
            }
        }

        if (merge_success) {
            break;
        }
        }

        if (!merge_success) {
        break;
        }
    }

    // Insert a new reconstruction manager for merged cluster.
    auto& reconstruction_manager = (*reconstruction_managers)[&cluster];
    for (const auto& reconstruction : reconstructions) {
        reconstruction_manager.Add();
        reconstruction_manager.Get(reconstruction_manager.Size() - 1) =
            *reconstruction;
    }

    // Delete all merged child cluster reconstruction managers.
    for (const auto& child_cluster : cluster.child_clusters) {
        reconstruction_managers->erase(&child_cluster);
    }
    }


    void MergeController::MergeClustersIntoLargestReconstructionRun(
    std::vector<Reconstruction*>& reconstructions,
    ReconstructionManager* reconstruction_manager,
    std::vector<std::string>& cluster_name)
    {
        PrintHeading2("Begin Combination");
        OptionManager options;
        options.AddBundleAdjustmentOptions();
         //options.AddRequiredOption("input_path", &input_path);
         //options.AddRequiredOption("output_path", &output_path);
         
        MergeClustersIntoLargestReconstruction(reconstructions,cluster_name);

        /////Debug
        if(reconstructions.size() != 1)
        {
            std::string debug_path = JoinPaths(mergeoptions_.output_path,"Debug");
            CreateDirIfNotExists(debug_path);
            for (size_t i = 0; i < reconstructions.size(); ++i) {
                const std::string reconstruction_path = JoinPaths(debug_path, 
                                                                    "partition_" + std::to_string(i));
                CreateDirIfNotExists(reconstruction_path);
                reconstructions[i]->Write(reconstruction_path);
                reconstructions[i]->WriteText(reconstruction_path);
            }
        }
        /////
        CHECK_EQ(reconstructions.size(),1);

        PrintHeading2("CombineReconstructions:adjust Global BA");
        //////////////////////////////////////////////////
        //            adjust Global BA                  //
        ////////////////////////////////////////////////// 
        // DACSFM::MyBundleAdjustmentController ba_controller(options, reconstructions[0]);
        // ba_controller.Start();
        // ba_controller.Wait();

        reconstruction_manager->Add();
        reconstruction_manager->Get(reconstruction_manager->Size() - 1) =*reconstructions[0]; 
        
         
    }

    void MergeController::MergeClustersIntoLargestReconstruction(std::vector<Reconstruction*>& reconstructions,
                                std::vector<std::string>& cluster_name)
    {
        const auto cmp = [](const Reconstruction* reconstruction1,const Reconstruction* reconstruction2){
            if(reconstruction1->NumRegImages() == reconstruction2->NumRegImages())
                return reconstruction1->NumPoints3D() > reconstruction2->NumPoints3D();
            return reconstruction1->NumRegImages()>reconstruction2->NumRegImages();  //??NumRegImages??NumPoints3D??
        };

        std::sort(reconstructions.begin(),reconstructions.end(),cmp);
        //merge other reconstructions into the largest reconstrcution
        CHECK_GE(reconstructions.size(),1);
        PrintHeading2("reconstructions sorted");


        const double kMaxReprojError = 8.0;
        Reconstruction* largest_reconstruction = reconstructions[0];
        int index = 0;
        while(reconstructions.size() > 1)
        {
            /*bool merge_success = false;
            auto it = reconstructions.begin();
            it++;
            
            printf("/////////\nThe %02d loops\n/////////",++index);
            for(size_t i=0;i<reconstructions.size();i++) 
                std::cout<<"Reconstruction "<<i<<".NumRegImages:"<<reconstructions[i]->NumRegImages()<<"\n"<<"      NumPoints3D:"<<reconstructions[i]->NumPoints3D()<<"\n";
            while(it!=reconstructions.end())
            {
                if (largest_reconstruction->Merge(*(*it), kMaxReprojError)) {
                    std::cout<<"merge_success!\n";
                    reconstructions.erase(it);
                    merge_success = true;
                }else{
                std::cout<<"merge_fail!\n";
                it++;
                }
            }

            if(!merge_success)
            {
                break;
            }*/

            bool merge_success = false;
            printf("/////////\nThe %02d loops\n/////////\n",++index);
            for(size_t i=0;i<reconstructions.size();i++) 
                std::cout<<"Reconstruction "<<i<<"("<<cluster_name[i]<<")"<<".NumRegImages:"<<reconstructions[i]->NumRegImages()<<"\n"<<"      NumPoints3D:"<<reconstructions[i]->NumPoints3D()<<"\n";
            for(size_t i = 1;i<reconstructions.size();i++)
            {
                if (largest_reconstruction->Merge(*reconstructions[i], kMaxReprojError)) {
                    std::cout<<"(0,"<<i<<")"<<"merge_success!\n";
                    std::cout<<"("<<cluster_name[0]<<","<<cluster_name[i]<<") merge_success!\n";
                    reconstructions.erase(reconstructions.begin()+i);
                    cluster_name.erase(cluster_name.begin()+i);
                    merge_success = true;
                    break;
                }
            }
            if(!merge_success)
            {
                std::cout<<"/////////adjust for merge.\n";
                for(size_t i = 1;i<reconstructions.size();i++)
                {
                    for(size_t j = i+1;j<reconstructions.size();j++)
                    {
                        if(reconstructions[i]->Merge(*reconstructions[j],kMaxReprojError)){
                            std::cout<<"("<<i<<","<<j<<")"<<"merge_success!\n";
                            std::cout<<"("<<cluster_name[i]<<","<<cluster_name[j]<<") merge_success!\n";
                            reconstructions.erase(reconstructions.begin()+j);
                            cluster_name.erase(cluster_name.begin()+j);
                            merge_success = true;
                            break;
                        }
                    }
                    if(merge_success) break;
                }
                if(!merge_success) break;

            }



        }

        


    }


    void MergeController::MergeClusterwithBARun(
    std::vector<Reconstruction*>& reconstructions,
    ReconstructionManager* reconstruction_manager,
    std::vector<std::string>& cluster_name)
    {
        
        PrintHeading2("Begin Combination");
        
        
        const auto cmp = [](const Reconstruction* reconstruction1,const Reconstruction* reconstruction2){
            if(reconstruction1->NumRegImages() == reconstruction2->NumRegImages())
                return reconstruction1->NumPoints3D() > reconstruction2->NumPoints3D();
            return reconstruction1->NumRegImages()>reconstruction2->NumRegImages();  //??NumRegImages??NumPoints3D??
        };


        std::sort(reconstructions.begin(),reconstructions.end(),cmp);
        PrintHeading2("BundleAdjustment for optimal SimilarityTransform");
        // invoke BA to obtain the optimal SimilarityTransform.
        const double kMaxReprojError = 8.0;
        std::unordered_map<Reconstruction*,const SimilarityTransform3> similaritytransform3_to_global;
        CallBundleAdjustment(reconstructions,kMaxReprojError,similaritytransform3_to_global);



        //merge other reconstructions into the largest reconstrcution
        CHECK_GE(reconstructions.size(),1);
        PrintHeading2("reconstructions sorted");


        //const double kMaxReprojError = 8.0;
        Reconstruction* global_reconstruction = reconstructions[0];
        int index = 0;

        
        while(reconstructions.size() > 1)
        {
            bool merge_success = false;
            printf("/////////\nThe %02d loops\n/////////\n",++index);
            for(size_t i=0;i<reconstructions.size();i++) 
                std::cout<<"Reconstruction "<<i<<"("<<cluster_name[i]<<")"<<".NumRegImages:"<<reconstructions[i]->NumRegImages()<<"\n"<<"      NumPoints3D:"<<reconstructions[i]->NumPoints3D()<<"\n";
            for(size_t i = 1;i<reconstructions.size();i++)
            {
                if(!similaritytransform3_to_global.count(reconstructions[i])){
                    std::cerr<<"Reconstruction "<<i<<"("<<cluster_name[i]<<") miss similaritytransform to global.\n";
                    continue;
                }
                if (MergeReconstructionsAmbiguously(*global_reconstruction,*reconstructions[i], kMaxReprojError,
                          similaritytransform3_to_global.at(reconstructions[i]))) {
                    std::cout<<"(0,"<<i<<")"<<"merge_success!\n";
                    std::cout<<"("<<cluster_name[0]<<","<<cluster_name[i]<<") merge_success!\n";
                    reconstructions.erase(reconstructions.begin()+i);
                    cluster_name.erase(cluster_name.begin()+i);
                    merge_success = true;
                    break;
                }
            }
            if(!merge_success) break;
        }

        

        /////Debug
        if(reconstructions.size() != 1)
        {
            std::string debug_path = JoinPaths(mergeoptions_.output_path,"Debug");
            CreateDirIfNotExists(debug_path);
            for (size_t i = 0; i < reconstructions.size(); ++i) {
                const std::string reconstruction_path = JoinPaths(debug_path, 
                                                                    "partition_" + std::to_string(i));
                CreateDirIfNotExists(reconstruction_path);
                reconstructions[i]->Write(reconstruction_path);
                reconstructions[i]->WriteText(reconstruction_path);
            }
        }
        /////
        CHECK_EQ(reconstructions.size(),1);
        

        reconstruction_manager->Add();
        reconstruction_manager->Get(reconstruction_manager->Size() - 1) =*reconstructions[0]; 
        
         
    }


    void MergeController::CallBundleAdjustment(std::vector<Reconstruction*>& reconstructions,const double kMaxReprojError,
                              std::unordered_map<Reconstruction*,const SimilarityTransform3>& similaritytransform3_to_global)
    {

        OptionManager options;
        options.AddBundleAdjustmentOptions();
         //options.AddRequiredOption("input_path", &input_path);
         //options.AddRequiredOption("output_path", &output_path);
        
        //global_reconstruction is in the first (default)
        Reconstruction* global_reconstruction = reconstructions[0];
        CHECK_NOTNULL(global_reconstruction);
        const double kMinInlierObservations = 0.3;
        std::unordered_map<Reconstruction*,SimilarityTransform*> srtform_to_global;
        PrintHeading2("Estimate initial SimilarityTransform");
        // obtain the initial SimilarityTransform.
        for(size_t i = 1; i <reconstructions.size();i++)
        {
                Eigen::Matrix3x4d alignment;
                CHECK_NOTNULL(reconstructions[i]);
                //compute global to reconstructions[i]
                if(ComputeAlignmentBetweenReconstructions(*reconstructions[i],*global_reconstruction,
                                                       kMinInlierObservations,kMaxReprojError,&alignment))
                {
                    std::cout<<"Reconstruction "<<i<<"'s srt estimation succeeded.\n";
                    if(!srtform_to_global.count(reconstructions[i]))
                    {
                        SimilarityTransform* srt = new SimilarityTransform(alignment);
                        srtform_to_global.insert(std::make_pair(reconstructions[i] ,srt));
                    }
                    else
                    {
                        std::cerr<<"Duplicate reconstruction in a array of reconstructions!\n";
                    }
                    
                }
                else
                {
                    std::cout<<"Reconstruction "<<i<<"'s srt estimation failed.\n";
                }
                std::cout<<"\n";
        }
        PrintHeading2("Begin Bundle Adjustment");
        //DACSFM::MyBundleAdjustmentController ba_controller(options, reconstructions,srtform_to_global);
        //ba_controller.Start();
        //ba_controller.Wait();
        
        //std::unordered_map<Reconstruction*,const SimilarityTransform3> similaritytransform3_to_global;
        for (auto srtform_item: srtform_to_global)
        {
            similaritytransform3_to_global.insert(std::make_pair(srtform_item.first,
            SimilarityTransform3(srtform_item.second->scale,srtform_item.second->Qvec,srtform_item.second->tvec)));
        }
        

    }

    // merge sub_reconstruction into main_reconstruction
    bool MergeController::MergeReconstructionsAmbiguously(Reconstruction& main_reconstruction,const Reconstruction& sub_reconstruction,
                           const double max_reproj_error,const SimilarityTransform3& tform) {

        // Find common and missing images in the two reconstructions.

        std::unordered_set<image_t> common_image_ids;
        common_image_ids.reserve(sub_reconstruction.NumRegImages());
        std::unordered_set<image_t> missing_image_ids;
        missing_image_ids.reserve(sub_reconstruction.NumRegImages());

        for (const auto& image_id : sub_reconstruction.RegImageIds()) {
            if (main_reconstruction.ExistsImage(image_id)) {   //main_reconstruction
            common_image_ids.insert(image_id);
            } else {
            missing_image_ids.insert(image_id);
            }
        }

        // Register the missing images in this sub_reconstruction.

        for (const auto image_id : missing_image_ids) {
            auto reg_image = sub_reconstruction.Image(image_id);
            reg_image.SetRegistered(false);
            //main_reconstruction
            main_reconstruction.AddImage(reg_image);
            main_reconstruction.RegisterImage(image_id);

            if (!main_reconstruction.ExistsCamera(reg_image.CameraId())) { //main_reconstruction
            main_reconstruction.AddCamera(sub_reconstruction.Camera(reg_image.CameraId()));//main_reconstruction
            }
            auto& image = main_reconstruction.Image(image_id); //main_reconstruction
            tform.TransformPose(&image.Qvec(), &image.Tvec());
        }

        // Merge the two point clouds using the following two rules:
        //    - copy points to this sub_reconstruction with non-conflicting tracks,
        //      i.e. points that do not have an already triangulated observation
        //      in this sub_reconstruction.
        //    - merge tracks that are unambiguous, i.e. only merge points in the two
        //      sub_reconstructions if they have a one-to-one mapping.
        // Note that in both cases no cheirality or reprojection test is performed.

        for (const auto& point3D : sub_reconstruction.Points3D()) {
            Track new_track;
            Track old_track;
            std::set<point3D_t> old_point3D_ids;
            for (const auto& track_el : point3D.second.Track().Elements()) {
            if (common_image_ids.count(track_el.image_id) > 0) {
                const auto& point2D =
                    main_reconstruction.Image(track_el.image_id).Point2D(track_el.point2D_idx);
                if (point2D.HasPoint3D()) {
                old_track.AddElement(track_el);
                old_point3D_ids.insert(point2D.Point3DId());
                } else {
                new_track.AddElement(track_el);
                }
            } else if (missing_image_ids.count(track_el.image_id) > 0) {
                main_reconstruction.Image(track_el.image_id).ResetPoint3DForPoint2D(track_el.point2D_idx);
                new_track.AddElement(track_el);
            }
            }

            const bool create_new_point = new_track.Length() >= 2;
            const bool merge_new_and_old_point =
                (new_track.Length() + old_track.Length()) >= 2 &&
                old_point3D_ids.size() == 1;
            if (create_new_point || merge_new_and_old_point) {
            Eigen::Vector3d xyz = point3D.second.XYZ();
            tform.TransformPoint(&xyz);
            const auto point3D_id =
                main_reconstruction.AddPoint3D(xyz, new_track, point3D.second.Color());//main_reconstruction
            if (old_point3D_ids.size() == 1) {
                main_reconstruction.MergePoints3D(point3D_id, *old_point3D_ids.begin());//main_reconstruction
            }
            }
        }
        //main_reconstruction
        FilterPoints3DWithLargeReprojectionError(main_reconstruction,max_reproj_error, main_reconstruction.Point3DIds());

        return true;
    }

    size_t MergeController::FilterPoints3DWithLargeReprojectionError(
    Reconstruction& reconstruction,
    const double max_reproj_error,
    const std::unordered_set<point3D_t>& point3D_ids) {
        const double max_squared_reproj_error = max_reproj_error * max_reproj_error;

        // Number of filtered points.
        size_t num_filtered = 0;

        for (const auto point3D_id : point3D_ids) {
            if (!reconstruction.ExistsPoint3D(point3D_id)) {
            continue;
            }

            class Point3D& point3D = reconstruction.Point3D(point3D_id);

            if (point3D.Track().Length() < 2) {
            reconstruction.DeletePoint3D(point3D_id);
            num_filtered += point3D.Track().Length();
            continue;
            }

            double reproj_error_sum = 0.0;

            std::vector<TrackElement> track_els_to_delete;

            for (const auto& track_el : point3D.Track().Elements()) {
            const class Image& image = reconstruction.Image(track_el.image_id);
            const class Camera& camera = reconstruction.Camera(image.CameraId());
            const Point2D& point2D = image.Point2D(track_el.point2D_idx);
            const double squared_reproj_error = CalculateSquaredReprojectionError(
                point2D.XY(), point3D.XYZ(), image.Qvec(), image.Tvec(), camera);
            if (squared_reproj_error > max_squared_reproj_error) {
                track_els_to_delete.push_back(track_el);
            } else {
                reproj_error_sum += std::sqrt(squared_reproj_error);
            }
            }

            if (track_els_to_delete.size() >= point3D.Track().Length() - 1) {
            num_filtered += point3D.Track().Length();
            reconstruction.DeletePoint3D(point3D_id);
            } else {
            num_filtered += track_els_to_delete.size();
            for (const auto& track_el : track_els_to_delete) {
                reconstruction.DeleteObservation(track_el.image_id, track_el.point2D_idx);
            }
            point3D.SetError(reproj_error_sum / point3D.Track().Length());
            }
        }

        return num_filtered;
    }
    
    


    //////////////////////////////////////////////////////
    //          Virtual Merge(SRT),BA,Merge
    //////////////////////////////////////////////////////

    void MergeController::MergeClustersintoMaximum(
    std::vector<Reconstruction*>& reconstructions,
    ReconstructionManager* reconstruction_manager,
    std::unordered_map<Reconstruction*,std::string>& cluster_name)
    {
        
        PrintHeading2("Begin Combination");
        

        const auto cmp = [](const Reconstruction* reconstruction1,const Reconstruction* reconstruction2){
            if(reconstruction1->NumRegImages() == reconstruction2->NumRegImages())
                return reconstruction1->NumPoints3D() > reconstruction2->NumPoints3D();
            return reconstruction1->NumRegImages()>reconstruction2->NumRegImages();  //??NumRegImages??NumPoints3D??
        };


        std::sort(reconstructions.begin(),reconstructions.end(),cmp);
        PrintHeading2("reconstructions sorted");
        
        
        for(size_t i=0;i<reconstructions.size();i++) 
        {
                std::cout<<"Reconstruction "<<i<<"("<<cluster_name.at(reconstructions[i])<<")"
                         <<".NumRegImages:"<<reconstructions[i]->NumRegImages()<<"\n"
                         <<"      NumPoints3D:"<<reconstructions[i]->NumPoints3D()<<"\n";
        }
        //////////////////////////////////////////////////////////////////////////
        //estimate  srt from all sub-reconstructions to the global reconstrcution
        //////////////////////////////////////////////////////////////////////////
        PrintHeading2("Virtual Merge");
        std::cout<<"PrintHeading2('Virtual Merge');"<<"\n";
        CHECK_GE(reconstructions.size(),1);
        std::cout<<"CHECK_GE(reconstructions.size(),1);"<<"\n";
        const double kMaxReprojError = 8.0;
        
        //Intialize the global reconstuction with largest reconstruction
        std::vector<bool> is_merged(reconstructions.size(),false);
        std::cout<<""<<"\n";
        
        Reconstruction global_reconstruction = *reconstructions[0];
        std::unordered_map<Reconstruction*,SimilarityTransform> srtform_to_maximum;
        std::cout<<"std::unordered_map<Reconstruction*,SimilarityTransform> srtform_to_maximum;"<<"\n";
        
        srtform_to_maximum.insert(std::make_pair(reconstructions[0],SimilarityTransform()));
        is_merged[0] = true;
        int index = 0;
        const int tot_recons = reconstructions.size();
        std::cout<<"const int tot_recons = reconstructions.size();"<<"\n";
        while(index < tot_recons)
        {
            bool merge_success = false;
            printf("/////////\nThe %02d loops\n/////////\n",++index);
            std::cout<<"global_reconstruction "<<".NumRegImages:"<<global_reconstruction.NumRegImages()<<"\n"<<"      NumPoints3D:"<<global_reconstruction.NumPoints3D()<<"\n";
            for(size_t i=0;i<reconstructions.size();i++) 
                if(!is_merged[i])std::cout<<"Reconstruction "<<i<<"("<<cluster_name.at(reconstructions[i])<<")"<<".NumRegImages:"<<reconstructions[i]->NumRegImages()<<"\n"<<"      NumPoints3D:"<<reconstructions[i]->NumPoints3D()<<"\n";
            for(size_t i = 1;i<reconstructions.size();i++)
            {
                if(is_merged[i]) continue;
                SimilarityTransform srt;
                if (MergeAllForSRT(global_reconstruction,*reconstructions[i], kMaxReprojError,srt)) {

                    CHECK_LE(srtform_to_maximum.count(reconstructions[i]),0);
                    srtform_to_maximum.insert(std::make_pair(reconstructions[i],srt));

                    std::cout<<"(0,"<<i<<")"<<"estimates srt successfully!\n";
                    std::cout<<"("<<cluster_name.at(reconstructions[0])<<","<<cluster_name.at(reconstructions[i])<<") estimates srt successfully!\n";
                    is_merged[i] = true;
                    merge_success = true;
                    break;
                }
            }
            if(!merge_success) break;
        }
        for(size_t i = 0 ; i < reconstructions.size();i++)
        {
            if(is_merged[i]==false)
            {
                std::cout<<"//////////////\n        Merge Fail\n//////////////\n";
                return;
            }
            std::cout<<"Recoonstruction "<<i<<"("<<cluster_name.at(reconstructions[i])<<"):\n"
                            <<"Qvec addresss:"<<srtform_to_maximum.at(reconstructions[i]).Qvec.data()<<"\n"
                            <<"Qvec:"<<srtform_to_maximum.at(reconstructions[i]).Qvec<<"\n"
                            <<"tvec:"<<srtform_to_maximum.at(reconstructions[i]).tvec<<"\n"
                            <<"scale:"<<srtform_to_maximum.at(reconstructions[i]).scale<<"\n";
        }

        //////////////////////////////////////////////////////
        ////////////////estimate the srt from local to global
        
        std::cout<<"////////////////estimate the srt from local to global\n";
        std::unordered_map<Reconstruction*,SimilarityTransform> srtform_to_global;
        for(size_t i = 0 ; i < reconstructions.size();i++)
        {
            SimilarityTransform srt;
            if (EstimateSRT(global_reconstruction,*reconstructions[i], kMaxReprojError,srt))
            {
                CHECK_LE(srtform_to_global.count(reconstructions[i]),0);
                    srtform_to_global.insert(std::make_pair(reconstructions[i],srt));
                std::cout<<"Recoonstruction "<<i<<"("<<cluster_name.at(reconstructions[i])<<"):\n"
                        <<"Qvec addresss:"<<srtform_to_global.at(reconstructions[i]).Qvec.data()<<"\n"
                        <<"Qvec:"<<srtform_to_global.at(reconstructions[i]).Qvec<<"\n"
                        <<"tvec:"<<srtform_to_global.at(reconstructions[i]).tvec<<"\n"
                        <<"scale:"<<srtform_to_global.at(reconstructions[i]).scale<<"\n";
            }
            else
            {
                std::cout<<"//////////////\n        Merge Fail\n//////////////\n";
                return;
            }
        }
        //////////////////////////////////////////////////////

        ///////////////////////////////
        //write for test
        std::string debug_path = JoinPaths(mergeoptions_.output_path,"0_withoutba");
        CreateDirIfNotExists(debug_path);
        global_reconstruction.Write(debug_path);
        global_reconstruction.WriteText(debug_path);
        WriteKittiPose("pose_withoutba.txt",global_reconstruction);
        //////////////////

        //////////////////////////////////////////////////////////////////////////
        //call BundleAdjustment
        //////////////////////////////////////////////////////////////////////////

        PrintHeading2("BundleAdjustment for optimal SimilarityTransform");
        std::cout<<"/////";
        std::cout<<"before adjustment MSE:"<<MeanSquaredErrorByRecons(global_reconstruction)<<"//////\n";
        

        OptionManager options;
        options.AddBundleAdjustmentOptions();
        
        
        if(mergeoptions_.useba)
        {
            DACSFM::MyBundleAdjustmentController::Solve(options,reconstructions ,srtform_to_maximum);
        }
        
        for(size_t i = 0 ; i < reconstructions.size();i++)
        {
            std::cout<<"Recoonstruction "<<i<<"("<<cluster_name.at(reconstructions[i])<<"):\n"
                            <<"Qvec:"<<srtform_to_maximum.at(reconstructions[i]).Qvec<<"\n"
                            <<"tvec:"<<srtform_to_maximum.at(reconstructions[i]).tvec<<"\n"
                            <<"scale:"<<srtform_to_maximum.at(reconstructions[i]).scale<<"\n";
        }

        //////////////////////////////////////////////////////////////////////////
        //Merge
        //////////////////////////////////////////////////////////////////////////
        PrintHeading2("Merge");
        Reconstruction* guided_reconstruction=reconstructions[0];
        while(reconstructions.size() > 1)
        {
            bool merge_success = false;
            printf("/////////\nThe %02d loops\n/////////\n",++index);
            for(size_t i=0;i<reconstructions.size();i++) 
                std::cout<<"Reconstruction "<<i<<"("<<cluster_name.at(reconstructions[i])<<")"<<".NumRegImages:"<<reconstructions[i]->NumRegImages()<<"\n"<<"      NumPoints3D:"<<reconstructions[i]->NumPoints3D()<<"\n";
            for(size_t i = 1;i<reconstructions.size();i++)
            {
                if(!srtform_to_maximum.count(reconstructions[i])){
                    std::cerr<<"Reconstruction "<<i<<"("<<cluster_name.at(reconstructions[i])<<") miss similaritytransform to global.\n";
                    continue;
                }
                auto& srt = srtform_to_maximum.at(reconstructions[i]);
                SimilarityTransform3 srt3(srt.scale,srt.Qvec,srt.tvec);
                if (MergeReconstructionsAmbiguously(*guided_reconstruction,*reconstructions[i], kMaxReprojError,
                          srt3)) {
                    std::cout<<"(0,"<<i<<")"<<"merge_success!\n";
                    std::cout<<"("<<cluster_name.at(reconstructions[0])<<","<<cluster_name.at(reconstructions[i])<<") merge_success!\n";
                    cluster_name.erase(reconstructions[i]);
                    reconstructions.erase(reconstructions.begin()+i);
                    
                    merge_success = true;
                    break;
                }
            }
            if(!merge_success) break;
        }
        // /////////////////////////////////
        // //compard with ogba 
        // /////////////////////////////////
        // DACSFM::MyBundleAdjustmentController ba_controller(options, &global_reconstruction);
        // //colmap::BundleAdjustmentController ba_controller(options, &global_reconstruction);
        // ba_controller.Start();
        // ba_controller.Wait();
        // std::cout<<"/////";
        // std::cout<<"after ogba MSE:"<<MeanSquaredErrorByRecons(global_reconstruction)<<"//////\n";
        // std::string ogba_path = JoinPaths(mergeoptions_.output_path,"ogba");
        // CreateDirIfNotExists(ogba_path);
        // global_reconstruction.Write(ogba_path);
        // global_reconstruction.WriteText(ogba_path);
        // WriteKittiPose("pose_ogb.txt",global_reconstruction);

        ////////////////////////////////////////////
        std::cout<<"/////";
        std::cout<<"after adjustment MSE:"<<MeanSquaredErrorByRecons(*reconstructions[0])<<"//////\n";
        
        /////////////////////////////
        /////Debug
        if(reconstructions.size() != 1)
        {
            std::string debug_path = JoinPaths(mergeoptions_.output_path,"Debug");
            CreateDirIfNotExists(debug_path);
            for (size_t i = 0; i < reconstructions.size(); ++i) {
                const std::string reconstruction_path = JoinPaths(debug_path, 
                                                                    "partition_" + std::to_string(i));
                CreateDirIfNotExists(reconstruction_path);
                reconstructions[i]->Write(reconstruction_path);
                reconstructions[i]->WriteText(reconstruction_path);
            }
        }
        /////////////////////////////
        CHECK_EQ(reconstructions.size(),1);
        //////////////////////////////
        ///////////write pose
        WriteKittiPose("pose.txt",*reconstructions[0]);
        /////////////////////////////

        reconstruction_manager->Add();
        reconstruction_manager->Get(reconstruction_manager->Size() - 1) =*reconstructions[0]; 
    }

    bool MergeController::MergeImagesToObtainSRT(Reconstruction& main_reconstruction,const Reconstruction& sub_reconstruction,
                           const double max_reproj_error,SimilarityTransform* srt) {
        std::cout<<"MergeImagesToObtainSRT\n";
        const double kMinInlierObservations = 0.3;

        Eigen::Matrix3x4d alignment;
        if (!ComputeAlignmentBetweenReconstructions(sub_reconstruction, main_reconstruction,
                                                    kMinInlierObservations,
                                                    max_reproj_error, &alignment)) {
            return false;
        }
    
        const SimilarityTransform3 tform(alignment);    //!!!

        // Find common and missing images in the two reconstructions.

        for (const auto& image_id : sub_reconstruction.RegImageIds()) {
            if (!main_reconstruction.ExistsImage(image_id)) {   //main_reconstruction
                auto reg_image = sub_reconstruction.Image(image_id);
                reg_image.SetRegistered(false);
                main_reconstruction.AddImage(reg_image);
                main_reconstruction.RegisterImage(image_id);
                auto& image = main_reconstruction.Image(image_id); //main_reconstruction
                tform.TransformPose(&image.Qvec(), &image.Tvec());
            } 
        }

        srt->Update(tform);    //!!!!

        return true;
    }


    bool MergeController::MergeReconstructionsUnambiguously(Reconstruction& main_reconstruction,const Reconstruction& sub_reconstruction,
                           const double max_reproj_error,const SimilarityTransform& tform) {

        // Find common and missing images in the two reconstructions.

        std::unordered_set<image_t> common_image_ids;
        common_image_ids.reserve(sub_reconstruction.NumRegImages());
        std::unordered_set<image_t> missing_image_ids;
        missing_image_ids.reserve(sub_reconstruction.NumRegImages());

        for (const auto& image_id : sub_reconstruction.RegImageIds()) {
            if (main_reconstruction.ExistsImage(image_id)) {   //main_reconstruction
            common_image_ids.insert(image_id);
            } else {
            missing_image_ids.insert(image_id);
            }
        }

        // Register the missing images in this sub_reconstruction.

        for (const auto image_id : missing_image_ids) {
            auto reg_image = sub_reconstruction.Image(image_id);
            reg_image.SetRegistered(false);
            //main_reconstruction
            main_reconstruction.AddImage(reg_image);
            main_reconstruction.RegisterImage(image_id);

            if (!main_reconstruction.ExistsCamera(reg_image.CameraId())) { //main_reconstruction
            main_reconstruction.AddCamera(sub_reconstruction.Camera(reg_image.CameraId()));//main_reconstruction
            }
            auto& image = main_reconstruction.Image(image_id); //main_reconstruction
            tform.TransformPose(&image.Qvec(), &image.Tvec());
        }

        // Merge the two point clouds using the following two rules:
        //    - copy points to this sub_reconstruction with non-conflicting tracks,
        //      i.e. points that do not have an already triangulated observation
        //      in this sub_reconstruction.
        //    - merge tracks that are unambiguous, i.e. only merge points in the two
        //      sub_reconstructions if they have a one-to-one mapping.
        // Note that in both cases no cheirality or reprojection test is performed.

        for (const auto& point3D : sub_reconstruction.Points3D()) {
            Track new_track;
            //Track old_track;
            //std::set<point3D_t> old_point3D_ids;
            for (const auto& track_el : point3D.second.Track().Elements()) {
            if (common_image_ids.count(track_el.image_id) > 0) {
                const auto& point2D =
                    main_reconstruction.Image(track_el.image_id).Point2D(track_el.point2D_idx);
                if (!point2D.HasPoint3D()) {
                    new_track.AddElement(track_el);
                }
            } else if (missing_image_ids.count(track_el.image_id) > 0) {
                main_reconstruction.Image(track_el.image_id).ResetPoint3DForPoint2D(track_el.point2D_idx);
                new_track.AddElement(track_el);
            }
            }

            const bool create_new_point = new_track.Length() >= 2;
            // const bool merge_new_and_old_point =
            //     (new_track.Length() + old_track.Length()) >= 2 &&
            //     old_point3D_ids.size() == 1;
            if (create_new_point) {
            Eigen::Vector3d xyz = point3D.second.XYZ();
            tform.TransformPoint(&xyz);
            const auto point3D_id =
                main_reconstruction.AddPoint3D(xyz, new_track, point3D.second.Color());//main_reconstruction
            // if (old_point3D_ids.size() == 1) {
            //     main_reconstruction.MergePoints3D(point3D_id, *old_point3D_ids.begin());//main_reconstruction
            // }
            }
        }
        //main_reconstruction
        FilterPoints3DWithLargeReprojectionError(main_reconstruction,max_reproj_error, main_reconstruction.Point3DIds());

        return true;
    }



    bool MergeController::MergeAllForSRT(Reconstruction& main_reconstruction,const Reconstruction& sub_reconstruction,
                           const double max_reproj_error,SimilarityTransform& srt) {
        const double kMinInlierObservations = 0.3;

        Eigen::Matrix3x4d alignment;
        if (!ComputeAlignmentBetweenReconstructions(sub_reconstruction, main_reconstruction,
                                                    kMinInlierObservations,
                                                    max_reproj_error, &alignment)) {
            return false;
        }

        const SimilarityTransform3 tform(alignment);
        srt.Update(tform);
        // Find common and missing images in the two reconstructions.

        std::unordered_set<image_t> common_image_ids;
        common_image_ids.reserve(sub_reconstruction.NumRegImages());
        std::unordered_set<image_t> missing_image_ids;
        missing_image_ids.reserve(sub_reconstruction.NumRegImages());

        for (const auto& image_id : sub_reconstruction.RegImageIds()) {
            if (main_reconstruction.ExistsImage(image_id)) {   //main_reconstruction
            common_image_ids.insert(image_id);
            } else {
            missing_image_ids.insert(image_id);
            }
        }

        // Register the missing images in this sub_reconstruction.

        for (const auto image_id : missing_image_ids) {
            auto reg_image = sub_reconstruction.Image(image_id);
            reg_image.SetRegistered(false);
            //main_reconstruction
            main_reconstruction.AddImage(reg_image);
            main_reconstruction.RegisterImage(image_id);

            if (!main_reconstruction.ExistsCamera(reg_image.CameraId())) { //main_reconstruction
            main_reconstruction.AddCamera(sub_reconstruction.Camera(reg_image.CameraId()));//main_reconstruction
            }
            auto& image = main_reconstruction.Image(image_id); //main_reconstruction
            tform.TransformPose(&image.Qvec(), &image.Tvec());
        }

        // Merge the two point clouds using the following two rules:
        //    - copy points to this sub_reconstruction with non-conflicting tracks,
        //      i.e. points that do not have an already triangulated observation
        //      in this sub_reconstruction.
        //    - merge tracks that are unambiguous, i.e. only merge points in the two
        //      sub_reconstructions if they have a one-to-one mapping.
        // Note that in both cases no cheirality or reprojection test is performed.

        for (const auto& point3D : sub_reconstruction.Points3D()) {
            Track new_track;
            Track old_track;
            std::set<point3D_t> old_point3D_ids;
            for (const auto& track_el : point3D.second.Track().Elements()) {
            if (common_image_ids.count(track_el.image_id) > 0) {
                const auto& point2D =
                    main_reconstruction.Image(track_el.image_id).Point2D(track_el.point2D_idx);
                if (point2D.HasPoint3D()) {
                old_track.AddElement(track_el);
                old_point3D_ids.insert(point2D.Point3DId());
                } else {
                new_track.AddElement(track_el);
                }
            } else if (missing_image_ids.count(track_el.image_id) > 0) {
                main_reconstruction.Image(track_el.image_id).ResetPoint3DForPoint2D(track_el.point2D_idx);
                new_track.AddElement(track_el);
            }
            }

            const bool create_new_point = new_track.Length() >= 2;
            const bool merge_new_and_old_point =
                (new_track.Length() + old_track.Length()) >= 2 &&
                old_point3D_ids.size() == 1;
            if (create_new_point || merge_new_and_old_point) {
            Eigen::Vector3d xyz = point3D.second.XYZ();
            tform.TransformPoint(&xyz);
            const auto point3D_id =
                main_reconstruction.AddPoint3D(xyz, new_track, point3D.second.Color());//main_reconstruction
            if (old_point3D_ids.size() == 1) {
                main_reconstruction.MergePoints3D(point3D_id, *old_point3D_ids.begin());//main_reconstruction
            }
            }
        }
        //main_reconstruction
        FilterPoints3DWithLargeReprojectionError(main_reconstruction,max_reproj_error, main_reconstruction.Point3DIds());

        return true;
    }

    double MergeController::MeanSquaredErrorByRecons(const Reconstruction& reconstruction)
    {
        double sum = 0.0;
        size_t count = 0;
        double sqrt_sum = 0.0;
        bool debug_flag = false;
        for(const auto& point3d:reconstruction.Points3D())
        {
            for(const auto& track_el:point3d.second.Track().Elements())
            {
                const Image& image = reconstruction.Image(track_el.image_id);
                const Point2D& point2d = image.Point2D(track_el.point2D_idx);
                const Camera& camera = reconstruction.Camera(image.CameraId());
                double squard_error = CalculateSquaredReprojectionError(point2d.XY(),point3d.second.XYZ(),image.Qvec(),image.Tvec(),camera);
                if(squard_error>1e18)
                {
                    std::cout<<"False Points\n";
                    continue;
                }
                sum += squard_error;
                
                sqrt_sum +=std::sqrt(squard_error);
                count ++;
                //debug_flag = true;
                if(debug_flag) break;
            }
            if(debug_flag) break;
        }
        //std::cout<<"sum:"<<sum<<"\n";
        //std::cout<<"count:"<<count<<"\n";
        //return sum/1.0*count;
        /*std::cout<<"sum/count = "<<sum/(1.0*count)<<"\n";
        double sum2 = 0.0;
        size_t count2 = 0;
        for(const image_t image_id : reconstruction.RegImageIds())
        {
            const Image& image = reconstruction.Image(image_id);
            const Camera& camera = reconstruction.Camera(image.CameraId());
            for(const Point2D& point2d:image.Points2D())
            {
                if(point2d.HasPoint3D())
                {
                    const Point3D& point3d = reconstruction.Point3D(point2d.Point3DId());
                    sum2 += CalculateSquaredReprojectionError(point2d.XY(),point3d.XYZ(),image.Qvec(),image.Tvec(),camera);
                    count2 ++;
                }
            }
        }

        std::cout<<"sum2:"<<sum2<<"\n";
        std::cout<<"count2:"<<count2<<"\n";
        
        std::cout<<"sum2/count2 = "<<sum2/(1.0*count2)<<"\n";*/
        std::cout<<"(Avg. Reproj. Error:"<<sqrt_sum/(1.0*count)<<")";
        return sum/(1.0*count);
    }

    void MergeController::MergeClustersintoGlobal(
    std::vector<Reconstruction*>& reconstructions,
    ReconstructionManager* reconstruction_manager,
    std::unordered_map<Reconstruction*,std::string>& cluster_name)
    {
        
        PrintHeading2("Begin Combination");
        

        const auto cmp = [](const Reconstruction* reconstruction1,const Reconstruction* reconstruction2){
            if(reconstruction1->NumRegImages() == reconstruction2->NumRegImages())
                return reconstruction1->NumPoints3D() > reconstruction2->NumPoints3D();
            return reconstruction1->NumRegImages()>reconstruction2->NumRegImages();  //??NumRegImages??NumPoints3D??
        };


        std::sort(reconstructions.begin(),reconstructions.end(),cmp);
        PrintHeading2("reconstructions sorted");
        /*size_t num_recons = reconstructions.size();
        size_t max_size = 2;
        size_t recons_erase_index = 0;
        for(size_t i = 0 ;i<num_recons;i++)
        {
            if(recons_erase_index >= max_size)
            {
                reconstructions.erase(reconstructions.begin()+recons_erase_index);
            }
            else
            {
                recons_erase_index++;
            }
            
        }*/
        
        for(size_t i=0;i<reconstructions.size();i++) 
        {
                std::cout<<"Reconstruction "<<i<<"("<<cluster_name.at(reconstructions[i])<<")"
                         <<".NumRegImages:"<<reconstructions[i]->NumRegImages()<<"\n"
                         <<"      NumPoints3D:"<<reconstructions[i]->NumPoints3D()<<"\n";
        }
        //////////////////////////////////////////////////////////////////////////
        //estimate  srt from all sub-reconstructions to the global reconstrcution
        //////////////////////////////////////////////////////////////////////////
        PrintHeading2("Virtual Merge");
        std::cout<<"PrintHeading2('Virtual Merge');"<<"\n";
        CHECK_GE(reconstructions.size(),1);
        std::cout<<"CHECK_GE(reconstructions.size(),1);"<<"\n";
        const double kMaxReprojError = 8.0;
        
        //Intialize the global reconstuction with largest reconstruction
        std::vector<bool> is_merged(reconstructions.size(),false);
        std::cout<<""<<"\n";
        
        Reconstruction global_reconstruction = *reconstructions[0];
        std::unordered_map<Reconstruction*,SimilarityTransform> srtform_to_maximum;
        
        srtform_to_maximum.insert(std::make_pair(reconstructions[0],SimilarityTransform()));
        is_merged[0] = true;
        int index = 0;
        const int tot_recons = reconstructions.size();
        std::cout<<"const int tot_recons = reconstructions.size();"<<"\n";
        while(index < tot_recons)
        {
            bool merge_success = false;
            printf("/////////\nThe %02d loops\n/////////\n",++index);
            std::cout<<"global_reconstruction "<<".NumRegImages:"<<global_reconstruction.NumRegImages()<<"\n"<<"      NumPoints3D:"<<global_reconstruction.NumPoints3D()<<"\n";
            for(size_t i=0;i<reconstructions.size();i++) 
                if(!is_merged[i])std::cout<<"Reconstruction "<<i<<"("<<cluster_name.at(reconstructions[i])<<")"<<".NumRegImages:"<<reconstructions[i]->NumRegImages()<<"\n"<<"      NumPoints3D:"<<reconstructions[i]->NumPoints3D()<<"\n";
            for(size_t i = 1;i<reconstructions.size();i++)
            {
                if(is_merged[i]) continue;
                SimilarityTransform srt;
                if (MergeAllForSRT(global_reconstruction,*reconstructions[i], kMaxReprojError,srt)) {

                    CHECK_LE(srtform_to_maximum.count(reconstructions[i]),0);
                    srtform_to_maximum.insert(std::make_pair(reconstructions[i],srt));

                    std::cout<<"(0,"<<i<<")"<<"estimates srt successfully!\n";
                    std::cout<<"("<<cluster_name.at(reconstructions[0])<<","<<cluster_name.at(reconstructions[i])<<") estimates srt successfully!\n";
                    is_merged[i] = true;
                    merge_success = true;
                    break;
                }
            }
            if(!merge_success) break;
        }
        for(size_t i = 0 ; i < reconstructions.size();i++)
        {
            if(is_merged[i]==false)
            {
                std::cout<<"//////////////\n        Merge Fail\n//////////////\n";
                return;
            }
            std::cout<<"Recoonstruction "<<i<<"("<<cluster_name.at(reconstructions[i])<<"):\n"
                            <<"Qvec addresss:"<<srtform_to_maximum.at(reconstructions[i]).Qvec.data()<<"\n"
                            <<"Qvec:"<<srtform_to_maximum.at(reconstructions[i]).Qvec<<"\n"
                            <<"tvec:"<<srtform_to_maximum.at(reconstructions[i]).tvec<<"\n"
                            <<"scale:"<<srtform_to_maximum.at(reconstructions[i]).scale<<"\n";
        }
        //////////////////////////////////////////////////////
        ////////////////estimate the srt from local to global
        //////////////////////////////////////////////////////
        std::cout<<"////////////////estimate the srt from local to global\n";
        std::unordered_map<Reconstruction*,SimilarityTransform> srtform_to_global;
        for(size_t i = 0 ; i < reconstructions.size();i++)
        {
            SimilarityTransform srt;
            if (EstimateSRT(global_reconstruction,*reconstructions[i], kMaxReprojError,srt))
            {
                CHECK_LE(srtform_to_global.count(reconstructions[i]),0);
                    srtform_to_global.insert(std::make_pair(reconstructions[i],srt));
                std::cout<<"Recoonstruction "<<i<<"("<<cluster_name.at(reconstructions[i])<<"):\n"
                        <<"Qvec addresss:"<<srtform_to_global.at(reconstructions[i]).Qvec.data()<<"\n"
                        <<"Qvec:"<<srtform_to_global.at(reconstructions[i]).Qvec<<"\n"
                        <<"tvec:"<<srtform_to_global.at(reconstructions[i]).tvec<<"\n"
                        <<"scale:"<<srtform_to_global.at(reconstructions[i]).scale<<"\n";
            }
            else
            {
                std::cout<<"//////////////\n        Merge Fail\n//////////////\n";
                return;
            }

            
        }



        //////////////////////////////////////////////////////////////////////////
        //call BundleAdjustment
        //////////////////////////////////////////////////////////////////////////

        PrintHeading2("BundleAdjustment for optimal SimilarityTransform");
        std::cout<<"/////";
        std::cout<<"before adjustment MSE:"<<MeanSquaredErrorByRecons(global_reconstruction)<<"//////\n";
        

        OptionManager options;
        options.AddBundleAdjustmentOptions();
        // invoke BA to obtain the optimal SimilarityTransform.
        /*{
            std::vector<Reconstruction*> tmp_reoncs;
            std::cout<<"///////////////exchange vector except last///////////////";
            size_t tmp_size = (reconstructions.size()-1);
            for(size_t i = 0 ; i < tmp_size;i++)
            {
                std::cout<<"reconstructions[i]: "<<reconstructions[i]<<"\n";
                tmp_reoncs.emplace_back(reconstructions[i]);
                std::cout<<"tmp_reoncs[i]: "<<tmp_reoncs[i]<<"\n";
            }
            DACSFM::MyBundleAdjustmentController::Solve(options,tmp_reoncs ,srtform_to_global);
        }*/
        if(mergeoptions_.useba)
        {
            DACSFM::MyBundleAdjustmentController::Solve(options,reconstructions ,srtform_to_global);
        }
        /*for(size_t i = 0 ; i < reconstructions.size();i++)
        {
            std::cout<<"BundleAdjustment for Reconstruction "<<i<<"////////////////////////\n";
            std::vector<Reconstruction*>tmp_recon{reconstructions[i]};
            DACSFM::MyBundleAdjustmentController::Solve(options,tmp_recon ,srtform_to_global);
            //BundleAdjustmentController ba_controller(options, reconstructions[i]);
            //ba_controller.Start();
            //ba_controller.Wait();
        }*/
        for(size_t i = 0 ; i < reconstructions.size();i++)
        {
            std::cout<<"Recoonstruction "<<i<<"("<<cluster_name.at(reconstructions[i])<<"):\n"
                            <<"Qvec:"<<srtform_to_global.at(reconstructions[i]).Qvec<<"\n"
                            <<"tvec:"<<srtform_to_global.at(reconstructions[i]).tvec<<"\n"
                            <<"scale:"<<srtform_to_global.at(reconstructions[i]).scale<<"\n";
        }

        //////////////////////////////////////////////////////////////////////////
        //Merge
        //////////////////////////////////////////////////////////////////////////
        PrintHeading2("Merge");
        Reconstruction* guided_reconstruction=reconstructions[0];
        while(reconstructions.size() > 1)
        {
            bool merge_success = false;
            printf("/////////\nThe %02d loops\n/////////\n",++index);
            for(size_t i=0;i<reconstructions.size();i++) 
                std::cout<<"Reconstruction "<<i<<"("<<cluster_name.at(reconstructions[i])<<")"<<".NumRegImages:"<<reconstructions[i]->NumRegImages()<<"\n"<<"      NumPoints3D:"<<reconstructions[i]->NumPoints3D()<<"\n";
            for(size_t i = 1;i<reconstructions.size();i++)
            {
                if(!srtform_to_global.count(reconstructions[i])){
                    std::cerr<<"Reconstruction "<<i<<"("<<cluster_name.at(reconstructions[i])<<") miss similaritytransform to global.\n";
                    continue;
                }
                auto& srt = srtform_to_global.at(reconstructions[i]);
                SimilarityTransform3 srt3(srt.scale,srt.Qvec,srt.tvec);
                if (MergeReconstructionsAmbiguously(*guided_reconstruction,*reconstructions[i], kMaxReprojError,
                          srt3)) {
                    std::cout<<"(0,"<<i<<")"<<"merge_success!\n";
                    std::cout<<"("<<cluster_name.at(reconstructions[0])<<","<<cluster_name.at(reconstructions[i])<<") merge_success!\n";
                    cluster_name.erase(reconstructions[i]);
                    reconstructions.erase(reconstructions.begin()+i);
                    
                    merge_success = true;
                    break;
                }
            }
            if(!merge_success) break;
        }
        std::cout<<"/////";
        std::cout<<"after adjustment MSE:"<<MeanSquaredErrorByRecons(*reconstructions[0])<<"//////\n";

        /////Debug
        if(reconstructions.size() != 1)
        {
            std::string debug_path = JoinPaths(mergeoptions_.output_path,"Debug");
            CreateDirIfNotExists(debug_path);
            for (size_t i = 0; i < reconstructions.size(); ++i) {
                const std::string reconstruction_path = JoinPaths(debug_path, 
                                                                    "partition_" + std::to_string(i));
                CreateDirIfNotExists(reconstruction_path);
                reconstructions[i]->Write(reconstruction_path);
                reconstructions[i]->WriteText(reconstruction_path);
            }
        }
        /////
        CHECK_EQ(reconstructions.size(),1);
        

        reconstruction_manager->Add();
        reconstruction_manager->Get(reconstruction_manager->Size() - 1) =*reconstructions[0]; 
    }


    
    bool MergeController::EstimateSRT(Reconstruction& main_reconstruction,const Reconstruction& sub_reconstruction,
                           const double max_reproj_error,SimilarityTransform& srt) {
        const double kMinInlierObservations = 0.3;

        Eigen::Matrix3x4d alignment;
        if (!ComputeAlignmentBetweenReconstructions(sub_reconstruction, main_reconstruction,
                                                    kMinInlierObservations,
                                                    max_reproj_error, &alignment)) {
            return false;
        }

        const SimilarityTransform3 tform(alignment);
        srt.Update(tform);
        // Find common and missing images in the two reconstructions.
        return true;

    }

    void MergeController::WriteWithoutBA(const Reconstruction& reconstruction) {
        std::string debug_path = JoinPaths(mergeoptions_.output_path,"0_withoutba");
        CreateDirIfNotExists(debug_path);
        reconstruction.Write(debug_path);
        reconstruction.WriteText(debug_path);

        WriteKittiPose("pose_withoutba.txt",reconstruction);

        std::cout<<"///write reconstruction and pose withou ba successfully !\n";

    }

    void MergeController::WriteKittiPose(const std::string filename,
                                         const Reconstruction& reconstruction)
    {
        //write pose
        std::cout<<"write Kitti pose\n";
        std::ofstream file(JoinPaths(mergeoptions_.output_path,filename), std::ios::trunc);
        std::vector<image_t> imageid_sorted=reconstruction.RegImageIds();
        std::sort(imageid_sorted.begin(),imageid_sorted.end());
        for(image_t image_id:imageid_sorted)
        {
          const Image& image = reconstruction.Image(image_id);
          Eigen::Matrix3x4d pose_w_c = image.InverseProjectionMatrix();
          file<<std::setprecision(9)<<pose_w_c(0,0)<<" "<<pose_w_c(0,1)<<" "<<pose_w_c(0,2)<<" "<<pose_w_c(0,3)<<" "
                                    <<pose_w_c(1,0)<<" "<<pose_w_c(1,1)<<" "<<pose_w_c(1,2)<<" "<<pose_w_c(1,3)<<" "
                                    <<pose_w_c(2,0)<<" "<<pose_w_c(2,1)<<" "<<pose_w_c(2,2)<<" "<<pose_w_c(2,3)<<"\n";
        }
        file.close();
        
    }

    void MergeController::WriteEuRoC(const std::string filename,
                                     const Reconstruction& reconstruction)
    {
        std::cout<<"write EuRoC pose\n";
        std::ofstream file(JoinPaths(mergeoptions_.output_path,filename), std::ios::trunc);
        for(image_t image_id:reconstruction.RegImageIds())
        {
            const Image& image = reconstruction.Image(image_id);
            std::string timestamp,ext;
            SplitFileExtension(image.Name(),&timestamp,&ext);
            Eigen::Matrix3x4d pose_w_c = image.InverseProjectionMatrix();
            Eigen::Vector4d q_w_c = NormalizeQuaternion(RotationMatrixToQuaternion(pose_w_c.block<3,3>(0,0)));  //qw qx qy qz
            file<<timestamp<<" "<<std::setprecision(9)<<pose_w_c(0,3)<<" "<<pose_w_c(1,3)<<" "<<pose_w_c(2,3)<<" "
                                     <<q_w_c(0)<<" "<<q_w_c(1)<<" "<<q_w_c(2)<<" "<<q_w_c(3)<<"\n";
        } 
        file.close();

    }

}
