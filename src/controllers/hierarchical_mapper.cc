// Copyright (c) 2018, ETH Zurich and UNC Chapel Hill.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//
//     * Neither the name of ETH Zurich and UNC Chapel Hill nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: Johannes L. Schoenberger (jsch-at-demuc-dot-de)

#include "controllers/hierarchical_mapper.h"
#include "merge/mergeclusters.h"

#include "base/scene_clustering.h"
#include "util/misc.h"
#include <iomanip>
using namespace DACSFM;
namespace colmap {
namespace {

void MergeClusters(
    const SceneClustering::Cluster& cluster,
    std::unordered_map<const SceneClustering::Cluster*, ReconstructionManager>*
        reconstruction_managers) {
  // Extract all reconstructions from all child clusters.
  std::vector<Reconstruction*> reconstructions;
  for (const auto& child_cluster : cluster.child_clusters) {
    if (!child_cluster.child_clusters.empty()) {
      MergeClusters(child_cluster, reconstruction_managers);
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

}  // namespace

bool HierarchicalMapperController::Options::Check() const {
  CHECK_OPTION_GT(init_num_trials, -1);
  CHECK_OPTION_GE(num_workers, -1);
  return true;
}

HierarchicalMapperController::HierarchicalMapperController(
    const Options& options, const SceneClustering::Options& clustering_options,
    const IncrementalMapperOptions& mapper_options,
    ReconstructionManager* reconstruction_manager)
    : options_(options),
      clustering_options_(clustering_options),
      mapper_options_(mapper_options),
      reconstruction_manager_(reconstruction_manager) {
  CHECK(options_.Check());
  CHECK(clustering_options_.Check());
  CHECK(mapper_options_.Check());
  //CHECK_EQ(clustering_options_.branching, 2);
}

void HierarchicalMapperController::Run() {
  PrintHeading1("Partitioning the scene");

  ////////////HierarchicalMapperController Test/////////////
  if(0){
    if(CombineReconswithExistingPartitions()) return;
    
  }
  //////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  // Cluster scene
  //////////////////////////////////////////////////////////////////////////////
  auto start = std::chrono::system_clock::now();
  auto start_sceneclustering = std::chrono::system_clock::now();

  SceneClustering scene_clustering(clustering_options_);

  std::unordered_map<image_t, std::string> image_id_to_name;

  {
    Database database(options_.database_path);

    std::cout << "Reading images..." << std::endl;
    const auto images = database.ReadAllImages();
    for (const auto& image : images) {
      image_id_to_name.emplace(image.ImageId(), image.Name());
    }

    std::cout << "Reading scene graph..." << std::endl;
    std::vector<std::pair<image_t, image_t>> image_pairs;
    std::vector<int> num_inliers;

    //////filter image_pairs by num_inliers///
    int min_inliers = 100;
    database.ReadTwoViewGeometryNumInliers(&image_pairs, &num_inliers);
    std::cout<<"image_pairs.size()"<<image_pairs.size()<<"\n";
    {
      size_t index_imagepair = 0;
      size_t num_imagepairs = image_pairs.size();
      for(size_t i=0;i<num_imagepairs;i++)
      {
        if(num_inliers[index_imagepair]<min_inliers)
        {
            image_pairs.erase(image_pairs.begin()+index_imagepair);
            num_inliers.erase(num_inliers.begin()+index_imagepair);
        }
        else
        {
          index_imagepair++;
        }
        
      }
    }
    std::cout<<"image_pairs.size()"<<image_pairs.size()<<"\n";
    //////////////////////////////
    
    std::cout << "Partitioning scene graph..." << std::endl;
    scene_clustering.Partition(image_pairs, num_inliers);
    WriteNormazliedGrpah(image_pairs,num_inliers);
    WriteIntraCluster(image_pairs,num_inliers);
  }

  auto leaf_clusters = scene_clustering.GetLeafClusters();
  
  size_t total_num_images = 0;
  for (size_t i = 0; i < leaf_clusters.size(); ++i) {
    total_num_images += leaf_clusters[i]->image_ids.size();
    std::cout << StringPrintf("  Cluster %d with %d images", i + 1,
                              leaf_clusters[i]->image_ids.size())
              << std::endl;
  }

  std::cout << StringPrintf("Clusters have %d images", total_num_images)
            << std::endl;
  auto end_sceneclustering = std::chrono::system_clock::now();
  ///////////////////////////////////////////////////////////////////////
  // Start reconstructing the bigger clusters first for resource usage.
  /////////////////////////////////////////////////////////////////////
  std::sort(leaf_clusters.begin(), leaf_clusters.end(),
            [](const SceneClustering::Cluster* cluster1,
               const SceneClustering::Cluster* cluster2) {
              return cluster1->image_ids.size() > cluster2->image_ids.size();
            });
  ////////////////////////////////////////////////////////////
  ///// push the images of `leaf_cluster` into `sub_clusters`
  ////////////////////////////////////////////////////////////
  std::vector<SceneClustering::Cluster*> sub_clusters;
  
  size_t num_sole_cluster = 302;
  for (size_t i = 0; i < leaf_clusters.size(); ++i) {
    //if(leaf_clusters[i]->image_ids.size()!=num_sole_cluster) continue;
    sub_clusters.push_back(new SceneClustering::Cluster());
    SceneClustering::Cluster* curr_cluster = *(sub_clusters.rbegin());
    if(leaf_clusters[i]->image_ids.size()==333)
    {
      for(auto image_id:leaf_clusters[i]->image_ids)
      {
          if(image_id>=1&&image_id<=101) continue;
          curr_cluster->image_ids.push_back(image_id);
      }
      continue;
    }
    curr_cluster->image_ids.insert(curr_cluster->image_ids.end(),
                                     leaf_clusters[i]->image_ids.begin(),leaf_clusters[i]->image_ids.end());
    if(curr_cluster->image_ids.size() == 312)
    {
       for(image_t k = 1 ; k<=101 ;k++)
       {
        curr_cluster->image_ids.push_back(k);
       }
    }
  }
  
  WriteLeafClusters("original_cluster.txt",sub_clusters,image_id_to_name);
  //////////////////////////////////////////////////////////////////////////////
  // Reconstruct clusters
  //////////////////////////////////////////////////////////////////////////////
  auto start_reconstruction = std::chrono::system_clock::now();
  PrintHeading1("Reconstructing clusters");

   ////////////HierarchicalMapperController Test/////////////
  if(1){
    if(CombineReconswithExistingPartitions()) return;
    
  }
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////
  // Determine the number of workers and threads per worker.
  //////////////////////////////////////////////////////////////
  const int kMaxNumThreads = -1;
  const int num_eff_threads = GetEffectiveNumThreads(kMaxNumThreads);
  const int kDefaultNumWorkers = 8;
  const int num_eff_workers =
      options_.num_workers < 1
          ? std::min(static_cast<int>(leaf_clusters.size()),
                     std::min(kDefaultNumWorkers, num_eff_threads))
          : options_.num_workers;
  const int num_threads_per_worker =
      std::max(1, num_eff_threads / num_eff_workers);
  
  // Function to reconstruct one cluster using incremental mapping. 
  auto ReconstructCluster = [&, this](
                                const SceneClustering::Cluster& cluster,
                                ReconstructionManager* reconstruction_manager) {
    if (cluster.image_ids.empty()) {
      return;
    }

    IncrementalMapperOptions custom_options = mapper_options_;
    //custom_options.max_model_overlap = 3;
    //custom_options.init_num_trials = options_.init_num_trials;
    custom_options.num_threads = num_threads_per_worker;
    //bool  filtered_cluster_flag1 = (cluster.image_ids.size()==270);
    for (const auto image_id : cluster.image_ids) {
        // if(filtered_cluster_flag1&&
        // ((185<=image_id&&image_id<=188)||(921<=image_id&&image_id<=923))) 
        // {
        //   continue;
        // }
      custom_options.image_names.insert(image_id_to_name.at(image_id));
    }
    //  if(cluster.image_ids.size() == 309){
    //    custom_options.init_image_id1 = 356;
    //    custom_options.init_image_id2 = 352;
    //  }

    IncrementalMapperController mapper(&custom_options, options_.image_path,
                                       options_.database_path,
                                       reconstruction_manager);
    mapper.Start();
    mapper.Wait();
  };
  

  /////////////////////////////////////////////////////
  // Start the reconstruction workers.
  /////////////////////////////////////////////////////
  std::unordered_map<SceneClustering::Cluster*, ReconstructionManager>
      reconstruction_managers;
  reconstruction_managers.reserve(leaf_clusters.size());

  
  size_t p_subclusters = 0;
  size_t num_turns = 0;
  
  //SceneClustering::Cluster tmp1= sub_clusters[0];
  //SceneClustering::Cluster tmp2=  sub_clusters[1];
  size_t standard_re_reconstruction = 35;
  
  while(p_subclusters < sub_clusters.size()){

    std::cout<<"/////////////////////////\n"                
             <<"     num_turns "<<num_turns<<"  \n"
             <<"/////////////////////////\n";  
    ThreadPool thread_pool(num_eff_workers);
    
    
      for (size_t i = p_subclusters; i <  sub_clusters.size(); i++) {
          thread_pool.AddTask(ReconstructCluster, *sub_clusters[i],
                            &reconstruction_managers[sub_clusters[i]]);
         
    }                  

    thread_pool.Wait();
    size_t pre_subcluster = p_subclusters;
    p_subclusters = sub_clusters.size();
    
    std::cout<<"pre_subcluster:"<<pre_subcluster<<"\n";
    std::cout<<"sub_clusters.size()"<<sub_clusters.size()<<"\n";
    
    //FilterImageRegisteredInCluster(pre_subcluster,standard_re_reconstruction,sub_clusters,reconstruction_managers);
    //sub_clusters.emplace_back(tmp1);
    //sub_clusters.emplace_back(tmp2);
    std::cout<<"p_subclusters:"<<p_subclusters<<"\n";
    std::cout<<"sub_clusters.size()"<<sub_clusters.size()<<"\n";
    num_turns++ ;
  }
  //WriteLeafClusters("re_cluster.txt",sub_clusters,image_id_to_name); //*//

  auto end_reconstruction = std::chrono::system_clock::now();
  ///////////////////////////////////////////////////
  ///////////////////////////////////////////////////
  
  /*************************write for debug**********************************/
  std::vector<Reconstruction*> reconstructions;
  std::unordered_map<Reconstruction*,std::string> clusters_name;
  WriteClustersPartition(sub_clusters,image_id_to_name,reconstruction_managers,reconstructions,clusters_name);

  /************************************************************************/

  //////////////////////////////////////////////////////////////////////////////
  // My Merge clusters
  //////////////////////////////////////////////////////////////////////////////

  PrintHeading1("My Merging clusters");


  ////////////HierarchicalMapperController Test/////////////
  if(0){
    if(CombineReconswithExistingPartitions()) return;
    
  }
  //////////////////////////////////////////////////////////
  auto start_merge = std::chrono::system_clock::now();
  MergeController mergecontroller(options_.output_path);
  mergecontroller.MergeClustersintoMaximum(reconstructions, reconstruction_manager_,clusters_name);
  auto end_merge = std::chrono::system_clock::now();
  auto end = std::chrono::system_clock::now();
   
  std::cout << std::endl;
  GetTimer().PrintMinutes();
  TearDown(sub_clusters,end-start,end_sceneclustering-start_sceneclustering,end_reconstruction-start_reconstruction,
           end_merge-start_merge);  //delete all new object and write the pose
}


void HierarchicalMapperController::TearDown(std::vector<SceneClustering::Cluster*>& sub_clusters,
                                            std::chrono::duration<double> duration,
                                            std::chrono::duration<double> duration_sceneclustering,
                                            std::chrono::duration<double> duration_reconstruction,
                                            std::chrono::duration<double> duration_merge)
{
    for(SceneClustering::Cluster* cluster:sub_clusters)
    {
      delete cluster;
    }
    sub_clusters.clear();

    std::cout<<"/////////////hierarchical duration/////////////\n"
             <<"sceneclustering: "<<duration_sceneclustering.count()<<" s\n"
             <<"reconstruction: "<<duration_reconstruction.count()<<" s\n"
             <<"merge: "<<duration_merge.count()<<" s\n"
             <<"holistic: "<<duration.count()<<" s\n"
             <<"/////////////hierarchical duration/////////////\n";

}



void HierarchicalMapperController::WriteLeafClusters(std::string filename,const std::vector<SceneClustering::Cluster*>& clusters,
                                                     const std::unordered_map<image_t, std::string>& image_id_to_name) //*//
{
    const std::string clusters_path = JoinPaths(options_.output_path, "clusters");
    CreateDirIfNotExists(clusters_path);
    std::ofstream file(JoinPaths(clusters_path,filename), std::ios::trunc);
    for(size_t i =0;i<clusters.size();i++)
    {
        file<<"Cluster "<<i<<"("<<clusters[i]->image_ids.size()<<")"<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
        std::vector<image_t> image_ids = clusters[i]->image_ids;
        std::sort(image_ids.begin(),image_ids.end());
        for(auto image_id:image_ids)
        {
          file<<"("<<image_id<<")"<<image_id_to_name.at(image_id)<<" ";
        }
        file<<"\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
    }
    file<<"\n\n";
    file.close();
}


bool  HierarchicalMapperController::CombineReconswithExistingPartitions()
{
    printf("Judge If reconstructions exist\n///////////");
    std::vector<std::string> cluster_list;
    
    std::string clusters_path_prio = JoinPaths(options_.output_path,"clusters");
    if(!ExistsDir(clusters_path_prio)) return false;
    cluster_list = GetRecursiveDirList(clusters_path_prio);
    if(!cluster_list.empty())
    {
        printf("Extract Reconstructions existed\n///////////");
        std::vector<Reconstruction*> reconstructions_test(cluster_list.size());
        std::unordered_map<Reconstruction*,std::string> cluster_name;
        for(size_t i = 0;i<cluster_list.size();i++)
        {
          reconstructions_test[i]=new Reconstruction();
          reconstructions_test[i]->Read(cluster_list[i]);
          cluster_name.insert(std::make_pair(reconstructions_test[i],GetPathBaseName(cluster_list[i])));
          std::cout<<cluster_name[reconstructions_test[i]]<<".NumRegImages():"<<reconstructions_test[i]->NumRegImages()<<"\n"
                   <<".NumImages():"<<reconstructions_test[i]->NumImages()<<"\n"
                   <<".NumPoints():"<<reconstructions_test[i]->NumPoints3D()<<"\n\n";
        }
        PrintHeading1("My Merging clusters");
        std::cout<<"reconstructions_test.size()"<<reconstructions_test.size()<<"\n";

        MergeControllerOptions merge_options(options_.output_path);
        merge_options.useba=options_.merge_useba;
        MergeController mergecontroller(merge_options);
        mergecontroller.MergeClustersintoMaximum(reconstructions_test, reconstruction_manager_,cluster_name);


        std::cout<<"reconstructions_test.size()"<<reconstructions_test.size()<<"\n";
        std::cout << std::endl;
        GetTimer().PrintMinutes();


        return true;
    }
    return  false;
}

void HierarchicalMapperController::WriteClustersPartition(std::vector<SceneClustering::Cluster*>& leaf_clusters,
                                                          const std::unordered_map<image_t, std::string>& image_id_to_name,
                                                          std::unordered_map<SceneClustering::Cluster*, ReconstructionManager>& 
                                                            reconstruction_managers,
                                                          std::vector<Reconstruction*>& reconstructions,
                                                          std::unordered_map<Reconstruction*,std::string>& clusters_name)
{
  
  std::vector<int> reconstructions_clusterid;
  
  const std::string clusters_path = JoinPaths(options_.output_path, "clusters");
  std::ofstream file(JoinPaths(clusters_path,"clustersrecons.txt"), std::ios::trunc);
  std::cout<<"leaf_clusters.size()"<<leaf_clusters.size()<<"\n";
  for (size_t i = 0 ; i < leaf_clusters.size();i++){
            auto& recon_manager = reconstruction_managers.at(leaf_clusters[i]);
            
            std::map<image_t,std::string> reg_image_id;
            
            for (size_t j = 0; j < recon_manager.Size(); j++) {
                auto& recons = recon_manager.Get(j);
                reconstructions.push_back(&recons);
                reconstructions_clusterid.push_back(static_cast<int>(i));
                auto& reg_images = recons.RegImageIds();
                std::cout<<"push Cluster "<<i<<",reconstruction "<<j<<" into clustersrecons"<<"\n";
                for(auto image_id:reg_images)
                {
                  if(reg_image_id.count(image_id)==0)
                  {
                    reg_image_id.emplace(std::make_pair(image_id,image_id_to_name.at(image_id)));
                    
                  }
                }
            }
            file<<"Cluster "<<i<<"("<<reg_image_id.size()<<")"<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
            for(auto& item:reg_image_id)
            {
               file<<"("<<item.first<<")"<<item.second<<" ";
            }
            file<<"\n\n";
            

  }
  file.close();
  std::cout<<"write the reconstructions successfully!"<<"\n";
  std::cout<<"reconstructions.size():"<<reconstructions.size()<<"\n";
  
  
  // Export un-transformed partial reconstructions for debugging.

  std::ofstream recons_file(JoinPaths(clusters_path,"reconstructions.txt"), std::ios::trunc);
  for (size_t i = 0; i < reconstructions.size(); ++i) {
      const std::string reconstruction_path = JoinPaths(clusters_path, 
                                                        std::to_string(reconstructions_clusterid[i])+"_partition_" + std::to_string(i));
      std::string name = std::to_string(reconstructions_clusterid[i])+"_partition_" + std::to_string(i);
      clusters_name.insert(std::make_pair(reconstructions[i],name));
      CreateDirIfNotExists(reconstruction_path);
      reconstructions[i]->Write(reconstruction_path);
      reconstructions[i]->WriteText(reconstruction_path);
      recons_file<<name<<"("<<reconstructions[i]->NumRegImages()<<")"<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
      std::set<std::string> image_names;
      std::map<image_t,std::string> reg_image_id;
      for(image_t image_id: reconstructions[i]->RegImageIds())
      {
        if(reg_image_id.count(image_id)==0)
        {
          reg_image_id.emplace(std::make_pair(image_id,image_id_to_name.at(image_id)));
        }
      }

      for(auto& item:reg_image_id)
      {
          recons_file<<"("<<item.first<<")"<<item.second<<" ";
      }
      recons_file<<"\n\n";
      
  }

  std::cout<<"clusters_name.size()"<<clusters_name.size()<<"\n";
}


void HierarchicalMapperController::ClusterSequentialScene(std::unordered_map<image_t, std::string>& image_id_to_name,
                                                          std::vector<const SceneClustering::Cluster*>& leaf_clusters)
{
  SceneClustering scene_clustering(clustering_options_);

  

  
  Database database(options_.database_path);

  std::cout << "Reading images..." << std::endl;
  auto images = database.ReadAllImages();
  for (const auto& image : images) {
    image_id_to_name.emplace(image.ImageId(), image.Name());
  }
  

  std::cout << "Partitioning sequential scene..." << std::endl;
  auto cmp =[](const Image& image1,const Image& image2){
    return (image1.Name().compare(image2.Name())) < 0;
  };

  std::sort(images.begin(),images.end(),cmp);
  

  size_t image_overlap = static_cast<size_t>(clustering_options_.image_overlap);
  size_t half_overlap = image_overlap/2;
  size_t num_images_per_cluster = static_cast<size_t>(clustering_options_.leaf_max_num_images) - image_overlap;
  size_t images_size = images.size();
  CHECK_GT(num_images_per_cluster,image_overlap);
  size_t num_clusters = images_size / num_images_per_cluster + 
                        (images_size % num_images_per_cluster==0?0:1);
  
  leaf_clusters.resize(num_clusters);
  size_t total_num_images = 0;


  for(size_t i = 0 ; i < num_clusters ;i++){
    SceneClustering::Cluster* cluster = new SceneClustering::Cluster();
    size_t cluster_start = (i==0?0:
                                (i==(num_clusters-1)?i * num_images_per_cluster - image_overlap
                                                    :i * num_images_per_cluster - half_overlap));
    size_t cluster_end = (i==(num_clusters-1)?(images_size):cluster_start+num_images_per_cluster+half_overlap);
    
    for(size_t j = cluster_start;j<cluster_end; j++)
    {
      cluster->image_ids.push_back(images[j].ImageId());

    }
    
    total_num_images += cluster->image_ids.size();
    std::cout << StringPrintf("  Cluster %d with %d images", i + 1,
                              cluster->image_ids.size())
              << std::endl;
    leaf_clusters[i] = cluster;

  }

  std::cout << StringPrintf("Clusters have %d images", total_num_images)
            << std::endl;
}

void HierarchicalMapperController::WriteNormazliedGrpah(std::vector<std::pair<image_t, image_t>> image_pairs,std::vector<int> num_inliers)
{
  //const std::string edges_path = JoinPaths(options_.output_path, "clusters");
  std::map<int,std::vector<std::pair<image_t, image_t>>> connected_edges;
  std::ofstream file(JoinPaths(options_.output_path,"NormazliedGrpah.txt"), std::ios::trunc);
  //image_t interval[]={172,308,909,993,1,29,1085,1135};
  std::vector<image_t> interval_first{1,53, 185,188, 374,474, 921,923, 1027,1135};
  std::vector<image_t> interval_second{ 1,3000};
  image_t max_image=1;
  
  for(size_t i = 0 ; i < image_pairs.size();i++)
  {
    max_image=std::max(std::max(image_pairs[i].first,image_pairs[i].second),max_image);
    connected_edges[num_inliers[i]].push_back(image_pairs[i]);
  }
  
  for(image_t w = 1;w<=max_image;w++)   //
  //for(size_t k=0;k<interval_first.size();k+=2)
  {
    image_t k = 0;  //
    interval_first[k]=w,interval_first[k+1]=w; //
    file<<"\n("<<interval_first[k]<<"-"<<interval_first[k+1]<<")\n";
   
    size_t num_obver = 0;
    int val_obver = 0;
  for(auto num_vector:connected_edges)
  {
    
    for(size_t i=0;i<num_vector.second.size();i++)
    {
      image_t id1 = num_vector.second[i].first;
      image_t id2 =num_vector.second[i].second;
      bool flag = false;
      if(id2>=interval_first[k]&&id2<=interval_first[k+1]) std::swap(id1,id2);
      if(id1>=interval_first[k]&&id1<=interval_first[k+1])
      {
        
        for(size_t j=0;j<interval_second.size();j+=2)
        {
          if(interval_second[j]==interval_first[k]&&interval_second[j+1]==interval_first[k+1]) continue;
          if(id2>=interval_second[j]&&id2<=interval_second[j+1])
          {
              file<<id1<<"-"<<id2<<"(in cluster):"<<num_vector.first<<"\n";
              num_obver ++;
              val_obver+=num_vector.first;
              flag = true;
              break;
          }
        }
        if(!flag&&id2>=interval_first[k]&&id2<=interval_first[k+1])
        {
          //file<<id1<<"-"<<id2<<":"<<num_vector.first<<"***\n";
          
        }
        else if(!flag)
        {
          file<<id1<<"-"<<id2<<":"<<num_vector.first<<"\n";
          num_obver++;
          val_obver+=num_vector.first;
        }
      }
    }
  }
  file<<"edges num :"<<num_obver<<"\n"
      <<"edge tot val:"<<val_obver<<"\n\n";
      
  }
  file.close();
}

void HierarchicalMapperController::WriteIntraCluster(std::vector<std::pair<image_t, image_t>> image_pairs,std::vector<int> num_inliers)
{
  //const std::string edges_path = JoinPaths(options_.output_path, "clusters");
  std::map<int,std::vector<std::pair<image_t, image_t>>> connected_edges;
  std::ofstream file(JoinPaths(options_.output_path,"connected_edges.txt"), std::ios::trunc);
  //image_t interval[]={172,308,909,993,1,29,1085,1135};
  std::vector<image_t> interval_first{ 51,110, 265,266, 599,730, 834,861, 1028,1086};
  std::vector<image_t> interval_second{51,110, 265,266, 599,730, 834,861, 1028,1086};
  image_t max_image=1;
  
  for(size_t i = 0 ; i < image_pairs.size();i++)
  {
    max_image=std::max(std::max(image_pairs[i].first,image_pairs[i].second),max_image);
    connected_edges[num_inliers[i]].push_back(image_pairs[i]);
  }
 
  for(size_t k=0;k<interval_first.size();k+=2)
  {

    file<<"\n("<<interval_first[k]<<"-"<<interval_first[k+1]<<")\n";
   
    size_t num_obver = 0;
    int val_obver = 0;
    for(auto num_vector:connected_edges)
    {
      
      for(size_t i=0;i<num_vector.second.size();i++)
      {
        image_t id1 = num_vector.second[i].first;
        image_t id2 =num_vector.second[i].second;
        bool flag = false;
        if(id2>=interval_first[k]&&id2<=interval_first[k+1]) std::swap(id1,id2);
        if(id1>=interval_first[k]&&id1<=interval_first[k+1])
        {
          
          for(size_t j=0;j<interval_second.size();j+=2)
          {
            if(interval_second[j]==interval_first[k]&&interval_second[j+1]==interval_first[k+1]) continue;
            if(id2>=interval_second[j]&&id2<=interval_second[j+1])
            {
                file<<id1<<"-"<<id2<<"(in cluster):"<<num_vector.first<<"\n";
                num_obver ++;
                val_obver+=num_vector.first;
                flag = true;
                break;
            }
          }
          if(!flag&&id2>=interval_first[k]&&id2<=interval_first[k+1])
          {
            //file<<id1<<"-"<<id2<<":"<<num_vector.first<<"***\n";
            
          }
          else if(!flag)
          {
            file<<id1<<"-"<<id2<<":"<<num_vector.first<<"\n";
            num_obver++;
            val_obver+=num_vector.first;
          }
        }
      }
    }
  file<<"edges num :"<<num_obver<<"\n"
      <<"edge tot val:"<<val_obver<<"\n\n";
      
  }
  file.close();
}

void HierarchicalMapperController::FilterImageRegisteredInCluster(size_t cluster_begin,size_t standard_re_reconstruction,
                                                        std::vector<SceneClustering::Cluster*>& clusters,
                                                        const std::unordered_map<SceneClustering::Cluster*, ReconstructionManager>&
                                                        reconstruction_managers)
{
    size_t old_num_clusters = clusters.size();
    
    for(size_t i = cluster_begin ; i < old_num_clusters ; i++){
      {
        std::cout<<"///////////////check address in Filter////////////////\n";
        for(size_t k = 0;k<clusters.size();k++)
        {
          std::cout<<"clu "<<k<<": "<<&clusters[k]<<"\n";
        }

        for(auto& recon_item:reconstruction_managers)
        {
          std::cout<<"recon.first:"<<recon_item.first<<"\n";
        }
        std::cout<<"////////////////////////////////////////////\n";

      }
      SceneClustering::Cluster* cluster = clusters[i];
      const ReconstructionManager& reconstruction_manager = reconstruction_managers.at(cluster);
      std::vector<image_t> non_registered_images;
      std::cout<<"clusters "<<i<<" :\n";
      for(auto it=cluster->image_ids.begin();it!=cluster->image_ids.end();)
      {
        image_t image_id = *it;
        bool b_registered = false;
        //bool b_registered = true;
        //std::cout<<"reconstruction_manager.Size():"<<reconstruction_manager.Size()<<"\n";
        std::cout<<"search for image "<<*it<<"\n";
        for(size_t j= 0 ;j<reconstruction_manager.Size();j++)
        {
          const auto& reconstruction = reconstruction_manager.Get(j);
          //std::cout<<"reconstruction.NumImages():"<<reconstruction.NumImages()<<"\n";
          if(reconstruction.ExistsImage(image_id))
          {
            //std::cout<<"Image "<<image_id<<" registered !\n";
            b_registered = reconstruction.IsImageRegistered(image_id);
          }
          else
          {
            //std::cout<<"Image "<<image_id<<" missing !\n";
          }
          if(b_registered) break;
        }
        if(!b_registered)
        {
            cluster->image_ids.erase(it);
            std::cout<<"Image "<<image_id<<" is not registered !!!!\n\n";
            non_registered_images.push_back(image_id);
        }
        else
        {
          it++;
        }
        
        
      }
      std::cout<<"clusters "<<i<<" complete:\n";
      std::cout<<"non_registered_images.size():"<<non_registered_images.size()<<"\n";
      if(non_registered_images.size() >= standard_re_reconstruction)
      {
          std::cout<<"push a new cluster "<<i<<":\n";
          clusters.push_back(new SceneClustering::Cluster());
          SceneClustering::Cluster* last_cluster = *(clusters.rbegin());
          last_cluster->image_ids.insert(last_cluster->image_ids.end(),non_registered_images.begin(),non_registered_images.end());
          std::cout<<"last_cluster.image_ids.size():"<<last_cluster->image_ids.size()<<"\n";
      }

    }

    std::cout<<"FilterImageRegisteredInCluster done "<<"\n";
    
}

void HierarchicalMapperController::ReconstructionFunction(const SceneClustering::Cluster& cluster,
                                ReconstructionManager* reconstruction_manager,size_t num_threads_per_worker,
                                std::unordered_map<image_t, std::string>& image_id_to_name)
{
   
    if (cluster.image_ids.empty()) {
      return;
    }

    IncrementalMapperOptions custom_options = mapper_options_;
    //custom_options.max_model_overlap = 3;
    //custom_options.init_num_trials = options_.init_num_trials;
    custom_options.num_threads = num_threads_per_worker;

    for (const auto image_id : cluster.image_ids) {
      custom_options.image_names.insert(image_id_to_name.at(image_id));
    }

    IncrementalMapperController mapper(&custom_options, options_.image_path,
                                       options_.database_path,
                                       reconstruction_manager);
    mapper.Start();
    mapper.Wait();
  
}

void HierarchicalMapperController::FilteringAndBFSForCluterScene(std::unordered_map<image_t, std::string>& image_id_to_name,
                                                  std::vector<SceneClustering>& scene_clusterings)
{
   
  //SceneClustering scene_clustering();
  //std::vector<SceneClustering> scene_clusterings;
  

  {
    Database database(options_.database_path);

    std::cout << "Reading images..." << std::endl;
    const auto images = database.ReadAllImages();
    for (const auto& image : images) {
      image_id_to_name.emplace(image.ImageId(), image.Name());
    }

    std::cout << "Reading scene graph..." << std::endl;
    std::vector<std::pair<image_t, image_t>> image_pairs;
    std::vector<int> num_inliers;

    //////filter image_pairs by num_inliers///
    int min_inliers = 46;
    database.ReadTwoViewGeometryNumInliers(&image_pairs, &num_inliers);
    std::unordered_map<image_t,std::unordered_map<image_t,int>> image_edges;
    for(auto image_item: image_id_to_name)
    {
      image_edges.emplace(image_item.first,std::unordered_map<image_t,int>());
    }
    {
      size_t index_imagepair = 0;
      size_t num_imagepairs = image_pairs.size();
      for(size_t i=0;i<num_imagepairs;i++)
      {
        if(num_inliers[index_imagepair]<min_inliers)
        {
            image_pairs.erase(image_pairs.begin()+index_imagepair);
            num_inliers.erase(num_inliers.begin()+index_imagepair);
        }
        else
        {
          image_edges.at(image_pairs[index_imagepair].first).emplace(image_pairs[index_imagepair].second,num_inliers[index_imagepair]);
          image_edges.at(image_pairs[index_imagepair].second).emplace(image_pairs[index_imagepair].first,num_inliers[index_imagepair]);
          index_imagepair++;
        }
        
      }
    }
    //////////////////////////////

    /////////BFS for segmentation coarsely//////////
    std::cout << "BFS in scene graph..." << std::endl;
    std::unordered_set<image_t> vis_image;  //flag 
    std::vector<std::vector<image_t>> cluster_images;  //images' id in a segmentation
    std::queue<image_t> queue_images;
    std::vector<std::vector<std::pair<image_t, image_t>>> cluster_image_pairs;  //images' pair in a segmentation
    std::vector<std::vector<int>> cluster_num_inliers;  ////images' pair in a segmentation
    for(auto image_item: image_id_to_name)
    {
      if(vis_image.count(image_item.first)) continue;
      //find a new independent segmentation and malloc a new space for it 
      cluster_images.push_back(std::vector<image_t>());
      cluster_image_pairs.push_back(std::vector<std::pair<image_t, image_t>>());
      cluster_num_inliers.push_back(std::vector<int>());
      std::vector<image_t>& one_cluster_images = *(cluster_images.rbegin());
      std::vector<std::pair<image_t, image_t>>& one_cluster_edges = *(cluster_image_pairs.rbegin());
      std::vector<int>& one_cluster_inliers = *(cluster_num_inliers.rbegin());

      //begin BFS
      vis_image.insert(image_item.first);
      queue_images.push(image_item.first);

      while(!queue_images.empty())
      {
        auto image_id = queue_images.front();
        queue_images.pop();
        for(auto edge_item:image_edges.at(image_id))
        {
            if(!vis_image.count(edge_item.first))
            {
              vis_image.insert(edge_item.first);
              queue_images.push(edge_item.first);
              one_cluster_images.push_back(edge_item.first);
              //one_cluster_edges.push_back(std::make_pair(image_id,edge_item.first));
              //one_cluster_inliers.push_back(edge_item.second);
            }
            
        }
      }

    }
    size_t all_num_pairs = 0;
    size_t all_images = 0;
    for(size_t i = 0;i<cluster_images.size();i++)
    {
      std::cout<<"Coarse cluster "<<i<<":\n"
               <<"               images: "<<cluster_images[i].size()<<"\n"
               <<"               images_pairs: "<<cluster_image_pairs[i].size()<<"\n";
      all_num_pairs+=cluster_image_pairs[i].size();
      all_images += cluster_images[i].size();
    }
    std::cout<<"\nimages are be divided into "<<cluster_images.size()<<"coarsely\n";
    std::cout<<"all_num_pairs:"<<all_num_pairs<<"\n";
    std::cout<<"image_pairs.size()"<<image_pairs.size()<<"\n";
    //////////////////////

    /////////use MLKKM for Subdivision////////////
    
    std::cout << "Partitioning scene graph..." << std::endl;
    
    for(size_t i = 0;i<cluster_images.size();i++)
    {
      scene_clusterings.emplace_back(SceneClustering(clustering_options_));
      scene_clusterings[i].Partition(cluster_image_pairs[i],cluster_num_inliers[i]);
    }
    
    
    WriteNormazliedGrpah(image_pairs,num_inliers);
    WriteIntraCluster(image_pairs,num_inliers);


    //////////////////////
  }

  
  
}

bool  HierarchicalMapperController::TestwithExistingPartitions(std::vector<SceneClustering::Cluster*>& clusters,
                                                        std::unordered_map<SceneClustering::Cluster*, ReconstructionManager>&
                                                        reconstruction_managers)
{
    printf("Judge If reconstructions exist\n///////////");
    std::vector<std::string> cluster_list;
    
    std::string clusters_path_prio = JoinPaths(options_.output_path,"clusters");
    cluster_list = GetRecursiveDirList(clusters_path_prio);
    std::unordered_map<std::string,SceneClustering::Cluster*> umap_cluster_name;
    std::cout<<"clusters[0]="<<clusters[0]->image_ids.size()<<"\n"
             <<"clusters[1]="<<clusters[1]->image_ids.size()<<"\n"
             <<"clusters[2]="<<clusters[2]->image_ids.size()<<"\n"
             <<"clusters[3]="<<clusters[3]->image_ids.size()<<"\n";
    std::cout<<"0_partition_0(343)\n1_partition_1(307)\n2_partition_2(220)\n3_partition_3(270)\n";
    umap_cluster_name["0_partition_0"] = clusters[0];
    umap_cluster_name["1_partition_1"] = clusters[1];
    umap_cluster_name["2_partition_2"] = clusters[2];
    umap_cluster_name["3_partition_3"] = clusters[3];
    if(!cluster_list.empty())
    {
        printf("Extract Reconstructions existed\n///////////");
        
      
        for(size_t i = 0;i<cluster_list.size();i++)
        {
          
          std::string clustername = GetPathBaseName(cluster_list[i]);
          reconstruction_managers[umap_cluster_name.at(clustername)].Read(cluster_list[i]);
          
        }
        //DACSFM::MergeClustersintoGlobal(reconstructions_test, reconstruction_manager_,cluster_name,options_.output_path);

        std::cout << std::endl;
        GetTimer().PrintMinutes();
        return true;
    }
    return  false;
}


}  // namespace colmap
