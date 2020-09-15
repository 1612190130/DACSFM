#ifndef SRC_CUT_SCENECUT_H_
#define SRC_CUT_SCENECUT_H_
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <queue>
#include <memory>
#include <iostream>
#include <utility>


#include "util/types.h"
#include "util/timer.h"
#include "cut/cluster.h"


using namespace colmap;

namespace DACSFM
{
    typedef int cluster_t;
    struct Edge
    {
        image_t src;
        image_t dst;
        float weight;

        Edge(image_t i, image_t j)
        {
            src = i;
            dst = j;
            weight = 0.0f;
        }

        Edge(image_t i, image_t j,float w)
        {
            src = i;
            dst = j;
            weight = w;
        }

    };
    struct ImageCluster
    {
        cluster_t cluster_id;

        std::vector<image_t> image_ids;

        std::unordered_map<std::pair<image_t,image_t>,int> edges;

        bool is_condition_satisfied = false;

        void ShowInfo()
        {
            LOG(INFO)<<image_ids.size()<<" nodes";
            for(auto image_id:image_ids)
            {
                std::cout<<image_id<<" ";
            }
            std::cout<<"\n";

            LOG(INFO)<<edges.size()<<"\n";
        }
    };
    class SceneCut
    {
    public:
        struct Options
        {
            // The maximum number of images in a cluster
            uint max_num_images = 100;

            uint image_overlap = 50;

            float overlap_ratio = 0.5;

            float relax_retio = 1.3;
            // The number of clusters in one cut
            int clusters_partitioned = 2;

            bool Check();
            
        };

        SceneCut(const Options& options,const ImageCluster& root_cluster);

        void Cut();
        void Expand();
        std::vector<ImageCluster> GetInterClusters();

    private:
        Options options_;
        Timer timer_;
        const ImageCluster root_cluster_;

        std::vector<ImageCluster> intra_clusters_;

        std::vector<ImageCluster> inter_clusters_;

        std::unordered_map<std::pair<cluster_t,cluster_t>,std::vector<Edge>> clusters_cut_edges_;

        void AddCutEdgesToExpandClusters(ImageCluster& cluster1,
                                           ImageCluster& cluster2,
                                           std::vector<Edge> cut_edges);
        bool IsSatisfyOverlapRatio(ImageCluster& cluster);
        int OverlapImagesBetweenClusters(const ImageCluster& cluster1,const ImageCluster& cluster2);
        std::unique_ptr<Cluster> CreateCluster();


        

    };

}



#endif