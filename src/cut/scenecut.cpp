#include "cut/scenecut.h"
#include "cut/ncut_cluster.h"

#include "util/misc.h"


#include <glog/logging.h>
#include <unordered_set>
namespace DACSFM{

bool SceneCut::Options::Check()
{
    CHECK_GT(max_num_images,0);
    CHECK_GT(image_overlap,0);
    CHECK_GT(relax_retio,1.0);
    CHECK_GT(clusters_partitioned,0);
    return true;

}

SceneCut::SceneCut(const Options& options,const ImageCluster& root_cluster):
    options_(options),root_cluster_(root_cluster){
        CHECK(options_.Check());
        
}
std::unique_ptr<Cluster> SceneCut::CreateCluster()
{
    return std::unique_ptr<Cluster>(new NCutCluster());
}


void SceneCut::Cut()
{
    timer_.Start();
    const uint num_clusters = root_cluster_.image_ids.size() / options_.max_num_images;

    CHECK_GE(num_clusters,1);

    std::vector<std::pair<int,int>> image_pairs;
    std::vector<int> weights;

    image_pairs.reserve(root_cluster_.edges.size());
    weights.reserve(root_cluster_.edges.size());

    for(auto& edge:root_cluster_.edges)
    {
        image_pairs.push_back(std::make_pair(edge.first.first,edge.first.second));
        weights.push_back(edge.second);
    }
    LOG(INFO)<<"Images Clustering Started";
    const std::unique_ptr<Cluster> cluster = this->CreateCluster();
    CHECK_NOTNULL(cluster.get());
    LOG(INFO)<<"cluster num:"<<num_clusters;

    std::unordered_map<cluster_t,cluster_t> labels=cluster->ComputeCluster(image_pairs,weights,num_clusters);

    LOG(INFO)<<"Cutting Complete,grouping images...";

    intra_clusters_.resize(num_clusters);

    for(const auto label:labels)
    {
        intra_clusters_[label.second].image_ids.push_back(label.first);
    }

    for(uint i=0;i<intra_clusters_.size();i++)
    {
        intra_clusters_[i].cluster_id= i;
    }


    LOG(INFO)<<"Grouping edges...";

    for(size_t k=0;k<image_pairs.size();k++)
    {
        const image_t i = image_pairs[k].first;
        const image_t j = image_pairs[k].second;
        const cluster_t cluster_id1 = labels[i];
        const cluster_t cluster_id2 = labels[j];

        if(cluster_id1 == cluster_id2)
        {
            intra_clusters_[cluster_id1].edges.insert(std::make_pair(std::make_pair(i,j),weights[k]));
        }
        else
        {
            const std::pair<cluster_t,cluster_t> cluster_pair = cluster_id1<cluster_id2?
                                                        (std::make_pair(cluster_id1,cluster_id2)):(std::make_pair(cluster_id2,cluster_id1));
            clusters_cut_edges_[cluster_pair].emplace_back(Edge(i,j,weights[k]));
        }
        

    }
    timer_.Pause();

    //summary

}

void SceneCut::Expand()
{
    LOG(INFO)<<"Expanding Images...";

    const uint num_clusters = intra_clusters_.size();

    inter_clusters_.reserve(num_clusters);  //resize?
    for(auto cluster:intra_clusters_)
    {
        inter_clusters_.emplace_back(cluster);
    }

    timer_.Start();

    if(num_clusters > 1)
    {
        for(auto it:clusters_cut_edges_)
        {
            const std::pair<cluster_t,cluster_t> cluster_pair = it.first;
            std::vector<Edge> cut_edges = it.second;
            ImageCluster& cluster1 = inter_clusters_[cluster_pair.first];
            ImageCluster& cluster2 = inter_clusters_[cluster_pair.second];

            AddCutEdgesToExpandClusters(cluster1,cluster2,cut_edges);
        }
    }

    timer_.Pause();
    //summary
    //analyzestatistic
}


void SceneCut::AddCutEdgesToExpandClusters(ImageCluster& cluster1,
                                           ImageCluster& cluster2,
                                           std::vector<Edge> cut_edges)
{
    if(IsSatisfyOverlapRatio(cluster1)&&IsSatisfyOverlapRatio(cluster2))
    {
        return;
    }

    const auto cmp = [](const Edge& edge1,const Edge& edge2){
        return edge1.weight > edge2.weight;
    };

    std::sort(cut_edges.begin(),cut_edges.end());

    for(int k=0 ; k< options_.image_overlap && k < cut_edges.size();k++)
    {
        const std::pair<image_t,image_t> image_pair = cut_edges[k].src < cut_edges[k].dst?
                                                      std::make_pair(cut_edges[k].src,cut_edges[k].dst):
                                                      std::make_pair(cut_edges[k].dst,cut_edges[k].src);
        int choice = cluster1.image_ids.size()>cluster2.image_ids.size()?0:1;

        if(choice)
        {
            std::unordered_set<image_t> images1(cluster1.image_ids.begin(),cluster1.image_ids.end());
            image_t addimage = images1.find(cut_edges[k].src)==images1.end()?cut_edges[k].src:cut_edges[k].dst;
            if(!IsSatisfyOverlapRatio(cluster1)&&images1.find(addimage)==images1.end())
            {
                cluster1.image_ids.push_back(addimage);
                cluster1.edges[image_pair] = cut_edges[k].weight;
            }
        }
        else
        {
            std::unordered_set<image_t> images2(cluster2.image_ids.begin(),cluster2.image_ids.end());
            image_t addimage = images2.find(cut_edges[k].src)==images2.end()?cut_edges[k].src:cut_edges[k].dst;
            if(!IsSatisfyOverlapRatio(cluster2)&&images2.find(addimage)==images2.end())
            {
                cluster2.image_ids.push_back(addimage);
                cluster2.edges[image_pair] = cut_edges[k].weight;
            }
        }
        
    }
}


bool SceneCut::IsSatisfyOverlapRatio(ImageCluster& cluster)
{
    if(cluster.is_condition_satisfied) return true;

    const int i = cluster.cluster_id;

    std::unordered_set<image_t> image_sets(cluster.image_ids.begin(),cluster.image_ids.end());

    uint overlap_images = 0;

    for(int j=0;j<inter_clusters_.size();j++)
    {
        if(i==j) continue;

        overlap_images += OverlapImagesBetweenClusters(inter_clusters_[i],inter_clusters_[j]);  //?

    }

    const float overlap_ratio = (float) overlap_images / (float)(inter_clusters_[i].image_ids.size());  //?

    if(overlap_ratio <= options_.overlap_ratio)
    {
        VLOG(2)<<"overlap ratio: "<<overlap_ratio;
        return false;
    }
    else
    {
        inter_clusters_[i].is_condition_satisfied = true;  //?
        return true;
    }
    

}

int SceneCut::OverlapImagesBetweenClusters(const ImageCluster& cluster1,const ImageCluster& cluster2)
{
    int overlap_images = 0;

    const std::unordered_set<image_t> images1(cluster1.image_ids.begin(),cluster1.image_ids.end());
    const std::unordered_set<image_t> images2(cluster2.image_ids.bgein(),cluster2.image_ids.end());
    for(auto it = images1.begin();it!=images1.end();it++)
    {
        if(images2.find(*it) != images2.end())
        {
            overlap_images++;
        }
    }

    return overlap_images;
}


std::vector<ImageCluster> SceneCut::GetInterClusters()
{
    return inter_clusters_;
}
}


