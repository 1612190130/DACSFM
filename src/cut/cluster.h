#ifndef SRC_CUT_CLUSTER_H
#define SRC_CUT_CLUSTER_H

#include <unordered_map>
#include <cstdio>
#include <glog/logging.h>

#include "util/types.h"

using namespace colmap;

namespace DACSFM{
    enum ClusterType {NCUT};

    class Cluster
    {
        protected:
        std::unordered_map<int,int> labels_;

        public:
        Cluster(){}
        ~Cluster(){}

        virtual std::unordered_map<int,int> ComputeCluster(
            const std::vector<std::pair<int,int>>& edges,
            const std::vector<int>& weights,
            const int num_partitions)=0;
        
    };
}


#endif