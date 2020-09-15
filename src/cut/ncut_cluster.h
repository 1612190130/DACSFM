#ifndef SRC_CUT_NCUT_CLUSTER_H_
#define SRC_CUT_NCUT_CLUSTER_H_

#include <vector>

#include "cut/cluster.h"

namespace DACSFM {

class NCutCluster : public Cluster
{
public:
    virtual std::unordered_map<int, int> ComputeCluster(
        const std::vector<std::pair<int, int>>& edges,
        const std::vector<int>& weights,
        const int num_partitions) override;
};

} 

#endif