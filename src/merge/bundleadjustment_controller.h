#ifndef DACSFM_SRC_MERGE_BUNDLEADJUSTMENTCONTROLLER_H_
#define DACSFM_SRC_MERGE_BUNDLEADJUSTMENTCONTROLLER_H_


#include "base/reconstruction.h"
#include "util/option_manager.h"
#include "util/threading.h"
#include "base/pose.h"
#include "merge/similaritytransform.h"

using namespace colmap;
namespace DACSFM
{


    class MyBundleAdjustmentController: public Thread
    {
        public :
        MyBundleAdjustmentController(const OptionManager& options,
                             Reconstruction* reconstruction);
        static void Solve(const OptionManager& options,std::vector<Reconstruction*>& reconstructions,
                          std::unordered_map<Reconstruction*,SimilarityTransform>& srtforms);

        private:
        void Run();

        const OptionManager options_;
        Reconstruction* reconstruction_;
        std::unordered_map<Reconstruction*,SimilarityTransform*> srtform_to_global;
    };
}



#endif