#ifndef DACSFM_SRC_MERGE_SIMILARITYTRANSFORM_H_
#define DACSFM_SRC_MERGE_SIMILARITYTRANSFORM_H_


#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "util/alignment.h"
#include "util/types.h"
#include "base/pose.h"
#include "base/similarity_transform.h"
#include "base/projection.h"

using namespace colmap;

namespace DACSFM{
class SimilarityTransform
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    double scale;
    Eigen::Vector4d Qvec;
    Eigen::Vector3d tvec;
    SimilarityTransform()
    {
        scale = 1.0;
        Qvec = Eigen::Vector4d(1, 0, 0, 0);
        tvec = Eigen::Vector3d(0, 0, 0);
    }

    SimilarityTransform(Eigen::Matrix3x4d srt_matrix);

    SimilarityTransform(SimilarityTransform3 similaritytransform3); 

    void Update(Eigen::Matrix3x4d srt_matrix);

    void Update(SimilarityTransform3 similaritytransform3); 

    Eigen::Matrix4d Inverse() const;

    void TransformPoint(Eigen::Vector3d* xyz) const;

    void TransformPose(Eigen::Vector4d* Pose_qvec,Eigen::Vector3d* Pose_tvec) const;

    

};

}

EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION_CUSTOM(DACSFM::SimilarityTransform)

#endif