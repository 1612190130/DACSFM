#include "merge/similaritytransform.h"


namespace  DACSFM{

SimilarityTransform::SimilarityTransform(Eigen::Matrix3x4d srt_matrix)
{
    scale = srt_matrix.block<1, 3>(0, 0).norm();
    Qvec = NormalizeQuaternion(RotationMatrixToQuaternion(srt_matrix.block<3, 3>(0, 0) / scale));
    tvec = srt_matrix.block<3, 1>(0, 3);
}

SimilarityTransform::SimilarityTransform(SimilarityTransform3 similaritytransform3)
{   
    scale = similaritytransform3.Scale();
    Qvec = NormalizeQuaternion(similaritytransform3.Rotation());
    tvec = similaritytransform3.Translation();
}

void SimilarityTransform::Update(Eigen::Matrix3x4d srt_matrix)
{
    scale = srt_matrix.block<1, 3>(0, 0).norm();
    Qvec = NormalizeQuaternion(RotationMatrixToQuaternion(srt_matrix.block<3, 3>(0, 0) / scale));
    tvec = srt_matrix.block<3, 1>(0, 3);
}

void SimilarityTransform::Update(SimilarityTransform3 similaritytransform3)
{   
    scale = similaritytransform3.Scale();
    Qvec = NormalizeQuaternion(similaritytransform3.Rotation());
    tvec = similaritytransform3.Translation();
}

Eigen::Matrix4d SimilarityTransform::Inverse() const
{
    double scale_inv = 1.0/scale;
    Eigen::Vector4d Qvec_inv{Qvec(0),-Qvec(1),-Qvec(2),-Qvec(3)};
    Eigen::Vector3d tvec_inv = QuaternionRotatePoint(Qvec_inv,tvec) * (-scale_inv);

    Eigen::Matrix4d matrix = Eigen::MatrixXd::Identity(4, 4);
    matrix.topLeftCorner<3, 4>() = ComposeProjectionMatrix(Qvec_inv, tvec_inv);
    matrix.block<3, 3>(0, 0) *= scale_inv;
    return matrix;
}

void SimilarityTransform::TransformPose(Eigen::Vector4d* Pose_qvec,Eigen::Vector3d* Pose_tvec) const {
  // Projection matrix P1 projects 3D object points to image plane and thus to
  // 2D image points in the source coordinate system:
  //    x' = P1 * X1
  // 3D object points can be transformed to the destination system by applying
  // the similarity transformation S:
  //    X2 = S * X1
  // To obtain the projection matrix P2 that transforms the object point in the
  // destination system to the 2D image points, which do not change:
  //    x' = P2 * X2 = P2 * S * X1 = P1 * S^-1 * S * X1 = P1 * I * X1
  // and thus:
  //    P2' = P1 * S^-1
  // Finally, undo the inverse scaling of the rotation matrix:
  //    P2 = s * P2'

  
  

  Eigen::Matrix4d src_matrix = Eigen::MatrixXd::Identity(4, 4);
  src_matrix.topLeftCorner<3, 4>() = ComposeProjectionMatrix(*Pose_qvec, *Pose_tvec);
  Eigen::Matrix4d dst_matrix =
      src_matrix.matrix() * Inverse();
  //dst_matrix *= Scale();

  *Pose_qvec = RotationMatrixToQuaternion(dst_matrix.block<3, 3>(0, 0));
  *Pose_tvec = dst_matrix.block<3, 1>(0, 3);
}

void SimilarityTransform::TransformPoint(Eigen::Vector3d* xyz) const
{
    *xyz = QuaternionRotatePoint(Qvec,*xyz) * scale + tvec ;
}


}