#ifndef DACSFM_SRC_MERGE_COSTFUNCTION_H_
#define DACSFM_SRC_MERGE_COSTFUNCTION_H_


#include <Eigen/Core>

#include <ceres/ceres.h>
#include <ceres/rotation.h>

using namespace colmap;


namespace DACSFM
{
template <typename CameraModel>
class  BundleAdjustmentConstantPoseCostFunction
{
public:
    BundleAdjustmentConstantPoseCostFunction(const Eigen::Vector4d& Qvec,
                                 const Eigen::Vector3d& tvec,const Eigen::Vector2d point2D)
    
    :qw_(Qvec(0)),qx_(Qvec(1)),qy_(Qvec(2)),qz_(Qvec(3)),tx_(tvec(0)),ty_(tvec(1)),tz_(tvec(2)),
    observed_x_(point2D(0)),observed_y_(point2D(1)){}

    

    template<typename T>
    bool operator() (const T* point3D,const T* camera_params,T* residuals)const
    {
        T projection[3];
        const T qvec[4]={T(qw_),T(qx_),T(qy_),T(qz_)};
        ceres::UnitQuaternionRotatePoint(qvec,point3D,projection);

        projection[0]+=T(tx_);
        projection[1]+=T(ty_);
        projection[2]+=T(tz_);

        projection[0]/=projection[2];
        projection[1]/=projection[2];

        CameraModel::WorldToImage(camera_params
        ,projection[0],projection[1],&residuals[0],&residuals[1]);

        residuals[0]-=T(observed_x_);
        residuals[1]-=T(observed_y_);
        return true;
        
    }

    static ceres::CostFunction* Create(const Eigen::Vector4d& Qvec,
                                 const Eigen::Vector3d& tvec,const Eigen::Vector2d point2D){

            return (new ceres::AutoDiffCostFunction<BundleAdjustmentConstantPoseCostFunction<CameraModel>,2,3,CameraModel::kNumParams>(
            new BundleAdjustmentConstantPoseCostFunction(Qvec,tvec,point2D)));
    }


private:
    const double qw_;
    const double qx_;
    const double qy_;
    const double qz_;
    const double tx_;
    const double ty_;
    const double tz_;

    const double observed_x_;
    const double observed_y_;

};




template <typename CameraModel>
class  BundleAdjustmentCostFunction
{
public:
    BundleAdjustmentCostFunction(const Eigen::Vector2d point2D):
        observed_x_(point2D(0)),observed_y_(point2D(1)){}
    

    template<typename T>
    bool operator() (const T* Qvec,const T* tvec,const T* point3D,const T* camera_params,T* residuals)const
    {
        T projection[3];
        
        ceres::UnitQuaternionRotatePoint(Qvec,point3D,projection);

        projection[0]+=tvec[0];
        projection[1]+=tvec[1];
        projection[2]+=tvec[2];

        projection[0]/=projection[2];
        projection[1]/=projection[2];

        CameraModel::WorldToImage(camera_params
        ,projection[0],projection[1],&residuals[0],&residuals[1]);

        residuals[0]-=T(observed_x_);
        residuals[1]-=T(observed_y_);
        return true;
        
    }

    static ceres::CostFunction* Create(const Eigen::Vector2d point2D){

            return (new ceres::AutoDiffCostFunction<BundleAdjustmentCostFunction<CameraModel>,2,4,3,3,CameraModel::kNumParams>(
            new BundleAdjustmentCostFunction(point2D)));
    }


private:
    const double observed_x_;
    const double observed_y_;

};

template <typename CameraModel>
class SRTBundleAdjustmentCostFunction
{
public:
    SRTBundleAdjustmentCostFunction(const Eigen::Vector2d point2D):
        observed_x_(point2D(0)),observed_y_(point2D(1)){}

    template<typename T>
    bool operator() (const T* srt_scale,const T* srt_Qvec,const T* srt_tvec,
                     const T* Qvec,const T* tvec,const T* point3D,const T* camera_params,T* residual)const
    {
        T simliaritytransform[3];
        
        ceres::UnitQuaternionRotatePoint(srt_Qvec,point3D,simliaritytransform);

        simliaritytransform[0]*=srt_scale[0];
        simliaritytransform[1]*=srt_scale[0];
        simliaritytransform[2]*=srt_scale[0];

        simliaritytransform[0]+=srt_tvec[0];
        simliaritytransform[1]+=srt_tvec[1];
        simliaritytransform[2]+=srt_tvec[2];

        T projection[3];

        ceres::UnitQuaternionRotatePoint(Qvec,simliaritytransform,projection);

        projection[0]+=tvec[0];
        projection[1]+=tvec[1];
        projection[2]+=tvec[2];

        projection[0]/=projection[2];
        projection[1]/=projection[2];

        CameraModel::WorldToImage(camera_params,
            projection[0],projection[1],&residual[0],&residual[1]);
        
        residual[0]-=T(observed_x_);
        residual[1]-=T(observed_y_);
        return true;

    }

    static ceres::CostFunction* Create(const Eigen::Vector2d point2D){
        return (new ceres::AutoDiffCostFunction<SRTBundleAdjustmentCostFunction<CameraModel>,2,1,4,3,4,3,3,CameraModel::kNumParams>
                (new SRTBundleAdjustmentCostFunction(point2D)));
    }

private:
    const double observed_x_;
    const double observed_y_;
};

template <typename CameraModel>
class DoubleSrtBundleAdjustmentCostFunction
{
public:
    DoubleSrtBundleAdjustmentCostFunction(const Eigen::Vector2d point2D):
        observed_x_(point2D(0)),observed_y_(point2D(1)){}

    template<typename T>
    bool operator() (const T* srt1_scale,const T* srt1_Qvec,const T* srt1_tvec,
                     const T* srt2_scale,const T* srt2_Qvec,const T* srt2_tvec,
                     const T* Qvec,const T* tvec,const T* point3D,const T* camera_params,T* residual)const
    {
        T simliaritytransform1[3];
        
        ceres::UnitQuaternionRotatePoint(srt1_Qvec,point3D,simliaritytransform1);

        simliaritytransform1[0]*=srt1_scale[0];
        simliaritytransform1[1]*=srt1_scale[0];
        simliaritytransform1[2]*=srt1_scale[0];

        simliaritytransform1[0]+=srt1_tvec[0];
        simliaritytransform1[1]+=srt1_tvec[1];
        simliaritytransform1[2]+=srt1_tvec[2];


        T simliaritytransform2[3];
        //the inverse of srt2:
        T srt2_scale_inv=1.0/srt2_scale[0];
        T srt2_Qvec_inv[4]={srt2_Qvec[0],-srt2_Qvec[1],-srt2_Qvec[2],-srt2_Qvec[3]};
        
        //the inverse of translation srt2_tvec_inv
        T srt2_tvec_inv[3];
        ceres::UnitQuaternionRotatePoint(srt2_Qvec_inv,srt2_tvec,srt2_tvec_inv);
        srt2_tvec_inv[0]*=(-srt2_scale_inv);
        srt2_tvec_inv[1]*=(-srt2_scale_inv);
        srt2_tvec_inv[2]*=(-srt2_scale_inv);
        
        //multiply srt2

        ceres::UnitQuaternionRotatePoint(srt2_Qvec_inv,simliaritytransform1,simliaritytransform2);

        simliaritytransform2[0]*=srt2_scale_inv;
        simliaritytransform2[1]*=srt2_scale_inv;
        simliaritytransform2[2]*=srt2_scale_inv;

        simliaritytransform2[0]+=srt2_tvec_inv[0];
        simliaritytransform2[1]+=srt2_tvec_inv[1];
        simliaritytransform2[2]+=srt2_tvec_inv[2];

        T projection[3];

        ceres::UnitQuaternionRotatePoint(Qvec,simliaritytransform2,projection);

        projection[0]+=tvec[0];
        projection[1]+=tvec[1];
        projection[2]+=tvec[2];

        projection[0]/=projection[2];
        projection[1]/=projection[2];

        CameraModel::WorldToImage(camera_params,
            projection[0],projection[1],&residual[0],&residual[1]);
        
        residual[0]-=T(observed_x_);
        residual[1]-=T(observed_y_);
        return true;

    }

    static ceres::CostFunction* Create(const Eigen::Vector2d point2D){
        return (new ceres::AutoDiffCostFunction<DoubleSrtBundleAdjustmentCostFunction<CameraModel>
                ,2,1,4,3,1,4,3,4,3,3,CameraModel::kNumParams>
                (new DoubleSrtBundleAdjustmentCostFunction(point2D)));
    }

private:
    const double observed_x_;
    const double observed_y_;
};
}




#endif