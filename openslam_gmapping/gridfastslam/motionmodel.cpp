#include "gmapping/gridfastslam/motionmodel.h"
#include <gmapping/utils/stat.h>
#include <iostream>

#define MotionModelConditioningLinearCovariance 0.01   /* 运动模型线性协方差 */
#define MotionModelConditioningAngularCovariance 0.001 /* 运动模型角度协方差 */

namespace GMapping
{
    /**
     * @brief 运动估算函数
     *
     * @param p 上一时刻运动估算位姿
     * @param linearMove 线性距离
     * @param angularMove 角度距离
     *
     * @return n 新位姿
     */
    OrientedPoint MotionModel::drawFromMotion(const OrientedPoint &p, double linearMove, double angularMove) const
    {
        OrientedPoint n(p);
        double lm = linearMove + fabs(linearMove) * sampleGaussian(srr) + fabs(angularMove) * sampleGaussian(str);
        double am = angularMove + fabs(linearMove) * sampleGaussian(srt) + fabs(angularMove) * sampleGaussian(stt);
        n.x += lm * cos(n.theta + .5 * am);
        n.y += lm * sin(n.theta + .5 * am);
        n.theta += am;
        n.theta = atan2(sin(n.theta), cos(n.theta));
        return n;
    }

    /**
     * @brief 运动估算函数
     *
     * @param p 上一时刻粒子群中的一个运动估算位姿
     * @param pnew 最新里程计坐标系下激光雷达传感器位姿（已依靠理论里程更新了位姿）
     * @param pold 上一时刻的系统位姿
     *
     * @return absoluteSum(p, noisypoint) 新位姿
     */
    OrientedPoint MotionModel::drawFromMotion(const OrientedPoint &p, const OrientedPoint &pnew, const OrientedPoint &pold) const
    {
        double sxy = 0.3 * srr;

        /* 通过最新理想里程计位姿与上一时刻的系统位姿做差得到控制量 */
        OrientedPoint delta = absoluteDifference(pnew, pold);

        /* 构建noisypoint，并通过辅助函数sampleGaussian为控制量添加上噪声项 */
        OrientedPoint noisypoint(delta);
        noisypoint.x += sampleGaussian(srr * fabs(delta.x) + str * fabs(delta.theta) + sxy * fabs(delta.y));
        noisypoint.y += sampleGaussian(srr * fabs(delta.y) + str * fabs(delta.theta) + sxy * fabs(delta.x));
        noisypoint.theta += sampleGaussian(stt * fabs(delta.theta) + srt * sqrt(delta.x * delta.x + delta.y * delta.y));
        noisypoint.theta = fmod(noisypoint.theta, 2 * M_PI);

        /* 角度转换 */
        if (noisypoint.theta > M_PI)
            noisypoint.theta -= 2 * M_PI;

        /* 将附加了噪声项的控制量添加到粒子的位姿向量上 */
        return absoluteSum(p, noisypoint);
    }

    /*
    OrientedPoint
    MotionModel::drawFromMotion(const OrientedPoint& p, const OrientedPoint& pnew, const OrientedPoint& pold) const{

        //compute the three stps needed for perfectly matching the two poses if the noise is absent

        OrientedPoint delta=pnew-pold;
        double aoffset=atan2(delta.y, delta.x);
        double alpha1=aoffset-pold.theta;
        alpha1=atan2(sin(alpha1), cos(alpha1));
        double rho=sqrt(delta*delta);
        double alpha2=pnew.theta-aoffset;
        alpha2=atan2(sin(alpha2), cos(alpha2));

        OrientedPoint pret=drawFromMotion(p, 0, alpha1);
        pret=drawFromMotion(pret, rho, 0);
        pret=drawFromMotion(pret, 0, alpha2);
        return pret;
    }
    */

    /**
     * @brief 高斯近似函数
     *
     * @param pnew 最新里程计坐标系下激光雷达传感器位姿（已依靠理论里程更新了位姿）
     * @param pold 上一时刻的系统位姿
     *
     * @return cov 协方差
     */
    Covariance3 MotionModel::gaussianApproximation(const OrientedPoint &pnew, const OrientedPoint &pold) const
    {
        OrientedPoint delta = absoluteDifference(pnew, pold);            /* 计算运动增量 */
        double linearMove = sqrt(delta.x * delta.x + delta.y * delta.y); /* 线性的增量 */
        double angularMove = fabs(delta.x);                              /* 角度的增量 */
        double s11 = srr * srr * linearMove * linearMove;
        double s22 = stt * stt * angularMove * angularMove;
        double s12 = str * angularMove * srt * linearMove;
        Covariance3 cov; /* 计算协方差 */
        double s = sin(pold.theta), c = cos(pold.theta);
        cov.xx = c * c * s11 + MotionModelConditioningLinearCovariance;
        cov.yy = s * s * s11 + MotionModelConditioningLinearCovariance;
        cov.tt = s22 + MotionModelConditioningAngularCovariance;
        cov.xy = s * c * s11;
        cov.xt = c * s12;
        cov.yt = s * s12;
        return cov;
    }

};
