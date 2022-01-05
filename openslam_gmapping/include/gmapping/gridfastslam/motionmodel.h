#ifndef MOTIONMODEL_H
#define MOTIONMODEL_H

#include <gmapping/utils/point.h>
#include <gmapping/utils/stat.h>
#include <gmapping/utils/macro_params.h>

namespace GMapping
{
    /* 运动模型结构体 */
    struct MotionModel
    {
        /* 运动估算函数 */
        OrientedPoint drawFromMotion(const OrientedPoint &p, double linearMove, double angularMove) const;
        OrientedPoint drawFromMotion(const OrientedPoint &p, const OrientedPoint &pnew, const OrientedPoint &pold) const;

        /* 高斯近似函数 */
        Covariance3 gaussianApproximation(const OrientedPoint &pnew, const OrientedPoint &pold) const;

        /* 平移时平移里程误差，平移时旋转里程误差，旋转时平移里程误差，旋转时旋转里程误差 */
        double srr, str, srt, stt;
    };
};

#endif
