#ifndef RANGEREADING_H
#define RANGEREADING_H

#include <vector>
#include <gmapping/sensor/sensor_base/sensorreading.h>
#include "gmapping/sensor/sensor_range/rangesensor.h"
#include <gmapping/sensor/sensor_range/sensor_range_export.h>

#ifdef _MSC_VER
namespace std
{
    extern template class __declspec(dllexport) vector<double>;
};
#endif

namespace GMapping
{
    /**
     * @brief RangeSensor 激光传感器数据类
     */
    class SENSOR_RANGE_EXPORT RangeReading : public SensorReading, public std::vector<double>
    {
    public:
        RangeReading(const RangeSensor *rs, double time = 0);                                        /* 激光传感器读取函数 */
        RangeReading(unsigned int n_beams, const double *d, const RangeSensor *rs, double time = 0); /* 激光扫描数据读取函数 */
        virtual ~RangeReading();
        inline const OrientedPoint &getPose() const { return m_pose; }    /* 获取位姿信息 */
        inline void setPose(const OrientedPoint &pose) { m_pose = pose; } /* 设置位姿信息 */
        unsigned int rawView(double *v, double density = 0.) const;       /* 原始视图函数 */
        std::vector<Point> cartesianForm(double maxRange = 1e6) const;    /* 激光扫描数据转笛卡尔坐标函数 */
        unsigned int activeBeams(double density = 0.) const;              /* 计算激活状态的扫描束函数 */
    protected:
        OrientedPoint m_pose; /* 读取传感器时的位姿 */
    };

};

#endif
