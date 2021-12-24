#include "gmapping/sensor/sensor_odometry/odometryreading.h"

namespace GMapping
{
    /**
     * @brief 里程计数据定义函数
     * 
     * @param odo 里程计传感器类
     * @param time 数据产生时间
     */
    OdometryReading::OdometryReading(const OdometrySensor *odo, double time) : SensorReading(odo, time) {}

};
