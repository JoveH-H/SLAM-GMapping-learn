#include "gmapping/sensor/sensor_odometry/odometrysensor.h"

namespace GMapping
{
    /**
     * @brief 里程计传感器定义函数
     * 
     * @param name 里程计传感器名称
     * @param ideal 标记里程计是否为理想情况下的传感器
     */
    OdometrySensor::OdometrySensor(const std::string &name, bool ideal) : Sensor(name) { m_ideal = ideal; }

};
