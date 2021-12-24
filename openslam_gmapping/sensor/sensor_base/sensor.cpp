#include "gmapping/sensor/sensor_base/sensor.h"

namespace GMapping
{
    /**
     * @brief 传感器名称定义函数
     * 
     * @param name 传感器名称
     */
    Sensor::Sensor(const std::string &name)
    {
        m_name = name;
    }

    Sensor::~Sensor()
    {
    }

}; // end namespace
