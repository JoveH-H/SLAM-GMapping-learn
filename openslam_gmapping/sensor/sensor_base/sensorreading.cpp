#include "gmapping/sensor/sensor_base/sensorreading.h"

namespace GMapping
{
    /**
     * @brief 传感器数据定义函数
     * 
     * @param s 传感器类
     * @param t 数据产生时间
     */
    SensorReading::SensorReading(const Sensor *s, double t)
    {
        m_sensor = s;
        m_time = t;
    }

    SensorReading::~SensorReading()
    {
    }

};
