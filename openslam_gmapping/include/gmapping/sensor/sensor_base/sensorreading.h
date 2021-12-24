#ifndef SENSORREADING_H
#define SENSORREADING_H

#include "gmapping/sensor/sensor_base/sensor.h"
#include <gmapping/sensor/sensor_base/sensor_base_export.h>

namespace GMapping
{

    /**
     * @brief SensorReading 传感器数据类
     * 
     * 描述传感器的对象和记录了数据产生时间
     */
    class SENSOR_BASE_EXPORT SensorReading
    {
    public:
        SensorReading(const Sensor *s = 0, double time = 0);
        virtual ~SensorReading();
        inline double getTime() const { return m_time; }
        inline void setTime(double t) { m_time = t; }
        inline const Sensor *getSensor() const { return m_sensor; }

    protected:
        double m_time;          /* 数据产生时间 */
        const Sensor *m_sensor; /* 传感器的对象 */
    };

}; //end namespace
#endif
