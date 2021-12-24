#ifndef ODOMETRYSENSOR_H
#define ODOMETRYSENSOR_H

#include <string>
#include <gmapping/sensor/sensor_base/sensor.h>
#include <gmapping/sensor/sensor_odometry/sensor_odometry_export.h>

namespace GMapping
{
    /**
     * @brief OdometrySensor 里程计传感器类
     * 
     * 描述里程计传感器名称和标记里程计是否为理想情况下的传感器
     */
    class SENSOR_ODOMETRY_EXPORT OdometrySensor : public Sensor
    {
    public:
        OdometrySensor(const std::string &name, bool ideal = false);
        inline bool isIdeal() const { return m_ideal; }

    protected:
        bool m_ideal;
    };

};

#endif
