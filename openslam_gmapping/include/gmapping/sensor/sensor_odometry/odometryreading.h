#ifndef ODOMETRYREADING_H
#define ODOMETRYREADING_H

#include <string.h>
#include <gmapping/sensor/sensor_base/sensorreading.h>
#include <gmapping/utils/point.h>
#include "gmapping/sensor/sensor_odometry/odometrysensor.h"
#include <gmapping/sensor/sensor_odometry/sensor_odometry_export.h>

namespace GMapping
{
    /**
     * @brief OdometryReading 里程计数据类
     * 
     * 用三个成员变量记录了机器人的位姿、速度和加速度，并分别提供访问和设置这些变量
     */
    class SENSOR_ODOMETRY_EXPORT OdometryReading : public SensorReading
    {
    public:
        OdometryReading(const OdometrySensor *odo, double time = 0);
        inline const OrientedPoint &getPose() const { return m_pose; }
        inline const OrientedPoint &getSpeed() const { return m_speed; }
        inline const OrientedPoint &getAcceleration() const { return m_acceleration; }
        inline void setPose(const OrientedPoint &pose) { m_pose = pose; }
        inline void setSpeed(const OrientedPoint &speed) { m_speed = speed; }
        inline void setAcceleration(const OrientedPoint &acceleration) { m_acceleration = acceleration; }

    protected:
        OrientedPoint m_pose;         /* 位姿 */
        OrientedPoint m_speed;        /* 速度 */
        OrientedPoint m_acceleration; /* 加速度 */
    };

};
#endif
