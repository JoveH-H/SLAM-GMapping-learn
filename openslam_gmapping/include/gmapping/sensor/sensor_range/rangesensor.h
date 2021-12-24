#ifndef RANGESENSOR_H
#define RANGESENSOR_H

#include <vector>
#include <gmapping/sensor/sensor_base/sensor.h>
#include <gmapping/utils/point.h>
#include <gmapping/sensor/sensor_range/sensor_range_export.h>

namespace GMapping
{
    /**
     * @brief RangeSensor 激光传感器类
     * 
     * 激光传感器的封装
     */
    class SENSOR_RANGE_EXPORT RangeSensor : public Sensor
    {
        friend class Configuration;
        friend class CarmenConfiguration;
        friend class CarmenWrapper;

    public:
        /* 定义结构体Beam，描述了扫描光束的物理特性 */
        struct Beam
        {
            /* 相对于传感器坐标系（传感器中心）的位姿 */
            OrientedPoint pose; //pose relative to the center of the sensor
            /* 宽度 0表示线状波束 */
            double span; //spam=0 indicates a line-like beam
            /* 最大范围 */
            double maxRange; //maximum range of the sensor
            /* 光束角的正余弦值 */
            double s, c; //sinus and cosinus of the beam (optimization);
        };
        /* 激光传感器名称定义函数 */
        RangeSensor(std::string name);

        /* 激光传感器数据设置函数 */
        RangeSensor(std::string name, unsigned int beams, double res, const OrientedPoint &position = OrientedPoint(0, 0, 0), double span = 0, double maxrange = 89.0);
        inline const std::vector<Beam> &beams() const { return m_beams; }
        inline std::vector<Beam> &beams() { return m_beams; }
        inline OrientedPoint getPose() const { return m_pose; }
        void updateBeamsLookup(); /* 扫描光束属性更新函数 */
        bool newFormat;

    protected:
        OrientedPoint m_pose;      /* 记录传感器的位姿 */
        std::vector<Beam> m_beams; /* 记录扫描光束的信息 */
    };

};

#endif
