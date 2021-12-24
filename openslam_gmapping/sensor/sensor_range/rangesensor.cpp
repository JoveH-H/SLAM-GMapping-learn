#include "gmapping/sensor/sensor_range/rangesensor.h"

namespace GMapping
{
    /**
     * @brief 激光传感器名称定义函数
     */
    RangeSensor::RangeSensor(std::string name) : Sensor(name) {}

    /**
     * @brief 激光传感器数据设置函数
     * 
     * @param name 传感器名称
     * @param beams_num 光束数量
     * @param res 角增量
     * @param position 位姿
     * @param span 宽度 0表示线状波束
     * @param maxrange 最大范围
     */
    RangeSensor::RangeSensor(std::string name, unsigned int beams_num, double res, const OrientedPoint &position, double span, double maxrange) : Sensor(name),
                                                                                                                                                  m_pose(position), m_beams(beams_num)
    {
        /* 计算起始角度 */
        double angle = -.5 * res * beams_num;
        for (unsigned int i = 0; i < beams_num; i++, angle += res)
        {
            RangeSensor::Beam &beam(m_beams[i]);
            beam.span = span;
            beam.pose.x = 0;
            beam.pose.y = 0;
            beam.pose.theta = angle;
            beam.maxRange = maxrange;
        }

        /* 新格式标志位置0并开始更新扫描光束正余弦值属性 */
        newFormat = 0;
        updateBeamsLookup();
    }

    /**
     * @brief 扫描光束属性更新函数
     * 
     * 根据光束的位置和方位角计算了正余弦值
     */
    void RangeSensor::updateBeamsLookup()
    {
        for (unsigned int i = 0; i < m_beams.size(); i++)
        {
            RangeSensor::Beam &beam(m_beams[i]);
            beam.s = sin(m_beams[i].pose.theta);
            beam.c = cos(m_beams[i].pose.theta);
        }
    }

};
