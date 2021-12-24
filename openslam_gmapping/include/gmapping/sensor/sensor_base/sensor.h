#ifndef SENSOR_H
#define SENSOR_H

#include <string>
#include <map>
#include <gmapping/sensor/sensor_base/sensor_base_export.h>

namespace GMapping
{
    /**
     * @brief Sensor 传感器类
     * 
     * 描述传感器的名称，获取或设置传感器的名称
     */
    class SENSOR_BASE_EXPORT Sensor
    {
    public:
        Sensor(const std::string &name = "");
        virtual ~Sensor();
        inline std::string getName() const { return m_name; }
        inline void setName(const std::string &name) { m_name = name; }

    protected:
        std::string m_name;
    };

    /* 声明定义SensorMap传感器关联容器 */
    typedef std::map<std::string, Sensor *> SensorMap;

}; //end namespace

#endif
