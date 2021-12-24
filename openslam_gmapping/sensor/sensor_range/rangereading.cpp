#include <limits>
#include <iostream>
#include <assert.h>
#include <sys/types.h>
#include <gmapping/utils/gvalues.h>
#include "gmapping/sensor/sensor_range/rangereading.h"

namespace GMapping
{

    using namespace std;
    /**
     * @brief 激光传感器数据定义函数
     *
     * @param rs 传感器对象
     * @param time 产生数据的时间
     */
    RangeReading::RangeReading(const RangeSensor *rs, double time) : SensorReading(rs, time) {}

    /**
     * @brief 激光扫描数据数据定义函数
     *
     * @param n_beams 光束数量
     * @param d 激光扫描数据
     * @param rs 传感器对象
     * @param time 产生数据的时间
     */
    RangeReading::RangeReading(unsigned int n_beams, const double *d, const RangeSensor *rs, double time) : SensorReading(rs, time)
    {
        assert(n_beams == rs->beams().size());
        resize(n_beams);
        for (unsigned int i = 0; i < size(); i++)
            (*this)[i] = d[i];
    }

    RangeReading::~RangeReading()
    {
        //    cerr << __func__ << ": CAZZZZZZZZZZZZZZZZZZZZOOOOOOOOOOO" << endl;
    }

    /**
     * @brief 原始视图函数
     *
     * 将激光扫描数据拷贝到一个数组
     *
     * @param v 保存数据的数组
     * @param density 过滤参数 density为0.0，过滤参数不起作用
     * @return static_cast<unsigned int>(size()) 数组大小
     */
    unsigned int RangeReading::rawView(double *v, double density) const
    {
        if (density == 0)
        {
            for (unsigned int i = 0; i < size(); i++)
                v[i] = (*this)[i];
        }
        else
        {
            /* 如果临近的几个激光扫描数据所探测的点距离小于该参数，就赋予一个很大的值。
           此时的扫描束被称为是抑制的(suppressed)。那些非抑制的扫描束称为激活的(actived)。 */
            Point lastPoint(0, 0);
            uint suppressed = 0;
            for (unsigned int i = 0; i < size(); i++)
            {
                const RangeSensor *rs = dynamic_cast<const RangeSensor *>(getSensor());
                assert(rs);
                Point lp(
                    cos(rs->beams()[i].pose.theta) * (*this)[i],
                    sin(rs->beams()[i].pose.theta) * (*this)[i]);
                Point dp = lastPoint - lp;
                double distance = sqrt(dp * dp);
                if (distance < density)
                {
                    //                v[i]=MAXDOUBLE;
                    v[i] = std::numeric_limits<double>::max();
                    suppressed++;
                }
                else
                {
                    lastPoint = lp;
                    v[i] = (*this)[i];
                }
                //std::cerr<< __func__ << std::endl;
                //std::cerr<< "suppressed " << suppressed <<"/"<<size() << std::endl;
            }
        }
        //    return size();
        return static_cast<unsigned int>(size());
    };

    /**
     * @brief 计算激活状态的扫描束函数
     * 
     * @param density 过滤参数
     * @return ab 激活的/非抑制的扫描束数量
     */
    unsigned int RangeReading::activeBeams(double density) const
    {
        if (density == 0.)
            return size();
        int ab = 0;
        Point lastPoint(0, 0);
        uint suppressed = 0;
        for (unsigned int i = 0; i < size(); i++)
        {
            const RangeSensor *rs = dynamic_cast<const RangeSensor *>(getSensor());
            assert(rs);
            Point lp(
                cos(rs->beams()[i].pose.theta) * (*this)[i],
                sin(rs->beams()[i].pose.theta) * (*this)[i]);
            Point dp = lastPoint - lp;
            double distance = sqrt(dp * dp);
            if (distance < density)
            {
                suppressed++;
            }
            else
            {
                lastPoint = lp;
                ab++;
            }
            //std::cerr<< __func__ << std::endl;
            //std::cerr<< "suppressed " << suppressed <<"/"<<size() << std::endl;
        }
        return ab;
    }

    /**
     * @brief 激光扫描数据转笛卡尔坐标函数
     * 
     * 极坐标 -> 笛卡尔坐标
     * 激光扫描数据都是用极坐标的形式描述的，使用到圆点距离以及相对于极轴的偏转角来确定空间中的一个点。
     * 
     * @param maxRange 限定测量距离的最大值
     * @return cartesianPoints 笛卡尔坐标点
     */
    std::vector<Point> RangeReading::cartesianForm(double maxRange) const
    {
        const RangeSensor *rangeSensor = dynamic_cast<const RangeSensor *>(getSensor());
        assert(rangeSensor && rangeSensor->beams().size());
        //    uint m_beams=rangeSensor->beams().size();
        uint m_beams = static_cast<unsigned int>(rangeSensor->beams().size());
        std::vector<Point> cartesianPoints(m_beams);
        double px, py, ps, pc;
        px = rangeSensor->getPose().x;
        py = rangeSensor->getPose().y;
        ps = sin(rangeSensor->getPose().theta);
        pc = cos(rangeSensor->getPose().theta);
        for (unsigned int i = 0; i < m_beams; i++)
        {
            const double &rho = (*this)[i];
            const double &s = rangeSensor->beams()[i].s;
            const double &c = rangeSensor->beams()[i].c;
            if (rho >= maxRange)
            {
                cartesianPoints[i] = Point(0, 0);
            }
            else
            {
                Point p = Point(rangeSensor->beams()[i].pose.x + c * rho, rangeSensor->beams()[i].pose.y + s * rho);
                cartesianPoints[i].x = px + pc * p.x - ps * p.y;
                cartesianPoints[i].y = py + ps * p.x + pc * p.y;
            }
        }
        return cartesianPoints;
    }

};
