#ifndef SMMAP_H
#define SMMAP_H
#include <gmapping/grid/map.h>
#include <gmapping/grid/harray2d.h>
#include <gmapping/utils/point.h>
#define SIGHT_INC 1

namespace GMapping
{

    struct PointAccumulator
    {
        /* 定义一个float类型的point模板 */
        typedef point<float> FloatPoint;
        /* before
        PointAccumulator(int i=-1): acc(0,0), n(0), visits(0){assert(i==-1);}
        */
        /*after begin*/
        PointAccumulator() : acc(0, 0), n(0), visits(0) {}
        PointAccumulator(int i) : acc(0, 0), n(0), visits(0) { assert(i == -1); }
        /*after end*/

        inline void update(bool value, const Point &p = Point(0, 0));      /* 点坐标累积更新函数 */
        inline Point mean() const { return 1. / n * Point(acc.x, acc.y); } /* 求均值 累积量除以计数总量 */
        inline operator double() const { return visits ? (double)n * SIGHT_INC / (double)visits : -1; }
        inline void add(const PointAccumulator &p)
        {
            acc = acc + p.acc;
            n += p.n;
            visits += p.visits;
        }

        static const PointAccumulator &Unknown();
        static PointAccumulator *unknown_ptr; /* 标记未知的指针 */
        FloatPoint acc;                       /* 累积点坐标 */
        int n, visits;                        /* 记录累积的次数（占用）和访问的次数（总数量） */
        inline double entropy() const;        /* 求熵 */
    };

    /**
     * @brief 点坐标累积更新函数
     *
     * @param value 是否为占用
     * @param p 点坐标
     */
    void PointAccumulator::update(bool value, const Point &p)
    {
        if (value) /* 点坐标占用 */
        {
            /* 累加占用点坐标 */
            acc.x += static_cast<float>(p.x);
            acc.y += static_cast<float>(p.y);
            /* 累加访问次数和占用计数 */
            n++;
            visits += SIGHT_INC;
        }
        else /* 点坐标空闲 */
        {
            visits++; /* 仅累加访问次数 */
        }
    }

    /**
     * @brief 求熵函数
     *
     * @return 熵
     */
    double PointAccumulator::entropy() const
    {
        if (!visits)
            return -log(.5);
        if (n == visits || n == 0)
            return 0;
        /* 二项分布的形式计算 */
        double x = (double)n * SIGHT_INC / (double)visits;

        /* x=1/2时熵最大，即一半数据认为是空闲的，一半认为是占有的，此时要付出的代价最大，可靠性最低 */
        return -(x * log(x) + (1 - x) * log(1 - x));
    }

    /* 声明定义ScanMatcherMap激光匹配地图类  原型Map openslam_gmapping\include\gmapping\grid\map.h */
    typedef Map<PointAccumulator, HierarchicalArray2D<PointAccumulator>> ScanMatcherMap;

};

#endif
