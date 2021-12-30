#ifndef MAP_H
#define MAP_H
#include <gmapping/utils/point.h>
#include <assert.h>
#include "gmapping/grid/accessstate.h"
#include "gmapping/grid/array2d.h"

namespace GMapping
{
    /**
    The cells have to define the special value Cell::Unknown to handle with the unallocated areas.
    The cells have to define (int) constructor;
    */
    /**
    单元格 必须定义特殊值 Cell::Unknown 来处理未分配的区域。
    单元格 必须定义(int)构造函数;
    */
    /* 声明定义DoubleArray2D双精度浮点二维数组类型 原型Array2D openslam_gmapping\include\gmapping\grid\array2d.h */
    typedef Array2D<double> DoubleArray2D;

    /**
     * @brief Map 地图模板类
     *
     * @param Cell 地图单元类型
     * @param Storage 存储结构
     * @param isClass TODO 功能未知
     */
    template <class Cell, class Storage, const bool isClass = true>
    class Map
    {
    public:
        /* Map 初始化函数 */
        Map(int mapSizeX, int mapSizeY, double delta);
        Map(const Point &center, double worldSizeX, double worldSizeY, double delta);
        Map(const Point &center, double xmin, double ymin, double xmax, double ymax, double delta);

        /* the standard implementation works filen in this case*/
        /* 在这种情况下，标准实现可以正常工作 */
        // Map(const Map& g);
        // Map& operator =(const Map& g);

        /* 栅格地图尺寸调整的相关函数 */
        void resize(double xmin, double ymin, double xmax, double ymax);
        void grow(double xmin, double ymin, double xmax, double ymax);

        /* 由于物理地图是一个连续的坐标空间，需要使用浮点数来描述，使用double类型来具现Point。
           而栅格地图的表示则是离散的，所以使用int类型来具现IntPoint */
        /* 点在物理地图与栅格地图的转换相关函数 */
        inline IntPoint world2map(const Point &p) const;
        inline Point map2world(const IntPoint &p) const;
        inline IntPoint world2map(double x, double y) const
        {
            return world2map(Point(x, y));
        }
        inline Point map2world(int x, int y) const
        {
            return map2world(IntPoint(x, y));
        }

        /* 获取基础参数的相关函数 */
        inline Point getCenter() const { return m_center; }
        inline double getWorldSizeX() const { return m_worldSizeX; }
        inline double getWorldSizeY() const { return m_worldSizeY; }
        inline int getMapSizeX() const { return m_mapSizeX; }
        inline int getMapSizeY() const { return m_mapSizeY; }
        inline double getDelta() const { return m_delta; }
        inline double getMapResolution() const { return m_delta; }
        inline double getResolution() const { return m_delta; }
        inline void getSize(double &xmin, double &ymin, double &xmax, double &ymax) const
        {
            Point min = map2world(0, 0), max = map2world(IntPoint(m_mapSizeX - 1, m_mapSizeY - 1));
            xmin = min.x, ymin = min.y, xmax = max.x, ymax = max.y;
        }

        /* 获取点的对应数值的相关函数 */
        inline Cell &cell(int x, int y)
        {
            return cell(IntPoint(x, y));
        }
        inline Cell &cell(const IntPoint &p);

        inline const Cell &cell(int x, int y) const
        {
            return cell(IntPoint(x, y));
        }
        inline const Cell &cell(const IntPoint &p) const;

        inline Cell &cell(double x, double y)
        {
            return cell(Point(x, y));
        }
        inline Cell &cell(const Point &p);

        inline const Cell &cell(double x, double y) const
        {
            return cell(Point(x, y));
        }

        /* 获取点的可访问性的相关函数 */
        inline bool isInside(int x, int y) const
        {
            return m_storage.cellState(IntPoint(x, y)) & Inside;
        }
        inline bool isInside(const IntPoint &p) const
        {
            return m_storage.cellState(p) & Inside;
        }

        inline bool isInside(double x, double y) const
        {
            return m_storage.cellState(world2map(x, y)) & Inside;
        }
        inline bool isInside(const Point &p) const
        {
            return m_storage.cellState(world2map(p)) & Inside;
        }

        inline const Cell &cell(const Point &p) const;

        /* 获取储存对象的相关函数 */
        inline Storage &storage() { return m_storage; }
        inline const Storage &storage() const { return m_storage; }

        /* 转换双精度浮点数数组和地图相关函数 */
        DoubleArray2D *toDoubleArray() const;
        Map<double, DoubleArray2D, false> *toDoubleMap() const;

    protected:
        Point m_center;                             /* 物理地图中心 */
        double m_worldSizeX, m_worldSizeY, m_delta; /* 物理地图尺寸和分辨率（缩放比例关系，一个像素大小(m)） */
        Storage m_storage;                          /* 用于存储地图数据的对象 */
        int m_mapSizeX, m_mapSizeY;                 /* 栅格地图的尺寸 */
        int m_sizeX2, m_sizeY2;                     /* 栅格地图的尺寸的一半，即是m_mapSizeX和m_mapSizeY的一半 */
        static const Cell m_unknown;                /* 未知值 */
    };

    /* 声明定义DoubleMap双精度浮点格式的地图 */
    typedef Map<double, DoubleArray2D, false> DoubleMap;

    /* 定义m_unknown未知值为-1 */
    template <class Cell, class Storage, const bool isClass>
    const Cell Map<Cell, Storage, isClass>::m_unknown = Cell(-1);

    /**
     * @brief Map初始化函数
     *
     * @param mapSizeX 栅格地图尺寸X
     * @param mapSizeY 栅格地图尺寸Y
     * @param delta 分辨率
     */
    template <class Cell, class Storage, const bool isClass>
    Map<Cell, Storage, isClass>::Map(int mapSizeX, int mapSizeY, double delta) : m_storage(mapSizeX, mapSizeY)
    {
        m_worldSizeX = mapSizeX * delta; /* 计算物理地图尺寸 = 栅格地图尺寸 * 分辨率 */
        m_worldSizeY = mapSizeY * delta;
        m_delta = delta;                                          /* 赋值分辨率 */
        m_center = Point(0.5 * m_worldSizeX, 0.5 * m_worldSizeY); /* 计算物理地图中心 */
        m_sizeX2 = m_mapSizeX >> 1;                               /* 计算栅格地图尺寸的一半，即栅格地图的中心 */
        m_sizeY2 = m_mapSizeY >> 1;
    }

    /**
     * @brief Map初始化函数
     *
     * @param center 物理地图中心
     * @param worldSizeX 物理地图尺寸X
     * @param worldSizeY 物理地图尺寸Y
     * @param delta 分辨率
     */
    template <class Cell, class Storage, const bool isClass>
    Map<Cell, Storage, isClass>::Map(const Point &center, double worldSizeX, double worldSizeY, double delta) : m_storage((int)ceil(worldSizeX / delta), (int)ceil(worldSizeY / delta))
    {
        m_center = center;
        m_worldSizeX = worldSizeX;
        m_worldSizeY = worldSizeY;
        m_delta = delta;

        /* 更新栅格地图尺寸 = 储存对象的XY放大间距倍数PatchSize，储存对象已在m_storage((int),(int))更新 */
        m_mapSizeX = m_storage.getXSize() << m_storage.getPatchSize();
        m_mapSizeY = m_storage.getYSize() << m_storage.getPatchSize();

        m_sizeX2 = m_mapSizeX >> 1;
        m_sizeY2 = m_mapSizeY >> 1;
    }

    /**
     * @brief Map初始化函数
     *
     * @param center 物理地图中心
     * @param xmin 物理地图x向初始最小尺寸
     * @param ymin 物理地图y向初始最小尺寸
     * @param xmax 物理地图x向初始最大尺寸
     * @param ymax 物理地图y向初始最大尺寸
     * @param delta 分辨率
     */
    template <class Cell, class Storage, const bool isClass>
    Map<Cell, Storage, isClass>::Map(const Point &center, double xmin, double ymin, double xmax, double ymax, double delta) : m_storage((int)ceil((xmax - xmin) / delta), (int)ceil((ymax - ymin) / delta))
    {
        m_center = center;
        m_worldSizeX = xmax - xmin; /* 计算物理地图尺寸 */
        m_worldSizeY = ymax - ymin;
        m_delta = delta;
        m_mapSizeX = m_storage.getXSize() << m_storage.getPatchSize();
        m_mapSizeY = m_storage.getYSize() << m_storage.getPatchSize();
        m_sizeX2 = (int)round((m_center.x - xmin) / m_delta);
        m_sizeY2 = (int)round((m_center.y - ymin) / m_delta);
    }

    /**
     * @brief 地图尺寸调整函数
     *
     * @param xmin 新物理地图x向初始最小尺寸
     * @param ymin 新物理地图y向初始最小尺寸
     * @param xmax 新物理地图x向初始最大尺寸
     * @param ymax 新物理地图y向初始最大尺寸
     */
    template <class Cell, class Storage, const bool isClass>
    void Map<Cell, Storage, isClass>::resize(double xmin, double ymin, double xmax, double ymax)
    {
        /* 转换获取新栅格地图的尺寸 */
        IntPoint imin = world2map(xmin, ymin);
        IntPoint imax = world2map(xmax, ymax);

        /* 计算新栅格地图的尺寸极值 */
        int pxmin, pymin, pxmax, pymax;
        pxmin = (int)floor((float)imin.x / (1 << m_storage.getPatchMagnitude())); /* 最小值向下取整 PatchMagnitude是个2的倍数缩放 PatchSize = 1 << PatchMagnitude */
        pxmax = (int)ceil((float)imax.x / (1 << m_storage.getPatchMagnitude()));  /* 最大值向上取整 */
        pymin = (int)floor((float)imin.y / (1 << m_storage.getPatchMagnitude()));
        pymax = (int)ceil((float)imax.y / (1 << m_storage.getPatchMagnitude()));

        /* 更新储存对象尺寸 */
        m_storage.resize(pxmin, pymin, pxmax, pymax);

        /* 更新基础参数 */
        m_mapSizeX = m_storage.getXSize() << m_storage.getPatchSize();
        m_mapSizeY = m_storage.getYSize() << m_storage.getPatchSize();
        m_worldSizeX = xmax - xmin;
        m_worldSizeY = ymax - ymin;
        m_sizeX2 -= pxmin * (1 << m_storage.getPatchMagnitude());
        m_sizeY2 -= pymin * (1 << m_storage.getPatchMagnitude());
    }

    /**
     * @brief 地图尺寸扩张函数
     *
     * 与地图尺寸调整函数resize相同效果
     *
     * @param xmin 新物理地图x向初始最小尺寸
     * @param ymin 新物理地图y向初始最小尺寸
     * @param xmax 新物理地图x向初始最大尺寸
     * @param ymax 新物理地图y向初始最大尺寸
     */
    template <class Cell, class Storage, const bool isClass>
    void Map<Cell, Storage, isClass>::grow(double xmin, double ymin, double xmax, double ymax)
    {
        IntPoint imin = world2map(xmin, ymin);
        IntPoint imax = world2map(xmax, ymax);
        if (isInside(imin) && isInside(imax))
            return;
        imin = min(imin, IntPoint(0, 0));
        imax = max(imax, IntPoint(m_mapSizeX - 1, m_mapSizeY - 1));
        int pxmin, pymin, pxmax, pymax;
        pxmin = (int)floor((float)imin.x / (1 << m_storage.getPatchMagnitude()));
        pxmax = (int)ceil((float)imax.x / (1 << m_storage.getPatchMagnitude()));
        pymin = (int)floor((float)imin.y / (1 << m_storage.getPatchMagnitude()));
        pymax = (int)ceil((float)imax.y / (1 << m_storage.getPatchMagnitude()));
        m_storage.resize(pxmin, pymin, pxmax, pymax);
        m_mapSizeX = m_storage.getXSize() << m_storage.getPatchSize();
        m_mapSizeY = m_storage.getYSize() << m_storage.getPatchSize();
        m_worldSizeX = xmax - xmin;
        m_worldSizeY = ymax - ymin;
        m_sizeX2 -= pxmin * (1 << m_storage.getPatchMagnitude());
        m_sizeY2 -= pymin * (1 << m_storage.getPatchMagnitude());
    }

    /**
     * @brief 将物理地图坐标点转换为栅格地图坐标点函数
     *
     * @param p 物理地图坐标点
     * @return 栅格地图坐标点
     */
    template <class Cell, class Storage, const bool isClass>
    IntPoint Map<Cell, Storage, isClass>::world2map(const Point &p) const
    {
        return IntPoint((int)round((p.x - m_center.x) / m_delta) + m_sizeX2, (int)round((p.y - m_center.y) / m_delta) + m_sizeY2);
    }

    /**
     * @brief 将栅格地图坐标点转换为物理地图坐标点函数
     *
     * @param p 栅格地图坐标点
     * @return 物理地图坐标点
     */
    template <class Cell, class Storage, const bool isClass>
    Point Map<Cell, Storage, isClass>::map2world(const IntPoint &p) const
    {
        return Point((p.x - m_sizeX2) * m_delta, (p.y - m_sizeY2) * m_delta) + m_center;
    }

    /* 访问地图数据相关函数 */

    /**
     * @brief 获取栅格地图对应点的值
     *
     * @param p 栅格地图坐标点
     * @return 对应点储存的栅格地图数值
     */
    template <class Cell, class Storage, const bool isClass>
    Cell &Map<Cell, Storage, isClass>::cell(const IntPoint &p)
    {
        /* 判断需要访问的点是否在范围内部 */
        AccessibilityState s = m_storage.cellState(p);
        if (!s & Inside)
            assert(0);
        // if (s&Allocated) return m_storage.cell(p); assert(0);

        // this will never happend. Just to satify the compiler..
        return m_storage.cell(p);
    }

    /**
     * @brief 获取物理地图点对应的栅格地图的值
     *
     * @param p 物理地图坐标点
     * @return 物理地图对应的栅格地图点储存数值
     */
    template <class Cell, class Storage, const bool isClass>
    Cell &Map<Cell, Storage, isClass>::cell(const Point &p)
    {
        IntPoint ip = world2map(p);
        AccessibilityState s = m_storage.cellState(ip);
        if (!s & Inside)
            assert(0);
        // if (s&Allocated) return m_storage.cell(ip); assert(0);

        // this will never happend. Just to satify the compiler..
        return m_storage.cell(ip);
    }

    /**
     * @brief 获取栅格地图对应点的值
     *
     * @param p 栅格地图坐标点
     * @return 对应的栅格地图点储存数值,未分配时为未知值,默认-1
     */
    template <class Cell, class Storage, const bool isClass>
    const Cell &Map<Cell, Storage, isClass>::cell(const IntPoint &p) const
    {
        AccessibilityState s = m_storage.cellState(p);
        // if (! s&Inside) assert(0);
        if (s & Allocated)
            return m_storage.cell(p);
        return m_unknown;
    }

    /**
     * @brief 获取物理地图点对应的栅格地图的值
     *
     * @param p 物理地图坐标点
     * @return 物理地图对应的栅格地图点储存数值,未分配时为未知值,默认-1
     */
    template <class Cell, class Storage, const bool isClass>
    const Cell &Map<Cell, Storage, isClass>::cell(const Point &p) const
    {
        IntPoint ip = world2map(p);
        AccessibilityState s = m_storage.cellState(ip);
        // if (! s&Inside) assert(0);
        if (s & Allocated)
            return m_storage.cell(ip);
        return m_unknown;
    }

    // FIXME check why the last line of the map is corrupted.
    /* FIXME 检查为什么地图的最后一行被损坏 */
    /**
     * @brief 转换双精度浮点数数组
     *
     * @return 双精度浮点数数组
     */
    template <class Cell, class Storage, const bool isClass>
    DoubleArray2D *Map<Cell, Storage, isClass>::toDoubleArray() const
    {
        DoubleArray2D *darr = new DoubleArray2D(getMapSizeX() - 1, getMapSizeY() - 1); /* 感觉这里不用-1 */
        for (int x = 0; x < getMapSizeX() - 1; x++)
            for (int y = 0; y < getMapSizeY() - 1; y++)
            {
                IntPoint p(x, y);
                darr->cell(p) = cell(p);
            }
        return darr;
    }

    /**
     * @brief 转换双精度浮点数数组地图
     *
     * @return 双精度浮点数数组地图
     */
    template <class Cell, class Storage, const bool isClass>
    Map<double, DoubleArray2D, false> *Map<Cell, Storage, isClass>::toDoubleMap() const
    {
        // FIXME size the map so that m_center will be setted accordingly
        /* FIXME 调整地图的大小，这样m_center就会被相应地设置 */
        Point pmin = map2world(IntPoint(0, 0));
        Point pmax = map2world(getMapSizeX() - 1, getMapSizeY() - 1);
        Point center = (pmax + pmin) * 0.5;

        Map<double, DoubleArray2D, false> *plainMap = new Map<double, DoubleArray2D, false>(center, (pmax - pmin).x, (pmax - pmin).y, getDelta());
        for (int x = 0; x < getMapSizeX() - 1; x++)
            for (int y = 0; y < getMapSizeY() - 1; y++)
            {
                IntPoint p(x, y);
                plainMap->cell(p) = cell(p);
            }
        return plainMap;
    }

};

#endif
