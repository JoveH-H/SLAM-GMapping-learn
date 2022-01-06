#ifndef SCANMATCHER_H
#define SCANMATCHER_H

#include "gmapping/scanmatcher/icp.h"
#include "gmapping/scanmatcher/smmap.h"
#include <gmapping/utils/macro_params.h>
#include <gmapping/utils/stat.h>
#include <iostream>
#include <gmapping/utils/gvalues.h>
#include <gmapping/scanmatcher/scanmatcher_export.h>
#define LASER_MAXBEAMS 2048

namespace GMapping
{

    class SCANMATCHER_EXPORT ScanMatcher
    {
    public:
        typedef Covariance3 CovarianceMatrix;

        ScanMatcher();
        ~ScanMatcher();
        double icpOptimize(OrientedPoint &pnew, const ScanMatcherMap &map, const OrientedPoint &p, const double *readings) const;

        /* 爬山优化器函数 */
        double optimize(OrientedPoint &pnew, const ScanMatcherMap &map, const OrientedPoint &p, const double *readings) const;
        double optimize(OrientedPoint &mean, CovarianceMatrix &cov, const ScanMatcherMap &map, const OrientedPoint &p, const double *readings) const;

        double registerScan(ScanMatcherMap &map, const OrientedPoint &p, const double *readings); /* 扫描数据注册函数 */
        void setLaserParameters(unsigned int beams, double *angles, const OrientedPoint &lpose);  /* 扫描匹配器激光数据设置函数 */

        /* 扫描匹配器参数设置函数 */
        void setMatchingParameters(double urange, double range, double sigma, int kernsize, double lopt, double aopt, int iterations, double likelihoodSigma = 1, unsigned int likelihoodSkip = 0);

        void invalidateActiveArea();                                                                 /* 有效面积未完成计算标志函数 */
        void computeActiveArea(ScanMatcherMap &map, const OrientedPoint &p, const double *readings); /* 有效面积计算函数 */

        inline double icpStep(OrientedPoint &pret, const ScanMatcherMap &map, const OrientedPoint &p, const double *readings) const;

        /* 获取粒子的匹配度函数 */
        inline double score(const ScanMatcherMap &map, const OrientedPoint &p, const double *readings) const;

        /* 获取粒子的似然度和匹配度函数 */
        inline unsigned int likelihoodAndScore(double &s, double &l, const ScanMatcherMap &map, const OrientedPoint &p, const double *readings) const;

        double likelihood(double &lmax, OrientedPoint &mean, CovarianceMatrix &cov, const ScanMatcherMap &map, const OrientedPoint &p, const double *readings);
        double likelihood(double &_lmax, OrientedPoint &_mean, CovarianceMatrix &_cov, const ScanMatcherMap &map, const OrientedPoint &p, Gaussian3 &odometry, const double *readings, double gain = 180.);
        inline const double *laserAngles() const { return m_laserAngles; }
        inline unsigned int laserBeams() const { return m_laserBeams; }

        static const double nullLikelihood;

    protected:
        // state of the matcher
        bool m_activeAreaComputed;

        /**laser parameters*/
        unsigned int m_laserBeams;
        double m_laserAngles[LASER_MAXBEAMS];
        // OrientedPoint m_laserPose;
        PARAM_SET_GET(OrientedPoint, laserPose, protected, public, public)
        PARAM_SET_GET(double, laserMaxRange, protected, public, public)
        /**scan_matcher parameters*/
        PARAM_SET_GET(double, usableRange, protected, public, public)
        PARAM_SET_GET(double, gaussianSigma, protected, public, public)
        PARAM_SET_GET(double, likelihoodSigma, protected, public, public)
        PARAM_SET_GET(int, kernelSize, protected, public, public)
        PARAM_SET_GET(double, optAngularDelta, protected, public, public)
        PARAM_SET_GET(double, optLinearDelta, protected, public, public)
        PARAM_SET_GET(unsigned int, optRecursiveIterations, protected, public, public)
        PARAM_SET_GET(unsigned int, likelihoodSkip, protected, public, public)
        PARAM_SET_GET(double, llsamplerange, protected, public, public)
        PARAM_SET_GET(double, llsamplestep, protected, public, public)
        PARAM_SET_GET(double, lasamplerange, protected, public, public)
        PARAM_SET_GET(double, lasamplestep, protected, public, public)
        PARAM_SET_GET(bool, generateMap, protected, public, public)
        PARAM_SET_GET(double, enlargeStep, protected, public, public)
        PARAM_SET_GET(double, fullnessThreshold, protected, public, public)
        PARAM_SET_GET(double, angularOdometryReliability, protected, public, public)
        PARAM_SET_GET(double, linearOdometryReliability, protected, public, public)
        PARAM_SET_GET(double, freeCellRatio, protected, public, public)
        PARAM_SET_GET(unsigned int, initialBeamsSkip, protected, public, public)

        // allocate this large array only once
        IntPoint *m_linePoints;
    };

    inline double ScanMatcher::icpStep(OrientedPoint &pret, const ScanMatcherMap &map, const OrientedPoint &p, const double *readings) const
    {
        const double *angle = m_laserAngles + m_initialBeamsSkip;
        OrientedPoint lp = p;
        lp.x += cos(p.theta) * m_laserPose.x - sin(p.theta) * m_laserPose.y;
        lp.y += sin(p.theta) * m_laserPose.x + cos(p.theta) * m_laserPose.y;
        lp.theta += m_laserPose.theta;
        unsigned int skip = 0;
        double freeDelta = map.getDelta() * m_freeCellRatio;
        std::list<PointPair> pairs;

        for (const double *r = readings + m_initialBeamsSkip; r < readings + m_laserBeams; r++, angle++)
        {
            skip++;
            skip = skip > m_likelihoodSkip ? 0 : skip;
            if (*r > m_usableRange || *r == 0.0)
                continue;
            if (skip)
                continue;
            Point phit = lp;
            phit.x += *r * cos(lp.theta + *angle);
            phit.y += *r * sin(lp.theta + *angle);
            IntPoint iphit = map.world2map(phit);
            Point pfree = lp;
            pfree.x += (*r - map.getDelta() * freeDelta) * cos(lp.theta + *angle);
            pfree.y += (*r - map.getDelta() * freeDelta) * sin(lp.theta + *angle);
            pfree = pfree - phit;
            IntPoint ipfree = map.world2map(pfree);
            bool found = false;
            Point bestMu(0., 0.);
            Point bestCell(0., 0.);
            for (int xx = -m_kernelSize; xx <= m_kernelSize; xx++)
                for (int yy = -m_kernelSize; yy <= m_kernelSize; yy++)
                {
                    IntPoint pr = iphit + IntPoint(xx, yy);
                    IntPoint pf = pr + ipfree;
                    // AccessibilityState s=map.storage().cellState(pr);
                    // if (s&Inside && s&Allocated){
                    const PointAccumulator &cell = map.cell(pr);
                    const PointAccumulator &fcell = map.cell(pf);
                    if (((double)cell) > m_fullnessThreshold && ((double)fcell) < m_fullnessThreshold)
                    {
                        Point mu = phit - cell.mean();
                        if (!found)
                        {
                            bestMu = mu;
                            bestCell = cell.mean();
                            found = true;
                        }
                        else if ((mu * mu) < (bestMu * bestMu))
                        {
                            bestMu = mu;
                            bestCell = cell.mean();
                        }
                    }
                    //}
                }
            if (found)
            {
                pairs.push_back(std::make_pair(phit, bestCell));
                // std::cerr << "(" << phit.x-bestCell.x << "," << phit.y-bestCell.y << ") ";
            }
            // std::cerr << std::endl;
        }

        OrientedPoint result(0, 0, 0);
        // double icpError=icpNonlinearStep(result,pairs);
        std::cerr << "result(" << pairs.size() << ")=" << result.x << " " << result.y << " " << result.theta << std::endl;
        pret.x = p.x + result.x;
        pret.y = p.y + result.y;
        pret.theta = p.theta + result.theta;
        pret.theta = atan2(sin(pret.theta), cos(pret.theta));
        return score(map, p, readings);
    }

    /**
     * @brief 获取粒子的匹配度函数
     *
     * @param map 粒子位姿地图
     * @param p 粒子位姿
     * @param readings 激光雷达传感器数据
     * @return s 匹配度
     */
    inline double ScanMatcher::score(const ScanMatcherMap &map, const OrientedPoint &p, const double *readings) const
    {
        double s = 0;                                             /* 初始化匹配度为0 */
        const double *angle = m_laserAngles + m_initialBeamsSkip; /* 使用指针angle获取传感器各个扫描数据所对应的转角 */

        /* 创建一个临时变量lp将粒子位姿p加上激光传感器的位置，即将粒子位姿p转换到地图坐标系下 */
        OrientedPoint lp = p;
        lp.x += cos(p.theta) * m_laserPose.x - sin(p.theta) * m_laserPose.y;
        lp.y += sin(p.theta) * m_laserPose.x + cos(p.theta) * m_laserPose.y;
        lp.theta += m_laserPose.theta;

        unsigned int skip = 0;                               /* 用于对传感器数据筛选的计数器 */
        double freeDelta = map.getDelta() * m_freeCellRatio; /* 空闲区域增量，用于确定空闲地图区域的参数，传感器数据减去该值作为空闲区域的判定标准 */

        /* 遍历激光数据 */
        for (const double *r = readings + m_initialBeamsSkip; r < readings + m_laserBeams; r++, angle++)
        {
            /* 累积筛选跳过忽略的计数 */
            skip++;

            /* 当达到跳跃步长时，将其清零，并进行一次匹配操作 */
            skip = skip > m_likelihoodSkip ? 0 : skip;

            /* 在跳过的范围内或超出设置的激光雷达最大使用量程时忽略进入下一个数据 */
            if (skip || *r > m_usableRange || *r == 0.0)
                continue;

            /* 根据扫描数据的距离读数和角度创建一个点 phit */
            Point phit = lp;
            phit.x += *r * cos(lp.theta + *angle);
            phit.y += *r * sin(lp.theta + *angle);

            /* 将扫描数据对应的点 phit 转换成地图数组中的索引点 iphit，该点在地图中应当是被占用的 */
            IntPoint iphit = map.world2map(phit);

            /* 根据扫描波束所指的方向构建一个点 pfree，该点在地图中应当是被空闲的 */
            Point pfree = lp;
            pfree.x += (*r - map.getDelta() * freeDelta) * cos(lp.theta + *angle); /* 由于传感器的读数存在噪声，所以减去一个合适的值作为空闲区域的判定标准 */
            pfree.y += (*r - map.getDelta() * freeDelta) * sin(lp.theta + *angle);

            /* 将扫描数据方向上应空闲的点 pfree与扫描数据对应的点 phit的差距 转换成地图数组上的差距索引 ipfree */
            pfree = pfree - phit;
            IntPoint ipfree = map.world2map(pfree);

            bool found = false;   /* 记录在扫描数据对应的点 phit 的一个邻域内有匹配点 */
            Point bestMu(0., 0.); /* 定义一个二维向量用于记录匹配点到邻域点之间的距离向量，该向量的长度越小说明两个点越接近，匹配度就越高 */

            /* 在扫描数据对应的点 phit 的邻域内搜索匹配点 */
            for (int xx = -m_kernelSize; xx <= m_kernelSize; xx++)
                for (int yy = -m_kernelSize; yy <= m_kernelSize; yy++)
                {
                    /* 在扫描数据对应点的一个邻域内搜索匹配点及其对应的空闲区域判定点 */
                    IntPoint pr = iphit + IntPoint(xx, yy);
                    IntPoint pf = pr + ipfree;

                    // AccessibilityState s=map.storage().cellState(pr);
                    // if (s&Inside && s&Allocated){

                    /* 查询地图获取邻域点和空闲区域判定点对应栅格的占用概率、均值和熵 */
                    const PointAccumulator &cell = map.cell(pr);
                    const PointAccumulator &fcell = map.cell(pf);

                    /* 当扫描数据对应邻域点被占用且对应的空闲区域判定点是空闲的，就认为匹配了一个点 */
                    if (((double)cell) > m_fullnessThreshold && ((double)fcell) < m_fullnessThreshold)
                    {
                        /* 计算匹配点到对应邻域点的矢量差 */
                        Point mu = phit - cell.mean();

                        /* 首次匹配成功时默认设置最优矢量差，并改变已匹配标志位 */
                        if (!found)
                        {
                            bestMu = mu;
                            found = true;
                        }
                        else
                        {
                            /* 更新最优矢量差， 匹配点到地图栅格中心的矢量差越小匹配程度就越高 */
                            bestMu = (mu * mu) < (bestMu * bestMu) ? mu : bestMu;
                        }
                    }
                    //}
                }

            /* 扫描数据对应点匹配成功 */
            if (found)
                s += exp(-1. / m_gaussianSigma * bestMu * bestMu); /* 累加扫描数据的匹配度 完全准确的话是 exp(-1/0.05 *0 *0)=1 */
        }
        return s;
    }

    /**
     * @brief 获取粒子的似然度和匹配度函数
     *
     * @param s 匹配度
     * @param l 似然度
     * @param map 粒子位姿地图
     * @param p 粒子位姿
     * @param readings 激光雷达传感器数据
     * @return c 匹配成功的数据点的数量
     */
    inline unsigned int ScanMatcher::likelihoodAndScore(double &s, double &l, const ScanMatcherMap &map, const OrientedPoint &p, const double *readings) const
    {
        using namespace std;

        /* 初始化似然度和匹配度为0 */
        l = 0;
        s = 0;

        /* 使用指针angle获取传感器各个扫描数据所对应的转角 */
        const double *angle = m_laserAngles + m_initialBeamsSkip;

        /* 创建一个临时变量lp将粒子位姿p加上激光传感器的位置，即将粒子位姿p转换到地图坐标系下 */
        OrientedPoint lp = p;
        lp.x += cos(p.theta) * m_laserPose.x - sin(p.theta) * m_laserPose.y;
        lp.y += sin(p.theta) * m_laserPose.x + cos(p.theta) * m_laserPose.y;
        lp.theta += m_laserPose.theta;

        double noHit = nullLikelihood / (m_likelihoodSigma); /* 未匹配的似然概率值 */
        unsigned int skip = 0;                               /* 用于对传感器数据筛选的计数器 */
        unsigned int c = 0;                                  /* 匹配计数器 */
        double freeDelta = map.getDelta() * m_freeCellRatio; /* 空闲区域增量，用于确定空闲地图区域的参数，传感器数据减去该值作为空闲区域的判定标准 */

        /* 遍历激光数据 */
        for (const double *r = readings + m_initialBeamsSkip; r < readings + m_laserBeams; r++, angle++)
        {
            /* 累积筛选跳过忽略的计数 */
            skip++;

            /* 当达到跳跃步长时，将其清零，并进行一次匹配操作 */
            skip = skip > m_likelihoodSkip ? 0 : skip;

            /* 超出设置的激光雷达最大使用量程时忽略进入下一个数据 */
            if (*r > m_usableRange)
                continue;

            /* 若在跳过的范围内则忽略进入下一个数据 */
            if (skip)
                continue;

            /* 根据扫描数据的距离读数和角度创建一个点 phit */
            Point phit = lp;
            phit.x += *r * cos(lp.theta + *angle);
            phit.y += *r * sin(lp.theta + *angle);

            /* 将扫描数据对应的点 phit 转换成地图数组中的索引点 iphit，该点在地图中应当是被占用的 */
            IntPoint iphit = map.world2map(phit);

            /* 根据扫描波束所指的方向构建一个点 pfree，该点在地图中应当是被空闲的 */
            Point pfree = lp;
            pfree.x += (*r - freeDelta) * cos(lp.theta + *angle); /* 由于传感器的读数存在噪声，所以减去一个合适的值作为空闲区域的判定标准 */
            pfree.y += (*r - freeDelta) * sin(lp.theta + *angle);

            /* 将扫描数据方向上应空闲的点 pfree与扫描数据对应的点 phit的差距 转换成地图数组上的差距索引 ipfree */
            pfree = pfree - phit;
            IntPoint ipfree = map.world2map(pfree);

            bool found = false;   /* 记录在扫描数据对应的点 phit 的一个邻域内有匹配点 */
            Point bestMu(0., 0.); /* 定义一个二维向量用于记录匹配点到邻域点之间的距离向量，该向量的长度越小说明两个点越接近，匹配度就越高 */

            /* 在扫描数据对应的点 phit 的邻域内搜索匹配点 */
            for (int xx = -m_kernelSize; xx <= m_kernelSize; xx++)
                for (int yy = -m_kernelSize; yy <= m_kernelSize; yy++)
                {
                    /* 在扫描数据对应点的一个邻域内搜索匹配点及其对应的空闲区域判定点 */
                    IntPoint pr = iphit + IntPoint(xx, yy);
                    IntPoint pf = pr + ipfree;

                    // AccessibilityState s=map.storage().cellState(pr);
                    // if (s&Inside && s&Allocated){

                    /* 查询地图获取邻域点和空闲区域判定点对应栅格的占用概率、均值和熵 */
                    const PointAccumulator &cell = map.cell(pr);
                    const PointAccumulator &fcell = map.cell(pf);

                    /* 当扫描数据对应邻域点被占用且对应的空闲区域判定点是空闲的，就认为匹配了一个点 */
                    if (((double)cell) > m_fullnessThreshold && ((double)fcell) < m_fullnessThreshold)
                    {
                        /* 计算匹配点到对应邻域点的矢量差 */
                        Point mu = phit - cell.mean();

                        /* 首次匹配成功时默认设置最优矢量差，并改变已匹配标志位 */
                        if (!found)
                        {
                            bestMu = mu;
                            found = true;
                        }
                        else
                        {
                            /* 更新最优矢量差， 匹配点到地图栅格中心的矢量差越小匹配程度就越高 */
                            bestMu = (mu * mu) < (bestMu * bestMu) ? mu : bestMu;
                        }
                    }
                    //}
                }

            /* 扫描数据对应点匹配成功 */
            if (found)
            {
                s += exp(-1. / m_gaussianSigma * bestMu * bestMu); /* 累加扫描数据的匹配度 完全准确的话是 exp(-1/0.05 *0 *0)=1 */
                c++;                                               /* 匹配计数器 */
            }

            if (!skip)
            {
                double f = (-1. / m_likelihoodSigma) * (bestMu * bestMu); /* 高斯分布的假设下计算该点的对数似然度 */
                l += (found) ? f : noHit;                                 /* 累加似然度，用于更新粒子的权重 */
            }
        }
        return c;
    }

};

#endif
