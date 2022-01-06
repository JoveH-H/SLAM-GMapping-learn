#include <cstring>
#include <limits>
#include <list>
#include <iostream>

#include "gmapping/scanmatcher/scanmatcher.h"
#include "gmapping/scanmatcher/gridlinetraversal.h"
//#define GENERATE_MAPS

namespace GMapping
{

    using namespace std;

    const double ScanMatcher::nullLikelihood = -.5; /* 空的似然 */

    ScanMatcher::ScanMatcher() : m_laserPose(0, 0, 0)
    {
        // m_laserAngles=0;
        m_laserBeams = 0;             /* 波束数量 */
        m_optRecursiveIterations = 3; /* 扫描匹配的迭代步数 */
        m_activeAreaComputed = false; /* 有效面积计算（地图已更新）标志位 */

        // This  are the dafault settings for a grid map of 5 cm
        /* 这是5厘米网格地图的默认设置 */

        m_llsamplerange = 0.01;
        m_llsamplestep = 0.01;
        m_lasamplerange = 0.005;
        m_lasamplestep = 0.005;

        m_enlargeStep = 10.;               /* 拓展地图时的边界膨胀值 */
        m_fullnessThreshold = 0.1;         /* 判定栅格被占用的阈值 */
        m_angularOdometryReliability = 0.; /* 里程计的角度置信度 */
        m_linearOdometryReliability = 0.;  /* 里程计的平移置信度 */
        m_freeCellRatio = sqrt(2.);        /* 空闲地图点（区域）的比例 */
        m_initialBeamsSkip = 0;            /* 初始光束跳过数量 */

        /*
            // This  are the dafault settings for a grid map of 10 cm
            m_llsamplerange=0.1;
            m_llsamplestep=0.1;
            m_lasamplerange=0.02;
            m_lasamplestep=0.01;
        */
        // This  are the dafault settings for a grid map of 20/25 cm
        /*
            m_llsamplerange=0.2;
            m_llsamplestep=0.1;
            m_lasamplerange=0.02;
            m_lasamplestep=0.01;
            m_generateMap=false;
        */

        m_linePoints = new IntPoint[20000];
    }

    ScanMatcher::~ScanMatcher()
    {
        delete[] m_linePoints;
    }

    /**
     * @brief 有效面积未完成计算标志函数
     *
     * 标记即将进入地图更新步骤
     */
    void ScanMatcher::invalidateActiveArea()
    {
        /* 有效面积计算标志位置false */
        m_activeAreaComputed = false;
    }

    /**
     * @brief 有效面积计算函数
     *
     * 地图有效区域更新，但未更新栅格地图上数组的值
     *
     * @param map 当前位姿的地图
     * @param p 当前位姿
     * @param readings 激光雷达传感器数据
     */
    void ScanMatcher::computeActiveArea(ScanMatcherMap &map, const OrientedPoint &p, const double *readings)
    {
        /* 检验是否已标注有效面积未完成计算 */
        if (m_activeAreaComputed)
            return;

        /* 创建一个临时变量lp将粒子位姿p加上激光传感器的位置，即将粒子位姿p转换到地图坐标系下 */
        OrientedPoint lp = p;
        lp.x += cos(p.theta) * m_laserPose.x - sin(p.theta) * m_laserPose.y;
        lp.y += sin(p.theta) * m_laserPose.x + cos(p.theta) * m_laserPose.y;
        lp.theta += m_laserPose.theta;

        /* 地图坐标系下粒子位姿lp 转换成地图数组中的索引点 p0 */
        IntPoint p0 = map.world2map(lp);

        /* 计算物理地图的边界尺寸 */
        Point min(map.map2world(0, 0));
        Point max(map.map2world(map.getMapSizeX() - 1, map.getMapSizeY() - 1));

        /* 根据粒子位姿更新物理地图边界数值 */
        if (lp.x < min.x)
            min.x = lp.x;
        if (lp.y < min.y)
            min.y = lp.y;
        if (lp.x > max.x)
            max.x = lp.x;
        if (lp.y > max.y)
            max.y = lp.y;

        /*determine the size of the area*/
        /* 确定区域的大小 */

        /* 使用指针angle获取传感器各个扫描数据所对应的转角 */
        const double *angle = m_laserAngles + m_initialBeamsSkip;

        /* 遍历激光数据 */
        for (const double *r = readings + m_initialBeamsSkip; r < readings + m_laserBeams; r++, angle++)
        {
            /* 筛选正常数据，除去大于激光雷达最大的量程或为0或非数值 */
            if (*r > m_laserMaxRange || *r == 0.0 || isnan(*r))
                continue;

            /* 激光雷达最大使用距离限幅 */
            double d = *r > m_usableRange ? m_usableRange : *r;

            /* 根据扫描数据的距离读数和角度创建一个点 phit */
            Point phit = lp;
            phit.x += d * cos(lp.theta + *angle);
            phit.y += d * sin(lp.theta + *angle);

            /* 根据激光数据更新物理地图边界数值 */
            if (phit.x < min.x)
                min.x = phit.x;
            if (phit.y < min.y)
                min.y = phit.y;
            if (phit.x > max.x)
                max.x = phit.x;
            if (phit.y > max.y)
                max.y = phit.y;
        }

        // min=min-Point(map.getDelta(),map.getDelta());
        // max=max+Point(map.getDelta(),map.getDelta());

        /* 粒子位姿或激光数据不在地图数组范围内，即需要拓展地图数组时 */
        if (!map.isInside(min) || !map.isInside(max))
        {
            /* 计算现在物理地图的边界尺寸 */
            Point lmin(map.map2world(0, 0));
            Point lmax(map.map2world(map.getMapSizeX() - 1, map.getMapSizeY() - 1));

            // cerr << "CURRENT MAP " << lmin.x << " " << lmin.y << " " << lmax.x << " " << lmax.y << endl;
            // cerr << "BOUNDARY OVERRIDE " << min.x << " " << min.y << " " << max.x << " " << max.y << endl;

            /* 计算需要调整的边界，并增加边缘膨胀m_enlargeStep */
            min.x = (min.x >= lmin.x) ? lmin.x : min.x - m_enlargeStep;
            max.x = (max.x <= lmax.x) ? lmax.x : max.x + m_enlargeStep;
            min.y = (min.y >= lmin.y) ? lmin.y : min.y - m_enlargeStep;
            max.y = (max.y <= lmax.y) ? lmax.y : max.y + m_enlargeStep;

            /* 地图数组尺寸调整 */
            map.resize(min.x, min.y, max.x, max.y);
            // cerr << "RESIZE " << min.x << " " << min.y << " " << max.x << " " << max.y << endl;
        }

        /* 定义类似金字塔的形式组织地图内存 */
        HierarchicalArray2D<PointAccumulator>::PointSet activeArea;

        /*allocate the active area*/
        /* 分配活动区域 */

        /* 更新传感器各个扫描数据所对应的转角，其实应该不用的，之前已经设置过了 */
        angle = m_laserAngles + m_initialBeamsSkip;

        /* 遍历激光数据 */
        for (const double *r = readings + m_initialBeamsSkip; r < readings + m_laserBeams; r++, angle++)

            /* 需要生成地图的情况 */
            if (m_generateMap)
            {
                /* 筛选正常数据，除去大于激光雷达最大量程或为0或非数值 */
                double d = *r;
                if (d > m_laserMaxRange || d == 0.0 || isnan(d))
                    continue;

                /* 激光雷达最大使用距离限幅 */
                if (d > m_usableRange)
                    d = m_usableRange;

                /* 计算粒子位姿（传感器位置）和激光扫描点（即光束的起点和终点）在物理地图的坐标 */
                Point phit = lp + Point(d * cos(lp.theta + *angle), d * sin(lp.theta + *angle));
                IntPoint p0 = map.world2map(lp);
                IntPoint p1 = map.world2map(phit);

                /* bresenham算法来计算激光起点到终点要经过的路径line */
                // IntPoint linePoints[20000] ;
                GridLineTraversalLine line;
                line.points = m_linePoints;
                GridLineTraversal::gridLine(p0, p1, &line);

                /* 获取画线算法计算光束路径上相关的点在栅格地图上的块索引，并增加有效区域 */
                for (int i = 0; i < line.num_points - 1; i++)
                {
                    assert(map.isInside(m_linePoints[i]));
                    activeArea.insert(map.storage().patchIndexes(m_linePoints[i]));
                    assert(m_linePoints[i].x >= 0 && m_linePoints[i].y >= 0);
                }

                /* 如果激光距离小于激光雷达最大使用距离 */
                if (d < m_usableRange)
                {
                    /* 获取扫描数据对应的点在栅格地图上的块索引，并增加有效区域 */
                    IntPoint cp = map.storage().patchIndexes(p1);
                    assert(cp.x >= 0 && cp.y >= 0);
                    activeArea.insert(cp);
                }
            }
            else /* 单独处理占用点的情况 */
            {
                /* 筛选正常数据，除去大于激光雷达最大量程或大于最大使用距离或为0或非数值 */
                if (*r > m_laserMaxRange || *r > m_usableRange || *r == 0.0 || isnan(*r))
                    continue;

                /* 获取扫描数据对应的点在地图的坐标，并转换到栅格地图的索引 */
                Point phit = lp;
                phit.x += *r * cos(lp.theta + *angle);
                phit.y += *r * sin(lp.theta + *angle);
                IntPoint p1 = map.world2map(phit);
                assert(p1.x >= 0 && p1.y >= 0);

                /* 获取扫描数据对应的点在栅格地图上的块索引，并增加有效区域 */
                IntPoint cp = map.storage().patchIndexes(p1);
                assert(cp.x >= 0 && cp.y >= 0);
                activeArea.insert(cp);
            }

        // this allocates the unallocated cells in the active area of the map
        // cout << "activeArea::size() " << activeArea.size() << endl;
        /*
            cerr << "ActiveArea=";
            for (HierarchicalArray2D<PointAccumulator>::PointSet::const_iterator it=activeArea.begin(); it!= activeArea.end(); it++){
                cerr << "(" << it->x <<"," << it->y << ") ";
            }
            cerr << endl;
        */

        /* 设置有效区域 */
        map.storage().setActiveArea(activeArea, true);

        /* 有效面积计算完成标志位 */
        m_activeAreaComputed = true;
    }

    /**
     * @brief 扫描数据注册函数
     *
     * 将激光扫描的占用信息注册到栅格地图上
     *
     * @param map 当前位姿的地图
     * @param p 当前位姿
     * @param readings 激光雷达传感器数据
     * @return esum 总体的熵的变化
     */
    double ScanMatcher::registerScan(ScanMatcherMap &map, const OrientedPoint &p, const double *readings)
    {
        /* 确保已更新有效区域 */
        if (!m_activeAreaComputed)
            computeActiveArea(map, p, readings);

        // this operation replicates the cells that will be changed in the registration operation
        /* 此操作将复制在注册操作中更改的单元格，即申请必要的存储空间 */
        map.storage().allocActiveArea();

        /* 创建一个临时变量lp将粒子位姿p加上激光传感器的位置，即将粒子位姿p转换到地图坐标系下 */
        OrientedPoint lp = p;
        lp.x += cos(p.theta) * m_laserPose.x - sin(p.theta) * m_laserPose.y;
        lp.y += sin(p.theta) * m_laserPose.x + cos(p.theta) * m_laserPose.y;
        lp.theta += m_laserPose.theta;

        /* 地图坐标系下粒子位姿lp 转换成地图数组中的索引点 p0 */
        IntPoint p0 = map.world2map(lp);

        /* 使用指针angle获取传感器各个扫描数据所对应的转角 */
        const double *angle = m_laserAngles + m_initialBeamsSkip;

        /* 总体的熵的变化 */
        double esum = 0;

        /* 遍历激光数据 */
        for (const double *r = readings + m_initialBeamsSkip; r < readings + m_laserBeams; r++, angle++)

            /* 需要生成地图的情况 */
            if (m_generateMap)
            {
                /* 筛选正常数据，除去大于激光雷达最大量程或为0或非数值 */
                double d = *r;
                if (d > m_laserMaxRange || d == 0.0 || isnan(d))
                    continue;

                /* 激光雷达最大使用距离限幅 */
                if (d > m_usableRange)
                    d = m_usableRange;

                /* 计算激光扫描点在物理地图的坐标 */
                Point phit = lp + Point(d * cos(lp.theta + *angle), d * sin(lp.theta + *angle));
                IntPoint p1 = map.world2map(phit);

                /* bresenham算法来计算激光起点到终点要经过的路径line */
                // IntPoint linePoints[20000] ;
                GridLineTraversalLine line;
                line.points = m_linePoints;
                GridLineTraversal::gridLine(p0, p1, &line);

                /* 把画线算法计算光束路径上相关的点都设置为空闲，并更新累积熵 */
                for (int i = 0; i < line.num_points - 1; i++)
                {
                    PointAccumulator &cell = map.cell(line.points[i]);
                    double e = -cell.entropy();
                    cell.update(false, Point(0, 0));
                    e += cell.entropy();
                    esum += e;
                }

                /* 如果激光距离小于激光雷达最大使用距离 */
                if (d < m_usableRange)
                {
                    /* 将扫描数据对应的点在栅格地图上设置为占用，并更新累积熵 */
                    double e = -map.cell(p1).entropy();
                    map.cell(p1).update(true, phit);
                    e += map.cell(p1).entropy();
                    esum += e;
                }
            }
            else /* 单独处理占用点的情况 */
            {
                /* 筛选正常数据，除去大于激光雷达最大量程或大于最大使用距离或为0或非数值 */
                if (*r > m_laserMaxRange || *r > m_usableRange || *r == 0.0 || isnan(*r))
                    continue;

                /* 获取扫描数据对应的点在地图的坐标，并转换到栅格地图的索引 */
                Point phit = lp;
                phit.x += *r * cos(lp.theta + *angle);
                phit.y += *r * sin(lp.theta + *angle);
                IntPoint p1 = map.world2map(phit);
                assert(p1.x >= 0 && p1.y >= 0);

                /* 设置占用格 */
                map.cell(p1).update(true, phit);
            }
        // cout  << "informationGain=" << -esum << endl;
        return esum;
    }

    double ScanMatcher::icpOptimize(OrientedPoint &pnew, const ScanMatcherMap &map, const OrientedPoint &init, const double *readings) const
    {
        double currentScore;
        double sc = score(map, init, readings);
        ;
        OrientedPoint start = init;
        pnew = init;
        int iterations = 0;
        do
        {
            currentScore = sc;
            sc = icpStep(pnew, map, start, readings);
            // cerr << "pstart=" << start.x << " " <<start.y << " " << start.theta << endl;
            // cerr << "pret=" << pnew.x << " " <<pnew.y << " " << pnew.theta << endl;
            start = pnew;
            iterations++;
        } while (sc > currentScore);
        cerr << "i=" << iterations << endl;
        return currentScore;
    }

    /**
     * @brief 爬山优化器函数
     *
     * @param pnew 矫正过后的新位姿
     * @param map 当前位姿的地图
     * @param init 当前位姿
     * @param readings 激光雷达传感器数据
     * @return bestScore 局部最优的位姿的匹配度
     */
    double ScanMatcher::optimize(OrientedPoint &pnew, const ScanMatcherMap &map, const OrientedPoint &init, const double *readings) const
    {
        double bestScore = -1;                                        /* 搜索过程中的最优匹配度 */
        OrientedPoint currentPose = init;                             /* 爬山算法当前考察节点 */
        double currentScore = score(map, currentPose, readings);      /* 当前考察节点的匹配度 */
        double adelta = m_optAngularDelta, ldelta = m_optLinearDelta; /* 爬山步长 */
        unsigned int refinement = 0;                                  /* 优化次数 */

        /* 定义枚举 爬山方向状态 */
        enum Move
        {
            Front,
            Back,
            Left,
            Right,
            TurnLeft,
            TurnRight,
            Done
        };
        /* cout << __func__<<  " readings: ";
            for (int i=0; i<m_laserBeams; i++){
                cout << readings[i] << " ";
            }
            cout << endl;
        */

        int c_iterations = 0; /* 迭代次数 */

        /* 开始爬山过程 */
        do
        {
            /* 当最优匹配度大于当前匹配度时，即优化无效果时 */
            if (bestScore >= currentScore)
            {
                refinement++; /* 优化次数+1 */
                adelta *= .5; /* 旋转优化步长减半，细化优化，可防止步长过大导致在最优值附近反复跳动 */
                ldelta *= .5; /* 平移优化步长减半 */
            }

            /* 设置最优匹配度为当前匹配度 */
            bestScore = currentScore;

            // cout <<"score="<< currentScore << " refinement=" << refinement;
            // cout <<  "pose=" << currentPose.x  << " " << currentPose.y << " " << currentPose.theta << endl;

            /* 定义记录当前最优的节点为当前原节点 */
            OrientedPoint bestLocalPose = currentPose;

            /* 定义本地节点，用于遍历当前节点周围邻接节点 */
            OrientedPoint localPose = currentPose;

            /* 定义爬山的方向 */
            Move move = Front;

            /* 遍历邻接节点 */
            do
            {
                /* 恢复当前本地节点，用于重新访问周围邻接节点 */
                localPose = currentPose;

                /* 选爬山方向并更改访问邻接节点 */
                switch (move)
                {
                case Front:
                    localPose.x += ldelta;
                    move = Back;
                    break;
                case Back:
                    localPose.x -= ldelta;
                    move = Left;
                    break;
                case Left:
                    localPose.y -= ldelta;
                    move = Right;
                    break;
                case Right:
                    localPose.y += ldelta;
                    move = TurnLeft;
                    break;
                case TurnLeft:
                    localPose.theta += adelta;
                    move = TurnRight;
                    break;
                case TurnRight:
                    localPose.theta -= adelta;
                    move = Done;
                    break;
                default:;
                }

                /* 里程计数据增益 */
                double odo_gain = 1;

                /* 根据里程计的角度置信度参数计算odo_gain, 默认为0，不会修改 */
                if (m_angularOdometryReliability > 0.)
                {
                    double dth = init.theta - localPose.theta;
                    dth = atan2(sin(dth), cos(dth));
                    dth *= dth;
                    odo_gain *= exp(-m_angularOdometryReliability * dth);
                }

                /* 根据里程计的平移置信度参数计算odo_gain, 默认为0，不会修改 */
                if (m_linearOdometryReliability > 0.)
                {
                    double dx = init.x - localPose.x;
                    double dy = init.y - localPose.y;
                    double drho = dx * dx + dy * dy;
                    odo_gain *= exp(-m_linearOdometryReliability * drho);
                }

                /* 计算当前邻接节点的匹配度 */
                double localScore = odo_gain * score(map, localPose, readings);

                /* 若当前邻接节点匹配度高于当前节点分数，则更新得到当前节点周边最优的节点 */
                if (localScore > currentScore)
                {
                    /* 更新最优分数和位姿 */
                    currentScore = localScore;
                    bestLocalPose = localPose;
                }

                c_iterations++; /* 迭代次数+1 */

            } while (move != Done);

            /* 更新当前位姿为遍历邻接节点后的最优位姿 */
            currentPose = bestLocalPose;

            // cout << "currentScore=" << currentScore<< endl;
            // here we look for the best move;

            /* 当前匹配度优化无效并且优化次数超出了指定的迭代次数时退出爬山寻优 */
        } while (currentScore > bestScore || refinement < m_optRecursiveIterations);

        // cout << __func__ << "bestScore=" << bestScore<< endl;
        // cout << __func__ << "iterations=" << c_iterations<< endl;

        /* 赋值传参pnew当前局部最优的位姿 */
        pnew = currentPose;

        /* 返回局部最优的位姿的匹配度 */
        return bestScore;
    }

    /* 移动得分结构体 */
    struct ScoredMove
    {
        OrientedPoint pose; /* 位姿 */
        double score;       /* 匹配度得分 */
        double likelihood;  /* 似然值 */
    };

    typedef std::list<ScoredMove> ScoredMoveList;

    double ScanMatcher::optimize(OrientedPoint &_mean, ScanMatcher::CovarianceMatrix &_cov, const ScanMatcherMap &map, const OrientedPoint &init, const double *readings) const
    {
        ScoredMoveList moveList;
        double bestScore = -1;
        OrientedPoint currentPose = init;
        ScoredMove sm = {currentPose, 0, 0};
        unsigned int matched = likelihoodAndScore(sm.score, sm.likelihood, map, currentPose, readings);
        double currentScore = sm.score;
        moveList.push_back(sm);
        double adelta = m_optAngularDelta, ldelta = m_optLinearDelta;
        unsigned int refinement = 0;
        int count = 0;
        enum Move
        {
            Front,
            Back,
            Left,
            Right,
            TurnLeft,
            TurnRight,
            Done
        };
        do
        {
            if (bestScore >= currentScore)
            {
                refinement++;
                adelta *= .5;
                ldelta *= .5;
            }
            bestScore = currentScore;
            // cout <<"score="<< currentScore << " refinement=" << refinement;
            // cout <<  "pose=" << currentPose.x  << " " << currentPose.y << " " << currentPose.theta << endl;
            OrientedPoint bestLocalPose = currentPose;
            OrientedPoint localPose = currentPose;

            Move move = Front;
            do
            {
                localPose = currentPose;
                switch (move)
                {
                case Front:
                    localPose.x += ldelta;
                    move = Back;
                    break;
                case Back:
                    localPose.x -= ldelta;
                    move = Left;
                    break;
                case Left:
                    localPose.y -= ldelta;
                    move = Right;
                    break;
                case Right:
                    localPose.y += ldelta;
                    move = TurnLeft;
                    break;
                case TurnLeft:
                    localPose.theta += adelta;
                    move = TurnRight;
                    break;
                case TurnRight:
                    localPose.theta -= adelta;
                    move = Done;
                    break;
                default:;
                }
                double localScore, localLikelihood;

                double odo_gain = 1;
                if (m_angularOdometryReliability > 0.)
                {
                    double dth = init.theta - localPose.theta;
                    dth = atan2(sin(dth), cos(dth));
                    dth *= dth;
                    odo_gain *= exp(-m_angularOdometryReliability * dth);
                }
                if (m_linearOdometryReliability > 0.)
                {
                    double dx = init.x - localPose.x;
                    double dy = init.y - localPose.y;
                    double drho = dx * dx + dy * dy;
                    odo_gain *= exp(-m_linearOdometryReliability * drho);
                }
                localScore = odo_gain * score(map, localPose, readings);
                // update the score
                count++;
                matched = likelihoodAndScore(localScore, localLikelihood, map, localPose, readings);
                if (localScore > currentScore)
                {
                    currentScore = localScore;
                    bestLocalPose = localPose;
                }
                sm.score = localScore;
                sm.likelihood = localLikelihood; //+log(odo_gain);
                sm.pose = localPose;
                moveList.push_back(sm);
                // update the move list
            } while (move != Done);
            currentPose = bestLocalPose;
            // cout << __func__ << "currentScore=" << currentScore<< endl;
            // here we look for the best move;
        } while (currentScore > bestScore || refinement < m_optRecursiveIterations);
        // cout << __func__ << "bestScore=" << bestScore<< endl;
        // cout << __func__ << "iterations=" << count<< endl;

        // normalize the likelihood
        double lmin = 1e9;
        double lmax = -1e9;
        for (ScoredMoveList::const_iterator it = moveList.begin(); it != moveList.end(); it++)
        {
            lmin = it->likelihood < lmin ? it->likelihood : lmin;
            lmax = it->likelihood > lmax ? it->likelihood : lmax;
        }
        // cout << "lmin=" << lmin << " lmax=" << lmax<< endl;
        for (ScoredMoveList::iterator it = moveList.begin(); it != moveList.end(); it++)
        {
            it->likelihood = exp(it->likelihood - lmax);
            // cout << "l=" << it->likelihood << endl;
        }
        // compute the mean
        OrientedPoint mean(0, 0, 0);
        double lacc = 0;
        for (ScoredMoveList::const_iterator it = moveList.begin(); it != moveList.end(); it++)
        {
            mean = mean + it->pose * it->likelihood;
            lacc += it->likelihood;
        }
        mean = mean * (1. / lacc);
        // OrientedPoint delta=mean-currentPose;
        // cout << "delta.x=" << delta.x << " delta.y=" << delta.y << " delta.theta=" << delta.theta << endl;
        CovarianceMatrix cov = {0., 0., 0., 0., 0., 0.};
        for (ScoredMoveList::const_iterator it = moveList.begin(); it != moveList.end(); it++)
        {
            OrientedPoint delta = it->pose - mean;
            delta.theta = atan2(sin(delta.theta), cos(delta.theta));
            cov.xx += delta.x * delta.x * it->likelihood;
            cov.yy += delta.y * delta.y * it->likelihood;
            cov.tt += delta.theta * delta.theta * it->likelihood;
            cov.xy += delta.x * delta.y * it->likelihood;
            cov.xt += delta.x * delta.theta * it->likelihood;
            cov.yt += delta.y * delta.theta * it->likelihood;
        }
        cov.xx /= lacc, cov.xy /= lacc, cov.xt /= lacc, cov.yy /= lacc, cov.yt /= lacc, cov.tt /= lacc;

        _mean = currentPose;
        _cov = cov;
        return bestScore;
    }

    /**
     * @brief 扫描匹配器激光数据设置函数
     *
     * @param beams 波束数量
     * @param angles 波束角度
     * @param lpose 传感器位置
     */
    void ScanMatcher::setLaserParameters(unsigned int beams, double *angles, const OrientedPoint &lpose)
    {
        /*if (m_laserAngles)
            delete [] m_laserAngles;
        */
        assert(beams < LASER_MAXBEAMS);
        m_laserPose = lpose;
        m_laserBeams = beams;
        // m_laserAngles=new double[beams];
        memcpy(m_laserAngles, angles, sizeof(double) * m_laserBeams);
    }

    double ScanMatcher::likelihood(double &_lmax, OrientedPoint &_mean, CovarianceMatrix &_cov, const ScanMatcherMap &map, const OrientedPoint &p, const double *readings)
    {
        ScoredMoveList moveList;

        for (double xx = -m_llsamplerange; xx <= m_llsamplerange; xx += m_llsamplestep)
            for (double yy = -m_llsamplerange; yy <= m_llsamplerange; yy += m_llsamplestep)
                for (double tt = -m_lasamplerange; tt <= m_lasamplerange; tt += m_lasamplestep)
                {

                    OrientedPoint rp = p;
                    rp.x += xx;
                    rp.y += yy;
                    rp.theta += tt;

                    ScoredMove sm;
                    sm.pose = rp;

                    likelihoodAndScore(sm.score, sm.likelihood, map, rp, readings);
                    moveList.push_back(sm);
                }

        // OrientedPoint delta=mean-currentPose;
        // cout << "delta.x=" << delta.x << " delta.y=" << delta.y << " delta.theta=" << delta.theta << endl;
        // normalize the likelihood
        double lmax = -1e9;
        double lcum = 0;
        for (ScoredMoveList::const_iterator it = moveList.begin(); it != moveList.end(); it++)
        {
            lmax = it->likelihood > lmax ? it->likelihood : lmax;
        }
        for (ScoredMoveList::iterator it = moveList.begin(); it != moveList.end(); it++)
        {
            // it->likelihood=exp(it->likelihood-lmax);
            lcum += exp(it->likelihood - lmax);
            it->likelihood = exp(it->likelihood - lmax);
            // cout << "l=" << it->likelihood << endl;
        }

        OrientedPoint mean(0, 0, 0);
        double s = 0, c = 0;
        for (ScoredMoveList::const_iterator it = moveList.begin(); it != moveList.end(); it++)
        {
            mean = mean + it->pose * it->likelihood;
            s += it->likelihood * sin(it->pose.theta);
            c += it->likelihood * cos(it->pose.theta);
        }
        mean = mean * (1. / lcum);
        s /= lcum;
        c /= lcum;
        mean.theta = atan2(s, c);

        CovarianceMatrix cov = {0., 0., 0., 0., 0., 0.};
        for (ScoredMoveList::const_iterator it = moveList.begin(); it != moveList.end(); it++)
        {
            OrientedPoint delta = it->pose - mean;
            delta.theta = atan2(sin(delta.theta), cos(delta.theta));
            cov.xx += delta.x * delta.x * it->likelihood;
            cov.yy += delta.y * delta.y * it->likelihood;
            cov.tt += delta.theta * delta.theta * it->likelihood;
            cov.xy += delta.x * delta.y * it->likelihood;
            cov.xt += delta.x * delta.theta * it->likelihood;
            cov.yt += delta.y * delta.theta * it->likelihood;
        }
        cov.xx /= lcum, cov.xy /= lcum, cov.xt /= lcum, cov.yy /= lcum, cov.yt /= lcum, cov.tt /= lcum;

        _mean = mean;
        _cov = cov;
        _lmax = lmax;
        return log(lcum) + lmax;
    }

    double ScanMatcher::likelihood(double &_lmax, OrientedPoint &_mean, CovarianceMatrix &_cov, const ScanMatcherMap &map, const OrientedPoint &p,
                                   Gaussian3 &odometry, const double *readings, double gain)
    {
        ScoredMoveList moveList;

        for (double xx = -m_llsamplerange; xx <= m_llsamplerange; xx += m_llsamplestep)
            for (double yy = -m_llsamplerange; yy <= m_llsamplerange; yy += m_llsamplestep)
                for (double tt = -m_lasamplerange; tt <= m_lasamplerange; tt += m_lasamplestep)
                {

                    OrientedPoint rp = p;
                    rp.x += xx;
                    rp.y += yy;
                    rp.theta += tt;

                    ScoredMove sm;
                    sm.pose = rp;

                    likelihoodAndScore(sm.score, sm.likelihood, map, rp, readings);
                    sm.likelihood += odometry.eval(rp) / gain;
                    assert(!isnan(sm.likelihood));
                    moveList.push_back(sm);
                }

        // OrientedPoint delta=mean-currentPose;
        // cout << "delta.x=" << delta.x << " delta.y=" << delta.y << " delta.theta=" << delta.theta << endl;
        // normalize the likelihood
        double lmax = -std::numeric_limits<double>::max();
        double lcum = 0;
        for (ScoredMoveList::const_iterator it = moveList.begin(); it != moveList.end(); it++)
        {
            lmax = it->likelihood > lmax ? it->likelihood : lmax;
        }
        for (ScoredMoveList::iterator it = moveList.begin(); it != moveList.end(); it++)
        {
            // it->likelihood=exp(it->likelihood-lmax);
            lcum += exp(it->likelihood - lmax);
            it->likelihood = exp(it->likelihood - lmax);
            // cout << "l=" << it->likelihood << endl;
        }

        OrientedPoint mean(0, 0, 0);
        double s = 0, c = 0;
        for (ScoredMoveList::const_iterator it = moveList.begin(); it != moveList.end(); it++)
        {
            mean = mean + it->pose * it->likelihood;
            s += it->likelihood * sin(it->pose.theta);
            c += it->likelihood * cos(it->pose.theta);
        }
        mean = mean * (1. / lcum);
        s /= lcum;
        c /= lcum;
        mean.theta = atan2(s, c);

        CovarianceMatrix cov = {0., 0., 0., 0., 0., 0.};
        for (ScoredMoveList::const_iterator it = moveList.begin(); it != moveList.end(); it++)
        {
            OrientedPoint delta = it->pose - mean;
            delta.theta = atan2(sin(delta.theta), cos(delta.theta));
            cov.xx += delta.x * delta.x * it->likelihood;
            cov.yy += delta.y * delta.y * it->likelihood;
            cov.tt += delta.theta * delta.theta * it->likelihood;
            cov.xy += delta.x * delta.y * it->likelihood;
            cov.xt += delta.x * delta.theta * it->likelihood;
            cov.yt += delta.y * delta.theta * it->likelihood;
        }
        cov.xx /= lcum, cov.xy /= lcum, cov.xt /= lcum, cov.yy /= lcum, cov.yt /= lcum, cov.tt /= lcum;

        _mean = mean;
        _cov = cov;
        _lmax = lmax;
        double v = log(lcum) + lmax;
        assert(!isnan(v));
        return v;
    }

    /**
     * @brief 扫描匹配器参数设置函数
     *
     * @param urange 激光雷达最大使用距离
     * @param range 激光雷达最大的量程
     * @param sigma 端点匹配的标准差
     * @param kernsize 内核中要查找对应关系的数量
     * @param lopt 平移优化步长
     * @param aopt 旋转优化步长
     * @param iterations 扫描匹配的迭代步数
     * @param likelihoodSigma 扫描匹配概率的激光标准差
     * @param likelihoodGain 评估可能性时使用的增益
     * @param likelihoodSkip 评估可能性时间隔跳过的激光束数量
     */
    void ScanMatcher::setMatchingParameters(double urange, double range, double sigma, int kernsize, double lopt, double aopt, int iterations, double likelihoodSigma, unsigned int likelihoodSkip)
    {
        m_usableRange = urange;
        m_laserMaxRange = range;
        m_kernelSize = kernsize;
        m_optLinearDelta = lopt;
        m_optAngularDelta = aopt;
        m_optRecursiveIterations = iterations;
        m_gaussianSigma = sigma;
        m_likelihoodSigma = likelihoodSigma;
        m_likelihoodSkip = likelihoodSkip;
    }

};
