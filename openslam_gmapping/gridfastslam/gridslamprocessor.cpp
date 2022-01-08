#include <string>
#include <deque>
#include <list>
#include <map>
#include <set>
#include <fstream>
#include <iomanip>
#include <gmapping/utils/stat.h>
#include "gmapping/gridfastslam/gridslamprocessor.h"

//#define MAP_CONSISTENCY_CHECK
//#define GENERATE_TRAJECTORIES

namespace GMapping
{
    /* 距离阀值 用于判断距离增量是否合理 */
    const double m_distanceThresholdCheck = 20;

    /* 定义标准库函数的标准命名空间 */
    using namespace std;

    /**
     * @brief 构造函数
     */
    GridSlamProcessor::GridSlamProcessor() : m_infoStream(cout)
    {
        period_ = 5.0;             /* 激光雷达数据扫描处理的间隔时间 */
        m_obsSigmaGain = 1;        /* 似然度的平滑因子 */
        m_resampleThreshold = 0.5; /* 重新采样阀值 */
        m_minimumScore = 0.;       /* 扫描匹配结果理想所需的最低分数线 */
    }

    /**
     * @brief 构造函数
     * 
     * @param gsp SLAM处理器对象
     */
    GridSlamProcessor::GridSlamProcessor(const GridSlamProcessor &gsp)
        : last_update_time_(0.0), m_particles(gsp.m_particles), m_infoStream(cout)
    {
        period_ = 5.0; /* 激光雷达的扫描周期 */

        m_obsSigmaGain = gsp.m_obsSigmaGain;
        m_resampleThreshold = gsp.m_resampleThreshold; /* 重采样阀值 */
        m_minimumScore = gsp.m_minimumScore;           /* 扫描匹配结果理想所需的最低分数线 */

        m_beams = gsp.m_beams;                         /* 激光波束数量 */
        m_indexes = gsp.m_indexes;                     /* 重采样后的粒子指标 */
        m_motionModel = gsp.m_motionModel;             /* 运动模型 */
        m_resampleThreshold = gsp.m_resampleThreshold; /* 重采样阀值（重复） */
        m_matcher = gsp.m_matcher;                     /* 扫描匹配器 */

        m_count = gsp.m_count;               /* 扫描处理迭代的计数器 */
        m_readingCount = gsp.m_readingCount; /* 读取激光雷达传感器数据次数 */
        m_lastPartPose = gsp.m_lastPartPose; /* 上一次的位姿 */
        m_pose = gsp.m_pose;
        m_odoPose = gsp.m_odoPose;                 /* 里程的位姿 */
        m_linearDistance = gsp.m_linearDistance;   /* 直线距离 */
        m_angularDistance = gsp.m_angularDistance; /* 角度距离 */
        m_neff = gsp.m_neff;                       /* 衡量粒子权重的相似性，Neff越大，粒子权重差距越小 */

        cerr << "FILTER COPY CONSTRUCTOR" << endl;
        cerr << "m_odoPose=" << m_odoPose.x << " " << m_odoPose.y << " " << m_odoPose.theta << endl;
        cerr << "m_lastPartPose=" << m_lastPartPose.x << " " << m_lastPartPose.y << " " << m_lastPartPose.theta << endl;
        cerr << "m_linearDistance=" << m_linearDistance << endl;
        cerr << "m_angularDistance=" << m_linearDistance << endl;

        m_xmin = gsp.m_xmin;   /* 物理地图x向初始最小尺寸 */
        m_ymin = gsp.m_ymin;   /* 物理地图y向初始最小尺寸 */
        m_xmax = gsp.m_xmax;   /* 物理地图x向初始最大尺寸 */
        m_ymax = gsp.m_ymax;   /* 物理地图y向初始最大尺寸 */
        m_delta = gsp.m_delta; /* 分辨率 */

        m_regScore = gsp.m_regScore;   /* 如果扫描分数高于该阈值，则在地图上注册 */
        m_critScore = gsp.m_critScore; /* 如果扫描分数低于该阈值，则上报扫描匹配失败 */
        m_maxMove = gsp.m_maxMove;     /* 在连续扫描之间允许的最大移动 */

        m_linearThresholdDistance = gsp.m_linearThresholdDistance;   /* 机器每平移该距离后处理一次激光扫描数据 */
        m_angularThresholdDistance = gsp.m_angularThresholdDistance; /* 机器每旋转该弧度后处理一次激光扫描数据 */
        m_obsSigmaGain = gsp.m_obsSigmaGain;                         /* 似然度的平滑因子 */

#ifdef MAP_CONSISTENCY_CHECK
        cerr << __func__ << ": trajectories copy.... ";
#endif
        TNodeVector v = gsp.getTrajectories();
        for (unsigned int i = 0; i < v.size(); i++)
        {
            m_particles[i].node = v[i];
        }
#ifdef MAP_CONSISTENCY_CHECK
        cerr << "end" << endl;
#endif

        cerr << "Tree: normalizing, resetting and propagating weights within copy construction/cloneing ...";
        updateTreeWeights(false);
        cerr << ".done!" << endl;
    }

    /**
     * @brief 构造函数
     *
     * @param infoS 控制台输出
     */
    GridSlamProcessor::GridSlamProcessor(std::ostream &infoS) : m_infoStream(infoS)
    {
        period_ = 5.0;             /* 激光雷达的扫描处理周期 */
        m_obsSigmaGain = 1;        /* 似然度的平滑因子 */
        m_resampleThreshold = 0.5; /* 重新采样阀值 */
        m_minimumScore = 0.;       /* 扫描匹配结果理想所需的最低分数线 */
    }

    GridSlamProcessor *GridSlamProcessor::clone() const
    {
#ifdef MAP_CONSISTENCY_CHECK
        cerr << __func__ << ": performing preclone_fit_test" << endl;
        typedef std::map<autoptr<Array2D<PointAccumulator>>::reference *const, int> PointerMap;
        PointerMap pmap;
        for (ParticleVector::const_iterator it = m_particles.begin(); it != m_particles.end(); it++)
        {
            const ScanMatcherMap &m1(it->map);
            const HierarchicalArray2D<PointAccumulator> &h1(m1.storage());
            for (int x = 0; x < h1.getXSize(); x++)
            {
                for (int y = 0; y < h1.getYSize(); y++)
                {
                    const autoptr<Array2D<PointAccumulator>> &a1(h1.m_cells[x][y]);
                    if (a1.m_reference)
                    {
                        PointerMap::iterator f = pmap.find(a1.m_reference);
                        if (f == pmap.end())
                            pmap.insert(make_pair(a1.m_reference, 1));
                        else
                            f->second++;
                    }
                }
            }
        }
        cerr << __func__ << ": Number of allocated chunks" << pmap.size() << endl;
        for (PointerMap::const_iterator it = pmap.begin(); it != pmap.end(); it++)
            assert(it->first->shares == (unsigned int)it->second);

        cerr << __func__ << ": SUCCESS, the error is somewhere else" << endl;
#endif
        GridSlamProcessor *cloned = new GridSlamProcessor(*this);

#ifdef MAP_CONSISTENCY_CHECK
        cerr << __func__ << ": trajectories end" << endl;
        cerr << __func__ << ": performing afterclone_fit_test" << endl;
        ParticleVector::const_iterator jt = cloned->m_particles.begin();
        for (ParticleVector::const_iterator it = m_particles.begin(); it != m_particles.end(); it++)
        {
            const ScanMatcherMap &m1(it->map);
            const ScanMatcherMap &m2(jt->map);
            const HierarchicalArray2D<PointAccumulator> &h1(m1.storage());
            const HierarchicalArray2D<PointAccumulator> &h2(m2.storage());
            jt++;
            for (int x = 0; x < h1.getXSize(); x++)
            {
                for (int y = 0; y < h1.getYSize(); y++)
                {
                    const autoptr<Array2D<PointAccumulator>> &a1(h1.m_cells[x][y]);
                    const autoptr<Array2D<PointAccumulator>> &a2(h2.m_cells[x][y]);
                    assert(a1.m_reference == a2.m_reference);
                    assert((!a1.m_reference) || !(a1.m_reference->shares % 2));
                }
            }
        }
        cerr << __func__ << ": SUCCESS, the error is somewhere else" << endl;
#endif
        return cloned;
    }

    GridSlamProcessor::~GridSlamProcessor()
    {
        cerr << __func__ << ": Start" << endl;
        cerr << __func__ << ": Deleting tree" << endl;
        for (std::vector<Particle>::iterator it = m_particles.begin(); it != m_particles.end(); it++)
        {
#ifdef TREE_CONSISTENCY_CHECK
            TNode *node = it->node;
            while (node)
                node = node->parent;
            cerr << "@" << endl;
#endif
            if (it->node)
                delete it->node;
            // cout << "l=" << it->weight<< endl;
        }

#ifdef MAP_CONSISTENCY_CHECK
        cerr << __func__ << ": performing predestruction_fit_test" << endl;
        typedef std::map<autoptr<Array2D<PointAccumulator>>::reference *const, int> PointerMap;
        PointerMap pmap;
        for (ParticleVector::const_iterator it = m_particles.begin(); it != m_particles.end(); it++)
        {
            const ScanMatcherMap &m1(it->map);
            const HierarchicalArray2D<PointAccumulator> &h1(m1.storage());
            for (int x = 0; x < h1.getXSize(); x++)
            {
                for (int y = 0; y < h1.getYSize(); y++)
                {
                    const autoptr<Array2D<PointAccumulator>> &a1(h1.m_cells[x][y]);
                    if (a1.m_reference)
                    {
                        PointerMap::iterator f = pmap.find(a1.m_reference);
                        if (f == pmap.end())
                            pmap.insert(make_pair(a1.m_reference, 1));
                        else
                            f->second++;
                    }
                }
            }
        }
        cerr << __func__ << ": Number of allocated chunks" << pmap.size() << endl;
        for (PointerMap::const_iterator it = pmap.begin(); it != pmap.end(); it++)
            assert(it->first->shares >= (unsigned int)it->second);
        cerr << __func__ << ": SUCCESS, the error is somewhere else" << endl;
#endif
    }

    /**
     * @brief 匹配参数设置函数
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
    void GridSlamProcessor::setMatchingParameters(double urange, double range, double sigma, int kernsize, double lopt, double aopt,
                                                  int iterations, double likelihoodSigma, double likelihoodGain, unsigned int likelihoodSkip)
    {
        m_obsSigmaGain = likelihoodGain;
        m_matcher.setMatchingParameters(urange, range, sigma, kernsize, lopt, aopt, iterations, likelihoodSigma, likelihoodSkip);
        if (m_infoStream)
            m_infoStream << " -maxUrange " << urange
                         << " -maxUrange " << range
                         << " -sigma     " << sigma
                         << " -kernelSize " << kernsize
                         << " -lstep " << lopt
                         << " -lobsGain " << m_obsSigmaGain
                         << " -astep " << aopt << endl;
    }

    /**
     * @brief 运动参数设置函数
     *
     * @param srr 平移时平移里程误差
     * @param srt 平移时旋转里程误差
     * @param str 旋转时平移里程误差
     * @param stt 旋转时旋转里程误差
     */
    void GridSlamProcessor::setMotionModelParameters(double srr, double srt, double str, double stt)
    {
        m_motionModel.srr = srr;
        m_motionModel.srt = srt;
        m_motionModel.str = str;
        m_motionModel.stt = stt;

        if (m_infoStream)
            m_infoStream << " -srr " << srr << " -srt " << srt
                         << " -str " << str << " -stt " << stt << endl;
    }

    /**
     * @brief 更新距离设置函数
     *
     * @param linear 机器每平移该距离后处理一次激光扫描数据
     * @param angular 机器每旋转该弧度后处理一次激光扫描数据
     * @param resampleThreshold 重新采样阀值 如果最新扫描处理的速度比更新速度慢，则处理扫描；小于零的值将关闭基于时间的更新
     */
    void GridSlamProcessor::setUpdateDistances(double linear, double angular, double resampleThreshold)
    {
        m_linearThresholdDistance = linear;
        m_angularThresholdDistance = angular;
        m_resampleThreshold = resampleThreshold;
        if (m_infoStream)
            m_infoStream << " -linearUpdate " << linear
                         << " -angularUpdate " << angular
                         << " -resampleThreshold " << m_resampleThreshold << endl;
    }

    // HERE STARTS THE BEEF

    /**
     * @brief 粒子构建函数
     *
     * @param m 激光地图对象
     */
    GridSlamProcessor::Particle::Particle(const ScanMatcherMap &m) : map(m), pose(0, 0, 0), weight(0), weightSum(0), gweight(0), previousIndex(0)
    {
        node = 0;
    }

    /**
     * @brief 激光传感器配置函数
     *
     * @param smap 传感器关联容器
     */
    void GridSlamProcessor::setSensorMap(const SensorMap &smap)
    {

        /*
          Construct the angle table for the sensor
          FIXME For now detect the readings of only the front laser, and assume its pose is in the center of the robot
        */
        /* 构建传感器角度表，假设"FLASER"激光雷达对机器人前面进行扫描，并且安装在机器人的中心 */

        /* 获取传感器"FLASER"，如果map中没有相应的传感器就报错，即只支持"FLASER"激光雷达 */
        SensorMap::const_iterator laser_it = smap.find(std::string("FLASER"));
        if (laser_it == smap.end())
        {
            cerr << "Attempting to load the new carmen log format" << endl;
            laser_it = smap.find(std::string("ROBOTLASER1"));
            assert(laser_it != smap.end());
        }

        /* 获取实际的激光雷达传感器对象（动态强制类型转换），并检查对象存在并且激光波束数量不为0 */
        const RangeSensor *rangeSensor = dynamic_cast<const RangeSensor *>((laser_it->second));
        assert(rangeSensor && rangeSensor->beams().size());

        /* 获取传感器的波束数量 */
        m_beams = static_cast<unsigned int>(rangeSensor->beams().size());

        /* 记录各个波束对应的角度 */
        double *angles = new double[rangeSensor->beams().size()];
        for (unsigned int i = 0; i < m_beams; i++)
        {
            angles[i] = rangeSensor->beams()[i].pose.theta;
        }

        /* 根据波束数量、波束角度和传感器位置配置扫描匹配器的激光参数 */
        m_matcher.setLaserParameters(m_beams, angles, rangeSensor->getPose());

        delete[] angles; /* 删除临时波束对应角度记录 */
    }

    /**
     * @brief SLAM处理器初始化函数
     *
     * @param size 初始化粒子个数
     * @param xmin 物理地图x向初始最小尺寸
     * @param ymin 物理地图y向初始最小尺寸
     * @param xmax 物理地图x向初始最大尺寸
     * @param ymax 物理地图y向初始最大尺寸
     * @param delta 分辨率
     * @param initialPose 建图初始位姿
     */
    void GridSlamProcessor::init(unsigned int size, double xmin, double ymin, double xmax, double ymax, double delta, OrientedPoint initialPose)
    {
        /* 记录下关于地图的配置信息 */
        m_xmin = xmin;
        m_ymin = ymin;
        m_xmax = xmax;
        m_ymax = ymax;
        m_delta = delta;

        /* 显示下关于粒子和地图的配置信息 */
        if (m_infoStream)
            m_infoStream
                << " -xmin " << m_xmin
                << " -xmax " << m_xmax
                << " -ymin " << m_ymin
                << " -ymax " << m_ymax
                << " -delta " << m_delta
                << " -particles " << size << endl;

        /* 清空粒子集合 */
        m_particles.clear();

        /* 创建了一个粒子运动轨迹对象node和地图对象lmap */
        TNode *node = new TNode(initialPose, 0, 0, 0);
        ScanMatcherMap lmap(Point(xmin + xmax, ymin + ymax) * .5, xmax - xmin, ymax - ymin, delta);

        /* 粒子初始化操作 */
        for (unsigned int i = 0; i < size; i++)
        {
            /* 为每个粒子赋予了刚刚定义的地图对象和初始位置，并设定粒子权重为0 */
            m_particles.push_back(Particle(lmap));
            m_particles.back().pose = initialPose;
            m_particles.back().previousPose = initialPose;
            m_particles.back().setWeight(0);
            m_particles.back().previousIndex = 0;

            // this is not needed
            // m_particles.back().node=new TNode(initialPose, 0, node, 0);

            // we use the root directly
            /* 我们直接使用根结点 */
            m_particles.back().node = node;
        }

        /* 初始化一些状态信息 */
        m_neff = (double)size;                    /* 衡量粒子权重的相似性，用于判断重采样 */
        m_count = 0;                              /* 扫描处理迭代的计数器 */
        m_readingCount = 0;                       /* 读取激光雷达传感器数据次数 */
        m_linearDistance = m_angularDistance = 0; /* 直线和角度距离 */
    }

    void GridSlamProcessor::processTruePos(const OdometryReading &o)
    {
        const OdometrySensor *os = dynamic_cast<const OdometrySensor *>(o.getSensor());
        if (os && os->isIdeal() && m_outputStream)
        {
            m_outputStream << setiosflags(ios::fixed) << setprecision(3);
            m_outputStream << "SIMULATOR_POS " << o.getPose().x << " " << o.getPose().y << " ";
            m_outputStream << setiosflags(ios::fixed) << setprecision(6) << o.getPose().theta << " " << o.getTime() << endl;
        }
    }

    /**
     * @brief 扫描数据处理函数
     *
     * 粒子滤波迭代过程
     *
     * @param reading 激光雷达和里程位姿数据
     * @param adaptParticles 在重采样过程中选用的粒子数量
     * @return processed 是否成功更新了建图引擎
     */
    bool GridSlamProcessor::processScan(const RangeReading &reading, int adaptParticles)
    {

        /**retireve the position from the reading, and compute the odometry*/
        /* 从读数中撤消该位置，并计算里程计   */

        /* 通过传感器读数的接口获取最新里程计坐标系下激光雷达传感器的位姿 此时已依靠理论里程更新了位姿 */
        OrientedPoint relPose = reading.getPose();

        /* 第一次扫描处理时，对上一次和里程计起始的位姿进行初始化 */
        if (!m_count)
        {
            m_lastPartPose = m_odoPose = relPose;
        }

        // write the state of the reading and update all the particles using the motion model
        /* 写入读取的状态，并使用运动模型更新所有的粒子 */

        /* 1.通过运动模型进行预测 */
        /* 遍历粒子群 */
        for (ParticleVector::iterator it = m_particles.begin(); it != m_particles.end(); it++)
        {
            /* 根据运动模型，预测更新最新时刻的粒子群的位姿 */
            OrientedPoint &pose(it->pose);
            pose = m_motionModel.drawFromMotion(it->pose, relPose, m_odoPose);
        }

        // update the output file
        /* 更新输出文件 */
        if (m_outputStream.is_open())
        {
            m_outputStream << setiosflags(ios::fixed) << setprecision(6);
            m_outputStream << "ODOM ";
            m_outputStream << setiosflags(ios::fixed) << setprecision(3) << m_odoPose.x << " " << m_odoPose.y << " ";
            m_outputStream << setiosflags(ios::fixed) << setprecision(6) << m_odoPose.theta << " ";
            m_outputStream << reading.getTime();
            m_outputStream << endl;
        }
        if (m_outputStream.is_open())
        {
            m_outputStream << setiosflags(ios::fixed) << setprecision(6);
            m_outputStream << "ODO_UPDATE " << m_particles.size() << " ";
            for (ParticleVector::iterator it = m_particles.begin(); it != m_particles.end(); it++)
            {
                OrientedPoint &pose(it->pose);
                m_outputStream << setiosflags(ios::fixed) << setprecision(3) << pose.x << " " << pose.y << " ";
                m_outputStream << setiosflags(ios::fixed) << setprecision(6) << pose.theta << " " << it->weight << " ";
            }
            m_outputStream << reading.getTime();
            m_outputStream << endl;
        }

        // invoke the callback
        /* 调用回调 */
        /* 处理里程计更新事件 实际上源码中这个函数什么也没有做，它存在的意义就是提供一种回调机制，方便日后的扩展 */
        onOdometryUpdate();

        // accumulate the robot translation and rotation
        /* 累积机器人的平移和旋转 */
        /* 通过对最新数据中的位置坐标与上次机器人的位置坐标做差，获得机器人的位移向量move */
        OrientedPoint move = relPose - m_odoPose;
        /* 通过三角函数运算，对转向做了调整， 使其处于区间(−π,π]，更新直线和转角距离 */
        move.theta = atan2(sin(move.theta), cos(move.theta));
        m_linearDistance += sqrt(move * move);
        m_angularDistance += fabs(move.theta);

        // if the robot jumps throw a warning
        /* 如果机器人位置暴增，就发出警告 */
        /* 机器人的速度是有限的，所以当机器人位置发生较大变化的时候，很可能系统出现了某些未知的错误 */
        if (m_linearDistance > m_distanceThresholdCheck)
        {
            cerr << "***********************************************************************" << endl;
            cerr << "********** Error: m_distanceThresholdCheck overridden!!!! *************" << endl;
            cerr << "m_distanceThresholdCheck=" << m_distanceThresholdCheck << endl;
            cerr << "Old Odometry Pose= " << m_odoPose.x << " " << m_odoPose.y
                 << " " << m_odoPose.theta << endl;
            cerr << "New Odometry Pose (reported from observation)= " << relPose.x << " " << relPose.y
                 << " " << relPose.theta << endl;
            cerr << "***********************************************************************" << endl;
            cerr << "** The Odometry has a big jump here. This is probably a bug in the   **" << endl;
            cerr << "** odometry/laser input. We continue now, but the result is probably **" << endl;
            cerr << "** crap or can lead to a core dump since the map doesn't fit.... C&G **" << endl;
            cerr << "***********************************************************************" << endl;
        }

        /* 更新里程计位置 */
        m_odoPose = relPose;

        /* 标记是否对当前的激光雷达传感器数据处理，将用作本函数的返回值 */
        bool processed = false;

        // process a scan only if the robot has traveled a given distance or a certain amount of time has elapsed
        /* 只有当机器人移动了给定的距离或经过了一定的时间后才进行扫描 */

        /* 初次扫描处理更新，或者机器人有显著的线距离和角距离变化，或者时间上有显著距离的时候，处理激光扫描数据 */
        if (!m_count || m_linearDistance >= m_linearThresholdDistance || m_angularDistance >= m_angularThresholdDistance || (period_ >= 0.0 && (reading.getTime() - last_update_time_) > period_))
        {
            /* 记录下当前传感器数据的获取时间 */
            last_update_time_ = reading.getTime();

            if (m_outputStream.is_open())
            {
                m_outputStream << setiosflags(ios::fixed) << setprecision(6);
                m_outputStream << "FRAME " << m_readingCount;
                m_outputStream << " " << m_linearDistance;
                m_outputStream << " " << m_angularDistance << endl;
            }

            if (m_infoStream)
                m_infoStream << "update frame " << m_readingCount << endl
                             << "update ld=" << m_linearDistance << " ad=" << m_angularDistance << endl;

            cerr << "Laser Pose= " << reading.getPose().x << " " << reading.getPose().y
                 << " " << reading.getPose().theta << endl;

            // this is for converting the reading in a scan-matcher feedable form
            /* 这是为了将读数转换成符合扫描匹配器的形式 */

            /* 断言传感器的激光波束数量与设置的一致 */
            assert(reading.size() == m_beams);

            /* 创建和赋值符合扫描匹配器运算的激光数据 */
            double *plainReading = new double[m_beams];
            for (unsigned int i = 0; i < m_beams; i++)
            {
                plainReading[i] = reading[i];
            }

            /* 显示当扫描处理前计数，从0开始 */
            m_infoStream << "m_count " << m_count << endl;

            /* 拷贝一份激光雷达传感器数据 */
            RangeReading *reading_copy =
                new RangeReading(reading.size(),
                                 &(reading[0]),
                                 static_cast<const RangeSensor *>(reading.getSensor()),
                                 reading.getTime());

            /* 非首次的情况 */
            if (m_count > 0)
            {
                /* 2.进行扫描匹配 */
                scanMatch(plainReading);

                /* 一些匹配后最优位姿日志输出 */
                if (m_outputStream.is_open())
                {
                    m_outputStream << "LASER_READING " << reading.size() << " ";
                    m_outputStream << setiosflags(ios::fixed) << setprecision(2);
                    for (RangeReading::const_iterator b = reading.begin(); b != reading.end(); b++)
                    {
                        m_outputStream << *b << " ";
                    }
                    OrientedPoint p = reading.getPose();
                    m_outputStream << setiosflags(ios::fixed) << setprecision(6);
                    m_outputStream << p.x << " " << p.y << " " << p.theta << " " << reading.getTime() << endl;
                    m_outputStream << "SM_UPDATE " << m_particles.size() << " ";
                    for (ParticleVector::const_iterator it = m_particles.begin(); it != m_particles.end(); it++)
                    {
                        const OrientedPoint &pose = it->pose;
                        m_outputStream << setiosflags(ios::fixed) << setprecision(3) << pose.x << " " << pose.y << " ";
                        m_outputStream << setiosflags(ios::fixed) << setprecision(6) << pose.theta << " " << it->weight << " ";
                    }
                    m_outputStream << endl;
                }

                /* 源码为空函数，提供了扩展的可能，和onOdometryUpdate一样 */
                onScanmatchUpdate();

                /* 3.更新粒子轨迹权重 */
                updateTreeWeights(false);

                /* 一些衡量粒子权重的相似性日志输出 */
                if (m_infoStream)
                {
                    m_infoStream << "neff= " << m_neff << endl;
                }
                if (m_outputStream.is_open())
                {
                    m_outputStream << setiosflags(ios::fixed) << setprecision(6);
                    m_outputStream << "NEFF " << m_neff << endl;
                }

                /* 4.重采样 */
                resample(plainReading, adaptParticles, reading_copy);
            }
            else
            {
                /* 首次更新 */
                m_infoStream << "Registering First Scan" << endl;

                /* 对每个粒子进行初始化 */
                for (ParticleVector::iterator it = m_particles.begin(); it != m_particles.end(); it++)
                {
                    /* 设置有效面积未完成计算标志 */
                    m_matcher.invalidateActiveArea();

                    /* 更新地图有效区域 */
                    m_matcher.computeActiveArea(it->map, it->pose, plainReading);

                    /* 将激光扫描的占用信息注册到栅格地图上 */
                    m_matcher.registerScan(it->map, it->pose, plainReading);

                    // cyr: not needed anymore, particles refer to the root in the beginning! 不再需要了，粒子指的是开头的词根!
                    TNode *node = new TNode(it->pose, 0., it->node, 0);
                    // node->reading=0;
                    node->reading = reading_copy;
                    it->node = node;
                }
            }
            // cerr  << "Tree: normalizing, resetting and propagating weights at the end..." ; 最后重新归一化，重置和传播权值
            /* 再次更新粒子权重 */
            updateTreeWeights(false);
            // cerr  << ".done!" <<endl;

            delete[] plainReading;      /* 删除临时用于扫描匹配器运算的激光数据 */
            m_lastPartPose = m_odoPose; // update the past pose for the next iteration 为下一次迭代更新上一次的位姿
            m_linearDistance = 0;       /* 累积直线距离清零 */
            m_angularDistance = 0;      /* 累积角度距离清零 */
            m_count++;                  /* 扫描处理前计数+1 */
            processed = true;           /* 标记当前当前的激光雷达传感器数据已处理 */

            // keep ready for the next step 为下一步做好准备
            /* 遍历粒子群 */
            for (ParticleVector::iterator it = m_particles.begin(); it != m_particles.end(); it++)
            {
                /* 先前的位姿等于当前位姿 */
                it->previousPose = it->pose;
            }
        }

        if (m_outputStream.is_open())
            m_outputStream << flush;

        m_readingCount++; /* 读取激光雷达传感器数据次数+1 */

        return processed; /* 返回是否处理数据标志，可能只读取不处理的 */
    }

    std::ofstream &GridSlamProcessor::outputStream()
    {
        return m_outputStream;
    }

    std::ostream &GridSlamProcessor::infoStream()
    {
        return m_infoStream;
    }

    /**
     * @brief 获取最优粒子的索引函数
     */
    int GridSlamProcessor::getBestParticleIndex() const
    {
        unsigned int bi = 0;
        double bw = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < m_particles.size(); i++)
            if (bw < m_particles[i].weightSum)
            {
                bw = m_particles[i].weightSum;
                bi = i;
            }
        return (int)bi;
    }

    /* 匹配器、重采样、里程更新扩展函数 */
    void GridSlamProcessor::onScanmatchUpdate() {}
    void GridSlamProcessor::onResampleUpdate() {}
    void GridSlamProcessor::onOdometryUpdate() {}

}; // end namespace
