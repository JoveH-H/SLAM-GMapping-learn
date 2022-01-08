#ifndef GRIDSLAMPROCESSOR_H
#define GRIDSLAMPROCESSOR_H

#include <climits>
#include <limits>
#include <fstream>
#include <vector>
#include <deque>
#include <gmapping/particlefilter/particlefilter.h>
#include <gmapping/utils/point.h>
#include <gmapping/utils/macro_params.h>
#include <gmapping/log/sensorlog.h>
#include <gmapping/sensor/sensor_range/rangesensor.h>
#include <gmapping/sensor/sensor_range/rangereading.h>
#include <gmapping/scanmatcher/scanmatcher.h>
#include "gmapping/gridfastslam/motionmodel.h"
#include <gmapping/gridfastslam/gridfastslam_export.h>

namespace GMapping
{

    /**This class defines the basic GridFastSLAM algorithm.  It
       implements a rao blackwellized particle filter. Each particle
       has its own map and robot pose.<br> This implementation works
       as follows: each time a new pair odometry/laser reading is
       received, the particle's robot pose is updated according to the
       motion model.  This pose is subsequently used for initalizing a
       scan matching algorithm.  The scanmatcher performs a local
       optimization for each particle.  It is initialized with the
       pose drawn from the motion model, and the pose is corrected
       according to the each particle map.<br>
       In order to avoid unnecessary computation the filter state is updated
       only when the robot moves more than a given threshold.
    */
    class GRIDFASTSLAM_EXPORT GridSlamProcessor
    {
    public:
        /**This class defines the the node of reversed tree in which the trajectories are stored.
           Each node of a tree has a pointer to its parent and a counter indicating the number of childs of a node.
           The tree is updated in a way consistent with the operation performed on the particles.
        */
        /* 定义了存储轨迹的反向树的节点。树的每个节点都有一个指向其父节点的指针和一个指示节点的子节点数量的计数器。 树方式更新操作与粒子的一致。 */
        struct TNode
        {
            /**Constructs a node of the trajectory tree.
             @param pose:      the pose of the robot in the trajectory
             @param weight:    the weight of the particle at that point in the trajectory
             @param accWeight: the cumulative weight of the particle
             @param parent:    the parent node in the tree
             @param childs:    the number of childs
            */
            /* 构造轨迹树的节点 */
            TNode(const OrientedPoint &pose, double weight, TNode *parent = 0, unsigned int childs = 0);

            /**Destroys a tree node, and consistently updates the tree. If a node whose parent has only one child is deleted,
             also the parent node is deleted. This because the parent will not be reacheable anymore in the trajectory tree.*/
            /* 销毁树节点，并持续更新树。 如果父节点只有一个子节点被删除，那么父节点也会被删除。这是因为父结点在轨迹树中不再存在。 */
            ~TNode();

            /**The pose of the robot*/
            /* 位姿 */
            OrientedPoint pose;

            /**The weight of the particle*/
            /* 权重 */
            double weight;

            /**The sum of all the particle weights in the previous part of the trajectory*/
            /* 所有粒子权重的总和是轨迹的前一节点的权重累加 */
            /* 在轨迹的某一节点的粒子权重总和 */
            double accWeight;

            double gweight;

            /**The parent*/
            /* 父节点 */
            TNode *parent;

            /**The range reading to which this node is associated*/
            /* 此节点关联的读取范围，即记录激光传感器数据  */
            const RangeReading *reading;

            /**The number of childs*/
            /* 子节点数量 */
            unsigned int childs;

            /**counter in visiting the node (internally used)*/
            /* 访问节点的计数器（内部使用） */
            mutable unsigned int visitCounter;

            /**visit flag (internally used)*/
            /* 访问标志（内部使用） */
            mutable bool flag;
        };

        typedef std::vector<GridSlamProcessor::TNode *> TNodeVector;
        typedef std::deque<GridSlamProcessor::TNode *> TNodeDeque;

        /**This class defines a particle of the filter. Each particle has a map, a pose, a weight and retains the current node in the trajectory tree*/
        /* 该类定义了过滤器的粒子。 每个粒子都有一个映射、一个位姿、一个权重，并在轨迹树中保留当前节点 */
        struct Particle
        {
            /** constructs a particle, given a map
                @param map: the particle map
            */
            /* 构造一个粒子，给定一个映射，参数map：粒子的地图 */
            Particle(const ScanMatcherMap &map);

            /** @returns the weight of a particle */
            /* 返回粒子的权重 */
            inline operator double() const { return weight; }

            /** @returns the pose of a particle */
            /* 返回粒子的位姿 */
            inline operator OrientedPoint() const { return pose; }

            /** sets the weight of a particle
                @param w the weight
            */
            /* 设置粒子的权重，参数w：权重 */
            inline void setWeight(double w) { weight = w; }

            /** The map */
            /* 地图 */
            ScanMatcherMap map;

            /** The pose of the robot */
            /* 位姿 */
            OrientedPoint pose;

            /** The pose of the robot at the previous time frame (used for computing thr odometry displacements) */
            /* 机器人在前一时间段的位姿(用于计算里程计位移)   */
            OrientedPoint previousPose;

            /** The weight of the particle */
            /* 权重 */
            double weight;

            /** The cumulative weight of the particle */
            /* 累积权重 */
            double weightSum;

            double gweight;

            /** The index of the previous particle in the trajectory tree */
            /* 轨迹树中前一个粒子的索引 */
            int previousIndex;

            /** Entry to the trajectory tree */
            /* 轨迹树的指针入口 */
            TNode *node;
        };

        typedef std::vector<Particle> ParticleVector;

        /** Constructs a GridSlamProcessor, initialized with the default parameters */
        GridSlamProcessor();

        /** Constructs a GridSlamProcessor, whose output is routed to a stream.
         @param infoStr: the output stream
        */
        GridSlamProcessor(std::ostream &infoStr);

        /** @returns  a deep copy of the grid slam processor with all the internal structures.
         */
        GridSlamProcessor *clone() const;

        /**Deleted the gridslamprocessor*/
        virtual ~GridSlamProcessor();

        // methods for accessing the parameters 访问参数的方法

        /* 激光传感器配置函数 */
        void setSensorMap(const SensorMap &smap);

        /* SLAM初始化函数 */
        void init(unsigned int size, double xmin, double ymin, double xmax, double ymax, double delta,
                  OrientedPoint initialPose = OrientedPoint(0, 0, 0));

        /* 配参数设置函数 */
        void setMatchingParameters(double urange, double range, double sigma, int kernsize, double lopt, double aopt,
                                   int iterations, double likelihoodSigma = 1, double likelihoodGain = 1, unsigned int likelihoodSkip = 0);

        void setMotionModelParameters(double srr, double srt, double str, double stt);    /* 运动参数设置函数 */
        void setUpdateDistances(double linear, double angular, double resampleThreshold); /* 更新距离设置函数 */
        void setUpdatePeriod(double p) { period_ = p; }                                   /* 设置激光雷达数据扫描处理的间隔时间 */

        // the "core" algorithm 核心算法
        void processTruePos(const OdometryReading &odometry);
        bool processScan(const RangeReading &reading, int adaptParticles = 0); /* 扫描数据处理函数 */

        /**This method copies the state of the filter in a tree.
         The tree is represented through reversed pointers (each node has a pointer to its parent).
         The leafs are stored in a vector, whose size is the same as the number of particles.
         @returns the leafs of the tree
        */
        TNodeVector getTrajectories() const;
        void integrateScanSequence(TNode *node);

        /**the scanmatcher algorithm*/
        ScanMatcher m_matcher;
        /**the stream used for writing the output of the algorithm*/
        std::ofstream &outputStream();
        /**the stream used for writing the info/debug messages*/
        std::ostream &infoStream();
        /**@returns the particles*/
        inline const ParticleVector &getParticles() const { return m_particles; }

        inline const std::vector<unsigned int> &getIndexes() const { return m_indexes; }
        int getBestParticleIndex() const; /* 获取最优粒子的索引函数 */

        // callbacks
        /* 里程、重采样、匹配器更新扩展函数 */
        virtual void onOdometryUpdate();
        virtual void onResampleUpdate();
        virtual void onScanmatchUpdate();

        // accessor methods 访问方式
        /**the maxrange of the laser to consider */
        MEMBER_PARAM_SET_GET(m_matcher, double, laserMaxRange, protected, public, public);

        /**the maximum usable range of the laser. A beam is cropped to this value. [scanmatcher]*/
        MEMBER_PARAM_SET_GET(m_matcher, double, usableRange, protected, public, public);

        /**The sigma used by the greedy endpoint matching. [scanmatcher]*/
        MEMBER_PARAM_SET_GET(m_matcher, double, gaussianSigma, protected, public, public);

        /**The sigma  of a beam used for likelihood computation [scanmatcher]*/
        MEMBER_PARAM_SET_GET(m_matcher, double, likelihoodSigma, protected, public, public);

        /**The kernel in which to look for a correspondence[scanmatcher]*/
        MEMBER_PARAM_SET_GET(m_matcher, int, kernelSize, protected, public, public);

        /**The optimization step in rotation [scanmatcher]*/
        MEMBER_PARAM_SET_GET(m_matcher, double, optAngularDelta, protected, public, public);

        /**The optimization step in translation [scanmatcher]*/
        MEMBER_PARAM_SET_GET(m_matcher, double, optLinearDelta, protected, public, public);

        /**The number of iterations of the scanmatcher [scanmatcher]*/
        MEMBER_PARAM_SET_GET(m_matcher, unsigned int, optRecursiveIterations, protected, public, public);

        /**the beams to skip for computing the likelihood (consider a beam every likelihoodSkip) [scanmatcher]*/
        MEMBER_PARAM_SET_GET(m_matcher, unsigned int, likelihoodSkip, protected, public, public);

        /**translational sampling range for the likelihood [scanmatcher]*/
        MEMBER_PARAM_SET_GET(m_matcher, double, llsamplerange, protected, public, public);

        /**angular sampling range for the likelihood [scanmatcher]*/
        MEMBER_PARAM_SET_GET(m_matcher, double, lasamplerange, protected, public, public);

        /**translational sampling range for the likelihood [scanmatcher]*/
        MEMBER_PARAM_SET_GET(m_matcher, double, llsamplestep, protected, public, public);

        /**angular sampling step for the likelihood [scanmatcher]*/
        MEMBER_PARAM_SET_GET(m_matcher, double, lasamplestep, protected, public, public);

        /**generate an accupancy grid map [scanmatcher]*/
        MEMBER_PARAM_SET_GET(m_matcher, bool, generateMap, protected, public, public);

        /**enlarge the map when the robot goes out of the boundaries [scanmatcher]*/
        MEMBER_PARAM_SET_GET(m_matcher, bool, enlargeStep, protected, public, public);

        /**pose of the laser wrt the robot [scanmatcher]*/
        MEMBER_PARAM_SET_GET(m_matcher, OrientedPoint, laserPose, protected, public, public);

        /**odometry error in translation as a function of translation (rho/rho) [motionmodel]*/
        STRUCT_PARAM_SET_GET(m_motionModel, double, srr, protected, public, public);

        /**odometry error in translation as a function of rotation (rho/theta) [motionmodel]*/
        STRUCT_PARAM_SET_GET(m_motionModel, double, srt, protected, public, public);

        /**odometry error in rotation as a function of translation (theta/rho) [motionmodel]*/
        STRUCT_PARAM_SET_GET(m_motionModel, double, str, protected, public, public);

        /**odometry error in  rotation as a function of rotation (theta/theta) [motionmodel]*/
        STRUCT_PARAM_SET_GET(m_motionModel, double, stt, protected, public, public);

        /**minimum score for considering the outcome of the scanmatching good*/
        PARAM_SET_GET(double, minimumScore, protected, public, public);

    protected:
        /**Copy constructor*/
        /* 构造函数 */
        GridSlamProcessor(const GridSlamProcessor &gsp);

        /**the laser beams*/
        unsigned int m_beams;     /* 激光光束数量 */
        double last_update_time_; /* 激光雷达的扫描最新更新时间 */
        double period_;           /* 激光雷达的扫描周期 */

        /**the particles*/
        ParticleVector m_particles; /* 粒子群 */

        /**the particle indexes after resampling (internally used)*/
        std::vector<unsigned int> m_indexes; /* 重采样后的粒子指标（内部用）  */

        /**the particle weights (internally used)*/
        std::vector<double> m_weights; /* 粒子权重（内部用） */

        /**the motion model*/
        MotionModel m_motionModel; /* 运动模型 */

        /**this sets the neff based resampling threshold*/
        /* 设置基于neff的重采样阈值 */
        PARAM_SET_GET(double, resampleThreshold, protected, public, public);

        // state 状态值
        int m_count, m_readingCount;  /* 处理激光雷达传感器数据次数，读取激光雷达传感器数据次数 */
        OrientedPoint m_lastPartPose; /* 上一次的位姿 */
        OrientedPoint m_odoPose;      /* 里程的位姿 */
        OrientedPoint m_pose;
        double m_linearDistance, m_angularDistance; /* 累积平移和角度距离 */

        PARAM_GET(double, neff, protected, public);

        // processing parameters (size of the map)
        PARAM_GET(double, xmin, protected, public);
        PARAM_GET(double, ymin, protected, public);
        PARAM_GET(double, xmax, protected, public);
        PARAM_GET(double, ymax, protected, public);
        // processing parameters (resolution of the map)
        PARAM_GET(double, delta, protected, public);

        // registration score (if a scan score is above this threshold it is registered in the map)
        PARAM_SET_GET(double, regScore, protected, public, public);
        // registration score (if a scan score is below this threshold a scan matching failure is reported)
        PARAM_SET_GET(double, critScore, protected, public, public);
        // registration score maximum move allowed between consecutive scans
        PARAM_SET_GET(double, maxMove, protected, public, public);

        // process a scan each time the robot translates of linearThresholdDistance
        PARAM_SET_GET(double, linearThresholdDistance, protected, public, public);

        // process a scan each time the robot rotates more than angularThresholdDistance
        PARAM_SET_GET(double, angularThresholdDistance, protected, public, public);

        // smoothing factor for the likelihood
        PARAM_SET_GET(double, obsSigmaGain, protected, public, public);

        // stream in which to write the gfs file
        std::ofstream m_outputStream;

        // stream in which to write the messages
        std::ostream &m_infoStream;

        // the functions below performs side effect on the internal structure,
        // should be called only inside the processScan method
    private:
        /**scanmatches all the particles*/
        /* 扫描匹配函数 */
        inline void scanMatch(const double *plainReading);

        /**normalizes the particle weights*/
        /* 标准化粒子重量 */
        inline void normalize();

        // return if a resampling occured or not
        /* 重采样函数 */
        inline bool resample(const double *plainReading, int adaptParticles,
                             const RangeReading *rr = 0);

        // tree utilities 树工具
        void updateTreeWeights(bool weightsAlreadyNormalized = false); /* 更新轨迹权重函数 */
        void resetTree();                                              /* 重置轨迹函数 */
        double propagateWeights();                                     /* 更新轨迹权重函数 */
    };

    typedef std::multimap<const GridSlamProcessor::TNode *, GridSlamProcessor::TNode *> TNodeMultimap;

#include "gmapping/gridfastslam/gridslamprocessor.hxx"

};

#endif
