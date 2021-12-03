/*
 * slam_gmapping
 * Copyright (c) 2008, Willow Garage, Inc.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *   * Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *   * Neither the names of Stanford University or Willow Garage, Inc. nor the names of its
 *     contributors may be used to endorse or promote products derived from
 *     this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */

/* Author: Brian Gerkey */
/* Modified by: Charles DuHadway */

/**

@mainpage slam_gmapping

@htmlinclude manifest.html

@b slam_gmapping is a wrapper around the GMapping SLAM library. It reads laser
scans and odometry and computes a map. This map can be
written to a file using e.g.

  "rosrun map_server map_saver static_map:=dynamic_map"

<hr>

@section topic ROS topics

Subscribes to (name/type):
- @b "scan"/<a href="../../sensor_msgs/html/classstd__msgs_1_1LaserScan.html">sensor_msgs/LaserScan</a> : data from a laser range scanner
- @b "/tf": odometry from the robot


Publishes to (name/type):
- @b "/tf"/tf/tfMessage: position relative to the map


@section services
 - @b "~dynamic_map" : returns the map


@section parameters ROS parameters

Reads the following parameters from the parameter server

Parameters used by our GMapping wrapper:

- @b "~throttle_scans": @b [int] throw away every nth laser scan
- @b "~base_frame": @b [string] the tf frame_id to use for the robot base pose
- @b "~map_frame": @b [string] the tf frame_id where the robot pose on the map is published
- @b "~odom_frame": @b [string] the tf frame_id from which odometry is read
- @b "~map_update_interval": @b [double] time in seconds between two recalculations of the map


Parameters used by GMapping itself:

Laser Parameters:
- @b "~/maxRange" @b [double] maximum range of the laser scans. Rays beyond this range get discarded completely. (default: maximum laser range minus 1 cm, as received in the the first LaserScan message)
- @b "~/maxUrange" @b [double] maximum range of the laser scanner that is used for map building (default: same as maxRange)
- @b "~/sigma" @b [double] standard deviation for the scan matching process (cell)
- @b "~/kernelSize" @b [int] search window for the scan matching process
- @b "~/lstep" @b [double] initial search step for scan matching (linear)
- @b "~/astep" @b [double] initial search step for scan matching (angular)
- @b "~/iterations" @b [int] number of refinement steps in the scan matching. The final "precision" for the match is lstep*2^(-iterations) or astep*2^(-iterations), respectively.
- @b "~/lsigma" @b [double] standard deviation for the scan matching process (single laser beam)
- @b "~/ogain" @b [double] gain for smoothing the likelihood
- @b "~/lskip" @b [int] take only every (n+1)th laser ray for computing a match (0 = take all rays)
- @b "~/minimumScore" @b [double] minimum score for considering the outcome of the scanmatching good. Can avoid 'jumping' pose estimates in large open spaces when using laser scanners with limited range (e.g. 5m). (0 = default. Scores go up to 600+, try 50 for example when experiencing 'jumping' estimate issues)

Motion Model Parameters (all standard deviations of a gaussian noise model)
- @b "~/srr" @b [double] linear noise component (x and y)
- @b "~/stt" @b [double] angular noise component (theta)
- @b "~/srt" @b [double] linear -> angular noise component
- @b "~/str" @b [double] angular -> linear noise component

Others:
- @b "~/linearUpdate" @b [double] the robot only processes new measurements if the robot has moved at least this many meters
- @b "~/angularUpdate" @b [double] the robot only processes new measurements if the robot has turned at least this many rads

- @b "~/resampleThreshold" @b [double] threshold at which the particles get resampled. Higher means more frequent resampling.
- @b "~/particles" @b [int] (fixed) number of particles. Each particle represents a possible trajectory that the robot has traveled

Likelihood sampling (used in scan matching)
- @b "~/llsamplerange" @b [double] linear range
- @b "~/lasamplerange" @b [double] linear step size
- @b "~/llsamplestep" @b [double] linear range
- @b "~/lasamplestep" @b [double] angular step size

Initial map dimensions and resolution:
- @b "~/xmin" @b [double] minimum x position in the map [m]
- @b "~/ymin" @b [double] minimum y position in the map [m]
- @b "~/xmax" @b [double] maximum x position in the map [m]
- @b "~/ymax" @b [double] maximum y position in the map [m]
- @b "~/delta" @b [double] size of one pixel [m]

*/

#include "slam_gmapping.h"

#include <iostream>

#include <time.h>

#include "ros/ros.h"
#include "ros/console.h"
#include "nav_msgs/MapMetaData.h"

#include "gmapping/sensor/sensor_range/rangesensor.h"
#include "gmapping/sensor/sensor_odometry/odometrysensor.h"

#include <rosbag/bag.h>
#include <rosbag/view.h>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// compute linear index for given map coords
/* 计算给定地图坐标的线性索引 */
#define MAP_IDX(sx, i, j) ((sx) * (j) + (i))

/**
 * @brief 构建函数
 */
SlamGMapping::SlamGMapping() : map_to_odom_(tf::Transform(tf::createQuaternionFromRPY(0, 0, 0), tf::Point(0, 0, 0))), /* 默认两个坐标系map和odom重合 */
                               laser_count_(0),
                               private_nh_("~"),
                               scan_filter_sub_(NULL),
                               scan_filter_(NULL),
                               transform_thread_(NULL)
{
    /* 设置高斯噪声的随机数种子和初始化函数 */
    seed_ = time(NULL);
    init();
}

/**
 * @brief 构建函数
 *
 * @param nh ros节点句柄
 * @param pnh 私有ros节点句柄
 */
SlamGMapping::SlamGMapping(ros::NodeHandle &nh, ros::NodeHandle &pnh) : map_to_odom_(tf::Transform(tf::createQuaternionFromRPY(0, 0, 0), tf::Point(0, 0, 0))),
                                                                        laser_count_(0),
                                                                        node_(nh),
                                                                        private_nh_(pnh),
                                                                        scan_filter_sub_(NULL),
                                                                        scan_filter_(NULL),
                                                                        transform_thread_(NULL)
{
    /* 设置高斯噪声的随机数种子和初始化函数 */
    seed_ = time(NULL);
    init();
}

/**
 * @brief 构建函数
 *
 * @param seed 高斯噪声的随机数种子
 * @param max_duration_buffer 最大持续时间缓冲
 */
SlamGMapping::SlamGMapping(long unsigned int seed, long unsigned int max_duration_buffer) : map_to_odom_(tf::Transform(tf::createQuaternionFromRPY(0, 0, 0), tf::Point(0, 0, 0))),
                                                                                            laser_count_(0),
                                                                                            private_nh_("~"),
                                                                                            scan_filter_sub_(NULL),
                                                                                            scan_filter_(NULL),
                                                                                            transform_thread_(NULL),
                                                                                            seed_(seed),
                                                                                            tf_(ros::Duration(max_duration_buffer))
{
    /* 设置初始化函数 */
    init();
}

/**
 * @brief SLAM初始化函数
 */
void SlamGMapping::init()
{
    // log4cxx::Logger::getLogger(ROSCONSOLE_DEFAULT_NAME)->setLevel(ros::console::g_level_lookup[ros::console::levels::Debug]);

    // The library is pretty chatty
    // gsp_ = new GMapping::GridSlamProcessor(std::cerr);
    /* 创建GridSlamSLAM的对象并断言检验 */
    gsp_ = new GMapping::GridSlamProcessor();
    ROS_ASSERT(gsp_);

    /* 创建 map_to_odom的tf变换广播并断言检验 */
    tfB_ = new tf::TransformBroadcaster();
    ROS_ASSERT(tfB_);

    /* 初始化GridSlamSLAM的激光和里程对象 */
    gsp_laser_ = NULL;
    gsp_odom_ = NULL;

    /* 初始化第一帧scan标志位和地图存在标志位 */
    got_first_scan_ = false;
    got_map_ = false;

    // Parameters used by our GMapping wrapper
    /* GMapping 参数设置 */
    if (!private_nh_.getParam("throttle_scans", throttle_scans_))
        throttle_scans_ = 1;
    if (!private_nh_.getParam("base_frame", base_frame_))
        base_frame_ = "base_link";
    if (!private_nh_.getParam("map_frame", map_frame_))
        map_frame_ = "map";
    if (!private_nh_.getParam("odom_frame", odom_frame_))
        odom_frame_ = "odom";

    private_nh_.param("transform_publish_period", transform_publish_period_, 0.05);

    double tmp;
    if (!private_nh_.getParam("map_update_interval", tmp))
        tmp = 5.0;
    map_update_interval_.fromSec(tmp);

    // Parameters used by GMapping itself
    /* 初始默认值，将在initMapper()中设置 */
    maxUrange_ = 0.0;
    maxRange_ = 0.0; // preliminary default, will be set in initMapper()

    if (!private_nh_.getParam("minimumScore", minimum_score_))
        minimum_score_ = 0;
    if (!private_nh_.getParam("sigma", sigma_))
        sigma_ = 0.05;
    if (!private_nh_.getParam("kernelSize", kernelSize_))
        kernelSize_ = 1;
    if (!private_nh_.getParam("lstep", lstep_))
        lstep_ = 0.05;
    if (!private_nh_.getParam("astep", astep_))
        astep_ = 0.05;
    if (!private_nh_.getParam("iterations", iterations_))
        iterations_ = 5;
    if (!private_nh_.getParam("lsigma", lsigma_))
        lsigma_ = 0.075;
    if (!private_nh_.getParam("ogain", ogain_))
        ogain_ = 3.0;
    if (!private_nh_.getParam("lskip", lskip_))
        lskip_ = 0;
    if (!private_nh_.getParam("srr", srr_))
        srr_ = 0.1;
    if (!private_nh_.getParam("srt", srt_))
        srt_ = 0.2;
    if (!private_nh_.getParam("str", str_))
        str_ = 0.1;
    if (!private_nh_.getParam("stt", stt_))
        stt_ = 0.2;
    if (!private_nh_.getParam("linearUpdate", linearUpdate_))
        linearUpdate_ = 1.0;
    if (!private_nh_.getParam("angularUpdate", angularUpdate_))
        angularUpdate_ = 0.5;
    if (!private_nh_.getParam("temporalUpdate", temporalUpdate_))
        temporalUpdate_ = -1.0;
    if (!private_nh_.getParam("resampleThreshold", resampleThreshold_))
        resampleThreshold_ = 0.5;
    if (!private_nh_.getParam("particles", particles_))
        particles_ = 30;
    if (!private_nh_.getParam("xmin", xmin_))
        xmin_ = -100.0;
    if (!private_nh_.getParam("ymin", ymin_))
        ymin_ = -100.0;
    if (!private_nh_.getParam("xmax", xmax_))
        xmax_ = 100.0;
    if (!private_nh_.getParam("ymax", ymax_))
        ymax_ = 100.0;
    if (!private_nh_.getParam("delta", delta_))
        delta_ = 0.05;
    if (!private_nh_.getParam("occ_thresh", occ_thresh_))
        occ_thresh_ = 0.25;
    if (!private_nh_.getParam("llsamplerange", llsamplerange_))
        llsamplerange_ = 0.01;
    if (!private_nh_.getParam("llsamplestep", llsamplestep_))
        llsamplestep_ = 0.01;
    if (!private_nh_.getParam("lasamplerange", lasamplerange_))
        lasamplerange_ = 0.005;
    if (!private_nh_.getParam("lasamplestep", lasamplestep_))
        lasamplestep_ = 0.005;

    if (!private_nh_.getParam("tf_delay", tf_delay_))
        tf_delay_ = transform_publish_period_;
}

/**
 * @brief SLAM启动函数
 */
void SlamGMapping::startLiveSlam()
{
    /* 创建发布者 entropy_publisher_，发布话题 entropy 熵 */
    entropy_publisher_ = private_nh_.advertise<std_msgs::Float64>("entropy", 1, true);

    /* 创建发布者 sst_，发布话题 OccupancyGrid 栅格地图 */
    sst_ = node_.advertise<nav_msgs::OccupancyGrid>("map", 1, true);

    /* 创建发布者 sstm_，发布话题 map_metadata 栅格地图的特征的基本信息 */
    sstm_ = node_.advertise<nav_msgs::MapMetaData>("map_metadata", 1, true);

    /* 创建服务端 ss_，创建服务 dynamic_map 动态地图，回调函数 mapCallback */
    ss_ = node_.advertiseService("dynamic_map", &SlamGMapping::mapCallback, this);

    /* 创建订阅者 scan_filter_sub_， 订阅 scan 消息过滤的激光数据 */
    scan_filter_sub_ = new message_filters::Subscriber<sensor_msgs::LaserScan>(node_, "scan", 5);

    /* tf::MessageFilter 可以订阅任何的ROS消息，然后将其缓存，直到这些消息可以转换到目标坐标系，然后进行相应的处理
       当 message_filters::Subscriber的消息能够由tf转换到目标坐标系时，调用回调函数 */

    /* 创建订阅者 scan_filter_，初始化设置 消息过滤的激光数据订阅器scan_filter_sub_ ，tf转换 tf_，目标坐标系 odom_frame_，等待时间5 */
    scan_filter_ = new tf::MessageFilter<sensor_msgs::LaserScan>(*scan_filter_sub_, tf_, odom_frame_, 5);

    /* 设置订阅者 scan_filter_ 的回调函数 laserCallback */
    scan_filter_->registerCallback(boost::bind(&SlamGMapping::laserCallback, this, _1));

    /* 通过boost::bind开启线程：map_to_odom的TF转换发布循环 publishLoop， 间隔 transform_publish_period_ */
    transform_thread_ = new boost::thread(boost::bind(&SlamGMapping::publishLoop, this, transform_publish_period_));
}

void SlamGMapping::startReplay(const std::string &bag_fname, std::string scan_topic)
{
    double transform_publish_period;
    ros::NodeHandle private_nh_("~");
    entropy_publisher_ = private_nh_.advertise<std_msgs::Float64>("entropy", 1, true);
    sst_ = node_.advertise<nav_msgs::OccupancyGrid>("map", 1, true);
    sstm_ = node_.advertise<nav_msgs::MapMetaData>("map_metadata", 1, true);
    ss_ = node_.advertiseService("dynamic_map", &SlamGMapping::mapCallback, this);

    rosbag::Bag bag;
    bag.open(bag_fname, rosbag::bagmode::Read);

    std::vector<std::string> topics;
    topics.push_back(std::string("/tf"));
    topics.push_back(scan_topic);
    rosbag::View viewall(bag, rosbag::TopicQuery(topics));

    // Store up to 5 messages and there error message (if they cannot be processed right away)
    std::queue<std::pair<sensor_msgs::LaserScan::ConstPtr, std::string>> s_queue;
    foreach (rosbag::MessageInstance const m, viewall)
    {
        tf::tfMessage::ConstPtr cur_tf = m.instantiate<tf::tfMessage>();
        if (cur_tf != NULL)
        {
            for (size_t i = 0; i < cur_tf->transforms.size(); ++i)
            {
                geometry_msgs::TransformStamped transformStamped;
                tf::StampedTransform stampedTf;
                transformStamped = cur_tf->transforms[i];
                tf::transformStampedMsgToTF(transformStamped, stampedTf);
                tf_.setTransform(stampedTf);
            }
        }

        sensor_msgs::LaserScan::ConstPtr s = m.instantiate<sensor_msgs::LaserScan>();
        if (s != NULL)
        {
            if (!(ros::Time(s->header.stamp)).is_zero())
            {
                s_queue.push(std::make_pair(s, ""));
            }
            // Just like in live processing, only process the latest 5 scans
            if (s_queue.size() > 5)
            {
                ROS_WARN_STREAM("Dropping old scan: " << s_queue.front().second);
                s_queue.pop();
            }
            // ignoring un-timestamped tf data
        }

        // Only process a scan if it has tf data
        while (!s_queue.empty())
        {
            try
            {
                tf::StampedTransform t;
                tf_.lookupTransform(s_queue.front().first->header.frame_id, odom_frame_, s_queue.front().first->header.stamp, t);
                this->laserCallback(s_queue.front().first);
                s_queue.pop();
            }
            // If tf does not have the data yet
            catch (tf2::TransformException &e)
            {
                // Store the error to display it if we cannot process the data after some time
                s_queue.front().second = std::string(e.what());
                break;
            }
        }
    }

    bag.close();
}

/**
 * @brief map_to_odom的TF转换发布循环函数
 *
 * @param transform_publish_period tf变换发布间隔
 */
void SlamGMapping::publishLoop(double transform_publish_period)
{
    if (transform_publish_period == 0)
        return;

    /* 计算循环频率 */
    ros::Rate r(1.0 / transform_publish_period);

    /* 循环发布map_to_odom的TF转换 */
    while (ros::ok())
    {
        publishTransform();
        r.sleep();
    }
}

/**
 * @brief 析构函数
 */
SlamGMapping::~SlamGMapping()
{
    /* 删除TF转换发布循环线程 */
    if (transform_thread_)
    {
        transform_thread_->join();
        delete transform_thread_;
    }

    /* 删除GridSlamSLAM相关的对象 */
    delete gsp_;
    if (gsp_laser_)
        delete gsp_laser_;
    if (gsp_odom_)
        delete gsp_odom_;

    /* 删除激光数据相关的对象 */
    if (scan_filter_)
        delete scan_filter_;
    if (scan_filter_sub_)
        delete scan_filter_sub_;
}

/**
 * @brief 获取当前里程计坐标系odom_frame下激光雷达的位姿函数
 *
 * @param gmap_pose 位姿
 * @param t 时间戳
 */
bool SlamGMapping::getOdomPose(GMapping::OrientedPoint &gmap_pose, const ros::Time &t)
{
    // Get the pose of the centered laser at the right time
    /* 更新时间戳 */
    centered_laser_pose_.stamp_ = t;

    // Get the laser's pose that is centered
    /* 让激光的姿势居中 */

    /* odom_pose存储输出的里程计坐标系下激光雷达tf::Stamped<tf::Pose>格式的位姿 */
    tf::Stamped<tf::Transform> odom_pose;
    try
    {
        tf_.transformPose(odom_frame_, centered_laser_pose_, odom_pose);
    }
    catch (tf::TransformException e)
    {
        ROS_WARN("Failed to compute odom pose, skipping scan (%s)", e.what());
        return false;
    }

    /* 获取偏移角 */
    double yaw = tf::getYaw(odom_pose.getRotation());

    /* 通过转化得到OrientedPoint格式的里程计坐标系下的位姿 */
    gmap_pose = GMapping::OrientedPoint(odom_pose.getOrigin().x(),
                                        odom_pose.getOrigin().y(),
                                        yaw);
    return true;
}

/**
 * @brief 初始化地图制作器函数
 *
 * @param scan 激光数据
 */
bool SlamGMapping::initMapper(const sensor_msgs::LaserScan &scan)
{
    /* 获取激光数据的坐标系 */
    laser_frame_ = scan.header.frame_id;

    // Get the laser's pose, relative to base.
    /* 得到激光的位置，相对于基座的位置 */

    /* 定义位姿 ident 变换 laser_pose */
    tf::Stamped<tf::Pose> ident;
    tf::Stamped<tf::Transform> laser_pose;

    /* 设置当前识别的参考坐标系和时间戳  */
    ident.setIdentity();
    ident.frame_id_ = laser_frame_;
    ident.stamp_ = scan.header.stamp;

    /* 进行位姿转换laser_frame -> base_frame，验证姿态 */
    try
    {
        tf_.transformPose(base_frame_, ident, laser_pose);
    }
    catch (tf::TransformException e)
    {
        ROS_WARN("Failed to compute laser pose, aborting initialization (%s)", e.what());
        return false;
    }

    // create a point 1m above the laser position and transform it into the laser-frame
    /* 在激光位置上方创建一个1米的点，并将其转换为 laser_frame， 判断laser方向 */
    tf::Vector3 v;
    v.setValue(0, 0, 1 + laser_pose.getOrigin().z());
    tf::Stamped<tf::Vector3> up(v, scan.header.stamp, base_frame_);
    try
    {
        tf_.transformPoint(laser_frame_, up, up);
        ROS_DEBUG("Z-Axis in sensor frame: %.3f", up.z());
    }
    catch (tf::TransformException &e)
    {
        ROS_WARN("Unable to determine orientation of laser: %s", e.what());
        return false;
    }

    // gmapping doesnt take roll or pitch into account. So check for correct sensor alignment.
    /* 不考虑滚动或倾斜。 所以检查正确的传感器对准 */
    if (fabs(fabs(up.z()) - 1) > 0.001)
    {
        /* 激光雷达必须安装在平面上! z坐标必须是1或-1，但现在:x */
        ROS_WARN("Laser has to be mounted planar! Z-coordinate has to be 1 or -1, but gave: %.5f", up.z());
        return false;
    }

    /* 获取激光束计数 */
    gsp_laser_beam_count_ = scan.ranges.size();

    /* 计算激光束中心角度 */
    double angle_center = (scan.angle_min + scan.angle_max) / 2;

    /* 判断扫描顺序方向do_reverse_range_和计算激光束中心位姿centered_laser_pose_ */
    if (up.z() > 0)
    {
        do_reverse_range_ = scan.angle_min > scan.angle_max;
        centered_laser_pose_ = tf::Stamped<tf::Pose>(tf::Transform(tf::createQuaternionFromRPY(0, 0, angle_center),
                                                                   tf::Vector3(0, 0, 0)),
                                                     ros::Time::now(), laser_frame_);
        ROS_INFO("Laser is mounted upwards.");
    }
    else
    {
        do_reverse_range_ = scan.angle_min < scan.angle_max;
        centered_laser_pose_ = tf::Stamped<tf::Pose>(tf::Transform(tf::createQuaternionFromRPY(M_PI, 0, -angle_center),
                                                                   tf::Vector3(0, 0, 0)),
                                                     ros::Time::now(), laser_frame_);
        ROS_INFO("Laser is mounted upside down.");
    }

    // Compute the angles of the laser from -x to x, basically symmetric and in increasing order
    /* 计算激光从-x到x的角度，基本对称且顺序递增 */

    /* 初始化存储每一个激光点的角度变量laser_angles_ */
    laser_angles_.resize(scan.ranges.size());

    // Make sure angles are started so that they are centered
    /* 确保开始的角度是居中的   */
    double theta = -std::fabs(scan.angle_min - scan.angle_max) / 2;

    /* 根据激光雷达数据中的angle_min,angle_increment等数据为每个激光束分配角度 */
    for (unsigned int i = 0; i < scan.ranges.size(); ++i)
    {
        laser_angles_[i] = theta;
        theta += std::fabs(scan.angle_increment);
    }

    /* 激光角度信息显示 */
    ROS_DEBUG("Laser angles in laser-frame: min: %.3f max: %.3f inc: %.3f", scan.angle_min, scan.angle_max, scan.angle_increment);
    ROS_DEBUG("Laser angles in top-down centered laser-frame: min: %.3f max: %.3f inc: %.3f", laser_angles_.front(), laser_angles_.back(), std::fabs(scan.angle_increment));

    /* 创建地图原点里程计位姿 */
    GMapping::OrientedPoint gmap_pose(0, 0, 0);

    // setting maxRange and maxUrange here so we can set a reasonable default
    /* 设置maxRange和maxUrange，这样我们就可以设置一个合理的默认值 */

    /* 通过真实的激光雷达数据，设置gmapping算法中激光的最大距离和最大使用距离，前面没有一起初始化 */
    ros::NodeHandle private_nh_("~");
    if (!private_nh_.getParam("maxRange", maxRange_))
        maxRange_ = scan.range_max - 0.01;
    if (!private_nh_.getParam("maxUrange", maxUrange_))
        maxUrange_ = maxRange_;

    // The laser must be called "FLASER".
    // We pass in the absolute value of the computed angle increment, on the
    // assumption that GMapping requires a positive angle increment.  If the
    // actual increment is negative, we'll swap the order of ranges before
    // feeding each scan to GMapping.
    /* 设置GridSlamSLAM的激光对象 */
    gsp_laser_ = new GMapping::RangeSensor("FLASER",
                                           gsp_laser_beam_count_,
                                           fabs(scan.angle_increment),
                                           gmap_pose,
                                           0.0,
                                           maxRange_);
    ROS_ASSERT(gsp_laser_);

    /* 定义GridSlamSLAM的SensorMap传感器关联容器 设置激光传感器 */
    GMapping::SensorMap smap;
    smap.insert(make_pair(gsp_laser_->getName(), gsp_laser_));
    gsp_->setSensorMap(smap);

    /* 定义GridSlamSLAM的里程 m_pose, m_speed和m_acceleration记录了机器人的位姿、速度和加速度 */
    gsp_odom_ = new GMapping::OdometrySensor(odom_frame_);
    ROS_ASSERT(gsp_odom_);

    /// @todo Expose setting an initial pose
    /* 定义里程计的初始位姿 */
    GMapping::OrientedPoint initialPose;

    /* 如果无法获取在odom坐标系的激光雷达坐标 则把初始位姿设置为(0,0,0)，也是建立地图的起始位置 */
    if (!getOdomPose(initialPose, scan.header.stamp))
    {
        /* 无法确定激光初始姿态! 起始点将被设置为零  */
        ROS_WARN("Unable to determine inital pose of laser! Starting point will be set to zero.");
        initialPose = GMapping::OrientedPoint(0.0, 0.0, 0.0);
    }

    /* 设置帧匹配参数 */
    gsp_->setMatchingParameters(maxUrange_, maxRange_, sigma_, kernelSize_, lstep_, astep_, iterations_, lsigma_, ogain_, lskip_);

    /* 设置运动模型参数 */
    gsp_->setMotionModelParameters(srr_, srt_, str_, stt_);

    /* 设置更新距离参数 */
    gsp_->setUpdateDistances(linearUpdate_, angularUpdate_, resampleThreshold_);

    /* 设置更新时间参数 */
    gsp_->setUpdatePeriod(temporalUpdate_);

    /* 设置不生成地图参数 */
    gsp_->setgenerateMap(false);

    /* 初始化粒子个数，地图尺寸，分辨率，建图初始位姿 */
    gsp_->GridSlamProcessor::init(particles_, xmin_, ymin_, xmax_, ymax_, delta_, initialPose);

    /* 设置似然计算的平移采样距离和步长 */
    gsp_->setllsamplerange(llsamplerange_);
    gsp_->setllsamplestep(llsamplestep_);
    /// @todo Check these calls; in the gmapping gui, they use
    /// llsamplestep and llsamplerange intead of lasamplestep and
    /// lasamplerange.  It was probably a typo, but who knows.
    /* 设置似然计算的角度采样距离和步长 */
    gsp_->setlasamplerange(lasamplerange_);
    gsp_->setlasamplestep(lasamplestep_);

    /* 设置最低分数线 */
    gsp_->setminimumScore(minimum_score_);

    // Call the sampling function once to set the seed.
    /* 设置高斯噪声的随机数种子 */
    GMapping::sampleGaussian(1, seed_);

    ROS_INFO("Initialization complete");

    return true;
}

/**
 * @brief 加载处理激光数据函数
 *
 * @param scan 激光数据
 * @param gmap_pose 里程位姿
 */
bool SlamGMapping::addScan(const sensor_msgs::LaserScan &scan, GMapping::OrientedPoint &gmap_pose)
{
    /* 转换获取在odom坐标系的里程坐标 */
    if (!getOdomPose(gmap_pose, scan.header.stamp))
        return false;

    /* 验证scan完整性 */
    if (scan.ranges.size() != gsp_laser_beam_count_)
        return false;

    // GMapping wants an array of doubles...
    /* GMapping需要一个双精度浮点数组 */
    double *ranges_double = new double[scan.ranges.size()];

    // If the angle increment is negative, we have to invert the order of the readings.
    /* 如果角度的增量是负的，必须颠倒读数的顺序   */
    if (do_reverse_range_)
    {
        ROS_DEBUG("Inverting scan");
        int num_ranges = scan.ranges.size();

        /* 对一帧激光数据的预处理，对距离小于range_min的数据，最大化 */
        for (int i = 0; i < num_ranges; i++)
        {
            // Must filter out short readings, because the mapper won't
            /* 必须转double */
            if (scan.ranges[num_ranges - i - 1] < scan.range_min)
                ranges_double[i] = (double)scan.range_max;
            else
                ranges_double[i] = (double)scan.ranges[num_ranges - i - 1];
        }
    }
    else
    {
        /* 对一帧激光数据的预处理，对距离小于range_min的数据，最大化 */
        for (unsigned int i = 0; i < scan.ranges.size(); i++)
        {
            // Must filter out short readings, because the mapper won't
            /* 必须转double */
            if (scan.ranges[i] < scan.range_min)
                ranges_double[i] = (double)scan.range_max;
            else
                ranges_double[i] = (double)scan.ranges[i];
        }
    }

    /* 把ROS的激光雷达数据信息转换为 GMapping算法要求的形式 */
    GMapping::RangeReading reading(scan.ranges.size(),
                                   ranges_double,
                                   gsp_laser_,
                                   scan.header.stamp.toSec());

    // ...but it deep copies them in RangeReading constructor, so we don't
    // need to keep our array around.
    /* 在RangeReading构造函数中深度复制它们，所以不需要保留数组 */
    delete[] ranges_double;

    /* 为每一个reading设置里程计位姿 */
    reading.setPose(gmap_pose);

    /* ROS_DEBUG("scanpose (%.3f): %.3f %.3f %.3f\n",
            scan.header.stamp.toSec(),
            gmap_pose.x,
            gmap_pose.y,
            gmap_pose.theta);
            */
    ROS_DEBUG("processing scan");

    /* 用gmapping算法进行处理 */
    return gsp_->processScan(reading);
}

/**
 * @brief 消息过滤的激光数据订阅并转换于odom_frame_坐标系回调函数
 *
 * @param scan 消息过滤的激光数据
 */
void SlamGMapping::laserCallback(const sensor_msgs::LaserScan::ConstPtr &scan)
{
    /* 更新接收到激光数据次数 */
    laser_count_++;

    /* 根据throttle_scans_节流设置判断是否需要处理 */
    if ((laser_count_ % throttle_scans_) != 0)
        return;

    /* 静态变量 存储上一次地图更新的时间 */
    static ros::Time last_map_update(0, 0);

    // We can't initialize the mapper until we've got the first scan
    /* 第一次在得到扫描数据时，需初始化mapper */
    if (!got_first_scan_)
    {
        if (!initMapper(*scan))
            return;
        got_first_scan_ = true;
    }

    /* 定义当前里程计坐标系下的激光雷达位姿的临时变量，带方向 */
    GMapping::OrientedPoint odom_pose;

    /* 加载处理激光数据 */
    if (addScan(*scan, odom_pose))
    {
        ROS_DEBUG("scan processed");

        /* 获取最优粒子在地图坐标系下的位姿 */
        GMapping::OrientedPoint mpose = gsp_->getParticles()[gsp_->getBestParticleIndex()].pose;

        /* 位姿信息显示 */
        ROS_DEBUG("new best pose: %.3f %.3f %.3f", mpose.x, mpose.y, mpose.theta);
        ROS_DEBUG("odom pose: %.3f %.3f %.3f", odom_pose.x, odom_pose.y, odom_pose.theta);
        ROS_DEBUG("correction: %.3f %.3f %.3f", mpose.x - odom_pose.x, mpose.y - odom_pose.y, mpose.theta - odom_pose.theta);

        /* 设置激光雷达和地图的转换、里程和激光雷达的转换 */
        tf::Transform laser_to_map = tf::Transform(tf::createQuaternionFromRPY(0, 0, mpose.theta), tf::Vector3(mpose.x, mpose.y, 0.0)).inverse();
        tf::Transform odom_to_laser = tf::Transform(tf::createQuaternionFromRPY(0, 0, odom_pose.theta), tf::Vector3(odom_pose.x, odom_pose.y, 0.0));

        /* 上锁计算地图和里程坐标系之间转换 */
        map_to_odom_mutex_.lock();
        map_to_odom_ = (odom_to_laser * laser_to_map).inverse();
        map_to_odom_mutex_.unlock();

        /* 若已存在地图 或 时间大于地图更新间隔 */
        if (!got_map_ || (scan->header.stamp - last_map_update) > map_update_interval_)
        {
            /* 更新地图及最新时间戳 */
            updateMap(*scan);
            last_map_update = scan->header.stamp;
            ROS_DEBUG("Updated the map");
        }
    }
    else
        ROS_DEBUG("cannot process scan");
}

/**
 * @brief 位姿熵计算函数
 *
 * @return -entropy 位姿熵
 */
double SlamGMapping::computePoseEntropy()
{
    /* 累计权重 */
    double weight_total = 0.0;
    for (std::vector<GMapping::GridSlamProcessor::Particle>::const_iterator it = gsp_->getParticles().begin(); it != gsp_->getParticles().end(); ++it)
    {
        weight_total += it->weight;
    }

    /* 熵计算 */
    double entropy = 0.0;
    for (std::vector<GMapping::GridSlamProcessor::Particle>::const_iterator it = gsp_->getParticles().begin(); it != gsp_->getParticles().end(); ++it)
    {
        /* 信息的大小跟随机事件的概率有关。越小概率的事情发生了产生的信息量越大，如地震了；越大概率的事情发生了产生的信息量越小，如太阳从东边升起来了
           信息量：h（x） = - log(p(x)) ，越小越接近无穷大，越大越接近0
           信息熵就是我要获得某些信息的代价，当信息的稀有度越高，得到这个信息需要付出的代价越高
           比如猜谜语，如果知道每个答案正确的概率，采用概率猜对求出的信息熵（log2）就是答对需要的次数 */
        if (it->weight / weight_total > 0.0)
            entropy += it->weight / weight_total * log(it->weight / weight_total);
    }
    return -entropy;
}

/**
 * @brief 地图更新函数
 *
 * @param scan 消息过滤的激光数据
 */
void SlamGMapping::updateMap(const sensor_msgs::LaserScan &scan)
{
    ROS_DEBUG("Update map");

    /* 设置区域锁，避免与mapCallback动态地图回调处理冲突 */
    boost::mutex::scoped_lock map_lock(map_mutex_);

    /* 定义匹配器 */
    GMapping::ScanMatcher matcher;

    /* 设置matcher的激光数据 */
    matcher.setLaserParameters(scan.ranges.size(), &(laser_angles_[0]), gsp_laser_->getPose());

    /* 设置matcher的激光雷达最大的量程和使用距离 */
    matcher.setlaserMaxRange(maxRange_);
    matcher.setusableRange(maxUrange_);

    /* 设置生成地图参数 */
    matcher.setgenerateMap(true);

    /* 得到权重最高的粒子 */
    GMapping::GridSlamProcessor::Particle best = gsp_->getParticles()[gsp_->getBestParticleIndex()];

    /* 计算位姿熵 */
    std_msgs::Float64 entropy;
    entropy.data = computePoseEntropy();
    if (entropy.data > 0.0)
        entropy_publisher_.publish(entropy);

    /* 如果还没有地图 则初始化一个地图信息（分辨率和原点） */
    if (!got_map_)
    {
        map_.map.info.resolution = delta_;
        map_.map.info.origin.position.x = 0.0;
        map_.map.info.origin.position.y = 0.0;
        map_.map.info.origin.position.z = 0.0;
        map_.map.info.origin.orientation.x = 0.0;
        map_.map.info.origin.orientation.y = 0.0;
        map_.map.info.origin.orientation.z = 0.0;
        map_.map.info.origin.orientation.w = 1.0;
    }

    /* 获取地图的中点 */
    GMapping::Point center;
    center.x = (xmin_ + xmax_) / 2.0;
    center.y = (ymin_ + ymax_) / 2.0;

    /* 初始化一个激光匹配地图，可视化用 */
    GMapping::ScanMatcherMap smap(center, xmin_, ymin_, xmax_, ymax_, delta_);

    ROS_DEBUG("Trajectory tree:");
    /* 遍历粒子的整条轨迹，按照轨迹上各个节点存储的信息来重新绘制一个地图 */
    for (GMapping::GridSlamProcessor::TNode *n = best.node; n; n = n->parent)
    {
        ROS_DEBUG("  %.3f %.3f %.3f", n->pose.x, n->pose.y, n->pose.theta);
        if (!n->reading)
        {
            ROS_DEBUG("Reading is NULL");
            continue;
        }

        /* 标记有效区域还没计算 */
        matcher.invalidateActiveArea();

        /* 拓展地图大小、找到地图的有效区域 */
        matcher.computeActiveArea(smap, n->pose, &((*n->reading)[0]));

        /* 更新每个被经过的栅格数据 */
        matcher.registerScan(smap, n->pose, &((*n->reading)[0]));
    }

    /* 判断地图是否已经扩展 */
    if (map_.map.info.width != (unsigned int)smap.getMapSizeX() || map_.map.info.height != (unsigned int)smap.getMapSizeY())
    {

        // NOTE: The results of ScanMatcherMap::getSize() are different from the parameters given to the constructor
        //       so we must obtain the bounding box in a different way

        /* 获取smap栅格地图尺寸并更新全局地图尺寸和地图服务消息 */
        GMapping::Point wmin = smap.map2world(GMapping::IntPoint(0, 0));
        GMapping::Point wmax = smap.map2world(GMapping::IntPoint(smap.getMapSizeX(), smap.getMapSizeY()));
        xmin_ = wmin.x;
        ymin_ = wmin.y;
        xmax_ = wmax.x;
        ymax_ = wmax.y;

        ROS_DEBUG("map size is now %dx%d pixels (%f,%f)-(%f, %f)", smap.getMapSizeX(), smap.getMapSizeY(), xmin_, ymin_, xmax_, ymax_);

        map_.map.info.width = smap.getMapSizeX();
        map_.map.info.height = smap.getMapSizeY();
        map_.map.info.origin.position.x = xmin_;
        map_.map.info.origin.position.y = ymin_;
        map_.map.data.resize(map_.map.info.width * map_.map.info.height);

        ROS_DEBUG("map origin: (%f, %f)", map_.map.info.origin.position.x, map_.map.info.origin.position.y);
    }

    /* 根据smap地图中存储的栅格数据，更新地图服务消息map_.map.data[]的数据 */
    for (int x = 0; x < smap.getMapSizeX(); x++)
    {
        for (int y = 0; y < smap.getMapSizeY(); y++)
        {
            /// @todo Sort out the unknown vs. free vs. obstacle thresholding
            /* 区分未知、自由和障碍阈值 */

            /* 获取栅格数据 */
            GMapping::IntPoint p(x, y);
            double occ = smap.cell(p);
            assert(occ <= 1.0);

            /* 未知区域 */
            if (occ < 0)
                map_.map.data[MAP_IDX(map_.map.info.width, x, y)] = -1;
            /* 障碍区域 */
            else if (occ > occ_thresh_)
            {
                // map_.map.data[MAP_IDX(map_.map.info.width, x, y)] = (int)round(occ*100.0);
                map_.map.data[MAP_IDX(map_.map.info.width, x, y)] = 100;
            }
            /* 自由区域 */
            else
                map_.map.data[MAP_IDX(map_.map.info.width, x, y)] = 0;
        }
    }

    /* 设置地图存在标志位 */
    got_map_ = true;

    // make sure to set the header information on the map
    /* 确保在地图上设置头部信息 */
    map_.map.header.stamp = ros::Time::now();
    map_.map.header.frame_id = tf_.resolve(map_frame_);

    /* 发布地图和地图信息话题 */
    sst_.publish(map_.map);
    sstm_.publish(map_.map.info);
}

/**
 * @brief 服务 dynamic_map 动态地图回调函数
 *
 * @param req 服务请求接收 无
 * @param res 服务返回 栅格地图
 * @return true/false请求状态
 */
bool SlamGMapping::mapCallback(nav_msgs::GetMap::Request &req,
                               nav_msgs::GetMap::Response &res)
{
    /* 设置区域锁，避免与updateMap更新地图冲突 */
    boost::mutex::scoped_lock map_lock(map_mutex_);

    /* 若已存在地图且长宽不为0，返回全局栅格地图 */
    if (got_map_ && map_.map.info.width && map_.map.info.height)
    {
        res = map_;
        return true;
    }
    else
        return false;
}

/**
 * @brief 发布map_to_odom的TF转换函数
 */
void SlamGMapping::publishTransform()
{
    /* 地图和里程坐标系之间转换操作上锁 */
    map_to_odom_mutex_.lock();

    /* 发布map_to_odom的TF转换 没搞懂为什么加这个缓冲时间 */
    ros::Time tf_expiration = ros::Time::now() + ros::Duration(tf_delay_);
    tfB_->sendTransform(tf::StampedTransform(map_to_odom_, tf_expiration, map_frame_, odom_frame_));

    /* 地图和里程坐标系之间转换操作解锁 */
    map_to_odom_mutex_.unlock();
}
