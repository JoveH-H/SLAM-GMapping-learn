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

#include "ros/ros.h"
#include "sensor_msgs/LaserScan.h"
#include "std_msgs/Float64.h"
#include "nav_msgs/GetMap.h"
#include "tf/transform_listener.h"
#include "tf/transform_broadcaster.h"
#include "message_filters/subscriber.h"
#include "tf/message_filter.h"

#include "gmapping/gridfastslam/gridslamprocessor.h"
#include "gmapping/sensor/sensor_base/sensor.h"

#include <boost/thread.hpp>

/**
 * @brief SlamGMapping类
 * SLAM GMapping的集成类
 */
class SlamGMapping
{
    /* 公共 */
public:
    SlamGMapping();
    SlamGMapping(ros::NodeHandle &nh, ros::NodeHandle &pnh);
    SlamGMapping(unsigned long int seed, unsigned long int max_duration_buffer);
    ~SlamGMapping(); /* 析构函数 */

    void init();                                                            /* SLAM初始化函数 */
    void startLiveSlam();                                                   /* SLAM启动函数 */
    void startReplay(const std::string &bag_fname, std::string scan_topic); /* bag回放启动函数 */
    void publishTransform();                                                /* tf坐标转换发布函数 */

    void laserCallback(const sensor_msgs::LaserScan::ConstPtr &scan);                  /* 激光数据话题回调函数 */
    bool mapCallback(nav_msgs::GetMap::Request &req, nav_msgs::GetMap::Response &res); /* 动态地图服务回调函数 */
    void publishLoop(double transform_publish_period);                                 /* map_to_odom的TF转换发布循环函数 */

    /* 私有 */
private:
    ros::NodeHandle node_;                                                 /* ros节点句柄 */
    ros::Publisher entropy_publisher_;                                     /* 话题 entropy 熵 的发布者 */
    ros::Publisher sst_;                                                   /* 话题 OccupancyGrid 栅格地图 的发布者 */
    ros::Publisher sstm_;                                                  /* 话题 map_metadata 栅格地图的特征的基本信息 的发布者 */
    ros::ServiceServer ss_;                                                /* 服务 dynamic_map 动态地图 的服务端  */
    tf::TransformListener tf_;                                             /* tf转换 */
    message_filters::Subscriber<sensor_msgs::LaserScan> *scan_filter_sub_; /* 话题 scan 消息过滤的激光数据 的订阅者 */
    tf::MessageFilter<sensor_msgs::LaserScan> *scan_filter_;               /* 话题 scan 消息过滤的激光数据 的tf消息筛选器 */
    tf::TransformBroadcaster *tfB_;                                        /* map_to_odom的tf变换广播 */

    GMapping::GridSlamProcessor *gsp_; /* GridSlamSLAM的对象 */
    GMapping::RangeSensor *gsp_laser_; /* GridSlamSLAM的激光对象 */
    // The angles in the laser, going from -x to x (adjustment is made to get the laser between
    // symmetrical bounds as that's what gmapping expects)
    std::vector<double> laser_angles_; /* 存储每一个激光点的角度 */
    // The pose, in the original laser frame, of the corresponding centered laser with z facing up
    tf::Stamped<tf::Pose> centered_laser_pose_; /* 激光束中心位姿 */
    // Depending on the order of the elements in the scan and the orientation of the scan frame,
    // We might need to change the order of the scan
    bool do_reverse_range_;              /* 反向扫描获取 */
    unsigned int gsp_laser_beam_count_;  /* 激光束计数 */
    GMapping::OdometrySensor *gsp_odom_; /* GridSlamSLAM的里程对象 */

    bool got_first_scan_; /* 获得第一帧scan标志位 */

    bool got_map_;                   /* 地图存在标志位 */
    nav_msgs::GetMap::Response map_; /* 接收地图服务消息 */

    ros::Duration map_update_interval_; /* 地图更新间隔 */
    tf::Transform map_to_odom_;         /* 地图和里程坐标系之间转换 */
    boost::mutex map_to_odom_mutex_;    /* 地图和里程坐标系之间转换操作区域互斥锁 */
    boost::mutex map_mutex_;            /* 地图操作区域互斥锁 */

    int laser_count_;    /* 接收到激光数据次数 */
    int throttle_scans_; /* 节流激光扫描次数 */

    boost::thread *transform_thread_; /* map_to_odom的TF转换发布循环线程 */

    std::string base_frame_;  /*机器人中心参考系*/
    std::string laser_frame_; /*激光雷达参考系*/
    std::string map_frame_;   /*地图参考系*/
    std::string odom_frame_;  /*里程参考系*/

    void updateMap(const sensor_msgs::LaserScan &scan);
    bool getOdomPose(GMapping::OrientedPoint &gmap_pose, const ros::Time &t);             /* 获取当前里程计坐标系odom_frame下激光雷达的位姿函数 */
    bool initMapper(const sensor_msgs::LaserScan &scan);                                  /* 初始化地图制作器函数 */
    bool addScan(const sensor_msgs::LaserScan &scan, GMapping::OrientedPoint &gmap_pose); /* 加载处理激光数据函数 */
    double computePoseEntropy();                                                          /* 位姿熵计算函数 */

    // Parameters used by GMapping
    double maxRange_;          /* 激光雷达最大的量程 */
    double maxUrange_;         /* 激光雷达最大使用距离 */
    double maxrange_;          /* 激光雷达最大的量程 无用*/
    double minimum_score_;     /* 最低分数线 */
    double sigma_;             /* 端点匹配的标准差 */
    int kernelSize_;           /* 内核中要查找对应关系的数量 */
    double lstep_;             /* 平移优化步长 */
    double astep_;             /* 旋转优化步长 */
    int iterations_;           /* 扫描匹配的迭代步数 */
    double lsigma_;            /* 扫描匹配概率的激光标准差 */
    double ogain_;             /* 评估可能性时使用的增益 */
    int lskip_;                /* 评估可能性时使用的增益 */
    double srr_;               /* 平移时平移里程误差 */
    double srt_;               /* 平移时旋转里程误差 */
    double str_;               /* 旋转时平移里程误差 */
    double stt_;               /* 旋转时旋转里程误差 */
    double linearUpdate_;      /* 机器每平移该距离后处理一次激光扫描数据 */
    double angularUpdate_;     /* 机器每旋转该弧度后处理一次激光扫描数据 */
    double temporalUpdate_;    /* 如果最新扫描处理的速度比更新速度慢，则处理扫描；小于零的值将关闭基于时间的更新 */
    double resampleThreshold_; /* 基于Neff的重采样门限 */
    int particles_;            /* 滤波器中粒子数目 */
    double xmin_;              /* 地图x向初始最小尺寸 */
    double ymin_;              /* 地图x向初始最大尺寸 */
    double xmax_;              /* 地图y向初始最小尺寸 */
    double ymax_;              /* 地图y向初始最大尺寸 */
    double delta_;             /* 地图分辨率 */
    double occ_thresh_;        /* 栅格地图占用值的阈值 */
    double llsamplerange_;     /* 似然计算的平移采样距离 */
    double llsamplestep_;      /* 似然计算的平移采样步长 */
    double lasamplerange_;     /* 似然计算的角度采样距离 */
    double lasamplestep_;      /* 似然计算的角度采样步长 */

    ros::NodeHandle private_nh_; /* ros节点句柄 */

    unsigned long int seed_; /* 高斯噪声的随机数种子 */

    double transform_publish_period_; /* tf变换发布间隔 */
    double tf_delay_;                 /* tf变换发布间隔 */
};
