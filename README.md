slam_gmapping

[![Build Status](https://travis-ci.com/ros-perception/slam_gmapping.svg?branch=melodic-devel)](https://travis-ci.org/ros-perception/slam_gmapping)
[![License](https://img.shields.io/badge/license-Apache%202-blue.svg)](./LICENSE)
[![Release](https://img.shields.io/badge/release-v1.0-blue)](https://github.com/JoveH-H/SLAM-GMapping-learn/releases/tag/v1.0)
[![Author](https://img.shields.io/badge/Author-Jove-%2300a8ff)](https://github.com/JoveH-H)
================================================================================================================================================================
## 版本

| 发布版本 | 发布模块 |
| --- | --- |
| `[v1.0]` | `[介绍和注释Gmapping的ROS封装和主要核心算法]` |

## 使用
```shell
$ mkdir -p ~/gmapping_ws/src
$ cd ~/gmapping_ws/src
$ git clone https://github.com/JoveH-H/SLAM-GMapping-learn.git
$ mv SLAM-GMapping-learn src
$ cd ~/gmapping_ws
$ catkin_make
$ source ~/gmapping_ws/devel/setup.sh
$ roslaunch gmapping slam_gmapping_pr2.launch
```

## 资料
#### 相关源码
| 博客 |
| --- |
| [《SLAM GMapping（1）ROS封装》](https://blog.csdn.net/qq_32618327/article/details/121706071) |
| [《SLAM GMapping（2）传感器》](https://blog.csdn.net/qq_32618327/article/details/122110185) |
| [《SLAM GMapping（3）地图结构》](https://blog.csdn.net/qq_32618327/article/details/122244370) |
| [《SLAM GMapping（4）SLAM处理器》](https://blog.csdn.net/qq_32618327/article/details/122318169) |
| [《SLAM GMapping（5）运动模型》](https://blog.csdn.net/qq_32618327/article/details/122330515) |
| [《SLAM GMapping（6）扫描匹配器》](https://blog.csdn.net/qq_32618327/article/details/122369165) |
| [《SLAM GMapping（7）粒子和轨迹》](https://blog.csdn.net/qq_32618327/article/details/122393554) |
| [《SLAM GMapping（8）重采样》](https://blog.csdn.net/qq_32618327/article/details/122393781) |

#### 相关应用
| 博客 |
| --- |
| [《ROS笔记（22） Gmapping》](https://joveh-h.blog.csdn.net/article/details/98219810) |

## 问题
欢迎以 [GitHub Issues](https://github.com/JoveH-H/SLAM-GMapping-learn/issues) 的形式提交问题和bug报告

## 声明
#### 个人声明
> CSDN博客专栏《OpenSLAM》讲述了一些开源SLAM算法，希望能帮助学者加快入门SLAM

#### 免责声明
> 以任何方式登录平台或直接、间接使用平台代码，均视为自愿接受免责声明

## 版权和许可
SLAM-GMapping-learn 由 [GNU General Public License v3.0 许可](https://github.com/JoveH-H/SLAM-GMapping-learn/blob/main/LICENSE) 提供

谢谢!