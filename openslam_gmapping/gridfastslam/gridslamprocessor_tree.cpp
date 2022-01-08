#include <string>
#include <deque>
#include <list>
#include <map>
#include <set>
#include <fstream>
//#include <gsl/gsl_blas.h>

#include <gmapping/utils/stat.h>
#include "gmapping/gridfastslam/gridslamprocessor.h"

namespace GMapping
{
    /* 定义标准库函数的标准命名空间 */
    using namespace std;

    /**
     * @brief 粒子运动轨迹节点构建函数
     *
     * @param p 粒子位姿
     * @param w 粒子权重
     * @param n 所要构造的节点的父节点
     * @param c 该节点的子节点的数量
     */
    GridSlamProcessor::TNode::TNode(const OrientedPoint &p, double w, TNode *n, unsigned int c)
    {
        pose = p;    /* 粒子位姿设置 */
        weight = w;  /* 粒子权重设置 */
        childs = c;  /* 子节点数量设置 */
        parent = n;  /* 父节点设置 */
        reading = 0; /* 记录激光传感器数据初始化 */
        gweight = 0;

        /* 增加父节点的子节点数量 */
        if (n)
        {
            n->childs++;
        }

        flag = 0;      /* 访问标志初始化 */
        accWeight = 0; /* 权重总和初始化 */
    }

    /**
     * @brief 粒子运动轨迹节点析构函数
     */
    GridSlamProcessor::TNode::~TNode()
    {
        /* 判定父节点是否存在，并减少父节点的子节点数量，释放资源 */
        if (parent && (--parent->childs) <= 0)
            delete parent;
        assert(!childs);
    }

    // BEGIN State Save/Restore

    GridSlamProcessor::TNodeVector GridSlamProcessor::getTrajectories() const
    {
        TNodeVector v;
        TNodeMultimap parentCache;
        TNodeDeque border;

        for (ParticleVector::const_iterator it = m_particles.begin(); it != m_particles.end(); it++)
        {
            TNode *node = it->node;
            while (node)
            {
                node->flag = false;
                node = node->parent;
            }
        }

        for (ParticleVector::const_iterator it = m_particles.begin(); it != m_particles.end(); it++)
        {
            TNode *newnode = new TNode(*(it->node));

            v.push_back(newnode);
            assert(newnode->childs == 0);
            if (newnode->parent)
            {
                parentCache.insert(make_pair(newnode->parent, newnode));
                // cerr << __func__ << ": node " << newnode->parent << " flag=" << newnode->parent->flag<< endl;
                if (!newnode->parent->flag)
                {
                    // cerr << __func__ << ": node " << newnode->parent << " flag=" << newnode->parent->flag<< endl;
                    newnode->parent->flag = true;
                    border.push_back(newnode->parent);
                }
            }
        }

        // cerr << __func__ << ": border.size(INITIAL)=" << border.size() << endl;
        // cerr << __func__ << ": parentCache.size()=" << parentCache.size() << endl;
        while (!border.empty())
        {
            // cerr << __func__ << ": border.size(PREPROCESS)=" << border.size() << endl;
            // cerr << __func__ << ": parentCache.size(PREPROCESS)=" << parentCache.size() << endl;
            const TNode *node = border.front();
            // cerr << __func__ << ": node " << node << endl;
            border.pop_front();
            if (!node)
                continue;

            TNode *newnode = new TNode(*node);
            node->flag = false;

            // update the parent of all of the referring childs
            pair<TNodeMultimap::iterator, TNodeMultimap::iterator> p = parentCache.equal_range(node);
            double childs = 0;
            for (TNodeMultimap::iterator it = p.first; it != p.second; it++)
            {
                assert(it->second->parent == it->first);
                (it->second)->parent = newnode;
                // cerr << "PS(" << it->first << ", "<< it->second << ")";
                childs++;
            }
            ////cerr << endl;
            parentCache.erase(p.first, p.second);
            // cerr << __func__ << ": parentCache.size(POSTERASE)=" << parentCache.size() << endl;
            assert(childs == newnode->childs);

            // unmark the node
            if (node->parent)
            {
                parentCache.insert(make_pair(node->parent, newnode));
                if (!node->parent->flag)
                {
                    border.push_back(node->parent);
                    node->parent->flag = true;
                }
            }
            // insert the parent in the cache
        }
        // cerr << __func__ << " : checking cloned trajectories" << endl;
        for (unsigned int i = 0; i < v.size(); i++)
        {
            TNode *node = v[i];
            while (node)
            {
                // cerr <<".";
                node = node->parent;
            }
            // cerr << endl;
        }

        return v;
    }

    void GridSlamProcessor::integrateScanSequence(GridSlamProcessor::TNode *node)
    {
        // reverse the list
        TNode *aux = node;
        TNode *reversed = 0;
        double count = 0;
        while (aux != 0)
        {
            TNode *newnode = new TNode(*aux);
            newnode->parent = reversed;
            reversed = newnode;
            aux = aux->parent;
            count++;
        }

        // attach the path to each particle and compute the map;
        if (m_infoStream)
            m_infoStream << "Restoring State Nodes=" << count << endl;

        aux = reversed;
        bool first = true;
        double oldWeight = 0;
        OrientedPoint oldPose;
        while (aux != 0)
        {
            if (first)
            {
                oldPose = aux->pose;
                first = false;
                oldWeight = aux->weight;
            }

            OrientedPoint dp = aux->pose - oldPose;
            double dw = aux->weight - oldWeight;
            oldPose = aux->pose;

            double *plainReading = new double[m_beams];
            for (unsigned int i = 0; i < m_beams; i++)
                plainReading[i] = (*(aux->reading))[i];

            for (ParticleVector::iterator it = m_particles.begin(); it != m_particles.end(); it++)
            {
                // compute the position relative to the path;
                double s = sin(oldPose.theta - it->pose.theta),
                       c = cos(oldPose.theta - it->pose.theta);

                it->pose.x += c * dp.x - s * dp.y;
                it->pose.y += s * dp.x + c * dp.y;
                it->pose.theta += dp.theta;
                it->pose.theta = atan2(sin(it->pose.theta), cos(it->pose.theta));

                // register the scan
                m_matcher.invalidateActiveArea();
                m_matcher.computeActiveArea(it->map, it->pose, plainReading);
                it->weight += dw;
                it->weightSum += dw;

                // this should not work, since it->weight is not the correct weight!
                // it->node=new TNode(it->pose, it->weight, it->node);
                it->node = new TNode(it->pose, 0.0, it->node);
                // update the weight
            }

            delete[] plainReading;
            aux = aux->parent;
        }

        // destroy the path
        aux = reversed;
        while (reversed)
        {
            aux = reversed;
            reversed = reversed->parent;
            delete aux;
        }
    }

    // END State Save/Restore

    // BEGIN

    /**
     * @brief 更新轨迹权重函数
     *
     * @param weightsAlreadyNormalized 权重已归一化
     */
    void GridSlamProcessor::updateTreeWeights(bool weightsAlreadyNormalized)
    {
        /* 确保所有权重归一化 */
        if (!weightsAlreadyNormalized)
        {
            normalize();
        }

        /* 重置所有轨迹树 */
        resetTree();

        /* 更新所有轨迹树权重 */
        propagateWeights();
    }

    /**
     * @brief 重置所有轨迹树函数
     *
     * 遍历了所有的粒子，各个轨迹节点的累积权重和访问计数
     */
    void GridSlamProcessor::resetTree()
    {
        // don't calls this function directly, use updateTreeWeights(..) !
        /* 不要直接调用这个函数，使用 updateTreeWeights */

        /* 遍历了所有的粒子 */
        for (ParticleVector::iterator it = m_particles.begin(); it != m_particles.end(); it++)
        {
            /* 沿着粒子的node节点向上追溯到根节点，清除遍历过程中各个节点的累积权重和访问计数 */
            TNode *n = it->node;
            while (n)
            {
                n->accWeight = 0;
                n->visitCounter = 0;
                n = n->parent;
            }
        }
    }

    /**
     * @brief 获取单轨迹树根节点权重函数
     *
     * 累积当前轨迹树的节点的访问计数和的权重，返回当前轨迹树的根节点权重
     *
     * @param n 轨迹树的节点
     * @param weight 当前轨迹树的累积权重
     * @return w 当前轨迹树的根节点权重
     */
    double propagateWeight(GridSlamProcessor::TNode *n, double weight)
    {
        if (!n)
            return weight; /* 返回当前轨迹树的节点权重 */

        double w = 0;
        n->visitCounter++;      /* 累积访问当前节点的计数 */
        n->accWeight += weight; /* 累积当前节点的权重 */

        /* 递归查询根节点权重 */
        if (n->visitCounter == n->childs)
        {
            w = propagateWeight(n->parent, n->accWeight);
        }

        assert(n->visitCounter <= n->childs);

        return w;
    }

    /**
     * @brief 更新所有轨迹树权重函数
     *
     * @return lastNodeWeight 粒子群的运动轨迹根节点的累积权重
     */
    double GridSlamProcessor::propagateWeights()
    {
        // don't calls this function directly, use updateTreeWeights(..) !

        // all nodes must be resetted to zero and weights normalized

        // the accumulated weight of the root
        /* 粒子群的运动轨迹根节点的累积权重 */
        double lastNodeWeight = 0;
        // sum of the weights in the leafs
        double aw = 0; /* 粒子群的运动轨迹叶子节点权重的和 */

        /* 定义所有（粒子群的最新位姿）权重的叶子节点 */
        std::vector<double>::iterator w = m_weights.begin();

        /* 遍历粒子群 */
        for (ParticleVector::iterator it = m_particles.begin(); it != m_particles.end(); it++)
        {
            /* 获取叶子节点权重并累加 */
            double weight = *w;
            aw += weight;

            TNode *n = it->node;                                        /* 获取当前粒子的新轨迹的叶子节点 */
            n->accWeight = weight;                                      /* 设置叶子节点的初始累积权重 */
            lastNodeWeight += propagateWeight(n->parent, n->accWeight); /* 从叶子节点寻轨迹的根节点，更新轨迹上节点的访问计数和累积权重，并获取根节点权重， */
            w++;                                                        /* 叶子节点（粒子群的最新位姿）权重入口随遍历粒子群变化 */
        }

        /* 所有粒子的轨迹的叶子和根节点的权重和都应该为1，因都属于某次数据更新的权重总和 */
        if (fabs(aw - 1.0) > 0.0001 || fabs(lastNodeWeight - 1.0) > 0.0001)
        {
            cerr << "ERROR: ";
            cerr << "root->accWeight=" << lastNodeWeight << "    sum_leaf_weights=" << aw << endl;
            assert(0);
        }

        /* 返回粒子群的运动轨迹根节点的累积权重 */
        return lastNodeWeight;
    }

};

// END
