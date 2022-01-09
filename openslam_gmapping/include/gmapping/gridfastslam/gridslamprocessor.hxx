#ifdef MACOSX
// This is to overcome a possible bug in Apple's GCC.
#define isnan(x) (x == FP_NAN)
#endif

/**Just scan match every single particle.
If the scan matching fails, the particle gets a default likelihood.*/
/* 只需扫描匹配每一个粒子。如果扫描匹配失败，粒子获得默认的似然值。*/

/**
 * @brief 扫描匹配函数
 *
 * @param plainReading 激光雷达传感器数据
 */
inline void GridSlamProcessor::scanMatch(const double *plainReading)
{
    // sample a new pose from each scan in the reference
    /* 在参考图中，从每次扫描中取样一个新的姿势 */

    double sumScore = 0; /* 初始化总成绩 */

    /* 遍历粒子群 */
    for (ParticleVector::iterator it = m_particles.begin(); it != m_particles.end(); it++)
    {
        OrientedPoint corrected; /* 定义矫正过的姿态的传参变量 */
        double score, l, s;      /* 定义匹配得分、对数似然度、匹配度 */

        /* 爬山优化器 */
        score = m_matcher.optimize(corrected, it->map, it->pose, plainReading);

        // it->pose=corrected;

        /* 匹配得分超过最低分数线，匹配成功 */
        if (score > m_minimumScore)
        {
            /* 更新使用矫正过的新姿态 */
            it->pose = corrected;
        }
        else
        {
            /* 匹配失败,仍然使用传统的运动模型下的预测状态作为更新样本,并输出信息 */
            if (m_infoStream)
            {
                m_infoStream << "Scan Matching Failed, using odometry. Likelihood=" << l << std::endl;
                m_infoStream << "lp:" << m_lastPartPose.x << " " << m_lastPartPose.y << " " << m_lastPartPose.theta << std::endl;
                m_infoStream << "op:" << m_odoPose.x << " " << m_odoPose.y << " " << m_odoPose.theta << std::endl;
            }
        }

        /* 获取粒子的分数（s匹配度）和可能性（对数似然度） */
        m_matcher.likelihoodAndScore(s, l, it->map, it->pose, plainReading);

        /* 更新粒子群总成绩和粒子的权重 */
        sumScore += score;
        it->weight += l;
        it->weightSum += l;

        // set up the selective copy of the active area
        // by detaching the areas that will be updated
        /* 设置激活区域的选择副本，通过分离将要被更新的区域 */

        /* 标记有效区域还没计算 */
        m_matcher.invalidateActiveArea();

        /* 更新有效区域（地图）*/
        m_matcher.computeActiveArea(it->map, it->pose, plainReading);

        /* 这里只更新有效区域，不更新地图单元内容 */
        /* 更新地图单元内容，放在了重采样之后进行，重采样过程会删除一些权重较小的粒子，此时再对所有粒子更新地图，避免计算资源浪费 */
    }

    /* 显示平均扫描匹配分数 */
    if (m_infoStream)
        m_infoStream << "Average Scan Matching Score=" << sumScore / m_particles.size() << std::endl;
}

/**
 * @brief 粒子权重归一化函数
 *
 * 归一化所有粒子权重，计算权重相似度
 */
inline void GridSlamProcessor::normalize()
{
    // normalize the log m_weights 归一化log权重

    /* 计算粒子的权重增益 = 1 /（似然度的平滑因子 * size）  */
    double gain = 1. / (m_obsSigmaGain * m_particles.size());

    /* 计算粒子集合中最大的权重 */
    double lmax = -std::numeric_limits<double>::max();
    for (ParticleVector::iterator it = m_particles.begin(); it != m_particles.end(); it++)
    {
        lmax = it->weight > lmax ? it->weight : lmax;
    }
    // cout << "!!!!!!!!!!! maxwaight= "<< lmax << endl;

    m_weights.clear(); /* 粒子权重初始化清除 */
    double wcum = 0;   /* 累积权重初始化 */
    m_neff = 0;        /* 相似性初始化 */

    /* 遍历粒子群，更新所有粒子权重范围 */
    for (std::vector<Particle>::iterator it = m_particles.begin(); it != m_particles.end(); it++)
    {
        m_weights.push_back(exp(gain * (it->weight - lmax))); /* 更新权重在范围[0,1] */
        wcum += m_weights.back();                             /* 累积求和权重 */
        // cout << "l=" << it->weight<< endl;
    }

    m_neff = 0; /* 相似性初始化 */
    /* 遍历粒子权重，更新所有粒子权重归一化 */
    for (std::vector<double>::iterator it = m_weights.begin(); it != m_weights.end(); it++)
    {
        *it = *it / wcum; /* 权重除以权重和，得到归一化后的粒子权重 */

        /* 累加权重相似度倒数 */
        double w = *it;
        m_neff += w * w; /*  m_neff会一直小于1，若所有w都比较接近，则m_neff会较小，倒数后就会较大，就说明权重差距较小 */
    }

    m_neff = 1. / m_neff; /* 计算权重相似度 */
}

/**
 * @brief 重采样函数
 *
 * @param plainReading 激光雷达传感器数据
 * @param adaptSize 在重采样过程中选用的粒子数量
 * @param reading 激光雷达和里程位姿数据
 * @return hasResampled 已进行重采样
 */
inline bool GridSlamProcessor::resample(const double *plainReading, int adaptSize, const RangeReading *reading)
{
    /* 初始化已进行重采样标志位 */
    bool hasResampled = false;

    /* 定义运动轨迹节点容器 记录上代的粒子状态 */
    TNodeVector oldGeneration;

    /* 遍历粒子群，赋值运动轨迹叶子节点入口 */
    for (unsigned int i = 0; i < m_particles.size(); i++)
    {
        oldGeneration.push_back(m_particles[i].node);
    }

    /* 当权重相似性小于设置的阈值时，重采样 */
    if (m_neff < m_resampleThreshold * m_particles.size())
    {
        /* 提示开始重采样 */
        if (m_infoStream)
            m_infoStream << "*************RESAMPLE***************" << std::endl;

        /* 定义的重采样器 */
        uniform_resampler<double, double> resampler;

        /* 获取采样后的粒子在原始粒子集合中的索引 */
        m_indexes = resampler.resampleIndexes(m_weights, adaptSize);

        /* 显示采样后的粒子 */
        if (m_outputStream.is_open())
        {
            m_outputStream << "RESAMPLE " << m_indexes.size() << " ";
            for (std::vector<unsigned int>::const_iterator it = m_indexes.begin(); it != m_indexes.end(); it++)
            {
                m_outputStream << *it << " ";
            }
            m_outputStream << std::endl;
        }

        /* 源码为空函数，提供了扩展的可能 */
        onResampleUpdate();

        // BEGIN: BUILDING TREE 开始:构建树

        /* 定义记录重采样后的粒子集合的变量 */
        ParticleVector temp;

        unsigned int j = 0; /* 原始粒子群索引，用于计算需要去除的粒子deletedParticles的辅助变量 */

        /* 这是为了删除重新采样的粒子，记录需要去除的粒子 */
        std::vector<unsigned int> deletedParticles; // this is for deleteing the particles which have been resampled away.

        // cerr << "Existing Nodes:" ;
        /* 遍历保留的粒子索引 */
        for (unsigned int i = 0; i < m_indexes.size(); i++)
        {
            // cerr << " " << m_indexes[i];

            /* 当前原始粒子群索引不在保留粒子索引中 */
            while (j < m_indexes[i])
            {
                deletedParticles.push_back(j); /* 记录下不在保留索引中的粒子 */
                j++;                           /* 增加原始粒子群索引 */
            }

            /* 原始粒子群索引在保留粒子索引中，只增加原始粒子群索引 */
            if (j == m_indexes[i])
                j++;

            /* 设置粒子更新最新位姿和传感器数据，并记录在轨迹树中 */
            Particle &p = m_particles[m_indexes[i]];
            TNode *node = 0;
            TNode *oldNode = oldGeneration[m_indexes[i]];
            // cerr << i << "->" << m_indexes[i] << "B("<<oldNode->childs <<") ";
            node = new TNode(p.pose, 0, oldNode, 0);
            // node->reading=0;
            node->reading = reading;
            // cerr << "A("<<node->parent->childs <<") " <<endl;

            /* 将选取的样本记录在临时粒子集合temp中 */
            temp.push_back(p);
            temp.back().node = node;
            temp.back().previousIndex = m_indexes[i];
        }

        /* 当前原始粒子群索引小于保留索引数量，即当前原始粒子群索引后续的粒子均被放弃 */
        while (j < m_indexes.size())
        {
            deletedParticles.push_back(j); /* 记录下不在保留索引中的粒子 */
            j++;                           /* 增加原始粒子群索引 */
        }

        // cerr << endl;
        /* 删除放弃的粒子运动轨迹节点 */
        std::cerr << "Deleting Nodes:";
        for (unsigned int i = 0; i < deletedParticles.size(); i++)
        {
            std::cerr << " " << deletedParticles[i];
            delete m_particles[deletedParticles[i]].node;
            m_particles[deletedParticles[i]].node = 0;
        }
        std::cerr << " Done" << std::endl;

        // END: BUILDING TREE 结束:构建树

        /* 删除粒子群 */
        std::cerr << "Deleting old particles...";
        m_particles.clear();
        std::cerr << "Done" << std::endl;

        std::cerr << "Copying Particles and  Registering  scans...";

        /* 遍历重采样保留的临时粒子集合 */
        for (ParticleVector::iterator it = temp.begin(); it != temp.end(); it++)
        {
            it->setWeight(0);

            /* 设置有效面积未完成计算标志 */
            m_matcher.invalidateActiveArea();

            /* 将激光扫描的占用信息注册到栅格地图上 */
            m_matcher.registerScan(it->map, it->pose, plainReading);

            /* 粒子群新增重采样保留的粒子 */
            m_particles.push_back(*it);
        }
        std::cerr << " Done" << std::endl;

        /* 标记已进行了重采样 */
        hasResampled = true;
    }
    else /* 未满足重采样条件 */
    {
        int index = 0;
        std::cerr << "Registering Scans:";

        /* 获取上代的粒子群节点入口 */
        TNodeVector::iterator node_it = oldGeneration.begin();

        /* 遍历粒子群，更新各个粒子的轨迹和地图 */
        for (ParticleVector::iterator it = m_particles.begin(); it != m_particles.end(); it++)
        {
            // create a new node in the particle tree and add it to the old tree
            /* 在粒子树中创建一个新节点，并将其添加到旧树中  */

            // BEGIN: BUILDING TREE 开始:构建树
            TNode *node = 0;
            node = new TNode(it->pose, 0.0, *node_it, 0);

            // node->reading=0;
            node->reading = reading;
            it->node = node;

            // END: BUILDING TREE 结束:构建树

            /* 设置有效面积未完成计算标志 */
            m_matcher.invalidateActiveArea();

            /* 将激光扫描的占用信息注册到栅格地图上 */
            m_matcher.registerScan(it->map, it->pose, plainReading);

            /* 设置轨迹树中前一个粒子的索引 */
            it->previousIndex = index;

            index++;   /* 前一个粒子的索引+1 */
            node_it++; /* 上代的粒子群节点入口+1 */
        }
        std::cerr << "Done" << std::endl;
    }
    // END: BUILDING TREE

    return hasResampled;
}
