/*
 * Author: Stacey <StaceyXC@foxmail.com>
 * Affiliation: Nanjing University
 * Date: 2022
 */

#ifndef PSO_MOBILITY_MODEL_H
#define PSO_MOBILITY_MODEL_H

#include <stdint.h>
#include <deque>
#include "ns3/mobility-model.h"
#include "ns3/vector.h"
#include "ns3/constant-velocity-helper.h"
#include "ns3/position-allocator.h"
#include "ns3/ptr.h"
#include "ns3/object.h"
#include "ns3/nstime.h"
#include "ns3/event-id.h"
#include "ns3/box.h"
#include "ns3/random-variable-stream.h"
#include <map>
#include <iterator>
#include <vector>

namespace ns3
{
    class AFSAMobilityModel : public  MobilityModel
    {
    public:
        static TypeId GetTypeId(void);
        virtual TypeId GetInstanceTypeId() const;
        AFSAMobilityModel();
        virtual ~AFSAMobilityModel();
        void Update(void); //循环函数
        void Swarm(void);  //群聚行为
        void Chase(void);  //追尾行为
        void find(void);//觅食行为,返回速度向量
        void UpdatePosition(void);//更新所有粒子位置
        double func(Vector position);
        Vector ChasingVelocity(Vector destination);//追逐时的速度
        void RandMove(void);//觅食行为的兜底
        void FindBest(void);


        static std::map<double,Vector> PositionStorage;//存储所有粒子的位置以及对应适应度
    private:
        std::pair<double , Vector> personalBest;
        double visual;       //视野
        double step;         //最大步长
        double delta;        //拥挤度因子
        int try_times;       //一条鱼的觅食迭代次数
        int totaliterations; //总迭代次数
        double fitness;//适应度
        double Rmax;
        double Rmin;
        double chaseFitness;//追尾行为的适应度
        double centrefitness;//群聚行为适应度
        double preyFitness;//觅食行为适应度
        Vector  optimalSolution;
        std::map<double ,Vector> personalStorage;//存储个人不同行为的适应度以及对应速度

        ConstantVelocityHelper m_helper;
        Vector personalBestPosition;
        Vector m_position;
        Vector m_velocity;
        Vector CentrePosition;//中心位置
        // EventId m_event;
        Box m_bounds;
        virtual Vector DoGetPosition(void) const;
        virtual void DoSetPosition(const Vector &position);
        virtual Vector DoGetVelocity(void) const;
        virtual void DoSetVelocity(const Vector &velocity);
    };

}

#endif