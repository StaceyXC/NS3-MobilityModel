/*

author : Xingchen Lee

*/
#ifndef ACO_MOBILITY_MODEL_H
#define ACO_MOBILITY_MODEL_H

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

    class ACOMobilityModel : public MobilityModel
    {
    public:
        static TypeId GetTypeId(void);
        virtual TypeId GetInstanceTypeId() const;
        ACOMobilityModel();
        virtual ~ACOMobilityModel();
        double function(Vector position); //返回10 * (x - optimal_solution sub x)^2 + 20 * (y - optimal_solution sub y)^2 + 30 * (z -
        // optimal_solution sub z)^2的值
        void update(void);          //每一次迭代函数
        void EndMobility(void); //结束模拟
        void DoSetProbability(void);
        void DoSetPc(void);          //赌盘轮转函数设置Pc
        void DoPartCalculate(void);  //局部搜索
        void DoWholeCalculate(void); //全局搜索,随机速度
        void Volatilization(void);   //信息素挥发和积累函数
    private:
        int iteration; //现迭代次数
        int n = 9;     //城市数量
        double tau;    //信息素浓度
        double Heu_F;
        double alpha; //信息素重要性因子
        double beta;  //启发因子重要性因子
        double Rho;   //信息素蒸发因子
        double Q;     //信息素增加强度系数
        double L;     //记录本代每只蚂蚁走的路程
        double P;     //转移 概率
        double Pc;    //轮盘赌算法后得到的Pc;
        double P0;    // Pc的范围界限。
        int totalIterations; //总迭代次数
        double PartVeloMax;
        double PartVeloMin;
        double VeloRNGMin ;
        double VeloRNGMax ;
        ConstantVelocityHelper m_helper;
        Vector m_position;
        Vector m_velocity;
        // EventId m_event;
        Box m_bounds;
        Vector optimalSolution; 
        virtual Vector DoGetPosition(void) const;
        virtual void DoSetPosition(const Vector &position);
        virtual Vector DoGetVelocity(void) const;
        virtual void DoSetVelocity(const Vector &velocity);

    public:
        static std::vector<double> P_Storage; //存储每一个蚂蚁的概率
    };
}

#endif