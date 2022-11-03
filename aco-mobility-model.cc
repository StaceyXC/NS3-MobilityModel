/*
 * Author: Stacey Lee <stacey10@163.com>
 * Affiliation: Nanjing University - Electronic Science and Engineering
 * Date: 2022
 */

#include <limits>
#include "ns3/abort.h"
#include "ns3/simulator.h"
#include "ns3/uinteger.h"
#include "ns3/log.h"
#include "ns3/boolean.h"
#include "ns3/config.h"
#include "ns3/test.h"
#include <cmath>
#include <vector>
#include "ns3/random-variable-stream.h"
#include "ns3/double.h"
#include "ns3/vector.h"
#include "aco-mobility-model.h"
#include <numeric>
#include <cstdlib>

namespace ns3
{
    std::vector<double> ACOMobilityModel::P_Storage;
    double step = 1.0; //步长

    NS_LOG_COMPONENT_DEFINE("ACOMobilityModel");
    NS_OBJECT_ENSURE_REGISTERED(ACOMobilityModel);

    TypeId
    ns3::ACOMobilityModel::GetTypeId(void)
    {
        static TypeId tid = TypeId("ns3::ACOMobilityModel")
                                .SetParent<MobilityModel>()
                                .SetGroupName("Mobility")
                                .AddConstructor<ACOMobilityModel>();
        return tid;
    }
    TypeId ns3::ACOMobilityModel::GetInstanceTypeId() const
    {
        return ACOMobilityModel::GetTypeId();
    }

    ns3::ACOMobilityModel::ACOMobilityModel()
    {
        n = 9;                      //城市数量
        tau = function(m_position); //初始化信息素浓度
        iteration = 0;
        Heu_F = 1 / tau;
        alpha = 1.0; //信息素重要性，过大，则蚂蚁过度依赖之前信息。
        beta = 5.0;  //过大容易陷入局部最优解。
        Rho = 0.1;
        Q = 100;
        P0 = 0.2;
        totalIterations = 100;
        optimalSolution.x = 55.0; //已经定好了最佳位置，所以搜索结束条件应该是是否找到这个点，与迭代次数无关，所以目的是模拟过程而不是找答案
        optimalSolution.y = 25.0;
        optimalSolution.z = 15.0;
        VeloRNGMin = -10;
        VeloRNGMax = 10; //速度范围
        PartVeloMax = 1;
        PartVeloMin = -1;
        Ptr<UniformRandomVariable> randNumGen = CreateObject<UniformRandomVariable>();
        randNumGen->SetAttribute("Min", DoubleValue(VeloRNGMin));
        randNumGen->SetAttribute("Max", DoubleValue(VeloRNGMax));
        m_velocity.x = randNumGen->GetValue();
        m_velocity.y = randNumGen->GetValue(); //初始化速度
        m_velocity.z = randNumGen->GetValue();
        for (int i = 0; i <= 200; i++)
        {
            Simulator::Schedule(Seconds(i), &ACOMobilityModel::update, this);
        }
    }

    ns3::ACOMobilityModel::~ACOMobilityModel()
    {
    }

    inline Vector
    ACOMobilityModel::DoGetPosition(void) const
    {
        return Vector(m_position.x, m_position.y, m_position.z);
    }

    void
    ACOMobilityModel::DoSetPosition(const Vector &position)
    {
        // X sub i+1 = X sub i + V sub i+1
        Vector newPosition = position + m_velocity;
        if (function(newPosition) < function(m_position))
        {
            m_position = newPosition;
            m_helper.SetPosition(m_position);
        }
        else
            m_helper.SetPosition(m_position);
        // m_event.Cancel ();
        // m_event = Simulator::ScheduleNow (&PSOMobilityModel::Start, this);
    }
    double //返回某位置函数的值
    ACOMobilityModel::function(Vector position)
    {
        double x = position.x;
        double y = position.y;
        double z = position.z;
        double optSolX = optimalSolution.x;
        double optSolY = optimalSolution.y;
        double optSolZ = optimalSolution.z;
        double base1 = x - optSolX;
        double base2 = y - optSolY;
        double base3 = z - optSolZ;
        double power = 2.0;
        return 10 * pow(base1, power) + 20 * pow(base2, power) + 30 * pow(base3, power);
    }

    void
    ACOMobilityModel::update(void)
    {
        if (iteration == totalIterations)
        {                  //原判断是判断现在位置是否等同于最佳位置，但一直报错，所以现在用迭代次数作为结束条件
            EndMobility(); //如果找到了最佳位置，则结束模拟
        }
        DoSetProbability();
        DoSetPc();    //赌盘轮转法求出PC。一般P越大，说明这个蚂蚁正确性更高，所以抽到他概率越大
        if (Pc >= P0) // p大于一定值，则周围出现极限值可能性大。
            DoPartCalculate();
        else
            DoWholeCalculate();
        DoSetPosition(DoGetPosition());
        std::cout << "pos:" << DoGetPosition() << std::endl;
        std::cout << "vel: " << DoGetVelocity() << std::endl;
        iteration = iteration + 1;
        Volatilization();
    }

    void //局部搜索
    ACOMobilityModel::DoPartCalculate(void)
    {
        Ptr<UniformRandomVariable> randNumGen = CreateObject<UniformRandomVariable>();
        randNumGen->SetAttribute("Min", DoubleValue(PartVeloMin));
        randNumGen->SetAttribute("Max", DoubleValue(PartVeloMax));
        double newVelocityVectorX = randNumGen->GetValue(); //一个 -1到1之间的随机数
        double newVelocityVectorY = randNumGen->GetValue();
        double newVelocityVectorZ = randNumGen->GetValue();
        Vector newVelocity(newVelocityVectorX, newVelocityVectorY, newVelocityVectorZ);
        DoSetVelocity(newVelocity);
    }
    void //全局搜索,由于无人机在全区域瞬移有点太扯，所以我弄成速度更大的随机移动
    ACOMobilityModel::DoWholeCalculate(void)
    {
        Ptr<UniformRandomVariable> randNumGen = CreateObject<UniformRandomVariable>();
        randNumGen->SetAttribute("Min", DoubleValue(VeloRNGMin - 5));
        randNumGen->SetAttribute("Max", DoubleValue(VeloRNGMax - 5));
        double newVelocityVectorX = randNumGen->GetValue(); //一个 -5到5之间的随机数
        double newVelocityVectorY = randNumGen->GetValue();
        double newVelocityVectorZ = randNumGen->GetValue();
        Vector newVelocity(newVelocityVectorX, newVelocityVectorY, newVelocityVectorZ);
        DoSetVelocity(newVelocity);
    }
    void
    ACOMobilityModel::DoSetVelocity(const Vector &velocity)
    {
        Vector newVelocity(velocity.x, velocity.y, velocity.z);
        m_velocity = newVelocity;
        m_helper.SetVelocity(m_velocity);
    }
    void
    ACOMobilityModel::DoSetPc(void)
    {
        double totalP; //所有粒子概率和
        totalP = accumulate(P_Storage.begin(), P_Storage.end(), 0.0);
        Pc = P / totalP;
    }

    void
    ACOMobilityModel::DoSetProbability(void)
    {
        P = P + pow(tau, alpha) * pow(Heu_F, beta); // tau 越小，这个范围内有极值点概率变大 P越小。
        P_Storage.insert(P_Storage.begin(), P);
    }
    void //信息素挥发和积累
    ACOMobilityModel::Volatilization(void)
    {
        tau = (1 - Rho) * tau + function(m_position);
    }
    inline Vector
    ACOMobilityModel::DoGetVelocity(void) const
    {
        return Vector(m_velocity.x, m_velocity.y, m_velocity.z);
    }
    void
    ACOMobilityModel::EndMobility(void)
    {
    }

}