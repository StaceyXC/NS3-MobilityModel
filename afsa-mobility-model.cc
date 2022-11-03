/*
 * Author:Stacey
 * Date: 2022
 * 试错了两天，除了一些细节问题外，最关键的问题是没有清空Personalstorage和positionstorage，这导致
 * 出错。
 * 
 */

#include "afsa-mobility-model.h"
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
#define MAX 100000; //定义适应度上限

namespace ns3
{
  std::map<double, Vector> AFSAMobilityModel::PositionStorage;
  NS_LOG_COMPONENT_DEFINE("AFSAMobilityModel");
  NS_OBJECT_ENSURE_REGISTERED(AFSAMobilityModel);

  TypeId
  ns3::AFSAMobilityModel::GetTypeId(void)
  {
    static TypeId tid = TypeId("ns3::AFSAMobilityModel")
                            .SetParent<MobilityModel>()
                            .SetGroupName("Mobility")
                            .AddConstructor<AFSAMobilityModel>();
    return tid;
  }
  TypeId ns3::AFSAMobilityModel::GetInstanceTypeId() const
  {
    return AFSAMobilityModel::GetTypeId();
  }

  ns3::AFSAMobilityModel::AFSAMobilityModel()
  {
    double VeloRNGMin = -10;
    double VeloRNGMax = 10; //速度范围
    Rmax = 1.0;
    Rmin = -1.0;
    step = 2;
    visual = 20.0;
    delta = 10;
    try_times = 50;
    totaliterations = 200;
    optimalSolution.x = 55.0; //已经定好了最佳位置，所以搜索结束条件应该是是否找到这个点，与迭代次数无关，所以目的是模拟过程而不是找答案
    optimalSolution.y = 25.0;
    optimalSolution.z = 15.0;
    Ptr<UniformRandomVariable> randNumGen = CreateObject<UniformRandomVariable>();
    randNumGen->SetAttribute("Min", DoubleValue(VeloRNGMin));
    randNumGen->SetAttribute("Max", DoubleValue(VeloRNGMax));
    m_velocity.x = randNumGen->GetValue();
    m_velocity.y = randNumGen->GetValue(); //应该是初始化速度
    m_velocity.z = randNumGen->GetValue();
    for (int i = 0; i <= 200; i++)
    {
      Simulator::Schedule(Seconds(i), &AFSAMobilityModel::Update, this);
    }
  }

  ns3::AFSAMobilityModel::~AFSAMobilityModel()
  {
  }

  inline Vector
  AFSAMobilityModel::DoGetPosition(void) const //其中变量值不变
  {
    return Vector(m_position.x, m_position.y, m_position.z);
  }

  void
  AFSAMobilityModel::DoSetPosition(const Vector &position)
  {
    // X sub i+1 = X sub i + V sub i+1
    Vector newPosition = position + m_velocity;
    m_position = newPosition;
    m_helper.SetPosition(m_position);
    // m_event.Cancel ();
    // m_event = Simulator::ScheduleNow (&PSOMobilityModel::Start, this);
  }

  inline Vector
  AFSAMobilityModel::DoGetVelocity(void) const
  {
    return Vector(m_velocity.x, m_velocity.y, m_velocity.z);
  }
  double
  AFSAMobilityModel::func(Vector position)
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
  AFSAMobilityModel::Update(void)
  {
    personalStorage.erase(personalStorage.begin(),personalStorage.end());
    PositionStorage.erase(PositionStorage.begin(),PositionStorage.end());
    fitness = func(DoGetPosition());
    if (fitness < 2)
      return;
    UpdatePosition();
    personalStorage.insert({fitness, DoGetPosition()});
    Swarm();
    Chase();
    //找到最优行为
    FindBest();
    if (personalBest.first == fitness)
      RandMove();
    else
      DoSetVelocity(personalBest.second);
    DoSetPosition(DoGetPosition());
    std::cout << "fitness:" << personalBest.first << std::endl;
    //std::cout << "position:" << DoGetPosition() << std::endl;
  }

  void
  AFSAMobilityModel::Swarm(void) //群聚行为
  {
    int n_swarm = 0; //初始化
    CentrePosition = Vector(0.0, 0.0, 0.0);
    //寻找视野范围内粒子个数
    for (std::map<double, Vector>::iterator it = PositionStorage.begin(); it != PositionStorage.end(); ++it)
    {
      if (CalculateDistance(m_position, it->second) < visual)
      {
        n_swarm += 1;
        CentrePosition = CentrePosition + it->second;
      }
    }
    CentrePosition = CentrePosition - m_position; //
    n_swarm -= 1;                                 //减掉自己
    CentrePosition.x = CentrePosition.x / n_swarm;
    CentrePosition.y = CentrePosition.y / n_swarm;
    CentrePosition.z = CentrePosition.z / n_swarm; //中心位置
    centrefitness = func(CentrePosition);          //得到中心位置适应度
    if (centrefitness < fitness)                   //判断成立则中心位置更优秀，同时不能太拥挤
    {
      personalStorage.insert({centrefitness, ChasingVelocity(CentrePosition)});
    }
    else
      find();
  }
  void
  AFSAMobilityModel::Chase(void) //追尾行为
  {
    int n_swarm = 0;
    Vector temp;
    std::pair<double, Vector> entryWithMinValueG = std::make_pair(5000000.0, Vector(0.0, 0.0, 0.0)); // for reference to find max
    for (std::map<double, Vector>::iterator it = PositionStorage.begin(); it != PositionStorage.end(); it++)
    {
      if (CalculateDistance(m_position, it->second) < visual)
      {
        n_swarm += 1;
        if (it->first < entryWithMinValueG.first)
          entryWithMinValueG = std::make_pair(it->first, it->second);
      }
    }
   // std::cout << "nswarm: " << n_swarm << std::endl;
    chaseFitness = entryWithMinValueG.first;
    if (entryWithMinValueG.first < fitness)
    {
      personalStorage.insert({chaseFitness, ChasingVelocity(entryWithMinValueG.second)});
    }
    else
      find();
  }

  void
  AFSAMobilityModel::find(void) //=prey觅食行为
  {
    preyFitness = MAX; //初始化
    Ptr<UniformRandomVariable> randNumGen = CreateObject<UniformRandomVariable>();
    randNumGen->SetAttribute("Min", DoubleValue(Rmin));
    randNumGen->SetAttribute("Max", DoubleValue(Rmax));
    Vector temp;
    for (int i = 0; i < try_times; i++)
    {
      temp.x = m_position.x + visual * randNumGen->GetValue();
      temp.y = m_position.y + visual * randNumGen->GetValue();
      temp.z = m_position.z + visual * randNumGen->GetValue();
      double tempfitness = func(temp);
      if (tempfitness < fitness)
      {
        personalStorage.insert({tempfitness, ChasingVelocity(temp)});
        preyFitness = tempfitness;
        break;
      }
    }
  }

  Vector
  AFSAMobilityModel::ChasingVelocity(Vector destination)
  {
    //计算单位方向向量
    Vector direction = (destination - m_position);
    double distance = CalculateDistance(m_position, destination);
    direction.x = direction.x / distance;
    direction.y = direction.y / distance;
    direction.z = direction.z / distance;
    Ptr<UniformRandomVariable> randNumGen = CreateObject<UniformRandomVariable>();
    randNumGen->SetAttribute("Min", DoubleValue(0));
    randNumGen->SetAttribute("Max", DoubleValue(Rmax * step));
    double speed = randNumGen->GetValue();
    Vector newVelocity(direction.x * speed, direction.y * speed, direction.z * speed);
    return newVelocity;
  }

  void
  AFSAMobilityModel::RandMove(void)
  {
    Ptr<UniformRandomVariable> randNumGen = CreateObject<UniformRandomVariable>();
    randNumGen->SetAttribute("Min", DoubleValue(Rmin * step));
    randNumGen->SetAttribute("Max", DoubleValue(Rmax * step));
    Vector newVelocity(randNumGen->GetValue(), randNumGen->GetValue(), randNumGen->GetValue());
    DoSetVelocity(newVelocity);
  }

  void
  AFSAMobilityModel::UpdatePosition(void)
  {
    PositionStorage.insert({fitness, m_position});
  }

  void
  AFSAMobilityModel::FindBest(void)
  {
    std::map<double, Vector>::iterator currentEntryG;
    std::pair<double, Vector> entryWithMinValueG = std::make_pair(5000000.0, Vector(0.0, 0.0, 0.0));
    for (currentEntryG = personalStorage.begin();
         currentEntryG != personalStorage.end();
         ++currentEntryG)
    {
      if (currentEntryG->first < entryWithMinValueG.first)
      {
        entryWithMinValueG = std::make_pair(
            currentEntryG->first,
            currentEntryG->second);
      }
    }
    personalBest = entryWithMinValueG;
  }

  void
  AFSAMobilityModel::DoSetVelocity(const Vector &velocity)
  {
    Vector newVelocity(velocity.x, velocity.y, velocity.z);
    m_velocity = newVelocity;
    m_helper.SetVelocity(m_velocity);
    // m_event.Cancel ();
    // m_event = Simulator::ScheduleNow (&PSOMobilityModel::Start, this);
  }

}
