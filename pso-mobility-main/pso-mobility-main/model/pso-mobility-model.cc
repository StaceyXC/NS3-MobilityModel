/*
* Author: Benjamin Westburg <benjamin.westburg@temple.edu>
* Affiliation: Temple University - College of Science and Technology
* Date: 2021
*/


#include "pso-mobility-model.h"
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

namespace ns3 {

Vector PSOMobilityModel::optimalSolution;
Vector PSOMobilityModel::groupBestPosition;
std::map<double,Vector> PSOMobilityModel::globalStorage;//MAP相当于PYTHON中的字典，一一对应。

NS_LOG_COMPONENT_DEFINE ("PSOMobilityModel");
NS_OBJECT_ENSURE_REGISTERED (PSOMobilityModel);

TypeId
ns3::PSOMobilityModel::GetTypeId (void)
{
  static TypeId tid = TypeId ("ns3::PSOMobilityModel")
    .SetParent<MobilityModel>()
    .SetGroupName ("Mobility")
    .AddConstructor<PSOMobilityModel>() 		   		    
   ;
   return tid;
}
TypeId ns3::PSOMobilityModel::GetInstanceTypeId() const
{
  return PSOMobilityModel::GetTypeId();
}

ns3::PSOMobilityModel::PSOMobilityModel ()
{
  totalIterations = 100;//迭代次数
  iteration = 0;//现迭代次数
  individualComponent = 2.0;//自我学习因子    
  groupComponent = 2.0;//群体学习因子
  randomComponent1 = 0.5;
  randomComponent2 = 0.5;
  optimalSolution.x = 55.0;//已经定好了最佳位置，所以搜索结束条件应该是是否找到这个点，与迭代次数无关，所以目的是模拟过程而不是找答案
  optimalSolution.y = 25.0;
  optimalSolution.z = 15.0;
  inertiaWeightDecrement = 0.5/totalIterations;//惯性质量，有一个问题是没有随着迭代次数增加而递减。
  double VeloRNGMin = -10;
  double VeloRNGMax = 10;
  Ptr<UniformRandomVariable> randNumGen = CreateObject<UniformRandomVariable> ();
  randNumGen->SetAttribute ("Min", DoubleValue (VeloRNGMin));
  randNumGen->SetAttribute ("Max", DoubleValue (VeloRNGMax));
  m_velocity.x = randNumGen->GetValue ();
  m_velocity.y = randNumGen->GetValue ();//应该是初始化速度
  m_velocity.z = randNumGen->GetValue ();
  for (int i = 0; i <= 200; i++) {
    Simulator::Schedule (Seconds (i), &PSOMobilityModel::Update, this);
  }
  //m_helper.Unpause ();
}

ns3::PSOMobilityModel::~PSOMobilityModel ()
{
}

void
PSOMobilityModel::Update (void)
{

  //每一次更新，粒子的和速度和位置不同，于是每一次计算得适应度不同，这也会使得PBEST可能改变。
  if (iteration == totalIterations) {//原判断是判断现在位置是否等同于最佳位置，但一直报错，所以现在用迭代次数作为结束条件
  	EndMobility();//如果找到了最佳位置，则结束模拟
  }
  DoSetFitnessValue();//每一次迭代，重新设置fitness value
  DoSetPersonalBestPosition();
  DoSetGroupBestPosition();
  CalculateVelocity();
  DoSetPosition(DoGetPosition());
  std::cout << "pos:" << DoGetPosition() << std::endl;
  std::cout << "fit val:" << DoGetFitnessValue() << std::endl;
  std::cout << "opt val:" << DoGetOptimalSolution() << std::endl;
  std::cout << "group best pos:" << DoGetGroupBestPosition() << std::endl;
  std::cout << "personalbest :" << DoGetPersonalBestPosition() << std::endl;
  iteration = iteration + 1;
  inertiaWeight = inertiaWeight -inertiaWeightDecrement;//随着迭代次数增加，惯性权重从0.9向0.4递减，因为越往后，集体最优解越正确。
  //Simulator::Schedule (Seconds (1), &PSOMobilityModel::Update, this);
}

inline Vector
PSOMobilityModel::DoGetPosition (void) const
{
  return Vector (m_position.x, m_position.y, m_position.z);
}

void
PSOMobilityModel::DoSetPosition (const Vector &position)
{
  // X sub i+1 = X sub i + V sub i+1
  Vector newPosition = position + m_velocity;
  m_position = newPosition;
  m_helper.SetPosition (m_position);
  //m_event.Cancel ();
  //m_event = Simulator::ScheduleNow (&PSOMobilityModel::Start, this);
}

inline Vector
PSOMobilityModel::DoGetVelocity (void) const
{
  return Vector (m_velocity.x, m_velocity.y, m_velocity.z);
}

void
PSOMobilityModel::CalculateVelocity (void)
{
  //V sub i+1 = w * V sub i + c1r1 (PBEST sub i - current position) + c2r2 (GBEST sub i - current position)
  
  //inertia component calculation
  double inertiaComponentVectorX = m_velocity.x * inertiaWeight;
  double inertiaComponentVectorY = m_velocity.y * inertiaWeight;
  double inertiaComponentVectorZ = m_velocity.z * inertiaWeight;
  
  //personal component calculation
  double personalComponentVectorX = personalBestPosition.x - m_position.x;
  double personalComponentVectorY = personalBestPosition.y - m_position.y;
  double personalComponentVectorZ = personalBestPosition.z - m_position.z;
  double personalComponentScalar = individualComponent * randomComponent1;//个人学习因子*随机因子
  personalComponentVectorX = personalComponentVectorX * personalComponentScalar;
  personalComponentVectorY = personalComponentVectorY * personalComponentScalar;
  personalComponentVectorZ = personalComponentVectorZ * personalComponentScalar;
  
  //group component calculation
  Vector GBP = DoGetGroupBestPosition();
  double groupComponentVectorX = GBP.x - m_position.x;
  double groupComponentVectorY = GBP.y - m_position.y;
  double groupComponentVectorZ = GBP.z - m_position.z;
  double groupComponentScalar = groupComponent * randomComponent2;
  groupComponentVectorX = groupComponentVectorX * groupComponentScalar;
  groupComponentVectorY = groupComponentVectorY * groupComponentScalar;
  groupComponentVectorZ = groupComponentVectorZ * groupComponentScalar;
  
  //new velocity vector calculation
  double newVelocityVectorX = inertiaComponentVectorX + personalComponentVectorX + groupComponentVectorX;
  double newVelocityVectorY = inertiaComponentVectorY + personalComponentVectorY + groupComponentVectorY;
  double newVelocityVectorZ = inertiaComponentVectorZ + personalComponentVectorZ + groupComponentVectorZ;
  
  Vector newVelocity (newVelocityVectorX,newVelocityVectorY,newVelocityVectorZ);
  DoSetVelocity(newVelocity);
}

void
PSOMobilityModel::DoSetVelocity (const Vector &velocity)
{
  Vector newVelocity (velocity.x,velocity.y,velocity.z);
  m_velocity = newVelocity;
  m_helper.SetVelocity (m_velocity);
  //m_event.Cancel ();
  //m_event = Simulator::ScheduleNow (&PSOMobilityModel::Start, this);
}

Vector
PSOMobilityModel::DoGetOptimalSolution (void) 
{
  return optimalSolution;
}

void 
PSOMobilityModel::DoSetOptimalSolution (void)
{
}

Vector
PSOMobilityModel::DoGetPersonalBestPosition (void)
{
  return personalBestPosition;
}

void
PSOMobilityModel::DoSetPersonalBestPosition (void)
{
  // find largest double entry in map, and set corresponding positional vector as PBEST (personalBestPosition).
  // calculation of PBEST
  std::map<double, Vector>::iterator currentEntryP; // iteration over personal storage
  std::pair<double, Vector> entryWithMinValueP = std::make_pair(5000000.0, Vector(0.0,0.0,0.0)); // for reference to find max
  for (currentEntryP = personalStorage.begin();
       currentEntryP != personalStorage.end();
       ++currentEntryP) {
       if (currentEntryP->first
           < entryWithMinValueP.first) {
           entryWithMinValueP
               = std::make_pair(
                   currentEntryP->first,
                   currentEntryP->second);
       }
   }//相当于遍历 personal storage 中的数据，找到最小的适应度。
   Vector entryWithMinValuePVector = entryWithMinValueP.second;//设置最佳位置。
   personalBestPosition = entryWithMinValuePVector;
}

Vector
PSOMobilityModel::DoGetGroupBestPosition (void)
{
  return groupBestPosition;
}

void 
PSOMobilityModel::DoSetGroupBestPosition (void)
{
  // find largest double entry in map, and set corresponding positional vector as GBEST (groupBestPosition).
  // calculation of GBEST 
  std::map<double, Vector>::iterator currentEntryG; // iteration over personal storage
  std::pair<double, Vector> entryWithMinValueG = std::make_pair(5000000.0, Vector(0.0,0.0,0.0)); // for reference to find max
  for (currentEntryG = globalStorage.begin();
       currentEntryG != globalStorage.end();
       ++currentEntryG) {
       if (currentEntryG->first
           < entryWithMinValueG.first) {
           entryWithMinValueG
               = std::make_pair(
                   currentEntryG->first,
                   currentEntryG->second);
       }
   }
   Vector entryWithMinValueGVector = entryWithMinValueG.second;
   groupBestPosition = entryWithMinValueGVector;
}

double
PSOMobilityModel::DoGetFitnessValue (void)
{
  return fitnessValue;
}

void
PSOMobilityModel::DoSetFitnessValue (void) 
{
  // fitness(optimization) function: 10 * (x - optimal_solution sub x)^2 + 20 * (y - optimal_solution sub y)^2 + 30 * (z -
  // optimal_solution sub z)^2
  // global optimum (solution) is x sub 1 = 55, x sub 2 = 25, x sub 3 = 15
  // goal is to attain minimum fitness value.//计算适应度的函数。 目的是寻找这个函数的最值，所以适应度即为函数值。模拟过程就是寻找函数最值
  double x = m_position.x;
  double y = m_position.y;
  double z = m_position.z;
  double optSolX = optimalSolution.x;
  double optSolY = optimalSolution.y;
  double optSolZ = optimalSolution.z;
  double base1 = x - optSolX;
  double base2 = y - optSolY;
  double base3 = z - optSolZ;
  double power = 2.0;
  fitnessValue = 10 * pow(base1,power) + 20 * pow(base2,power) + 30 * pow(base3,power);
  personalStorage.insert({ fitnessValue, m_position });
  globalStorage.insert({ fitnessValue, m_position });
}

void
PSOMobilityModel::EndMobility (void)
{
}

}


    		   

