#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/assignment.hpp>
#include "bayesopt/bayesopt.hpp"
#include "param_loader.hpp"

//g++ bayesopttest1.cpp -o bayesopttest1 -I /usr/include/eigen3  -lceres -lglog -lpthread -fopenmp -lcholmod -lblas -llapack -lbayesopt -lnlopt -lm
//把源文件重新编译了之后要想显示出效果这个文件也得再编译一次,并sudo make install
class BraninNormalized: public bayesopt::ContinuousModel
{
public:
  BraninNormalized(bayesopt::Parameters par):
    ContinuousModel(2,par) {}

  double evaluateSample( const vectord& xin)
  {
     if (xin.size() != 2)
      {
	std::cout << "WARNING: This only works for 2D inputs." << std::endl
		  << "WARNING: Using only first two components." << std::endl;
      }
    //std::cout<<"xin(0) is.... "<<xin(0)<<std::endl;
    //std::cout<<"xin(1) is.... "<<xin(1)<<std::endl;
    double x = xin(0) * 15 - 5;
    double y = xin(1) * 15;
    
    return branin(x,y); //in practice use of bayesian optimization, we don't know what this function looks like. But we need know every input, what kind of output it has.
  }                     // in my use, we need know reprojection error at any position

  double branin(double x, double y)
  {
    const double pi = boost::math::constants::pi<double>();
    const double rpi = pi*pi;
    return sqr(y-(5.1/(4*rpi))*sqr(x)
	       +5*x/pi-6)+10*(1-1/(8*pi))*cos(x)+10;
  };

  bool checkReachability(const vectord &query)
  {
   bool flag=true;
   //if (abs(query(0))<0.6) {flag=false;}  // if (abs(query(0))<0.15|abs(query(1))<0.5) {flag=false;} 
   return flag;
  };

  inline double sqr( double x ){ return x*x; };

  void printOptimal()
  {
    vectord sv(2);  
    sv(0) = 0.1238938; sv(1) = 0.818333;
    std::cout << "Solutions: " << sv << "->" 
	      << evaluateSample(sv) << std::endl;
    sv(0) = 0.5427728; sv(1) = 0.151667;
    std::cout << "Solutions: " << sv << "->" 
	      << evaluateSample(sv) << std::endl;
    sv(0) = 0.961652; sv(1) = 0.1650;
    std::cout << "Solutions: " << sv << "->" 
	      << evaluateSample(sv) << std::endl;
  }

};


int main(int nargs, char *args[])
{
bayesopt::Parameters par;
 if(nargs > 1){
    if(!bayesopt::utils::ParamLoader::load(args[1], par)){
        std::cout << "ERROR: provided file \"" << args[1] << "\" does not exist" << std::endl;
        return -1;
    }
  }
  else{
    par = initialize_parameters_to_default();
    par.n_iterations = 20;
    par.random_seed = 0;//if this is non negative, every time the program will produce the same initial x points. otherwise it will be different.
    par.verbose_level = 1;
    par.noise = 1e-10;//about there parameters see https://rmcantin.bitbucket.io/html/usemanual.html#basicparams usage should be object'name.n_iterations and so on

    par.SLAM_Fail_Threshold = 100; //I add this parameter into parameter class. Sometimes if the parameter for SLAM is totally wrong. There may be just a few frame detected or there might be mo key frame. Then the error output is 0, which of course will be think as minimum. This should not happen. What I do is if the keyframe number is smaller than like 10, I'll set the error output from slam to be  a pretty big number. Then when bayesopt receive this number, it just jump from next iteration, it won't take this point into Gaussian Process. The cons is iteration time will be one less.
      }
    BraninNormalized branin(par);
  vectord result(2);
  //branin.checkReachability(result); system will check reachability itself. But the rsult is strange,, it gives all the same result na matter waht kind of constraint I add
  boost::numeric::ublas::vector<double> lowerBound(2);
  boost::numeric::ublas::vector<double> upperBound(2);
  lowerBound(0)=0.1;lowerBound(1)=0.15;  upperBound(0)=0.6;upperBound(1)=0.6;//Use this to add upper bound and lower bound for variable. default is 0 to 1
  branin.setBoundingBox(lowerBound,upperBound);
  branin.optimize(result);
 
  std::cout << "Result: " << result << "->" 
	    << branin.evaluateSample(result) << std::endl;
  branin.printOptimal();
return 0;
}
