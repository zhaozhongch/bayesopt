#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <unistd.h>
#include <fstream>
#include <string>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/assignment.hpp>
#include "bayesopt/bayesopt.hpp"
#include "param_loader.hpp"
#include "/home/zhaozhong/bayesopt/utils/displaygp.hpp"

//g++ 2D_EKF_Fusion_NIS.cpp -o 2D_EKF_Fusion_NIS -I /usr/include/eigen3  -lceres -lglog -lpthread -fopenmp -lcholmod -lblas -llapack -lbayesopt -lnlopt -lm

using namespace std;

class Fusion: public bayesopt::ContinuousModel
{
public:
  Fusion(bayesopt::Parameters par):
    ContinuousModel(2,par) {}

  double evaluateSample(const vectord& xin)
  {
    if (xin.size() != 2)
      {
	std::cout << "WARNING: This only works for 2D inputs." << std::endl
		  << "WARNING: Using only first component." << std::endl;
      }

    ofstream wT("/home/zhaozhong/ekf_bayesopt_noupload/applications/kf1/data/pn_on.csv",ios::out);
    //wT<<xin(0)<<","<<xin(1);
    wT<<xin(0)<<","<<xin(1);
    wT.close();
    //cout<<"finish writing process noise and observation noise"<<endl;
    //sleep(6); //relaunch.cpp keep reading if there is "relaunch" so wait for several sencods can make sure it has realuched

    system("./../../ekf_bayesopt_noupload/applications/kf1/2D_bayes_NIS");
    
    //sleep(3);//bag may stop about 6 seconds before the launch file ends. So wait for several seonds. Depends on when I choose to stop
    ifstream rT("/home/zhaozhong/ekf_bayesopt_noupload/applications/kf1/data/NIS.csv");
    if(!rT) {cout<<"no such file"<<endl;}
    string NIS;
    getline(rT,NIS);
    rT.close();
   
    double re = (double)atof(NIS.c_str());
    //cout<<"NEES is "<<re<<endl;

    return re;
  };

  bool checkReachability(const vectord &query)
  {return true;};

  void printOptimal()
  {
    //std::cout << "Optimal:" << 0.23719 << std::endl;
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
    par.n_init_samples = 20;
    par.n_iterations = 130;
    par.random_seed = -1;
    par.verbose_level = 1;
    par.surr_name = "sStudentTProcessNIG";
    //par.epsilon = 0.2;
    //par.RANSAC_Noise = 2;
    par.noise = 1e-7;
    par.force_jump = 5;
    //about noise see the source code: in bayesoptbase about line 123 we can see if (std::pow(mYPrev - yNext,2) < mParameters.noise)
    //so if the yNext is similar to mYPrev, we say it stucks. In SLAM as we use RANSAC, even for the same input, Y output will be different. See cost&k1_His2,3 when xin(0) is about -0.2817. So I can assume the noise is 1(Normally it is small)
    //Also see the code following, when we satisfy std::pow(mYPrev - yNext,2) < mParameters.noise, mCounterStuck++. Then if (mCounterStuck > mParameters.force_jump) we have "Forced random query". Default  mParameters.force_jump is 20. we set it to 10.
    //about there parameters see https://rmcantin.bitbucket.io/html/usemanual.html#basicparams usage should be object'name.n_iterations and so on
      } //par.l_type=L_MCMC; //par.crit_name="cLCB"; }//for this simple problem, 200 iteration. L_MCMC use 42.8 seconds while L_empiracal use 4.5s .But as the website says the former one may leads to higher accuracy(for this simple problem of course not). The two methods are learning methods for kernal hyperparameters.
// the l_type is not a string variable, so not need "" but the cLCB is. You can see these in parameters.hpp                          
  Fusion D1(par);
  vectord result(2);
  //branin.checkReachability(result); system will check reachability itself. But the rsult is strange,, it gives all the same result na matter waht kind of constraint I add
  boost::numeric::ublas::vector<double> lowerBound(2);
  boost::numeric::ublas::vector<double> upperBound(2);
  lowerBound(0)=0.001;upperBound(0)=10.0;
  lowerBound(1)=0.001;upperBound(1)=10.0;
  
  D1.setBoundingBox(lowerBound,upperBound);
  D1.optimize(result);
 
  std::cout << "Result: " << result << "->" 
	    << D1.evaluateSample(result) << std::endl;
  D1.printOptimal();
return 0;
}
