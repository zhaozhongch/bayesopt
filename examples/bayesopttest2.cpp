#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/assignment.hpp>
#include "bayesopt/bayesopt.hpp"
//#include "bayesopt/parameters.hpp"
#include "param_loader.hpp"//if wants to change the parameters. Must include parameters.hpp, but the param_loader.hpp has already include it. Para_loader is not in the system file /usr/... So if wants to include it you must use something -I to where it is or just copy the file to where the source file is
//#include "criteria_lcb.hpp"
//#include "criteria_functors.hpp"

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

    double x = xin(0) * 10;
    double y = xin(1) * 10;
    
    return x*x+y*y;              // branin(x,y); //in  practice use of bayesian optimization, we don't know what this function looks like. But we need know every input, what kind of output it has.
  }                     // in my use, we need know reprojection error at any position
/*********************
  double branin(double x, double y)
  {
    const double pi = boost::math::constants::pi<double>();
    const double rpi = pi*pi;
    return sqr(y-(5.1/(4*rpi))*sqr(x)
	       +5*x/pi-6)+10*(1-1/(8*pi))*cos(x)+10;
  };
*********************/
  bool checkReachability(const vectord &query)
  {
   bool flag=true;
   //if (abs(query(0))<0.6) {flag=false;}  // if (abs(query(0))<0.15|abs(query(1))<0.5) {flag=false;} 
   return flag;
  };

 // inline double sqr( double x ){ return x*x; };

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
    par.n_iterations = 190;
    par.random_seed = 0;
    par.verbose_level = 1;
    par.noise = 1e-10;//about there parameters see https://rmcantin.bitbucket.io/html/usemanual.html#basicparams usage should be object'name.n_iterations and so on
    par.crit_name="cLCB";  } //par.l_type=L_MCMC;  }//for this simple problem, 200 iteration. L_MCMC use 42.8 seconds while L_empiracal use 4.5s .But as the website says the former one may leads to higher accuracy(for this simple problem of course not). The two methods are learning methods for kernal hyperparameters.
// the l_type is not a string variable, so not need "" but the cLCB is. You can see these in parameters.hpp                          
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
