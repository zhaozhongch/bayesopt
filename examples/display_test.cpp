#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <unistd.h>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/assignment.hpp>
#include "bayesopt/bayesopt.hpp"
#include "param_loader.hpp"
#include "/home/zhaozhong/bayesopt/utils/displaygp.hpp"


//g++ -std=c++11 display_test.cpp /home/zhaozhong/bayesopt/matplotpp/matplotpp.cc /home/zhaozhong/bayesopt/matplotpp/gl2ps.c /home/zhaozhong/bayesopt/utils/displaygp.cpp -o displaytest -L /usr/local/include/bayesopt/ -l bayesopt -lnlopt -lm -lGLU -lGL -lglut

// True function as red line, estimation is blue line, Data points as black star, uncertainty is green line

class ExampleOneD: public bayesopt::ContinuousModel
{
public:
  ExampleOneD(bayesopt::Parameters par):
    ContinuousModel(1,par) {}

  double evaluateSample(const vectord& xin)
  {
    int i = 0;
    int j = 0;
    if (xin.size() != 1)
      {
	std::cout << "WARNING: This only works for 1D inputs." << std::endl
		  << "WARNING: Using only first component." << std::endl;
      }
    //usleep(0.5); 
    //for(i = 0;i<100000; i++){;}//delay
      //for(j = 0; j <1000; j++){;}   
    //std::cout<<"xin is............................ "<<xin(0)<<std::endl;//很奇怪。。。这个函数一直在被调用或者检查？？都不知道为什么
    //从源代码看evaluateSampleInternal调用一次evaluateSample就调用一次，但是我cout显示的是就算evaluateSampleInternal没调用时这个也在调用
    //不画图就不会出现这种情况，发现utils/displaygp.cpp　里有evaluateSample并且似乎注释有个重复一千次什么的，问题就在那儿了
    double x = xin(0);
    return (x-0.3)*(x-0.3) + sin(20*x)*0.2;
  };

  bool checkReachability(const vectord &query)
  {return true;};

  void printOptimal()
  {
    std::cout << "Optimal:" << 0.23719 << std::endl;
  }
};

bayesopt::utils::DisplayProblem1D GLOBAL_MATPLOT;

void display( void ){ GLOBAL_MATPLOT.display(); }
void reshape( int w,int h ){ GLOBAL_MATPLOT.reshape(w,h); }
void idle( void ) { glutPostRedisplay(); } 

void mouse(int button, int state, int x, int y ){ GLOBAL_MATPLOT.mouse(button,state,x,y); }
void motion(int x, int y ){ GLOBAL_MATPLOT.motion(x,y); }
void passive(int x, int y ){ GLOBAL_MATPLOT.passivemotion(x,y); }

void keyboard(unsigned char key, int x, int y)
{
    GLOBAL_MATPLOT.keyboard(key, x, y); 
    if(key=='r')   //Toogle run/stop
      { 
	GLOBAL_MATPLOT.toogleRUN();
      }
    if(key=='s')   //Activate one step
      { 
	GLOBAL_MATPLOT.setSTEP();
      }
};

int menu()
{
  std::string input;
  int option = 0;
  while ((option < 1) || (option > 5))
    {
      std::cout << "Please select an option for the parameters:\n\n"
		<< "  1- Default parameters.\n"
		<< "  2- Student t process model.\n"
		<< "  3- Combined kernel.\n"
		<< "  4- Lower Confidence Bound.\n"
		<< "  5- A-optimality criteria.\n\n"
		<< "Select [1-5]>";

      std::cin >> input;
      std::istringstream is(input);
      is >> option;
    }
  return option;
}

int main(int nargs, char *args[])
{
  bopt_params parameters = initialize_parameters_to_default();
  parameters.n_init_samples = 7;
  parameters.n_iterations = 20;
  parameters.verbose_level = 2;
  //parameters.kernel.name = "kMaternARD3";

  switch( menu() )
    {
    case 1: break;
    case 2: 
      {
	set_surrogate(&parameters,"sStudentTProcessNIG"); 
	parameters.n_iter_relearn = 5;
	break;
      }
    case 3:   
      { 
	set_kernel(&parameters,"kSum(kPoly3,kRQISO)");
	double mean[128] = {1, 1, 1, 1};
	double std[128] = {5, 5, 5, 5};
	size_t nhp = 4;
	memcpy(parameters.kernel.hp_mean, mean, nhp * sizeof(double));
	memcpy(parameters.kernel.hp_std,std, nhp * sizeof(double));
	parameters.kernel.n_hp = nhp;
	break;
      }
    case 4:
      set_criteria(&parameters,"cLCB");
      parameters.crit_params[0] = 5;
      parameters.n_crit_params = 1;
      break;      
    case 5:
      set_criteria(&parameters,"cAopt");
      parameters.n_crit_params = 0;
      break;
    default:
      break;
    };

  boost::numeric::ublas::vector<double> lowerBound(1);
  boost::numeric::ublas::vector<double> upperBound(1);
  //lowerBound(0)=0;upperBound(0)=0.25;//这里设置的范围是出的图的横轴哈，是x不是y.　另外横轴是归一化(normalized)了的。所以图示结果始终在０到１之间
  //scoped_ptr是一个类似标准库中的auto_ptr智能指针，它包装了new操作符在椎上分配的动态对象，能够保证动态创建的对象在任何时候都能够被正确的删除。 看不懂算求
  bayesopt::Parameters parameters_class(parameters);
  boost::scoped_ptr<ExampleOneD> opt(new ExampleOneD(parameters_class));
  //opt->setBoundingBox(lowerBound,upperBound);//scoped_ptr把Example定义的对象赋值到了指针中，所有原来带有的函数由指针对象引用就是了。即原来ExampleOneD test(par); test.SetBoundingBox() 现在opt->setBoundingBox(lowerBound,upperBound);
  GLOBAL_MATPLOT.init(opt.get(),1);

  glutInit(&nargs, args);
  glutCreateWindow(50,50,800,650);
  glutDisplayFunc( display );
  glutReshapeFunc( reshape );
  glutIdleFunc( idle );
  glutMotionFunc( motion );
  glutMouseFunc( mouse );
  glutPassiveMotionFunc(passive);    
  glutKeyboardFunc( keyboard );        
  glutMainLoop();    

  return 0;
}
