/*
-------------------------------------------------------------------------
   This file is part of BayesOpt, an efficient C++ library for 
   Bayesian optimization.

   Copyright (C) 2011-2015 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
   BayesOpt is free software: you can redistribute it and/or modify it 
   under the terms of the GNU Affero General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   BayesOpt is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Affero General Public License for more details.

   You should have received a copy of the GNU Affero General Public License
   along with BayesOpt.  If not, see <http://www.gnu.org/licenses/>.
------------------------------------------------------------------------
*/

#include <ctime>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <string>
//#include <chrono>
#include <sys/time.h>
#include "bayesopt/bayesoptbase.hpp"
#include "bayesopt/parameters.hpp"

#include "log.hpp"
#include "posteriormodel.hpp"
#include "specialtypes.hpp"
#include "bopt_state.hpp"

namespace bayesopt
{
  BayesOptBase::BayesOptBase(size_t dim, Parameters parameters):
    mParameters(parameters), mDims(dim)
  {
    // Random seed
    if (mParameters.random_seed < 0) mParameters.random_seed = std::time(0); 
    mEngine.seed(mParameters.random_seed);
    
    // Setting verbose stuff (files, levels, etc.)
    int verbose = mParameters.verbose_level;
    if (verbose>=3)
      {
        FILE* log_fd = fopen( mParameters.log_filename.c_str() , "w" );
        Output2FILE::Stream() = log_fd; 
        verbose -= 3;
      }

    switch(verbose)
      {
      case 0: FILELog::ReportingLevel() = logWARNING; break;
      case 1: FILELog::ReportingLevel() = logINFO; break;
      case 2: FILELog::ReportingLevel() = logDEBUG4; break;
      default:
	FILELog::ReportingLevel() = logERROR; break;
      }
  }

  BayesOptBase::~BayesOptBase()
  { } // Default destructor


  // OPTIMIZATION INTERFACE
  void BayesOptBase::optimize(vectord &bestPoint)
  {
   // auto start  = std::chrono::system_clock::now();
/*
    struct timeval tp;
    gettimeofday(&tp, NULL);
    long int ms_ini = tp.tv_sec * 1000 + tp.tv_usec / 1000;
*/
    assert(mDims == bestPoint.size());
    
    // Restore state from file
    if(mParameters.load_save_flag == 1 || mParameters.load_save_flag == 3)
      {
        BOptState state;
        bool load_succeed = state.loadFromFile(mParameters.load_filename, 
					       mParameters);
        if(load_succeed)
	  {
            restoreOptimization(state);
            FILE_LOG(logINFO) << "State succesfully restored from file \"" 
			      << mParameters.load_filename << "\"";
	  }
        else
	  {
	    // If load not succeed, print info message and then
	    // initialize a new optimization
            FILE_LOG(logINFO) << "File \"" << mParameters.load_filename 
			      << "\" does not exist,"
			      << " starting a new optimization";
            initializeOptimization();
	  }
      }
    else
      {
	// Initialize a new state
        initializeOptimization();
      }
    
    for (size_t ii = mCurrentIter; ii < mParameters.n_iterations; ++ii)
      {      
        stepOptimization();
        //std::cout<<"noise is.........................."<< mParameters.noise<<std::endl;
        //FILE_LOG(logINFO)<<"noise is "<<mParameters.noise;
      }
   
    bestPoint = getFinalResult();
   // auto end  = std::chrono::system_clock::now();
   // std::chrono::duration<double> diff = end-start;
   // std::cout<<"spend "<<diff.count()<<" seconds for bayesopt"<<std::endl;
/*
   gettimeofday(&tp, NULL);
   long int ms_end = tp.tv_sec * 1000 + tp.tv_usec / 1000; 
   long int diff = ms_end - ms_ini;
   std::cout<<"we spend "<<diff<<" msecs"<<std::endl;
*/
  } // optimize


  void BayesOptBase::stepOptimization()
  {
    // Find what is the next point.
    vectord xNext = nextPoint(); 
    //std::cout<<"xNext is \n" << xNext<<std::endl;
    //exit(0);
    double yNext = evaluateSampleInternal(xNext);
    

    //clear the mean, std and so on data here so that we can reget the new iterations value. In evaluateSample function we'll get those data.
    x.clear(); y.clear(); su.clear(); sl.clear(); c.clear();
    xx.clear(); yy.clear(); susu.clear(); slsl.clear(); cc.clear();
    
    if(yNext>mParameters.SLAM_Fail_Threshold)
    {
      //std::cout<<"new step......"<<std::endl;
      //std::cout<<"SLAM failed because the parameters are too bad, ignore this iteration and find another point"<<std::endl;
      //std::cout<<"Set yNext as a large number, this will change the Gaussian Process"<<std::endl;
      //yNext = mParameters.SLAM_Fail_Threshold;
      //return;
    }
    
    // If we are stuck in the same point for several iterations, try a random jump!
    if (mParameters.force_jump)
      {
        //std::cout<<"program is in force jump..................."<<std::endl;
        if (std::pow(mYPrev - yNext,2) < mParameters.noise+mParameters.RANSAC_Noise)//RANSAC造成的noise如果算到kernel里就影响太大
          {
            mCounterStuck++;
            FILE_LOG(logDEBUG) << "Stuck for "<< mCounterStuck << " steps";
          }
        else
          {
            mCounterStuck = 0;
          }
        mYPrev = yNext;

        if (mCounterStuck > mParameters.force_jump)
          {
            FILE_LOG(logINFO) << "Forced random query!";
            xNext = samplePoint();
            yNext = evaluateSampleInternal(xNext);
            mCounterStuck = 0;
          }
      }
    
    
    mModel->addSample(xNext,yNext);

    // int n_sample = mModel->getData()->getNSamples();
    // if(n_sample >= 2){
    //   vectord lastx = mModel->getData()->getSampleX(n_sample-2);
    //   double lasty = mModel->getData()->getSampleY(n_sample-2);
    //   vectord currentx = mModel->getData()->getSampleX(n_sample-1);
    //   double currenty = mModel->getData()->getSampleY(n_sample-1);
    //   double norm_diff_x = 0,norm_diff_y = 0;
    //   //std::cout<<"last x and current x ";
    //   for(int i = 0; i<lastx.size(); i++){
    //     norm_diff_x += (lastx[i] - currentx[i]) * (lastx[i] - currentx[i]);
    //     std::cout<<lastx[i]<<","<<currentx[i]<<"    ";
    //   }
    //   norm_diff_y = sqrt((lasty - currenty)*(lasty-currenty));
    //   if(norm_diff_y <= 0.005 && norm_diff_x<=0.01)
    //     exit(0);
    //   //std::cout<<"last y and current y "<<lasty<<","<<currenty<<", final x norm "<<sqrt(norm_diff_x)<<", final y norm "<<norm_diff_y<<std::endl;
    // }

    // Update surrogate model
    bool retrain = ((mParameters.n_iter_relearn > 0) && 
		    ((mCurrentIter + 1) % mParameters.n_iter_relearn == 0));

    if (retrain)  // Full update
      {
        mModel->updateHyperParameters();
        mModel->fitSurrogateModel();
      }
    else          // Incremental update
      {
        mModel->updateSurrogateModel();
      } 

    plotStepData(mCurrentIter,xNext,yNext);
    mModel->updateCriteria(xNext);
    mCurrentIter++;
    
    // Save state if required
    if(mParameters.load_save_flag == 2 || mParameters.load_save_flag == 3)
      {
        BOptState state;
        saveOptimization(state);
        state.saveToFile(mParameters.save_filename);
      }
    //store the data
   //std::cout<<"program is before store the data"<<std::endl;
   if(mDims == 1)
   {

    //std::cout<<"program is in store the data"<<std::endl;
    for(size_t j=0;j<mDims;j++)
    {
      vectord tmp_x(1);
      vectord unnormalize_x(1);
      tmp_x(0) = xNext[j];
      unnormalize_x = remapPoint(tmp_x);
      AllX.push_back(unnormalize_x(0));
    }
      AllY.push_back(yNext);
    static int plot_flag = 0;
    int n=1000;//we'll plot 1000 points

    //private member now
    //std::vector<double> x,y,z,su,sl,c;
    for(size_t m = 0; m<n; m++)
    {
      x.push_back(m*0.001);
    }
    y = x; z = x; su = x; sl = x; c = x;
    
    const std::string filename_first = "/home/zhaozhong/bayesopt/bayesplot/iteration";
    const std::string filename_last  = ".csv";
    std::ostringstream oss;
    oss<<filename_first<<plot_flag<<filename_last;   //don not add <<std::endl, or file name will be "A.csv "(cannot be read) instead of "A.csv"        
    std::ofstream wr(oss.str().c_str(),std::ios::out|std::ios::app);//.str() to string, .str().c_str() to char*
    //comment the txt output
    vectord q(1);
    for(size_t i=0; i<n; ++i)//n==1000, we'll sample 1000 points
    {
        q(0) = x[i];        // Query, x is from 0.000,0.001,0.002 to 1.000 because the domain(定义域) of getPrediction is normalised
        
        ProbabilityDistribution* pd = mModel->getPrediction(q);
        y[i] = pd->getMean();                                //Expected value
        su[i] = y[i] + 2*pd->getStd();                       //Upper bound (95 %)
        sl[i] = y[i] - 2*pd->getStd();                       //Lower bound (95 %)
        c[i] = -mModel->evaluateCriteria(q);             //Criteria value
        //z[i] = mModel->evaluateSample(q);                //Target function true value
        
        /*unormalize x so that the plot's x axis we got can be lower boung to upper bound instead of 0 to 1*/
        vectord tmp_x(1);
        vectord unnormalize_x(1);
        tmp_x(0) = x[i];
        unnormalize_x = remapPoint(tmp_x);
        x[i] = unnormalize_x(0);
        wr<<x[i]<<","<<y[i]<<","<<su[i]<<","<<sl[i]<<","<<c[i]<<std::endl; //comment the txt output
     }
     
     for(size_t k=0;k<(plot_flag+mParameters.n_init_samples);k++)
     {
        //std::cout<<"sample x is "<<AllX[k]<<", sample y is "<<AllY[k]<<std::endl;
        wr<<AllX[k]<<","<<AllY[k]<<std::endl;//so last several rows is sample point, //comment the txt output
     }
     
      plot_flag++;
   }

   if(mDims == 2)
   {
        vectord tmp_x(2);
        vectord unnormalize_x(2);
        tmp_x(0) = xNext[0];
        tmp_x(1) = xNext[1];
        unnormalize_x = remapPoint(tmp_x);
        for(size_t l=0;l<mDims;l++)
        {
            AllX.push_back(unnormalize_x[l]);//every two number for one Y!;
        }
        AllY.push_back(yNext);

       static int plot_flag = 0;
       
       const std::string filename_first = "/home/zhaozhong/bayesopt/bayesplot/iteration";
       const std::string filename_last  = ".csv";
       std::ostringstream oss;
       oss<<filename_first<<plot_flag<<filename_last;   //don not add <<std::endl, or file name will be "A.csv "(cannot be read) instead of "A.csv"        
       std::ofstream wr(oss.str().c_str(),std::ios::out|std::ios::app);//.str() to string, .str().c_str() to char*
        //comment the txt output
       vectord q(2);
       q(0) = 0.0; q(1) = 0.0;
       /*
       struct timeval tp;
       gettimeofday(&tp, NULL);
       long int ms_ini = tp.tv_sec * 1000 + tp.tv_usec / 1000;
       */
       size_t n=50;
       //std::vector<std::pair<double,double> > xx;

       //private member
       //std::vector<double> x;
       //std::vector<std::vector<double> >xx,yy,zz,susu,slsl,cc;
       for(size_t k = 0; k<n; k++)
       {
           x.push_back(k*0.02);
       }
       std::vector<std::vector<double> > set_size;
       for(size_t k = 0; k<n; k++){
          set_size.push_back(x);
       }
       yy = set_size; zz = set_size; susu = set_size; slsl = set_size; cc = set_size;
       
       for(int i=0; i<n; i++)
       {
          vectord tmp_x0(2);
          vectord unnormalize_x0(2);
          tmp_x0(0) = x[i];
          tmp_x0(1) = 0;
          unnormalize_x0 = remapPoint(tmp_x0);
          xx.push_back(unnormalize_x0(0));  

          for(int j=0; j<n; j++)
            {
              q(0) = x[i];
              q(1) = x[j];
              //std::cout<<"q(0) q(1)  is "<<q(0)<<","<<q(1)<<std::endl;
              ProbabilityDistribution* pd = mModel->getPrediction(q);
              yy[i][j]   = pd->getMean();
              susu[i][j] = yy[i][j] + 2*pd->getStd();                       //Upper bound (95 %)
              slsl[i][j] = yy[i][j] - 2*pd->getStd();                       //Lower bound (95 %)
              cc[i][j]    = -mModel->evaluateCriteria(q);             //Criteria value
      //z[i] = mModel->evaluateSample(q);                //Target function true value
              
              //std::cout<<"before mapping the q is "<<q(0)<<","<<q(1)<<std::endl;
              q = remapPoint(q);
              //std::cout<<"after mapping the q is "<<q(0)<<","<<q(1)<<std::endl;
              wr<<q(0)<<","<<q(1)<<","<<yy[i][j]<<","<<susu[i][j]<<","<<slsl[i][j]<<","<<cc[i][j]<<std::endl;//comment the txt output

              /*unormalize x so that the plot's x axis we got can be lower bound to upper bound instead of 0 to 1*/
              if(i==n-1){
                vectord tmp_x1(2);
                vectord unnormalize_x1(2);
                tmp_x1(0) = x[i];
                tmp_x1(1) = x[j];
                unnormalize_x1 = remapPoint(tmp_x1);
                xx.push_back(unnormalize_x1(1));
              }
            }//x[i],x[j] is the corrdinate because it is two dimensional
       }
       
       for(size_t k=0;k<(plot_flag+mParameters.n_init_samples);k++)
       {
            wr<<AllX[2*k]<<","<<AllX[2*k+1]<<","<<AllY[k]<<std::endl;//so last several rows is sample point
       }
       //comment the txt output
       plot_flag++;
       /*
       gettimeofday(&tp, NULL);
       long int ms_end = tp.tv_sec * 1000 + tp.tv_usec / 1000; 
       long int diff = ms_end - ms_ini;
       std::cout<<"we spend "<<diff<<" msecs for 10000 mean"<<std::endl;
       */
      //std::cout<<"we are in 2D and x size, yysize, yy[0] size "<<x.size()<<" "<<yy.size()<<" "<<yy[0].size()<<std::endl;
   }
   //exit(0); 
  }
  

  void BayesOptBase::initializeOptimization()
  {
    // Posterior surrogate model
    mModel.reset(PosteriorModel::create(mDims,mParameters,mEngine));
    
    // Configure iteration parameters
    if (mParameters.n_init_samples <= 0)
      {
        mParameters.n_init_samples = 
          static_cast<size_t>(ceil(0.1*mParameters.n_iterations));	
      }
    
    size_t nSamples = mParameters.n_init_samples;
    size_t count = 0;
    size_t i     = 0;
    // Generate xPoints for initial sampling
    matrixd xPoints(nSamples,mDims);
    vectord yPoints(nSamples,0);
    //std::cout<<"initial x point is \n "<<xPoints<<std::endl;
    // Save generated xPoints before its evaluation
    //if(mParameters.his_address.c_str() == "generate initial points by running program")
    //{
        std::cout<<"generate initial points by runnning program"<<std::endl;
        generateInitialPoints(xPoints);//each time you use this you'll get different points
        // Save on each evaluation for safety reasons
        //generate points from 0~1, these points are also saved in the system. When use evaluateSampleInternal(row(xPoints,i))
        //points will be mapped to the bounding box range and then forward to evaluateSample() function
        //the points I saved is from evaluateSample(), they are all in the bounding box range not from 0~1 
        //if we want read those data as history directly, we need unmap it first, which means we need normalize the range
        //that's why we have unremapPoint(row(xPoints,count)) below.
        //unremapPoint(row(xPoints,count)) is added by myself based on the normalize function in bounding_box.hpp
        if(mDims == 1){
          for(int i = 0; i<nSamples; i++){
            vectord tmpx(1);
            tmpx = remapPoint(row(xPoints,i));
            AllX.push_back(tmpx(0));
          }
        }
        else if(mDims==2){
          for(int i = 0; i<nSamples; i++){
            vectord tmpx(2);
            tmpx = remapPoint(row(xPoints,i));
            AllX.push_back(tmpx(0));
            AllX.push_back(tmpx(1));
          }
        }
        for( i=0; i<yPoints.size(); i++)
        {
            yPoints[i] = evaluateSampleInternal(row(xPoints,i));
            //std::cout<<"y point is "<<yPoints[i] <<std::endl;
    //We clear the vector in the first iteration
            saveResponse(yPoints[i], i==0);
            AllY.push_back(yPoints[i]);
        } //most of the time we don't need plot the initial random sample points
          // count++;
    //}
    
    //else
    if(0)
    {
        std::cout<<"generate initial points by history"<<std::endl;
        std::ifstream read(mParameters.his_address.c_str());
        if(!read)
        {
            std::cout<<"the address is "<<mParameters.his_address.c_str()<<std::endl;
            std::cout<<"not a right address"<<std::endl;
            exit(0);
        }
        else
        {
            std::string st;
            while(getline(read,st))
            {
                 if(count==nSamples)
                     break;
                 std::istringstream sin(st);
                 std::vector<double> vc;
                 std::string st2;
                 while(getline(sin, st2, ','))
                 {
                     vc.push_back((double)atof(st2.c_str()));//be careful about the data you input in
                 }
                 yPoints[count] = vc[0];
                 saveResponse(yPoints[count], count==0);
                 //if(count<10)
                     //std::cout<<yPoints[count]<<",";
                 for(i=1;i<mDims+1;i++)//first one is SLAM output
                 {
                     xPoints(count,i-1) = vc[i];
                 }
                 row(xPoints,count) = unremapPoint(row(xPoints,count));//be careful
                 /*
                 if(count<10)
                 {
                   std::cout<<"point before remap is "<< row(xPoints,count)<<std::endl;
                   remapPoint(row(xPoints,count));
                   std::cout<<"point after remap is "<< row(xPoints,count)<<std::endl;                   
                 } 
                 */
                 count++;                  
            }
         read.close();
        }
    
    } 
    
    //std::cout<<"initial x point is after generation \n "<<xPoints<<std::endl;
    //std::cout<<"initial y point is \n"<<yPoints<<std::endl;
    //exit(0);
    //saveInitialSamples(xPoints);
    //mModel->setSamples(xPoints); //原程序在这儿save并set了初始sample,但是这样下面在for循环中重新设置会出错，所以在for循环后，确定了有效的采样再设置
    
    
    // Put samples into model
    saveInitialSamples(xPoints);
    mModel->setSamples(xPoints);
    mModel->setSamples(yPoints);
 
    if(mParameters.verbose_level > 0)
      {
        mModel->plotDataset(logDEBUG);
      }
    
    mModel->updateHyperParameters();    
    mModel->fitSurrogateModel();
    mCurrentIter = 0;

    mCounterStuck = 0;
    mYPrev = 0.0;

    //exit(0);
    
  }

  vectord BayesOptBase::getFinalResult()
  {
    return remapPoint(getPointAtMinimum());
  }


  // SAVE-RESTORE INTERFACE
  void BayesOptBase::saveOptimization(BOptState &state)
  {   
    // BayesOptBase members
    state.mCurrentIter = mCurrentIter;
    state.mCounterStuck = mCounterStuck;
    state.mYPrev = mYPrev;

    state.mParameters = mParameters;

    // Samples
    state.mX = mModel->getData()->mX;
    state.mY = mModel->getData()->mY;
  }

  void BayesOptBase::restoreOptimization(BOptState state)
  {
    
    std::cout<<"restore optimization............."<<std::endl;
    // Restore parameters
    mParameters = state.mParameters; 
    
    // Posterior surrogate model
    mModel.reset(PosteriorModel::create(mDims, mParameters, mEngine));
    
    // Load samples, putting mX vecOfvec into a matrixd
    matrixd xPoints(state.mX.size(),state.mX[0].size());
    vectord yPoints(state.mX.size(),0);
    for(size_t i=0; i<state.mX.size(); i++)
      {
        row(xPoints, i) = state.mX[i];
        if(i < state.mY.size())
	  {
            yPoints[i] = state.mY[i];
	  }
	else
	  {
	    // Generate remaining initial samples saving in each evaluation	    
	    yPoints[i] = evaluateSampleInternal(row(xPoints,i));
	    saveResponse(yPoints[i], false);
	  }
      }
    
    // Set loaded and generated samples
    mModel->setSamples(xPoints,yPoints);
        
    if(mParameters.verbose_level > 0)
    {
        mModel->plotDataset(logDEBUG);
    }
    
    // Calculate the posterior model
    mModel->updateHyperParameters();
    mModel->fitSurrogateModel();
    
    mCurrentIter = state.mCurrentIter;
    mCounterStuck = state.mCounterStuck;
    mYPrev = state.mYPrev;
    
    // Check if optimization has already finished
    if(mCurrentIter >= mParameters.n_iterations)
      {
        FILE_LOG(logINFO) << "Optimization has already finished, delete \"" 
			  << mParameters.load_filename 
			  << "\" or give more n_iterations in parameters."; 
      }
  }

  
  // GETTERS AND SETTERS
  // Potential inline functions. Moved here to simplify API and header
  // structure.
  ProbabilityDistribution* BayesOptBase::getPrediction(const vectord& query)
  { return mModel->getPrediction(query); };
  
  const Dataset* BayesOptBase::getData()
  { return mModel->getData(); };

  Parameters* BayesOptBase::getParameters() 
  {return &mParameters;};

  double BayesOptBase::getValueAtMinimum()
  { return mModel->getValueAtMinimum(); };

  double BayesOptBase::evaluateCriteria(const vectord& query)
  {
    if (checkReachability(query)) return mModel->evaluateCriteria(query);
    else return 0.0;
  }

  size_t BayesOptBase::getCurrentIter()
  {return mCurrentIter;};

  std::vector<double>& BayesOptBase::get1Dx(){
    return x;
  }
  std::vector<double>& BayesOptBase::get1Dmean(){
    return y;
  };
  std::vector<double>& BayesOptBase::get1DStdLowerBound(){
    return sl;
  };
  std::vector<double>& BayesOptBase::get1DStdUpperBound(){
    return su;
  };
  std::vector<double>& BayesOptBase::get1DAcquisitionFunctionValue(){
    return c;
  };

  std::vector<double>& BayesOptBase::get2Dx(){
    return xx;
  };
  std::vector<std::vector<double> >& BayesOptBase::get2Dmean(){
    return yy;
  };
  std::vector<std::vector<double> >& BayesOptBase::get2DStdLowerBound(){
    return slsl;
  };
  std::vector<std::vector<double> >& BayesOptBase::get2DStdUpperBound(){
    return susu;
  };
  std::vector<std::vector<double> >& BayesOptBase::get2DAcquisitionFunctionValue(){
    return cc;
  };
  std::vector<double>& BayesOptBase::getSampleX(){
    return AllX;
  };
  std::vector<double>& BayesOptBase::getSampleY(){
    return AllY;
  }
  

  // PROTECTED
  vectord BayesOptBase::getPointAtMinimum() 
  { return mModel->getPointAtMinimum(); };

  double BayesOptBase::evaluateSampleInternal( const vectord &query )
  { 
    //FILE_LOG(logDEBUG)<<"the point is "<<query[0];
    const double yNext = evaluateSample(remapPoint(query)); 
    //cout<<"yNext is "<<yNext<<endl;
    if (yNext == HUGE_VAL)
      {
	throw std::runtime_error("Function evaluation out of range");
      }
    return yNext;
  }; 



  
  
  void BayesOptBase::plotStepData(size_t iteration, const vectord& xNext,
				     double yNext)
  {
    if(mParameters.verbose_level >0)
      { 
	FILE_LOG(logINFO) << "Iteration: " << iteration+1 << " of " 
			  << mParameters.n_iterations << " | Total samples: " 
			  << iteration+1+mParameters.n_init_samples ;
	FILE_LOG(logINFO) << "Query: "         << remapPoint(xNext); ;
	FILE_LOG(logINFO) << "Query outcome: " << yNext ;
	FILE_LOG(logINFO) << "Best query: "    << getFinalResult(); 
	FILE_LOG(logINFO) << "Best outcome: "  << getValueAtMinimum();
      }
  } //plotStepData


  void BayesOptBase::saveInitialSamples(matrixd xPoints)
  {
    // Save state if required
    if(mParameters.load_save_flag == 2 || mParameters.load_save_flag == 3)
      {
        BOptState state;
        saveOptimization(state);
        
        // Overwrite the state with initial samples so far
        state.mX.clear();
        std::cout<<"overwrite the initial sample"<<std::endl;
        for(size_t i=0; i<xPoints.size1(); i++)
	  {
            state.mX.push_back(row(xPoints,i));
	  }
        state.saveToFile(mParameters.save_filename);
      }
  }


  void BayesOptBase::saveResponse(double yPoint, bool clear)
  {
    // Save state if required
    if(mParameters.load_save_flag == 2 || mParameters.load_save_flag == 3)
      {
        BOptState state;
        saveOptimization(state);
	if (clear)
	  {
	    state.mY.clear();
	  }
	utils::append(state.mY,yPoint);
        state.saveToFile(mParameters.save_filename);
      }
  }






  // PRIVATE MEMBERS
  vectord BayesOptBase::nextPoint()
  {
    //Epsilon-Greedy exploration (see Bull 2011)
    if ((mParameters.epsilon > 0.0) && (mParameters.epsilon < 1.0))
      {
	randFloat drawSample(mEngine,realUniformDist(0,1));
	double result = drawSample();
	FILE_LOG(logINFO) << "Trying random jump with prob:" << result;
	if (mParameters.epsilon > result)
	  {
	    FILE_LOG(logINFO) << "Epsilon-greedy random query!";
	    return samplePoint();
	  }
      }

    vectord Xnext(mDims);    
    // GP-Hedge and related algorithms
    if (mModel->criteriaRequiresComparison())
      {
	bool changed = true;

	mModel->setFirstCriterium();
	while (changed)
	  {
	    findOptimal(Xnext);
	    changed = mModel->setNextCriterium(Xnext);
	  }
	std::string name = mModel->getBestCriteria(Xnext);
	FILE_LOG(logINFO) << name << " was selected.";
      }
    else  // Standard "Bayesian optimization"
      {
	FILE_LOG(logDEBUG) << "------ Optimizing criteria ------";
	findOptimal(Xnext);
      }
    return Xnext;
  }


} //namespace bayesopt

