/*-------------------------------------------------------------------------*/
/*
  This is the interface to BONMIN for solving simple MIQPs.
  
*/

/*-------------------------------------------------------------------------*/

#include "Bonmin.hpp"
#include <R.h>	

extern "C" {

void R_bonmin_solve(double *c, double *b, double *A, double *x, double *opt)
{

  using namespace Ipopt;
  using namespace Bonmin;
  SmartPtr<MyTMINLP> tminlp = new MyTMINLP;
  

  BonminSetup bonmin;
  bonmin.initializeOptionsAndJournalist();
  //Register an additional option
  bonmin.roptions()->AddStringOption2("print_solution","Do we print the solution or not?",
                                 "yes",
                                 "no", "No, we don't.",
                                 "yes", "Yes, we do.",
                                 "A longer comment can be put here");
  
  
  
  // Here we can change the default value of some Bonmin or Ipopt option
  bonmin.options()->SetNumericValue("bonmin.time_limit", 5); //changes bonmin's time limit
  bonmin.options()->SetStringValue("mu_oracle","loqo");
  
  //Here we read several option files
  //bonmin.readOptionsFile("Mybonmin.opt");
  //bonmin.readOptionsFile();// This reads the default file "bonmin.opt"
  
  // Options can also be set by using a string with a format similar to the bonmin.opt file
  bonmin.readOptionsString("bonmin.algorithm B-BB\n");
  
  // Now we can obtain the value of the new option
  /*int printSolution;
  bonmin.options()->GetEnumValue("print_solution", printSolution,"");
  if(printSolution == 1){
    tminlp->printSolutionAtEndOfAlgorithm();
  }
  */
  
  //Now initialize from tminlp
  bonmin.initialize(GetRawPtr(tminlp));



  //Set up done, now let's branch and bound
  double time1 = CoinCpuTime();
  try {
    Bab bb;
    bb(bonmin);//process parameter file using Ipopt and do branch and bound using Cbc


  }
  catch(TNLPSolver::UnsolvedError *E) {
    //There has been a failure to solve a problem with Ipopt.
    std::cerr<<"Ipopt has failed to solve a problem"<<std::endl;
  }
  catch(OsiTMINLPInterface::SimpleError &E) {
    std::cerr<<E.className()<<"::"<<E.methodName()
	     <<std::endl
	     <<E.message()<<std::endl;
  }
  catch(CoinError &E) {
    std::cerr<<E.className()<<"::"<<E.methodName()
	     <<std::endl
	     <<E.message()<<std::endl;
  }


}


}
