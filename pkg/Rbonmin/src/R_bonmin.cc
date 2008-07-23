/*-------------------------------------------------------------------------*/
/*
  This is the interface to BONMIN for solving simple MIQPs.
  
*/

/*-------------------------------------------------------------------------*/

#include "Bonmin.hpp"
#include <R.h>	

extern "C" {

void R_bonmin_solve(int *test)
{

  using namespace Ipopt;
  SmartPtr<TMINLP> tminlp = new MyTMINLP;
  IpoptInterface nlpSolver(tminlp);
  
  //Option can be set here directly either to bonmin or ipopt
  nlpSolver.retrieve_options()->SetNumericValue("bonmin.time_limit", 1); //changes bonmin's time limit
  nlpSolver.retrieve_options()->SetStringValue("mu_oracle","loqo");

  // we can also try and read an option file (can eventually change options set before, option file always have priority)
  nlpSolver.readOptionFile("My_bonmin.opt");

  //Set up done, now let's branch and bound
  double time1 = CoinCpuTime();
  try {
    BonminCbcParam par;
    BonminBB bb;
    par(nlpSolver);

    bb(nlpSolver, par);//process parameter file using Ipopt and do branch and bound

    std::cout.precision(10);

    std::string message;
    if(bb.mipStatus()==BonminBB::FeasibleOptimal) {
      std::cout<<"\t\"Finished\"\t";
      message = "\nbonmin: Optimal";
    }
    else if(bb.mipStatus()==BonminBB::ProvenInfeasible) {
      std::cout<<"\t\"Finished\"\t";
      message = "\nbonmin: Infeasible problem";
    }
    else if(bb.mipStatus()==BonminBB::Feasible) {
      std::cout<<"\t\"Not finished\""<<"\t";
      message = "\n Optimization not finished.";
    }
    else if(bb.mipStatus()==BonminBB::NoSolutionKnown) {
      std::cout<<"\t\"Not finished\""<<"\t";
      message = "\n Optimization not finished.";
    }
    std::cout<<CoinCpuTime()-time1<<"\t"
    <<bb.bestObj()<<"\t"
    <<bb.numNodes()<<"\t"
    <<bb.iterationCount()<<"\t"
    <<std::endl;

  }
  catch(IpoptInterface::UnsolvedError &E) {
    //There has been a failure to solve a problem with Ipopt.
    std::cerr<<"Ipopt has failed to solve a problem"<<std::endl;
  }
  catch(IpoptInterface::SimpleError &E) {
    std::cerr<<E.className()<<"::"<<E.methodName()
	     <<std::endl
	     <<E.message()<<std::endl;
  }
  catch(CoinError &E) {
    std::cerr<<E.className()<<"::"<<E.methodName()
	     <<std::endl
	     <<E.message()<<std::endl;
  }


  /* from SYMPHONY
   int i;
   

   sym_environment *env = sym_open_environment();
   sym_set_int_param(env, "verbosity", -2);


   char * int_vars = (char *) malloc (sizeof(char) * (*n_cols));
   for(i=0; i < (*n_cols); i++)
       if(is_int[i] == 1)
	   int_vars[i] = TRUE;
       else
	   int_vars[i] = FALSE;

   sym_explicit_load_problem(env, *n_cols, *n_rows, start, index, value,
			     col_lb, col_ub, int_vars, objective, NULL,
			     *row_sense, row_rhs, row_range, TRUE);
 			     			     

   sym_solve(env);


   double * solution = (double *) malloc (sizeof(double) * (*n_cols));
   double objective_value = 0.0;

   sym_get_col_solution(env, solution);
   sym_get_obj_val(env, &objective_value);
   
   *obj_final = objective_value;
   for(i=0; i < (*n_cols); i++)	
     sol_final[i] = solution[i];

   *solve_status = sym_get_status(env);

   sym_close_environment(env);

*/

};

}
