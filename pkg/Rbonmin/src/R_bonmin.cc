/*-------------------------------------------------------------------------*/
/*
  This is the interface to BONMIN for solving simple MIQPs.
  
*/

/*-------------------------------------------------------------------------*/

#include "bonmin.h"
#include "R.h"	

extern "C" {

void R_bonmin_solve(int *n_cols, int *n_rows, int *start, int *index, 
		      double *value, double *col_lb, double *col_ub,
		      int* is_int, double *objective, double *obj2,
		      char **row_sense, double *row_rhs,
		      double *row_range, double *obj_final,
		      double *sol_final, int *solve_status)
{



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
