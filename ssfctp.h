/********************************************************************************

FILE     : ssfctp.h
VERSION  : 7.52
DATE     : June 14, 2007
LANGUAGE : c, header file
AUTHOR   : Andreas Klose
SUBJECT  : header file of module ssfctp.c: 
           Solution  methods for the single-sink, fixed-charge transportation 
	   problem.

********************************************************************************/

#ifndef __ssfctp_H
#define __ssfctp_H
#endif

#ifdef CPLEX
#ifndef __CPXDEFS_H
#include <cplex.h>
#endif
#else
#define CPXMIP_TIME_LIM_FEAS 30212
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*-----------------------------------------------------------------------------*/
void ssfctp_param ( int param, double user_value );
/*-------------------------------------------------------------------------------
Purpose: Set parameter "param" to value user_value

Parameters:
- param : integer specifying the paramter to change
          param=0 -> change relative optimality tolerance;
                     a solution of value U will be deemed optimal if 
		     L + U*epgap > U, where epgap is the percentage 
		     optimality tolerance and L the lower bound
		     (default value is 1.0E-6)
          param=1 -> set/change the time limit in seconds for procedures
	             ssfctp_cpx and ssfctp_enu; default is 1.0E75
                     
                     REMARK: For activating this option the file
                     sscftp.c must be compiled with the additional
                     option -D__CHKTIM. For checking computation time,
                     the FORTRAN routine "setreftim()" provided by
                     the file rtime.for is used. To this end, the
                     object files must be linked with the g77 library
                     
- user_value: a double specifying the new value of the parameter		       		     
       	  
-------------------------------------------------------------------------------*/	 

/*-----------------------------------------------------------------------------*/
void ssfctp_gr ( int n, int demand, int* cap, double* fcost, double* tcost, 
                 int adapt, int* x, double* upbnd ); 
/*-------------------------------------------------------------------------------
Purpose: Applies the greedy or adaptive greedy method to the SSFCTP instance

Parameters:
- n      : number of suppliers
- demand : demand of the sink
- cap    : array of integers of length of at least n. 
           cap[j] is the capacity of the arc linking supplier j and the sink.
- fcost  : array of doubles of length of at least n. 
           fcost[j] is the fixed cost of arc j
- tcost  : array of doubles of length of at least n. 
           tcost[j] is the unit cost of flow on arc j.
- adapt  : adapt=0 -> just do greedy; otherwise do "adaptive greedy"
- x      : array of integers of length of at least n. 
           On output, x[j] is the computed flow on arc j. 
- upbnd  : objective function value of the computed solution       	  
-------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------*/
void ssfctp_gl ( int n, int demand, int* cap, double* fcost, double* tcost, 
                 int preproc, int* x, double* upbnd ); 
/*-------------------------------------------------------------------------------
Purpose: Adjusts the heuristic of Gens and Levner for the min-knapsack problem
         to the case of the SSFCTP

Parameters:
- n      : number of suppliers
- demand : demand of the sink
- cap    : array of integers of length of at least n. 
           cap[j] is the capacity of the arc linking supplier j and the sink.
- fcost  : array of doubles of length of at least n. 
           fcost[j] is the fixed cost of arc j
- tcost  : array of doubles of length of at least n. 
           tcost[j] is the unit cost of flow on arc j.
- preproc: if > 0, the reduction test is performed before executing the heuristic
- x      : array of integers of length of at least n. 
           On output, x[j] is the computed flow on arc j. 
- upbnd  : objective function value of the computed solution       	  
-------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------*/
void ssfctp_gli( int n, int demand, int* cap, double* fcost, double* tcost, 
                 int preproc, int* x, double* upbnd ); 
/*-------------------------------------------------------------------------------
Purpose: Improved Gens/Levner heuristic adjusted to the case of the SSFCTP.

Parameters:
- n      : number of suppliers
- demand : demand of the sink
- cap    : array of integers of length of at least n. 
           cap[j] is the capacity of the arc linking supplier j and the sink.
- fcost  : array of doubles of length of at least n. 
           fcost[j] is the fixed cost of arc j
- tcost  : array of doubles of length of at least n. 
           tcost[j] is the unit cost of flow on arc j.
- preproc: if > 0, the reduction test is performed before actually executing
           the heuristic
- x      : array of integers of length of at least n. 
           On output, x[j] is the computed flow on arc j. 
- upbnd  : objective function value of the computed solution       	  
-------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------*/
void ssfctp_gli2( int n, int demand, int* cap, double* fcost, double* tcost, 
                  int preproc, int* x, double* upbnd ); 
/*-------------------------------------------------------------------------------
Purpose: Second possible way of improving the Gens-Levner type heuristic 

Parameters:
- n      : number of suppliers
- demand : demand of the sink
- cap    : array of integers of length of at least n. 
           cap[j] is the capacity of the arc linking supplier j and the sink.
- fcost  : array of doubles of length of at least n. 
           fcost[j] is the fixed cost of arc j
- tcost  : array of doubles of length of at least n. 
           tcost[j] is the unit cost of flow on arc j.
- preproc: if > 0, the reduction test is performed before executing the heuristic.           
- x      : array of integers of length of at least n. 
           On output, x[j] is the computed flow on arc j. 
- upbnd  : objective function value of the computed solution       	  
-------------------------------------------------------------------------------*/
 
/*-----------------------------------------------------------------------------*/
int ssfctp_dp( int n, int demand, int* cap, double* fcost, double* tcost,
               int preproc, int* x, double* objval );
/*-------------------------------------------------------------------------------
PURPOSE: 
Dynamic programming procedure of Alidaee and Kochenberger to solve the 
single-sink, fixed-charge transportation problem.

PARAMETERS:
- n      : number of suppliers
- demand : demand of the sink
- cap    : array of integers of length of at least n. 
           cap[j] is the capacity of the arc linking supplier j and the sink.
- fcost  : array of doubles of length of at least n. 
           fcost[j] is the fixed cost of arc j
- tcost  : array of doubles of length of at least n. 
           tcost[j] is the unit cost of flow on arc j.
- preproc: preproc = 0 -> no preprocessing is applied
           preproc > 0 -> problem reduction based on reduced cost is applied
	   preproc > 1 -> additionally, bounds on L based on reduced cost
	                  information are computed.	   
- x      : array of integers of length of at least n. 
           On output, x[j] is the optimal flow on arc j. 
- objval : pointer to a double that on output will contain the cost of 
           an optimal solution.

RETURN VALUE: 
0 if no error occured and 1001 in case of memory problems.

LITERATURE:
Alidaee B, Kochenberger GA (2005) A note on a simple dynamic programming 
approach to the single-sink, fixed-charge transportation problem. 
Transp. Sci 39(1), 140-143.
-------------------------------------------------------------------------------*/

#ifdef CPLEX
/*-----------------------------------------------------------------------------*/
int ssfctp_cpx( CPXENVptr Env, int n, int demand, int* cap, double* fcost, 
                double* tcost, int* x, double* objval );
/*-------------------------------------------------------------------------------
PURPOSE: 
Solves the single-sink, fixed-charge transportation problem by means of CPLEX's
MIP solver. Before calling the solver, the LP relaxation is solved, heuristic
solutions are computed, and the problem is reduced.

PARAMETERS:
- Env    : pointer to the CPLEX environment as returned by CPXopenCPLEX
           ( may be NULL, in that case this routine tries to open CPLEX )
- n      : number of suppliers
- demand : demand of the sink
- cap    : array of integers of length of at least n. 
           cap[j] is the capacity of the arc linking supplier j and the sink.
- fcost  : array of doubles of length of at least n. 
           fcost[j] is the fixed cost of arc j
  tcost  : array of doubles of length of at least n. 
           tcost[j] is the unit cost of flow on arc j.
  x      : array of integers of length of at least n. 
           On output, x[j] is the optimal flow on arc j.
  objval : pointer to a double that on output will contain the cost of 
           an optimal solution.

RETURN VALUE: 
0 if no error occured and a positive integer otherwise
-------------------------------------------------------------------------------*/
#endif		

/*-----------------------------------------------------------------------------*/
int ssfctp_enu( int n, int demand, int* cap, double* fcost, double* tcost, 
                int preproc, int* x, double* objval );
/*-------------------------------------------------------------------------------
PURPOSE: 
Implicit enumeration algorithm of Herer et al. to solve the single-sink, 
fixed-charge transportation problem. 

PARAMETERS:
- n      : number of suppliers
- demand : demand of the sink
- cap    : array of integers of length of at least n. 
           cap[j] is the capacity of the arc linking supplier j and the sink.
- fcost  : array of doubles of length of at least n. 
           fcost[j] is the fixed cost of arc j
- tcost  : array of doubles of length of at least n. 
           tcost[j] is the unit cost of flow on arc j.
- preproc: integer paramter:
           - preproc=0 -> no preprocessing measures are taken
	   - preproc>0 -> the LP solution and a heuristic solution based on
	                  this LP solution is computed
           - preproc>1 -> the problem is reduced by reduction tests
- x      : array of integers of length of at least n. 
           On output, x[j] is the optimal flow on arc j. 
- objval : pointer to a double that on output will contain the cost of an 
           optimal solution.

RETURN VALUE: 
0 if no error occured, 107 if time limit reached, and 1001 in case of 
memory problems.

LITERATURE:
Herer YT, Rosenblatt MJ, Hefter I (1996). Fast algorithms for single-sink
fixed charge transportation problems with applications to manufacturing
and transportation. Transp. Sci., 30, 276-290.
-------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------*/
int ssfctp_mt1( int n, int demand, int* cap, double* fcost, double* tcost, 
                int preproc, int* x, double* objval );
/*-------------------------------------------------------------------------------
PURPOSE: 
Implicit enumeration algorithm for the single-sink fixed-charge transportation 
problem. This enumeration schemes relies on a branch-and-bound algorithm of
Martello & Toth for the binary knapsack problem.

PARAMETERS:
- n      : number of suppliers
- demand : demand of the sink
- cap    : array of integers of length of at least n. 
           cap[j] is the capacity of the arc linking supplier j and the sink.
- fcost  : array of doubles of length of at least n. 
           fcost[j] is the fixed cost of arc j
- tcost  : array of doubles of length of at least n. 
           tcost[j] is the unit cost of flow on arc j.
- preproc: integer paramter:
           - preproc=0 -> no preprocessing measures are taken
	   - preproc>0 -> the LP solution and a heuristic solution based on
	                  this LP solution is computed
           - preproc>1 -> the problem is reduced by reduction tests
- x      : array of integers of length of at least n. 
           On output, x[j] is the optimal flow on arc j. 
- objval : pointer to a double that on output will contain the cost of an 
           optimal solution.

RETURN VALUE: 
0 if no error occured and 107 if time limit reached.

-------------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif
