/********************************************************************************

FILE    : ssfctp.c

VERSION : 7.54
DATE    : June 26, 2008
LANGUAGE: C
AUTHOR  : Andreas Klose
SUBJECT : Algorithms for solving the Single-Sink, Fixed Charge Transportation
          Problem

EXTERNALS (in addition to routines provided by compiler libraries)

CHANGE_LOG
----------
 5th Nov 2005: Project start
27th Nov 2005: G[i][L] not explicitely set to infty in case that L is outside
               the range [Lmin,Lmax]. Macro introduced that returns infty
	              in case that L is outside this range
28th Nov 2005: Routine that solves the LP relaxation and applies reduction tests
 2nd Dec 2005: Routine ssfctp_red: reduction tests based on reduced costs
 3rd Dec 2005: Routine ssfctp_bol: bounds on L
 8th Dec 2005: Some bugs removed; second heuristic solution procedure implemented
 9th Dec 2005: Routine ssfctp_cpx implemented
15th Dec 2005: Routine ssfctp_enu implemented 
16th Dec 2005: Domination rule 1 implemented (part of routine ssfctp_enu)
 3rd Jan 2006: Domination rule 21 implemented (part of routine ssfctp_enu)
 6th Jan 2006: Implementation of LP lower bound for given partial solutions
               (part of routine ssfctp_enu)
19th Jan 2006: Some bugs in ssfctp_enu removed	       
23th Jan 2006: Started with implementation of procedure ssfct_mt1
24th Jan 2006: Small changes to procedure ssfctp_enu (only one dominance matrix
               instead of two upper diagonal matrices) and implementation of
	              procedure ssfctp_flow that determines without sorting 
	              the optimal flows for a given set of selected arcs/suppliers
26th Jan 2006: Procedure ssfctp_mt1: Domination matrix implemented 	       
27th Jan 2006: Procedure ssfctp_mt1: Domination matrix removed again since
               the domination tests based on these matrix are usually useless
	              if variables are sorted according to nondecreasing linearized
	              costs ( in this case if i < j and j strictly dominates i, than
	              it follows that cap[i]=cap[j] and lpcst[i]=lpcst[j], something
	              which should happen very rarely)
30th Jan 2006: Procedure ssfctp_dp: Application of preprocessing (reduced cost
               test, bounds on L) made optional. To this end, procecure ssfctp_nbol
	              included. ssfctp_nbol just computes the "natural" bounds on L.
31th Jan 2006: Procedure ssfctp_enu and ssfctp_mt1: integration of
	              information on minimum required flow on arcs by means of 
	              adjusting the problem data accordingly.
 1st Feb 2006: Procedure ssfctp_red: reduction test changed. If it is recognized
               that some flow variables must be positive, all these variables
	              except the one with largest unit cost must be at the upper bound.
	              Problem is reduced accordingly, all procedures changed
	              accordingly.	              
 2nd Feb 2006: Removed some bugs that were created on 1st Feb 2006	       
 5th Feb 2006: Relative/percentage optimality tolerance introduced.
 6th Feb 2006: Proc ssfctp_mt1: bug concerning uncancelled old implications
               removed
 7th Feb 2006: Proc ssfctp_cpx and ssfctp_enu: time limit for computation 
               introduced	       
 9th Feb 2006: Proc ssfctp_mt1: small bug eliminated (goto to an false point
               in case that upbnd equals lobnd0) 	       
14th Feb 2006: Proc ssfctp_red: additionally parameter included that equals 1
               if after reduction maxflows contradict the flows in the solution
	              (which means that this solution must be optimal)	       
28th Feb 2006: Optional time limit included for procedures ssfctp_enu and
               ssfctp_cpx
31th Mar 2006: Optional time limit included for procedure ssfctp_mt1	       	       
5th  Nov 2006: Heuristic of Gens and Levner for the min-knapsack adjusted
               for the SSFCTP
25th Nov 2006: Second version of improved Gens-Levner type heuristic added	       
 5th Jun 2007: Bug removed in ssfctp_mt1():
               "if ( fc[ii] + tc[ii]*mflow[ii] + ZERO < tc[jj]*mflow[jj] )
	              impl[j] = bvar+1;  (NOTE x[ii]=0 implies x[jj] < mflow[jj] ) "
	              is to be replaced by	 
	              "if ( fc[ii] + tc[ii]*mflow[ii] + ZERO < tc[jj]*mflow[ii] )
	                 impl[j] = bvar+1;" !!!
14th Jun 2007: General bug corrected: If maxflow[j] is reduced for some j
               of relative cost larger than those of the critical supplier,
	              we have to reorder the suppliers according to the new
	              relative cost e_j=c_j + f_j/maxflow[j]. This has otherwise
	              impact on the lower on L (ssfctp_bol) and may lead to wrong
	              values of Lmin and thus wrong solutions produced by the
	              dynamic programming procedure.
14th Jun 2007: In procedure ssfctp_lp, instruction upbnd2=ssfctp_flow(...):
               q replaced with n_arcs!	       
26th Jun 2008: Some caution added for handling trival cases.	       
********************************************************************************/

#include <stdlib.h>
#include <math.h>
#ifdef __CHKTIM
#include <rtime.h>
#endif

#ifdef CPLEX
#include <string.h>
#include <cplex.h>
#else
#define CPXMIP_TIME_LIM_FEAS 107
#endif

#define INFTY 1.0E31
#define ZERO 1.0E-05
#define ONE  0.99999
#define IDIV( a, b ) (  (a) % (b) == 0 ? ( (a)/(b)) : ( (a)/(b)+1 ) )
#define MIN( x, y ) ( (x) < (y) ? (x) : (y) )
#define MAX( x, y ) ( (x) > (y) ? (x) : (y) )
#define GFUNC( i, L ) ( ((L >= Lmin[i])&&(L<=Lmax[i])) ? (G_prev[L]) : (INFTY) )

#define TILIM 1.0E75 /* default value for the time limit */
#define EPGAP 1.0E-6 /* default value for the percentage optimality tolerance */
#define LBgUB( L, U ) ( ( (L)+(U*epgap) > U ) ? (1) : (0) )

typedef struct HList { /* List of H-function values mentioned above */
  int r;               /* argument of the function H */
  double H;            /* function value H_i(r)      */
  struct HList* pred;  /* pointer to previous pair */
  struct HList* succ;  /* pointer to next pair */
} HList;

typedef struct Tnode { /* structure to the data of a supplier node      */
  double fc, tc;       /* fixed cost and unit transporation cost        */
  int    lflow,uflow;  /* lower and upper bound on the flow on this arc */
  int    idx;          /* index of the node                             */
} Tnode;

double* SSFCTP_ptr = NULL; /* pointer to array of doubles thas has to be sorted*/
int     SSFCTP_sgn = 1;    /* sorting direction: 1->ascending, -1->descending  */

static double epgap=EPGAP;  /* percentage optimality tolerance value           */ 

#ifdef __CHKTIM
static double tilim=TILIM; /* time limit                                      */
#endif

/*-----------------------------------------------------------------------------*/

void ssfctp_param ( int param, double user_value ) {
/* Set a parameter "param" to the value "user_value" */

  if ( param == 0 ) {
    if ( ( user_value > 1.0E-11 ) && ( user_value < 1.0 ) ) epgap = user_value;
  } 
#ifdef __CHKTIM
  else if ( user_value > 0.0 ) tilim = user_value;  
#endif  
}

/*-----------------------------------------------------------------------------*/

int SSFCTP_compare( const int *order1, const int *order2 ) {
/* Comparison function */

  if ( SSFCTP_ptr[*order1] > SSFCTP_ptr[*order2] ) 
    return( SSFCTP_sgn );
  else if ( SSFCTP_ptr[*order1] < SSFCTP_ptr[*order2] ) 
    return( -SSFCTP_sgn );
  else
    return( 0 );  

}

/*-----------------------------------------------------------------------------*/

static 
double ssfctp_dta ( int n, int demand, double* tcost, double* fcost, double** tc, 
                    double** fc ) {
/*
  Description:
  ------------
  Transforms the data of an instance of the SSFCTP in such a way that all cost
  data are nonnegative.
  
  Parameters:
  ----------
  - n     : number of suppliers
  - demand: the demand of the sink 
  - tcost : array of doubles of length of at least n containing unit 
            transportation cost
  - fcost : array of doubles of lenght of at least n containing fixed-charges
  - tc    : pointer to an array of doubles. If tcost[j]>=0 for j=1,...n,
            tc points to tcost on output. Otherwise new memory for nonnegativ
	           transportation costs is created and tc points to this new array.  	    
  - fc    : pointer to an array of doubles. If fcost[j]>=0 for j=1,...n,
            fc points to tcost on output. Otherwise new memory for nonnegative
	    fixed costs is created and fc points to this new array.  	    
	    
  Return value: constant term in objective function in case that the data
  are transformed (otherwise 0 )	       	    
*/
  int    j, newfc;
  double cmin, zcnst=0.0;
  
  for ( j=newfc=0, cmin=tcost[0]; j < n; j++ ) { 
    cmin = MIN( cmin, tcost[j] );
    if ( fcost[j] < 0.0 ) newfc = 1;
  } 
  *fc = fcost; 
  if ( newfc ) {
    *fc = (double*) calloc ( n, sizeof(double) );
    for ( j=0; j < n; j++ ) {
      if ( fcost[j] > 0.0 ) 
        (*fc)[j] = fcost[j];
      else { 
        (*fc)[j] = 0.0; zcnst += fcost[j];
      }
    }  
  } 
  *tc = tcost;
  if ( cmin < 0.0 ) {
    zcnst += cmin*demand;
    *tc = (double*) calloc( n, sizeof(double) );
    for ( j=0; j < n; j++ ) (*tc)[j] = tcost[j] - cmin;
  }
  
  return( zcnst );

}

/*-----------------------------------------------------------------------------*/

static
int ssfctp_trivial( int n, int demand, int* cap, double* fcost, double* tcost, 
                    int* x, double* objval ) {
/* Description:
   ------------
   Checks if the SSFCTP instance is trivially solvable by supplying the sink's
   demand from node i where i has minimum cost c*=f_i + c_i D and maximum
   capacity under all nodes with minimum supply cost c*.
  
   Scope: Intern
   ------
   
   Parameters:
   ----------
   n      : number of suppliers
   demand : the sink's demand
   cap    : array of integers of size of at least n. cap[j] is the capacity the 
            arc linking supplier j with the sink.
   fcost  : array of doubles. fcost[j] is the fixed cost of arc j
   tcost  : array of doubles. tcost[j] is the unit cost of flow on arc j
   x      : array of integers. On output x[j] is the optimal flow on arc j
   objval : pointer to a double containing the cost of an optimal solution
   
   Return value: 1 if trivially solvable and 0 otherwise
   ------------		    
*/
   int    i=0, j, maxk=cap[0], lucky;
   double cost, mincost=fcost[0]+tcost[0]*demand;
   
   cap[0] = MIN( cap[0], demand );
   for ( j=1; j < n; j++ ) {
     cap[j] = MIN( cap[j], demand );
     cost   = fcost[j] + tcost[j]*demand;
     if ( cost < mincost ) 
       mincost = cost, maxk=cap[j], i=j; 
     else if ( ( !(cost > mincost) ) && (cap[j] > maxk) )
       mincost = cost, maxk=cap[j], i=j;
   }
   lucky = ( maxk >= demand );
   if ( lucky ) {
     if ( objval ) *objval = mincost;
     for ( j=0; j < n; j++ ) x[j]=0;
     x[i]=demand;
   }
   return( lucky );
}		    

/*-----------------------------------------------------------------------------*/

static 
double ssfctp_flow( int n, int demand, int* cap, double* fcost, double* tcost,
                    int narcs, int* sel_arc, int* flow ) {
/* Description:
   ------------
   Determines the optimal flows if the suppliers/arcs order[0 ... narcs-1] are
   selected.
   
   Parameters:
   ----------
   - n      : number of suppliers
   - demand : the sink's demand
   - cap    : array of integers of length of at least n containing arc capacities
   - fcost  : array of doubles of length of at least n containing fixed-charges
   - tcost  : array of doubles of length of at least n containint unit trans.cost
   - narcs  : number of selected suppliers/arcs
   - sel_arc: ordered set of arcs/suppliers such that arcs sel_arc[0 ... narcs-1]
              are the arcs where flow can be positive
   - flow   : array of integers of length of at least n that contains on output
              the computed flows on the arcs

   Return value: objective function value of the solution
   
   Remark   : set of selected arcs must have enough capacity to meet the demand!
*/
  int    excess, j, jj, jmax, arc[n];
  double cost, tcmax, fcmax;
  
  for ( j=0; j < n; j++ ) flow[j]=0;
  for ( j=0, excess=-demand; j < narcs; j++ ) 
    jj = sel_arc[j], arc[j] = jj, flow[jj] = cap[jj], excess += cap[jj];
  if ( excess < 0 ) return ( INFTY );
  while ( excess > 0 ) {
    for ( j=1,jmax=0,jj=arc[0],tcmax=tcost[jj],fcmax=fcost[jj]; j < narcs; j++ ) {
      jj = arc[j];
      if ( (tcost[jj] > tcmax) || ( (!(tcost[jj] < tcmax)) && (fcost[jj] > fcmax) ) )
        tcmax = tcost[jj], fcmax = fcost[jj], jmax=j;
    }  
    jj = arc[jmax];
    if ( excess < cap[jj] ) {
      flow[jj]-= excess;
      excess   = 0; 
    } else {
      arc[jmax] = arc[--narcs];
      arc[narcs]= jj; 
      flow[jj] -= cap[jj];
      excess   -= cap[jj];
    }  
  }
  for ( j=0, cost=0.0; j < narcs; j++ ) 
    jj=arc[j], cost += fcost[jj]+tcost[jj]*flow[jj];
  return( cost );
}		 
   
/*-----------------------------------------------------------------------------*/

static
int ssfctp_lp( int n, int demand, int* cap, double* fcost, double* tcost, 
               double* lpcst, int* order, int* critind, double* lpbnd, int* xlp,
               int* x, double* upbnd ) {
/* Description:
   ------------
   Solves the LP relaxation of the SSFCTP. The LP relaxation is to solve
   
       min \sum_j ( c_j + f_j/k_j ) x_j
     s.t.:
           \sum_j x_j = D
	   0 <= x_j <= k_j for j=1,...,n 
	   
   Assume that c_1+f_1/k_1 <= ... <= c_n+f_n/k_n and let q be such that
   \sum_j k_j < D and \sum_j k_j >= D. Then x_j=k_j for j=1,...,q-1,
   x_q=D-\sum_{j=1}^{q-1} k_j, x_j=0 for j=q+1,...,n is an optimal
   solution to this LP. Let \sigma and \eta_j denote the dual multipliers
   of the above constraints. Then \sigma=c_q + f_q/k_q and 
   \eta_j=\max\{ 0, \sigma - c_j - f_j/k_j \} is an optimal dual solution.
   
   Scope: intern
   -------------
   
   Parameters:
   -----------
   n       : number of suppliers
   demand  : the sink's demand
   cap     : array of integers of length of at least n, where cap[j] is
             the capacity of the arc linking node j and the sink
   fcost   : array of doubles of length of at least n containing the
             fixed costs f_j of the arcs
   tcost   : array of doubles of length of at least n containing the
             unit transportation costs c_j on the arcs
   lpcst   : array of doubles of length of at least n that on output
             contains the linearized costs tcost[j] + fcost[j]/cap[j]
	            of each arc
   order   : an array of integers of length of at least n. On output
             the array orders the suppliers according to non-decreasing
	            linearized cost lpcst[j]
   critind : index of the "critical supplier" q, i.e. q=order[critind]
   lpbnd   : objective function value of the LP solution
   x       : array of integers of length of at least n that on output
             contains the transportation quantities in a feasible
             solution constructed from the LP solution (in case that the
	            LP is integer-valued this solution is an optimal one).
	            x may be NULL. In that case the solution is not returned.
   upbnd   : pointer to a double that on output contains the objective
             function value of the solution x. 
	     
  Return value:
  ------------
  1 if the LP solution is integer or lower bound equals upper bound;
  0 otherwise.
   
*/		
  double lobnd=0, upbnd2, minval, value;
  int    i, j, narcs, next, q, supply, sel_arc[n], x2[n];
  
  for ( j=0; j < n; j++ ) {
    order[j] = j;
    xlp[j]   = 0;
    if ( cap[j] > 0 ) 
      lpcst[j] = tcost[j] + fcost[j]/(double)cap[j];
    else
      lpcst[j] = INFTY;
  }
  SSFCTP_ptr = lpcst;
  SSFCTP_sgn = 1;
  qsort( order, n, sizeof(int), (void*) SSFCTP_compare );
  for ( q=supply=0, lobnd=0.0; supply < demand; q++ ) {
    i = order[q];
    xlp[i]  = cap[i];
    supply += cap[i];
    lobnd  += lpcst[i]*cap[i];
  }  
  if ( supply == demand ) { /* LP solution is integer */
    *upbnd = lobnd;
    for ( j=0; j < n; j++ ) x[j]=xlp[j];
    return( 1 );
  }
  xlp[i]-= supply-demand;
  lobnd -= (supply-demand)*lpcst[i];
  *lpbnd = lobnd;
  *critind = q-1;

  /* Determine a feasible solution by selecting those arcs that have positive 
     flow in the LP solution */
  *upbnd = ssfctp_flow( n, demand, cap, fcost, tcost, q, order, x );
  if  LBgUB( lobnd, *upbnd ) return( 1 );

  /* Determine a second feasible solution in the following way: first select the 
     arcs that have flow equal to the capacity in the LP solution. Let d* be the 
     remaining demand. As long as d* is greater than zero select an additional 
     arc j with minimum value of tcost[j]+fcost[j]/min(d*,cap[j]) */
  for ( j=0; j < n; j++ ) sel_arc[j] = order[j];
  narcs = *critind;
  for ( j=supply=0; j < narcs; j++ ) supply += cap[order[j]];
  while ( supply < demand ) {
    for ( j=narcs, minval=INFTY; j < n; j++ ) {
      i = sel_arc[j]; 
      value = tcost[i] + fcost[i]/MIN( cap[i], demand-supply );
      if ( value < minval ) minval = value, next = j;
    }
    i = sel_arc[next];
    supply += cap[i];
    sel_arc[next] = sel_arc[narcs];
    sel_arc[narcs++] = i;
  }
  /* Determine the flow on the arcs selected above */
  upbnd2 = ssfctp_flow( n, demand, cap, fcost, tcost, narcs, sel_arc, x2 );
  if ( upbnd2 < *upbnd ) {
    *upbnd = upbnd2;
    for ( j=0; j < n; j++ ) x[j] = x2[j];
    if  LBgUB( lobnd, *upbnd ) return( 1 );
  }

  return ( 0 );   
}

/*-----------------------------------------------------------------------------*/

void ssfctp_gr( int n, int demand, int* cap, double* fcost, double* tcost, 
 	        int adapt, int* x, double* upbnd ) {
/* Description:
   ------------
   Computes a feasible solution for the SSFCTP using either the greedy heuristic
   (adapt=0) or the adaptive greedy approach (adapt=1)      
*/		
  int    i, j, next, order[n], s, supply;
  double minval, value, relcst[n];

  for ( j=0; j < n; j++ ) order[j]=j, relcst[j]=tcost[j]+fcost[j]/(double)cap[j];
  SSFCTP_ptr = relcst;
  SSFCTP_sgn = 1;
  qsort( order, n, sizeof(int), (void*) SSFCTP_compare );            
  
  for ( j=0; j < n; j++ ) x[j] = 0;
  for ( s=supply=0, *upbnd=0.0; supply < demand; s++ ) {
    i       = order[s];
    x[i]    = cap[i];
    supply += cap[i];
    *upbnd += fcost[i]+cap[i]*tcost[i];
  }  
  if ( supply > demand ) {
    if ( adapt == 0 ) {
       x[i] -= supply-demand;
      *upbnd = ssfctp_flow( n, demand, cap, fcost, tcost, s, order, x );
    } 
    else { /* Apply adaptive greedy method */ 
      x[i]    = 0;
      supply -= cap[i];
      s--;
      while ( supply < demand ) {
        for ( j=s, minval=INFTY; j < n; j++ ) {
          i = order[j];
          value = tcost[i] + fcost[i]/MIN( cap[i], demand-supply );
          if ( value < minval ) minval = value, next = j;
        }
        supply += cap[i=order[next]];
        order[next] = order[s];
        order[s++] = i;
      }
      *upbnd = ssfctp_flow( n, demand, cap, fcost, tcost, s, order, x );
    }
  }
}

/*-----------------------------------------------------------------------------*/

static
void ssfctp_glo ( int n_dim, int n, int demand, int skip, int* cap, double* fcost, 
                  double* tcost, int* order, int* n_bigs, int* bigs, int* x, 
                  double* upbnd ) {
/*
  Heuristic of Gens and Levner for the min-knapsack adjust to the case of
  the SSFCTP. It is assumed that the array order already gives an ordering
  of the suppliers according to increasing relative costs. The original problem
  size is assumed to equal n_dim suppliers; the suppliers that may provide
  positive supply is given by the set order[0 .. n-1] where n can be smaller than
  n_dim. The procedure returns the big suppliers' positions in the ordering 
  order in the array bigs, if on input bigs is not NULL.
*/		
  int    cur_x[n_dim], i, ii, j, jj, li, n_small, small[n_dim], supply, s_supply;
  double cur_obj, s_cst; 
	
  /* Construct the trial solutions by completing the full supply of conscecutive
     small suppliers with a part of the supply of a big supplier */		
  if ( bigs ) *n_bigs=0;   
  for ( j=0; j < n; j++ ) cur_x[order[j]] = 0;     
  for ( j=0, cur_obj=0.0, *upbnd=INFTY, supply=n_small=0; j < n; j++ ) if ( (jj=order[j])!=skip ){     
    supply += cap[jj=order[j]];
    if ( supply < demand ) {  /* this is a small supplier */
      cur_x[jj] = cap[jj];
      cur_obj  += fcost[jj] + tcost[jj]*cur_x[jj];  
      small[n_small++] = jj;
    } 
    else { /* Big supplier found: complete the trial solution */
      /* Check if the solution can be improved by removing some small suppliers */
      if ( bigs ) bigs[ (*n_bigs)++] = j;
      s_supply  = supply-cap[jj];
      s_cst     = cur_obj;
      cur_x[jj] = cap[jj];
      cur_obj  += fcost[jj] + tcost[jj]*cap[jj];
      /* Improve candidate solution by removing small suppliers no more required */
      for ( i=n_small-1, li=-1; i >= 0; i-- ) {
        ii = small[i];
        if ( (supply - cap[ii]) < demand ) break;
          supply   -= cur_x[ii];
	         cur_obj  -= fcost[ii] + tcost[ii]*cur_x[ii];
	         cur_x[ii] = 0;
	         li        = i;  /* record the small supplies reset to zero */
        } 
      cur_x[jj]-= (supply - demand); 
      cur_obj  -= ( cap[jj] - cur_x[jj] )*tcost[jj];
      if ( cur_obj < *upbnd ) { /* improved solution found */
        *upbnd = cur_obj;
	       for ( i=0; i < n; i++ ) ii=order[i], x[ii] = cur_x[ii];	 
      }
      /* Remove the big supplier from the current solution and try the next one */
      cur_obj  = s_cst;
      supply   = s_supply;
      cur_x[jj]= 0;
      /* Restore the supplies of small suppliers that were reset to zero */
      if ( li >= 0 )
        for ( i=li; i < n_small; i++ ) ii=small[i], cur_x[ii] = cap[ii];
    }
  }
      
}		

/*-----------------------------------------------------------------------------*/

static
Tnode *ssfctp_red( int n, int* demand, int* cap, double* fcost, double* tcost, 
                   double* lpbnd, double* lpcst, int* order, int* q, int* xlp,
                   double* upbnd, int* x, int* nn, double* zfix, int* maxflow,
	                  int* mincflow ) {
/* Description:
   ------------
   Applies a reduced cost based reduction test to the SSFCTP instance. If
   possible, the procedure excludes arcs/suppliers from consideration (capacity
   of these suppliers is then zero) and tries to fix the supply of other 
   suppliers to the upper bound. Based on reduced cost information, the
   procedure might detect a subset Jpos of suppliers that have to supply a
   positive amount in a solution that improves the upper bound. Since there
   is an optimal solution where at most one selected supplier supplies less
   than it's capacity, all flow variables x[j], j \in Jpos, except one 
   can be fixed to the upper bound. The flow variable j_w with largest
   unit transportation cost from this set cannot be fixed to the upper bound,
   but a lower bound l on the flow x[j_w] might be established. The problem
   is then reduced by eliminating all fixed variables and replacing the 
   variable x[j_w] be x_[j_w]-l. The data, that is demand, solution values
   and solution are the adjusted accordingly.
   
   Parameters:
   -----------
   - n      : number of suppliers
   - demand : pointer to an integer containing the sink's demand. In case
              that some flow variables can be fixed to its upper bound, the
	             demand is reduced accordingly by this amount of flow. Furthermore,
	             there can be one flow variable that could not be fixed to the
	             upper bound, but which must be less than or equal to a certain
	             value l. The demand is also reduced by this value and the
	             corresponding flow variable x substituted by x-l.
   - cap    : array of integers of length of at least n containing arc 
              capacities
   - fcost  : array of doubles of length of at least n containing fixed costs
              (fixed cost of a supplier j_w are reset to zero)
   - tcost  : array of doubles of length of at least n containing unit costs
   - lpbnd  : pointer to a double that containts the objective function value 
              of the LP solution. The value might be reduced on output
	             in such a way that it matches to the reduced problem.
   - lpcst  : array of doubles of length of at least n containing the linearized
              cost tcost[j] + fcost[j]/cap[j] for j=0,...,n-1
	             The linearized cost of a supplier j_w are reset to tcost[j_w]
   - order  : array of integers of length of at least n such that on input
              lpcst[order[j]] <= lpcst[order[j+1]] for j=0,...,n-2
	             and
	             lpcst[order[j]] <= lpcst[order[j+1]] for j=0,...,nn-2 on output.
	             Furthermore, on output order[nn ... n-1] are the index numbers of 
	             arcs where flow is fixed to zero and the upper bound, resp.
   - q      : pointer to an integer that contains on input and on output the 
              index of the "critical supplier".
   - xlp    : array of doubles containing LP solution.
   - upbnd  : pointer to a double containing the objective function value of 
              a known feasible solution. If the problem is reduced, this value
	             is adjusted accordingly on output.
   - x      : feasible solution corresponding to upbnd.
   - nn     : number of remaining variables whose value could not be fixed.
   - zfix   : pointer to a double that on output contains the part of the
              objective function value that is already determined.   
   - maxflow: array of integers of length of at least n. On output maxflow[j]
              is the maximum amount of flow on this arc. For arcs 
	             j from order[nn ... n-1] maxflow[j] equals zero if the flow
	             variable was fixed to zero or equal to cap[j] if the flow
	             variable could be fixed to its upper bound.
	             bound
   - mincflow:array of integers of length of at least n. On output mincflow[j]
              is the minimum required amount of flow on this arc provided
	             that this arc is selected (y[j]=1)
	      
  Return value:
  -------------
  In case that a set Jpos as described above is detected, the return value
  points to the data of the above mentioned supplier j_w. If no such set
  Jpos could be found, NULL is returned.  	      
*/		 

  int    forced[n], ord[n], i, j, jjw, jw, jw_min, jw_max, k, lflow, num, 
         nforced, qq=*q, rdem;
  double eta, cstq, cmax, gap=*upbnd-*lpbnd, pen, tmp;
  char   isforced[n];
  Tnode  *wnode=NULL;

  *zfix = 0;  
  cstq  = lpcst[order[qq]];
  for ( j=nforced=0, cmax=-1.0; j <= qq; j++ ) {
    i           = order[j];
    ord[j]      = i;
    eta         = cstq-lpcst[i];
    lflow       = 0;
    mincflow[i] = 0;
    maxflow[i]  = cap[i];
    isforced[i] = 0;
    if ( eta > ZERO ) {
      tmp = (double)cap[i] - gap/eta;
      if ( tmp > ZERO ) { 
        lflow = (int) (tmp+ONE); /* minimum amount of flow on this arc */     
        if ( lflow > 0 ) {
          if (  (tcost[i] > cmax) || ( (!( tcost[i] < cmax)) 
              && ( lflow + jw_max > jw_min + maxflow[i] ) ) )
            jjw=j, jw=nforced, cmax=tcost[i], jw_min=lflow, jw_max=maxflow[i];
          isforced[i] = 1;
          forced[nforced++] = i;
        }     
      }	
    }
  }
  
  for ( j=qq+1; j < n; j++ ) {
    i = order[j];
    ord[j] = i;
    eta    = lpcst[i] - cstq;
    mincflow[i]= 0;
    maxflow[i] = cap[i];
    isforced[i]= 0;
    if ( tcost[i] < cstq ) 
      pen = *lpbnd + fcost[i] + ( tcost[i] - cstq )*cap[i];    
    else
      pen = *lpbnd + fcost[i] + tcost[i] - cstq;  
    if ( pen + ZERO > *upbnd ) 
      maxflow[i] = 0;
    else { 
      if ( eta > ZERO ) {
        tmp = gap/eta;
        maxflow[i]=MIN(cap[i], (int) (tmp+ZERO) );   
      }	
      if ( tcost[i] > cstq + ZERO ) {
        tmp = (gap - fcost[i])/(tcost[i]-cstq);
	       maxflow[i] = MIN( maxflow[i], (int) (tmp + ZERO));
      }
      else if ( cstq > tcost[i] + ZERO ) {
        tmp = (gap - fcost[i])/(tcost[i]-cstq);
	       mincflow[i] = MAX( 0, (int) (tmp+ONE) );
      }
    }  
  }
  
  if ( nforced > 0 ) {
    i = forced[jw]; 
    if ( jw_min < maxflow[i] ) {
      forced[jw]      = forced[--nforced]; 
      forced[nforced] = i;
      wnode = (Tnode*) calloc(1, sizeof(Tnode) );
      wnode->idx = i, wnode->lflow = jw_min, wnode->uflow = maxflow[i];
      wnode->fc  = fcost[i], wnode->tc = tcost[i];
      maxflow[i]-= jw_min;
      *zfix      = jw_min*tcost[i] + fcost[i];
      *demand   -= jw_min;
      x[i]      -= jw_min;
      xlp[i]     = 0;
      fcost[i]   = 0.0;
      lpcst[i]   = tcost[i];    
      isforced[i]= 0;
      order[0]   = wnode->idx;
    } else 
      jjw = -1;
    for ( j=0; j < nforced; j++ ) {
      i = forced[j];
      *zfix   += fcost[i] + tcost[i]*maxflow[i];
      *demand-= maxflow[i];
      x[i]   -= maxflow[i];
    }  
    *upbnd -= *zfix;
    num     = ( wnode ) ? 0 : -1;
    for ( j=0, k=n; j < n; j++ ) if ( j != jjw ) {
      xlp[i=ord[j]] = 0;
      if ( ( maxflow[i] == 0 ) || ( isforced[i] ) )
        order[--k] = i;
      else if ( (wnode) && (lpcst[i] < wnode->tc) )
        order[num]=i, order[++num] = wnode->idx;	
      else order[++num] = i;
    }
    for( num++, qq=0, rdem=*demand, *lpbnd=0.0; qq < num; qq++ ) {
      i = order[qq];
      xlp[i] = MIN( rdem, maxflow[i] );
      rdem  -= xlp[i];
      *lpbnd+= xlp[i]*lpcst[i];
      if ( rdem==0 ) break;
    }
    *q=qq;
  } 
  else for ( j=num=0, k=n; j < n; j++ ) {
    i = ord[j];
    if ( maxflow[i] > 0 ) order[num++] = i; else order[--k] = i;
  }
  *nn = num;
  return ( wnode );
		 
}		 
 
/*-----------------------------------------------------------------------------*/

void ssfctp_gl ( int n, int demand, int* cap, double* fcost, double* tcost, 
                 int preproc, int* x, double* upbnd ) {
/*
  First orders the suppliers according to increasing relative costs and then
  calls the heuristic of Gens and Levner. If preproc > 0, then also the
  complete reduction test is carried out before calling the Gens and Levner
  type heuristic.
*/		
  int    order[n], maxflow[n], mincflow[n], xlp[n], first_x[n];
  int    i, j, nn=n, q, *mflow=cap;
  double *fc, *tc, lpcst[n];  
  double first_ub=INFTY, lpbnd, zcnst=0.0, zfix=0.0;
  Tnode  *jw=NULL;	 
  
  /* Transform data such that fixed-charges and transp. cost are nonnegative */
  zcnst = ssfctp_dta ( n, demand, tcost, fcost, &tc, &fc ); 

  /* Check if problem is trivially solvable by supplying the sink's demand
     from the node that minimizes the cost */
  if ( ssfctp_trivial( n, demand, cap, fc, tc, x, upbnd ) ) goto RETURN;
                   
  /* If preproc > 0 then first solve the LP relaxation and then apply the
     reduction test */
  *upbnd = INFTY;
  if ( preproc > 0 ) {
    mflow = maxflow;
    if ( ssfctp_lp( n, demand, cap, fc, tc, lpcst, order, &q, &lpbnd, xlp,
                    x, upbnd ) ) goto RETURN;
    /* The reduction test ask if a penalty is larger or equal than the
       current upper bound. If the upper bound is optimal it might thus
       in some cases happen that, e.g., the maximum flow on arc i is
       reduced to zero also x[i]>0 must hold in an optimal solution.
       Hence, save the solution found so far */
    first_ub = *upbnd;
    for ( j=0; j < n; j++ ) first_x[j]=x[j];
    jw = ssfctp_red( n, &demand, cap, fc, tc, &lpbnd, lpcst, order, &q, xlp,
                     upbnd, x, &nn, &zfix, mflow, mincflow );        
    /* Some arc capacities may have been reduced. Recompute therefore
       the relative costs and reorder */
    for ( j=0; j < nn; j++ ) i=order[j], lpcst[i]=tc[i]+fc[i]/mflow[i];
    qsort( order, nn, sizeof(int), (void*) SSFCTP_compare ); 
  } 
  else { 
    for ( j=0, nn=n; j < nn; j++ ) { 
      order[j] = j;
      lpcst[j] = tcost[j] + fcost[j]/mflow[j];
    }  
    SSFCTP_ptr = lpcst;
    SSFCTP_sgn = 1;
    qsort( order, nn, sizeof(int), (void*) SSFCTP_compare );   
  }  
                                    
  /* Perform then Gens & Levner type heuristic */
  if ( nn > 0 )
    ssfctp_glo( n, nn, demand, n, mflow, fc, tc, order, NULL, NULL, x, upbnd );
   
RETURN:
  *upbnd += zfix;
  for ( j=nn; j < n; j++ ) i=order[j], x[i]=mflow[i];
  if ( jw ) {
    fc[jw->idx] = jw->fc;
    x[jw->idx] += jw->lflow;
    free ( jw );
  }
  if ( *upbnd > first_ub - ZERO ) {
    *upbnd = first_ub;
    for ( j=0; j < n; j++ ) x[j]=first_x[j];
  }
  if ( zcnst < 0.0 ) {
    *upbnd += zcnst;
    if ( tc != tcost ) free( tc );
    if ( fc != fcost ) free( fc );
  }
   
}		

/*-----------------------------------------------------------------------------*/

void ssfctp_gli( int n, int demand, int* cap, double* fcost, double* tcost, 
                 int preproc, int* x, double* upbnd ) {
/* 
  Improvement over the heuristic of Gens and Levner as proposed by Csirik et al.
  for the min-knapsack problem.
*/
  int    big[n], maxflow[n], mincflow[n], ord[n], order[n], xlp[n], first_x[n];
  int    i, ii, j, jb, jo, n_bigs=0, nn, q, *cur_x=xlp, *mflow=cap;
  double *fc, *tc, lpcst[n];  
  double cur_obj, fcstj, first_ub=INFTY, lpbnd, zcnst=0.0, zfix=0.0;
  Tnode  *jw=NULL;	 
  
  /* Transform data such that fixed-charges and transp. cost are nonnegative */
  zcnst = ssfctp_dta ( n, demand, tcost, fcost, &tc, &fc ); 

  /* Check if problem is trivially solvable by supplying the sink's demand
     from the node that minimizes the cost */
  if ( ssfctp_trivial( n, demand, cap, fc, tc, x, upbnd ) ) goto RETURN;
                   
  /* If preproc > 0 then first solve the LP relaxation and then apply the
     reduction test */
  *upbnd = INFTY;
  if ( preproc > 0 ) {
    mflow = maxflow;
    if ( ssfctp_lp( n, demand, cap, fc, tc, lpcst, order, &q, &lpbnd, xlp,
                    x, upbnd ) ) goto RETURN;
    /* Save first solution */
    first_ub = *upbnd;
    for ( j=0; j < n; j++ ) first_x[j]=x[j];    
    jw = ssfctp_red( n, &demand, cap, fc, tc, &lpbnd, lpcst, order, &q, xlp,
                     upbnd, x, &nn, &zfix, mflow, mincflow );     
    /* Some arc capacities may have been reduced. Recompute therefore
       the relative costs and reorder */
    for ( j=0; j < nn; j++ ) i=order[j], lpcst[i]=tc[i]+fc[i]/mflow[i];
    qsort( order, nn, sizeof(int), (void*) SSFCTP_compare ); 
  } 
  else { 
    for ( j=0, nn=n; j < nn; j++ ) { 
      order[j] = j;
      lpcst[j] = tcost[j] + fcost[j]/mflow[j];
    }  
    SSFCTP_ptr = lpcst;
    SSFCTP_sgn = 1;
    qsort( order, nn, sizeof(int), (void*) SSFCTP_compare );   
  }  
                    
  /* Apply Gens/Levner heuristic */
  if ( nn > 0 )
    ssfctp_glo( n, nn, demand, n, mflow, fcost, tcost, order, &n_bigs, big, x, upbnd );
  
  /* For every big supplier j_b do the following in turn: Reset its fixed cost
     to zero and reapply the Gens/Levner heuristic to the reduced problem */
  for ( j=0; j < n_bigs; j++ ) {
    jo    = big[j];     /* position of the j-th big supplier in ordering "order" */
    jb    = order[jo];  /* index of this supplier */
    fcstj = fcost[jb], fcost[jb] = 0.0;
    for ( i=jo+1; i < nn; i++ ) ord[i]=order[i];
    for ( i=0; i < jo; i++ )
      if ( lpcst[order[i]] < tcost[jb] ) ord[i] = order[i]; else break;
    ord[i] = jb;  
    for ( ii=i+1; ii <= jo; ii++ ) ord[ii] = order[ii-1];
    ssfctp_glo( n, nn, demand, n, mflow, fcost, tcost, ord, NULL, NULL, cur_x, &cur_obj );    
    fcost[jb] = fcstj;
    if ( cur_x[jb] > 0 ) cur_obj += fcstj;
    if ( cur_obj < *upbnd ) 
      for ( i=0, *upbnd=cur_obj; i < nn; i++ ) ii=order[i], x[ii] = cur_x[ii];
  }

RETURN:
  *upbnd += zfix;
  for ( j=nn; j < n; j++ ) i=order[j], x[i]=mflow[i];
  if ( jw ) {
    fc[jw->idx] = jw->fc;
    x[jw->idx] += jw->lflow;
    free ( jw );
  }     
  if ( *upbnd > first_ub - ZERO ) {
    *upbnd = first_ub;
    for ( j=0; j < n; j++ ) x[j]=first_x[j];
  }
  if ( zcnst < 0.0 ) {
    *upbnd += zcnst;
    if ( tc != tcost ) free( tc );
    if ( fc != fcost ) free( fc );
  }

}		 

/*-----------------------------------------------------------------------------*/

void ssfctp_gli2( int n, int demand, int* cap, double* fcost, double* tcost, 
                  int preproc, int* x, double* upbnd ) {
/* 
  Second way of improving the Gens-Levner type heuristic for the SSFCTP
*/
  int    big[n], order[n], maxflow[n], mincflow[n], xlp[n], first_x[n];
  int    i, ii, j, n_bigs=0, nn, q, *cur_x=xlp, *mflow=cap;
  double *fc, *tc, lpcst[n];  
  double cur_obj, first_ub=INFTY, lpbnd, zcnst=0.0, zfix=0.0;
  Tnode  *jw=NULL;	 
  
  /* Transform data such that fixed-charges and transp. cost are nonnegative */
  zcnst = ssfctp_dta ( n, demand, tcost, fcost, &tc, &fc ); 

  /* Check if problem is trivially solvable by supplying the sink's demand
     from the node that minimizes the cost */
  if ( ssfctp_trivial( n, demand, cap, fc, tc, x, upbnd ) ) goto RETURN;
                   
  /* If preproc > 0 then first solve the LP relaxation and then apply the
     reduction test */
  *upbnd = INFTY;
  if ( preproc > 0 ) {
    mflow = maxflow;
    if ( ssfctp_lp( n, demand, cap, fc, tc, lpcst, order, &q, &lpbnd, xlp,
                    x, upbnd ) ) goto RETURN;
    /* Save first solution */
    first_ub = *upbnd;
    for ( j=0; j < n; j++ ) first_x[j]=x[j];
    jw = ssfctp_red( n, &demand, cap, fc, tc, &lpbnd, lpcst, order, &q, xlp,
                     upbnd, x, &nn, &zfix, mflow, mincflow );     
    /* Some arc capacities may have been reduced. Recompute therefore
       the relative costs and reorder */
    for ( j=0; j < nn; j++ ) i=order[j], lpcst[i]=tc[i]+fc[i]/mflow[i];
    qsort( order, nn, sizeof(int), (void*) SSFCTP_compare ); 
  } 
  else { 
    for ( j=0, nn=n; j < nn; j++ ) { 
      order[j] = j;
      lpcst[j] = tcost[j] + fcost[j]/mflow[j];
    }  
    SSFCTP_ptr = lpcst;
    SSFCTP_sgn = 1;
    qsort( order, nn, sizeof(int), (void*) SSFCTP_compare );   
  }  

  /* Apply Gens/Levner heuristic */
  if ( nn > 0 )
    ssfctp_glo( n, nn, demand, n, mflow, fcost, tcost, order, &n_bigs, big, x, upbnd );
  
  /* For every big supplier j_b do the following in turn: Set his supply
     to the upper bound and solve the remaining reduced problem by
     the Gens-Levner type heuristic */
  for ( j=0; j < n_bigs; j++ ) {
    i = order[big[j]];  /* index of current big supplier */
    ssfctp_glo( n, nn, demand-mflow[i], i, mflow, fcost, tcost, order, 
                NULL, NULL, cur_x, &cur_obj );    
    cur_x[i] = mflow[i];
    cur_obj += fcost[i] + tcost[i]*mflow[i];
    if ( cur_obj < *upbnd ) 
      for ( i=0, *upbnd=cur_obj; i < nn; i++ ) ii=order[i], x[ii] = cur_x[ii];
  }

RETURN:
  *upbnd += zfix;
  for ( j=nn; j < n; j++ ) i=order[j], x[i]=mflow[i];
  if ( jw ) {
    fc[jw->idx] = jw->fc;
    x[jw->idx] += jw->lflow;
    free ( jw );
  }     
  if ( *upbnd > first_ub - ZERO ) {
    *upbnd = first_ub;
    for ( j=0; j < n; j++ ) x[j]=first_x[j];
  }
  if ( zcnst < 0.0 ) {
    *upbnd += zcnst;
    if ( tc != tcost ) free( tc );
    if ( fc != fcost ) free( fc );
  }

}		 

/*-----------------------------------------------------------------------------*/

static
void ssfctp_bol( int n, int demand, int q, double gap, int* maxflow, 
                 double* fcost, double* tcost, int* order, int* Lmin, int* Lmax ) {
/* Description:
   ------------
   For a given ordering of the nodes and for every node i, bounds on the amount
   L that has to be supplied by nodes j=1,...,i-1 in feasible/optimal solutions
   are computed. Beside demand and capacities, reduced cost values are employed
   for computing these bounds.
   
   Parameters:
   -----------
   - n      : number of suppliers
   - demand : demand of the sink
   - q      : index of the "critical supplier" as returned by ssfctp_lp
   - gap    : duality gap (difference between upper bound and LP solution)
   - maxflow: array of integers of length of at least n where maxflow[j]
              is maximum possible flow on arc j
   - fcost  : array of doubles of length of at least n containing fixed costs
   - tcost  : array of doubles of length of at least n containing unit costs
   - order  : array of integers of length of at least n. 
              On input: cost[order[j]] <= cost[order[j+1]] for j=0,...,n-2,
	             where cost[i]=tcost[i]+fcost[i]/cap[i]. 
   - Lmin   : array of integers of length of at least n. Lmin[i] is a lower
              bound on L
   - Lmax   : array of integers of length of at least n. Lmin[i] is an upper
              bound on L   	      
*/		 
  int    capsum_i, i_capsum, i, j, minL, prev;
  double cst, cstq, eta, tmp;
  
  for ( j=capsum_i=0; j < n; j++ ) capsum_i += maxflow[order[j]];
  i    = order[q];
  cstq = tcost[i] + fcost[i]/maxflow[i];
  i    = order[0];
  cst  = tcost[i] + fcost[i]/maxflow[i];
  Lmin[i]  = 0;
  Lmax[i]  = 0;
  prev     = 0;
  i_capsum = maxflow[i];
  capsum_i-= maxflow[i];
  
  for ( j=1; j <= q; j++ ) {
    i = order[j];
    Lmax[i] = MIN( demand, i_capsum );
    Lmin[i] = MAX( demand - capsum_i, 0 );
    eta     = cstq - cst;
    if ( Lmin[i] < prev ) Lmin[i] = prev;
    if ( eta > ZERO ) {
      tmp  = (double)i_capsum - gap/eta;
      minL = ( int ) ( tmp + ONE );
      if ( minL  > Lmin[i] ) Lmin[i] = minL;
    }  
    i_capsum += maxflow[i];
    capsum_i -= maxflow[i];
    cst  = tcost[i] + fcost[i]/maxflow[i];
    prev = Lmin[i];
  }
  
  for ( j=q+1; j < n; j++ ) {
    i   = order[j];
    cst = tcost[i] + fcost[i]/maxflow[i];
    Lmax[i] = MIN(demand, i_capsum );
    Lmin[i] = MAX( demand - capsum_i, 0 );
    eta     = cst - cstq;
    if ( Lmin[i] < prev ) Lmin[i] = prev;
    if ( eta > ZERO ) {
      tmp = (double)demand - gap/eta;
      minL= ( int ) ( tmp + ONE );
      if ( minL > Lmin[i] ) Lmin[i] = minL;
    }
    i_capsum += maxflow[i];
    capsum_i -= maxflow[i];
    prev = Lmin[i];
  }
  /* We must have that Lmin[i] >= Lmin[i+1]-k_i */
  for ( j=n-2; j >= 0; j-- ) {
    i    = order[j];
    minL = Lmin[order[j+1]] - maxflow[i];
    if ( Lmin[i] < minL ) Lmin[i] = minL;
  }
  
}  

/*-----------------------------------------------------------------------------*/

static
void ssfctp_nbol( int n, int demand, int* maxflow, int* order, int* Lmin, 
                  int* Lmax ) {
/* Description:
   ------------
   For a given ordering of the nodes and for every node i, bounds on the amount
   L that has to be supplied by nodes j=1,...,i-1 in feasible/optimal solutions
   are computed. In contrast to procedure ssfctp_bol, these bounds do now only
   depend on the arc capacities and the sink's demand.
   
   Parameters:
   -----------
   - n      : number of suppliers
   - demand : demand of the sink
   - minflow: array of integers of length of at least n where minflow[j]
              is minimum flow required on arc j
   - maxflow: array of integers of length of at least n where maxflow[j]
              is maximum possible flow on arc j
   - order  : array of integers of length of at least n. 
              On input: cost[order[j]] <= cost[order[j+1]] for j=0,...,n-2,
	      where cost[i]=tcost[i]+fcost[i]/cap[i]. 
   - Lmin   : array of integers of length of at least n. Lmin[i] is a lower
              bound on L
   - Lmax   : array of integers of length of at least n. Lmin[i] is an upper
              bound on L   	      
*/		 
  int    capsum_i, i_capsum, i, j;
  
  for ( j=capsum_i=0; j < n; j++ ) capsum_i += maxflow[order[j]];  
  for ( i=order[0], i_capsum=Lmin[i]=Lmax[i]=0, j=1; j < n ; j++ ) {
    i_capsum += maxflow[i];
    capsum_i -= maxflow[i];
    i = order[j];
    Lmax[i] = MIN( demand, i_capsum );
    Lmin[i] = MAX( demand - capsum_i, 0 );
  }
  
}  

/*-----------------------------------------------------------------------------*/

static int ssfctp_gcd( int n, int demand, int* cap ) {
/* Compute greatest common divisor of demand and capacitites */

  int c, gcd, j, tmp;
 
  for ( j=0, gcd=demand; (j < n) && (gcd > 1); j++ ) {
    c  = cap[j];
    if ( gcd > c ) tmp=gcd, gcd=c, c=tmp;
    do {
      c = c % gcd;
      if ( c > 0 ) tmp=gcd, gcd=c, c=tmp;
    } while ( c > 0 ); 
  }
  
  return( gcd );
  
}

/*-----------------------------------------------------------------------------*/

static 
void addtoHList( HList** firstH, HList** lastH, int r, double hval ) {
/* Insert the element (r,hval) into the sorted list of function values
   (r_q, H_i(r_q) ) while removing all elements (r',H_i(r') ) with
   H_i(r') >= hval from the list (r' < r has to be ensured)
*/   
  HList *hptr=*firstH, *nexth;

  while ( hptr != NULL ) {
    nexth = hptr->succ;
    if ( hptr->H >= hval ) { /* delete the element to which hptr points to */
      if ( hptr->pred ) 
        (hptr->pred)->succ = hptr->succ;
      else
        *firstH = hptr->succ;	
      if ( hptr->succ ) 
        (hptr->succ)->pred = hptr->pred;
      else
        *lastH = hptr->pred;	
      free( hptr );
    }
    hptr = nexth;
  }
  hptr       = calloc( 1, sizeof( HList ) );
  hptr->r    = r;
  hptr->H    = hval;
  hptr->succ = NULL;
  hptr->pred = *lastH;
  if ( *lastH ) (*lastH)->succ = hptr;
  *lastH     = hptr;
  if ( *firstH == NULL ) *firstH=hptr;      
}

/*-----------------------------------------------------------------------------*/

static
void rearrangeHList( HList** firstH ) {
/* Remove the first element of the list {(r,H_i(r)}. Then rearrange the list
   such that the first element points to (r*, H*), where H* = min_r H_i(r) 
   and r*= max { r : H_i(r)=H* }. Furthermore, remove all elements 
   (r,H_i(r)) with r < r* and H_i(r) >= H*
*/
  HList *hptr, *mptr=NULL, *nexth, *fptr=*firstH;
   
  for ( hptr=fptr->succ, mptr=hptr; hptr != NULL; hptr=hptr->succ ) {
    if ( hptr->H < mptr->H ) mptr=hptr;
    else if ( ( !(hptr->H > mptr->H) ) && (hptr->r > mptr->r) ) mptr=hptr;
  }
  fptr->r=mptr->r;
  fptr->H=mptr->H;
  mptr->r = -1;
  hptr = fptr->succ;
  while ( hptr != NULL ) {
    nexth = hptr->succ;
    if ( hptr->r < fptr->r ) {
      if ( hptr->pred ) (hptr->pred)->succ = hptr->succ;
      if ( hptr->succ ) (hptr->succ)->pred = hptr->pred;
      free( hptr );
    }
    hptr = nexth;
  } 
  *firstH=fptr;
}

/*-----------------------------------------------------------------------------*/

int ssfctp_dp( int n, int demand, int* cap, double* fcost, double* tcost,
               int preproc, int* x, double* objval ) {
/*
  Description:	       
  ------------
  Computes an optimal solution to the single-sink fixed charge transportation
  problem (SSFCTP). SSFCTP is the mixed-integer program
  
    Z = min \sum_{j=1}^n ( c_j x_j + f_j y_j )
    
      s.t.: \sum_{j=1}^n x_j = D
             0 \le x_j \le k_j y_j for j=1, ...,n,
	     y_j = 0,1 for j=1, ..., n,
	     
   D > 0 is the demand of a single customer. Demand is met by product flows
   x_j that are send from n suppliers. Capacities of the suppliers (arcs from
   suppliers to the customers) are given by k_1, ..., k_n. Flow costs involve
   costs c_j per unit of flow as well as a fixed charge f_j. The binary variable
   y_j equals 1 if and only if x_j > 0.	It is assumed w.l.o.g. that k_j \le D
   for each j.     
   
   An optimal solution to the problem is computed by means of the dynamic programming
   algorithm proposed in [1]. Define the function
   
   G_i(L) = min  \sum_{j=i}^n ( c_j x_j + f_j y_j )
           s.t.: \sum_{j=i}^n x_j = D-L,
	          0 \le x_j \le k_j y_j for j=i, ..., n
		  y_j = 0,1 for j=i, ..., n
		  
   and observe that
   G_i(D) = 0
   G_i(L) = \infty for L < Lmin or L > Lmax, where  		  
            L_min = D - \sum_{j=i}^n k_j
	    L_max = \min\{ D, \sum_{j=1}^{i-1} k_j \}
   and for L = Lmin, ..., Lmax:
   
   G_i(L) = \min\{ G_{i+1}(L), \min_{x_i=1,...,k_i}\{f_i + c_i x_i + G_{i+1}(L+x_i)\}\}	    
   
   The second term in brackets can be rewritten as
   
   f_i + c_i x_i + G_{i+1}(L+x_i) 
     =  f_i - c_i L + c_i ( L+x_i ) + G_{i+1}( L+x_i )  
     =  f_i - c_i L + H_i( r_i )
     
   where r_i := L+x_i and H_i( r_i ) = c_i r_i + G_{i+1}( r_i )  
   
   Assume that for L=Lmin the function H_i(r_i) has been computed for each
   possible value r_i=L+1, ..., L+k_i and the result is stored in a list
   of pairs HList={ (r_1, H_i(r_1)), ..., (r_m, H_i(r_m)) } where
   H_i( r_1 ) < ... < H_i( r_m) and r_q is the largest value r < r_{q+1} with
   H_i(r) < H_i( r_{q+1} ). The function G_i( L ) is then given by
   
   G_i(L) = \min\{ G_{i+1}(L) + f_i - c_i L + H_1( r_1 )
   
   When moving from L to L+1, KList has to be updated. First H_i(L+1+k_i) has to 
   be computed and to be inserted in the list. Doing this all elements r with
   H_i(r) >= H_i(L+1+k_i) can be removed. 
   
   
   [1] Alidaee B, Kochenberger GA (2005) A note on a simple dynamic programming
       approach to the single-sink, fixed-charge transportation problem.
       Transp. Sci 39(1):140-143.
       
         
  Scope: To be exported
  ------
  
  Parameters:
  -----------       
  n      : number of suppliers
  demand : the sink's demand
  cap    : array of integers of size of at least n. cap[j] is the capacity the 
           arc linking supplier j with the sink.
  fcost  : array of doubles. fcost[j] is the fixed cost of arc j
  tcost  : array of doubles. tcost[j] is the unit cost of flow on arc j
  preproc: preproc=0 -> do no preprocessing; preproc=1 -> apply problem reduction test;
           preproc=2 -> apply reduction test and compute bounds on L
  x      : array of integers. On output x[j] is the optimal flow on arc j
  objval : pointer to a double containing the cost of an optimal solution
  
  The function returns 0 if no error is encountered and a positive integer otherwise.
*/

  double *G_cur=NULL, *G_prev=NULL; /*G_cur[L]:=G[i][L], G_prev[L]:=G[i+1][L]*/
  double *fc, first_ub=INFTY, hval, gval, lpbnd, lpcst[n], *dptr, *tc, upbnd, 
          value, zcnst=0.0, zfix=0.0; 
  HList  *firstH=NULL, *lastH=NULL, *hptr;
  int    **X=NULL, first_x[n], Lmin[n], Lmax[n], maxflow[n], mincflow[n], order[n], 
         xlp[n];
  int    err=0, firstr, flow, gcd=1, i, iprev, j, k, L, lastr, nn=n, q, r, xdim;
  Tnode  *jw=NULL;
 
  /* Transform data such that all costs are non-negative */
  zcnst = ssfctp_dta ( n, demand, tcost, fcost, &tc, &fc );
    
  /* Check if problem is trivially solvable by supplying the sink's demand
     from the node that minimizes the cost */
  if ( ssfctp_trivial( n, demand, cap, fc, tc, x, objval ) ) goto RETURN;

  /* Scale problem using greatest common divisor of capacities and demand */
  gcd = ssfctp_gcd( n, demand, cap );
  if ( gcd > 1 ) {
    demand /= gcd;
    for ( j=0; j < n; j++ ) tc[j] *= gcd, cap[j] /= gcd;
  }   
  
  /* Solve the LP relaxation and return if the solution is integral */
  if ( ssfctp_lp( n, demand, cap, fc, tc, lpcst, order, &q, &lpbnd, xlp, 
                  x, &upbnd ) ) {
    *objval = upbnd;
    goto RETURN;
  }  

  /* Try to reduce the problem */
  if ( preproc > 0 ) {
    /* Save solution obtained so far */
    first_ub = upbnd;
    for ( j=0; j < n; j++ ) first_x[j]=x[j];
    jw = ssfctp_red( n, &demand, cap, fc, tc, &lpbnd, lpcst, order, &q, xlp,
           &upbnd, x, &nn, &zfix, maxflow, mincflow );     
    if ( nn==0 ) {
      *objval = upbnd;
      goto RETURN;
    }	   
    /* Some arc capacities may have been reduced. Recompute therefore
       the relative costs and reorder */
    for ( j=0; j < nn; j++ ) k=order[j], lpcst[k]=tc[k]+fc[k]/maxflow[k];
    qsort( order, nn, sizeof(int), (void*) SSFCTP_compare );   
  }	   
  else 
    for ( j=0, nn=n; j < nn; j++ ) mincflow[j]=0, maxflow[j]=cap[j];

  /* Compute bounds on L */
  if ( preproc > 1 ) 
    ssfctp_bol( nn, demand, q, upbnd-lpbnd, maxflow, fc, tc, order, Lmin, Lmax );
  else
    ssfctp_nbol( nn, demand, maxflow, order, Lmin, Lmax );		
	        
  /* Allocate memory for the table G_i(L) and the corresponding solution values,
     say, X_i(L) */
  err   = 1001;   
  G_cur = (double*) calloc(demand+1,sizeof(double)); if (G_cur==NULL) goto RETURN;       
  G_prev= (double*) calloc(demand+1,sizeof(double)); if (G_cur==NULL) goto RETURN;       
  for ( i=xdim=0; i < nn; i++ ) xdim += Lmax[order[i]]-Lmin[order[i]]+1;
  X     = (int**) calloc( nn, sizeof(int*) ); if ( X==NULL ) goto RETURN;
  X[0]  = (int*) calloc( xdim, sizeof(int) ); if ( X[0]==NULL ) goto RETURN;
  for ( i=0; i < nn-1; i++ ) X[i+1] = X[i]+(Lmax[order[i]]-Lmin[order[i]]+1);
  err = 0;
             
  /* At the last stage, i=n, the function G_n(L) is easily computed */
  j = nn-1;
  i = order[j];
  for ( L=Lmin[i]; L < Lmax[i]; L++ ) {
    flow = demand-L;
    G_cur[L] = fc[i] + tc[i]*flow;  
    X[j][L-Lmin[i]] = flow;
  }
  L    = Lmax[i];
  flow = demand-L;
  G_cur[L] = (flow > 0) ? fc[i]+tc[i]*flow : 0.0;
  X[j][L-Lmin[i]] = flow;
  
  /* Backward recursion: for i=n-2, ...,0, compute the function G_i(L) and
     the corresponding solution values X_i(L) */
  for ( j=nn-2; j >= 0; j-- ) {
    iprev  = i;
    i      = order[j];
    dptr   = G_prev;
    G_prev = G_cur;
    G_cur  = dptr;   
    /* Compute the function values H_i(r) for r=L+1, ..., L+k_i */
    lastr  = MIN( Lmin[i]+maxflow[i], Lmax[iprev] );
    firstr = Lmin[i]+MAX( 1, mincflow[i] );
    firstr = MAX( firstr, Lmin[iprev] );
    for ( r=firstr; r <= lastr; r++ ) {
      hval = tc[i]*r + GFUNC(iprev,r);
      addtoHList( &firstH, &lastH, r, hval );
    }
    /* The function G_i(L) can now be evaluated for L=Lmin,...,Lmax */
    if ( firstr > lastr ) {
      /* firstr > lastr occurs if Lmin[i]+mincflow[i] > demand. Thus
         the variable x[i] must be equal to zero. This results in
         G_i(L)=G_{i+1}(L) and Lmin[i]=Lmin[i+1], Lmax[i]=Lmax[i+1].
         Furthermore, we must obviously have in any case that
         Lmin[k] >= Lmin[k+1]-cap[k] */
      dptr = G_prev, G_prev = G_cur, G_cur = dptr;
      Lmin[i] = Lmin[iprev], Lmax[i]=Lmax[iprev];
      for ( k=j-1; k >=0; k-- ) 
        Lmin[order[k]] = MAX( Lmin[order[k]], Lmin[order[k+1]]-maxflow[order[k]] );
    }
    else for ( L=Lmin[i]; L <= Lmax[i]; L++ ) {
      if ( L==demand ) {
        G_cur[L] = 0.0;
	       X[j][L-Lmin[i]]  = 0;
        break;
      }
      value = firstH->H + fc[i] - tc[i]*L;
      gval  = GFUNC(iprev,L);
      if ( value < gval ){
        G_cur[L] = value;
        X[j][L-Lmin[i]]  = (firstH->r)-L;
      } else {
        G_cur[L] = gval;
        X[j][L-Lmin[i]]  = 0;
      }  
      /* When moving from L to L+1 the range of the function H_i(r)
         is shifted from [L+1,min{L+k_i,D}] to [L+2,min{L+1+k_i,D}].
	        Compute H_i(r) for the largest new feasible value of r
	        (if not already done) and add it to the list. Remove the 
	        left most element and let the first element in the list
	        point to the minimum */
      r = MIN( lastr+1, demand );
      if ( r > lastr ) { /* insert new element in the list */
        if ( firstH->r == L+1 ) firstH->H = INFTY;
        hval = tc[i]*r+GFUNC(iprev,r);	 
        addtoHList( &firstH, &lastH, r, hval );      
      } else if ( firstH->r == L+1 ) { 
        /* remove the first element from the list and let the new first
	         element point to the minimum value */
	       if ( firstH->succ ) rearrangeHList( &firstH );   
      }
      lastr = r;
    }
    /* Free the linked list Hlist */
    while ( firstH != NULL ) {
      hptr = firstH->succ;
      free( firstH );
      firstH = hptr;
    }
    lastH=firstH;
  }   
  
  /* Reconstruct the optimal solution */
  *objval = MIN( *G_cur, upbnd )+zfix;
  if ( *G_cur < upbnd ) {
    for ( i=L=0;  i < nn; i++ ) {
      j = order[i];
      x[j] = X[i][L-Lmin[j]];
      L += x[j];  
    }
  }

RETURN: 
  for ( j=nn; j < n; j++ ) i = order[j], x[i] = maxflow[i];
  if ( jw ) {
    fc[jw->idx] = jw->fc;
    x[jw->idx] += jw->lflow;
    free( jw );
  }
  if ( *objval > first_ub - ZERO ) {
    *objval = first_ub;
    for ( j=0; j < n; j++ ) x[j]=first_x[j];
  }  
  if ( gcd > 1 ) {
    for ( j=0; j < n; j++ ) cap[j] *= gcd;
    if ( tc == tcost ) for ( j=0; j < n; j++ ) tc[j] /= gcd;
    for ( j=0; j < n; j++ ) x[j] *= gcd;
  }
  if ( zcnst < 0.0 ) {
    *objval += zcnst;
    if ( tc != tcost ) free( tc );
    if ( fc != fcost ) free( fc );
  }
  if ( G_cur ) free( G_cur );
  if ( G_prev ) free ( G_prev ); 
  if ( X ) {
    if ( X[0] ) free ( X[0] );
    free( X );
  }
  return( err );  
	       
}	        

/*-----------------------------------------------------------------------------*/

#ifdef CPLEX
int ssfctp_cpx( CPXENVptr Env, int n, int demand, int* cap, double* fcost, 
                double* tcost, int* x, double* objval ) {
/*
  Description:
  ------------
  Applies CPLEX's MIP solver to solve the SSFCTP. Before calling this procedure
  the problem is reduced by means of some reduction tests.

  Scope: To be exported
  ------
  
  Parameters:
  -----------       
  Env    : Pointer to CPLEX's enviroment (may be NULL)
  n      : number of suppliers
  demand : the sink's demand
  cap    : array of integers of size of at least n. cap[j] is the capacity the 
           arc linking supplier j with the sink.
  fcost  : array of doubles. fcost[j] is the fixed cost of arc j
  tcost  : array of doubles. tcost[j] is the unit cost of flow on arc j
  x      : array of integers. On output x[j] is the optimal flow on arc j
  objval : pointer to a double containing the cost of an optimal solution
  
  The function returns 0 if no error is encountered and a positive integer otherwise.
*/

  double cpxobj, first_ub=INFTY, lpbnd, upbnd, lpcst[n];
  int    maxflow[n], mincflow[n], order[n], xlp[n], first_x[n];
  int    col, err=0, i, j, ncols, nrows, nvlb, numnz, nn, nz, q, row, status;
  int    *matbeg=NULL, *matcnt=NULL, *matind=NULL;
  char   *pname=NULL, *sense=NULL, *ctype=NULL;
  double *lb=NULL, *matval=NULL, *obj=NULL, *rhs=NULL, *ub=NULL, *xsol = NULL;
  double *fc, *tc, zcnst=0.0, zfix=0.0;
  Tnode  *jw=NULL;
  
  CPXENVptr myEnv=NULL;
  CPXLPptr  lp=NULL;
  
  /* Transform data such that fixed-charges and transp. cost are nonnegative */
  zcnst = ssfctp_dta ( n, demand, tcost, fcost, &tc, &fc ); 

  /* Check if problem is trivially solvable by supplying the sink's demand
     from the node that minimizes the cost */
  if ( ssfctp_trivial( n, demand, cap, fc, tc, x, objval ) ) goto RETURN;
  
  /* Solve the LP relaxation and return if the solution is integral */
  if ( ssfctp_lp( n, demand, cap, fc, tc, lpcst, order, &q, &lpbnd, xlp, 
                  x, &upbnd ) ) {
    *objval = upbnd;
    goto RETURN;
  }  

  /* Save solution obtained so far */
  first_ub = upbnd;
  for ( j=0; j < n; j++ ) first_x[j]=x[j];
  /* Try to reduce the problem */
  jw = ssfctp_red( n, &demand, cap, fc, tc, &lpbnd, lpcst, order, &q, xlp,
         &upbnd, x, &nn, &zfix, maxflow, mincflow );     
  if ( nn==0 ) {
    *objval = upbnd;
    goto RETURN;
  }	 
  
  /* Create LP problem object */  
  pname = strdup("SSFCTP");
  if ( Env == NULL ) {
    Env = CPXopenCPLEX( &err );  
    myEnv = Env;
  }
  if ( err==0 ) err = CPXsetintparam( Env, CPX_PARAM_SCRIND, CPX_OFF );
  if ( err==0 ) lp = CPXcreateprob( Env, &err, pname );
  if ( err ) goto RETURN;

  /* Allocate required memory */
  for ( j=nvlb=0; j < nn; j++ ) if ( mincflow[order[j]] > 0 ) nvlb++;
  ncols = 2*nn;
  nrows = nn + 1 + nvlb;
  numnz = 3*nn + 2*nvlb;
  ctype  = (char*) calloc( ncols, sizeof(char) );
  obj    = ( double*) calloc( ncols, sizeof(double) );
  rhs    = ( double*) calloc( nrows, sizeof(double) );  
  sense  = (char*) calloc(nrows,sizeof(char));
  matbeg = (int*) calloc(ncols,sizeof(int));
  matcnt = (int*) calloc(ncols,sizeof(int));
  matind = (int*) calloc(numnz,sizeof(int));
  matval = (double*) calloc(numnz,sizeof(double));
  lb     = (double*) calloc(ncols, sizeof(double));
  ub     = (double*) calloc(ncols, sizeof(double));
  xsol   = (double*) calloc(ncols, sizeof(double));
  
  sense[0] = 'E',  rhs[0] = (double) demand;
  for ( row=1; row < nrows; row++ ) sense[row]='G', rhs[row]=0.0;
  
  /* columns of the flow variables */
  for ( col=nz=nvlb=0; col < nn; col++ ) {
    i = order[col];
    lb[col] = 0.0;
    ub[col] = maxflow[i];
    obj[col]= tc[i];
    ctype[col] = 'C';
    matbeg[col] = nz;
    matcnt[col] = 2;
    matval[nz]  = 1.0;
    matind[nz++]= 0;
    matval[nz]  = -1.0; 
    matind[nz++]= col+1;
    if ( mincflow[i] > 0 ) {
      nvlb++;
      matcnt[col]++;
      matval[nz] = 1.0;
      matind[nz++] = nn+nvlb;
    }
  }
  /* columns of the fixed-charge variables */
  for ( j=nvlb=0; j < nn; j++,col++ ) {
    i = order[j];
    lb[col]     = 0.0;
    ub[col]     = 1.0;
    obj[col]    = fc[i];
    ctype[col]  = 'B';
    matcnt[col] = 1;
    matbeg[col] = nz;
    matval[nz]  = (double) cap[i];
    matind[nz++]= j+1;
    if ( mincflow[i] > 0 ) {
      nvlb++;
      matcnt[col]++;
      matval[nz] = (double)(-mincflow[i]);
      matind[nz++] = nn+nvlb;
    }
  }

  err = CPXcopylp ( Env, lp, ncols, nrows, CPX_MIN, obj, rhs, sense, matbeg,
                    matcnt, matind, matval, lb, ub, NULL );    
  if ( err == 0 ) err = CPXsetdblparam( Env, CPX_PARAM_CUTUP, upbnd );		    
  if ( err == 0 ) err = CPXsetdblparam( Env, CPX_PARAM_EPGAP, epgap );
#ifdef __CHKTIM
  if ( err == 0 ) err = CPXsetdblparam( Env, CPX_PARAM_TILIM, tilim );  
#endif  
  if ( err == 0 ) err = CPXcopyctype( Env, lp, ctype );
  if ( err == 0 ) err = CPXmipopt( Env, lp );
  if ( err == 0 ) { 
    *objval = upbnd;
    status = CPXgetstat( Env, lp );
    if ( ( status == CPXMIP_OPTIMAL )     || 
         ( status == CPXMIP_OPTIMAL_TOL ) ||
	        ( status == CPXMIP_TIME_LIM_FEAS ) ) {
      err = CPXgetmipobjval( Env, lp, &cpxobj );
      if ( (err==0) && (cpxobj < upbnd) ) {
        *objval = cpxobj;
        err = CPXgetmipx( Env, lp, xsol, 0, nn-1 );
        if ( err == 0 ) 
	         for ( j=0; j < nn; j++ ) x[order[j]] = (int)( xsol[j]+ZERO); 
      }	
    } 
    if ( ( status==CPXMIP_TIME_LIM_FEAS ) || ( status==CPXMIP_TIME_LIM_INFEAS ) ) 
      err=status;
  }

RETURN:  
  *objval += zfix;
  for ( j=nn; j < n; j++ ) i=order[j], x[i] = maxflow[i];  
  if ( jw ) {
    x[jw->idx] += jw->lflow;
    fc[jw->idx] = jw->fc;
    free( jw );
  }
  if ( *objval > first_ub - ZERO ) {
    *objval = first_ub;
    for ( j=0; j < n; j++ ) x[j]=first_x[j];
  }  
  if ( zcnst < 0.0 ) {
    *objval += zcnst;
    if ( tc != tcost ) free( tc );
    if ( fc != fcost ) free( fc );
  }
  if ( pname ) free( pname );
  if ( ctype ) free( ctype );
  if ( obj ) free( obj );
  if ( rhs ) free( rhs );
  if ( sense ) free( sense );
  if ( matbeg ) free( matbeg );
  if ( matcnt ) free( matcnt );
  if ( matind ) free( matind );
  if ( matval ) free( matval );  
  if ( lb ) free( lb );
  if ( ub ) free( ub );
  if ( xsol ) free( xsol );  
  if ( lp ) CPXfreeprob( Env, &lp );
  if ( myEnv ) CPXcloseCPLEX( &myEnv );	       
  return( err );
  
}	       
#endif

/*-----------------------------------------------------------------------------*/

static
char Dom1( int j, int n, int bvar, int* order, int* xcur, char** dom ) {
/* Checks if there is a supplier i <= bvar with zero flow in xcur that dominates
   supplier j > bvar */

  int  i;  
  for ( i=0; i <= bvar; i++ ) if ( xcur[order[i]]==0 ) {
    if ( dom[i][j] ) return( 1 );
  }
  return( 0 );

}

/*-----------------------------------------------------------------------------*/

static
char Dom2( int j, int n, int* order, int* xcur, char** dom ) {
/* Checks if there is a supplier i < j with positive flow that is dominated
   by supplier j */

  int  i;  
  for ( i=0; i < j; i++ ) if ( xcur[order[i]] > 0 ) {
    if ( dom[j][i] ) return( 1 );
  }
  return( 0 );

}

/*-----------------------------------------------------------------------------*/
static
double bndpart ( int nn, int require, int* cap, double* lpcst, int* xcur, 
                 char* isfree, int* lpord, char* isdom ){
/*
  Description:
  ------------
  Computes the LP bound given a partial solution x[0...bvar]. Furthermore, the 
  following constraint is partially taken into account. Let F_D denote the set
  of free variables (suppliers) j for which there exists a fixed variable 
  (supplier) i with x[i]=0 that dominates supplier j. Then x[j] < cap[j] must
  hold in any improved completed solution. Since there is an optimal solution x
  where 0 < x[l] < cap[l] holds for at most one l, the constraint $\sum_{j\in
  F_D} x_j/k_j \le 1$ can be added to the LP relaxation. In order to avoid the
  larger effort, this constraint is again relaxed in the following manner.
  Assume that the set F_D is given by {1,...,q}. Define lpcst[j] =  tcost[j] +
  fcost[j]/cap[j]. Let the free variables be sorted  according to  nondecreasing
  values of lpcst. Let the variables j=1,..i, be already set to the value
  cap[j]. Set cap[l]= max{ cap[j] : j=1,...i-1 }. If cap[i] <= cap[l], then
  x[i]=0 must hold in an optimal LP solution. Otherwise  consider the (implied)
  constraint x[l] + x[i] <= cap[i]. Thus the next supplier from the set F_D can
  only be allocated an amount of cap[i]-cap[l].  		
*/		 

  double lobnd=0.0;
  int    capl=0, j, jj, xlp;

  for ( j=0, lobnd=0.0, capl=0; (j < nn) && (require > 0); j++ ) {
    jj = lpord[j];
    if ( isfree[jj] ) {
      if ( isdom[jj] ) {
        if ( cap[jj] > capl ) {
	  xlp = MIN( require, cap[jj] - capl );
	  capl= cap[jj];
	} 
	else 
	  xlp = 0;
      } 
      else
        xlp = MIN( require, cap[jj] );
      require -= xlp;
      lobnd   += xlp*lpcst[jj];	
    }
  }
  if ( require > 0 ) lobnd = INFTY;
  return( lobnd );  

}

/*-----------------------------------------------------------------------------*/

int ssfctp_enu( int n, int demand, int* cap, double* fcost, double* tcost, 
                int preproc, int* x, double* objval ) {
/*
  Description:
  ------------
  Applies an the implicit enumeration algorithm of Herer et al (1996) to the
  SSFCTP.
  
  Herer YT, Rosenblatt MJ, Hefter I (1996). Fast algorithms for single-sink
  fixed charge transportation problems with applications to manufacturing
  and transportation. Transp. Sci., 30, 276-290.

  Scope: To be exported
  ------
  
  Parameters:
  -----------       
  n      : number of suppliers
  demand : the sink's demand
  cap    : array of integers of size of at least n. cap[j] is the capacity the 
           arc linking supplier j with the sink.
  fcost  : array of doubles. fcost[j] is the fixed cost of arc j
  tcost  : array of doubles. tcost[j] is the unit cost of flow on arc j
  x      : array of integers. On output x[j] is the optimal flow on arc j
  objval : pointer to a double containing the cost of an optimal solution
  
  The function returns 0 if no error is encountered and a positive integer otherwise.
*/

  int    first_x[n], forced[n], lpord[n], maxflow[n], mincflow[n], order[n], xcur[n];
  int    bvar, err=0, i, ii, itmp, j, jj, kmin, last, *mflow=cap, nn=n, 
         nforced, q, require;
  char   ok, **dom=NULL, isdom[n], isfree[n];
  double csti, cstj, cmax, *fc, first_ub=INFTY, lpbnd, lpcst[n], *tc, tmp, upbnd, 
         zcnst=0.0, zcur, zfix=0.0; 
  Tnode  *jw=NULL;	 
#ifdef __CHKTIM
  char   checktim = (tilim < TILIM);
  int    iter=0;
  if ( checktim ) setreftime();
#endif
  
  /* Transform data such that fixed-charges and transp. cost are nonnegative */
  zcnst = ssfctp_dta ( n, demand, tcost, fcost, &tc, &fc ); 

  /* Check if problem is trivially solvable by supplying the sink's demand
     from the node that minimizes the cost */
  if ( ssfctp_trivial( n, demand, cap, fc, tc, x, objval ) ) goto RETURN;
  
  /* Solve the LP relaxation and return if the solution is integral */
  if ( preproc > 0 ) {
    if ( ssfctp_lp( n, demand, cap, fc, tc, lpcst, lpord, &q, &lpbnd, xcur,
                    x, &upbnd ) ) {
      *objval = upbnd;
      goto RETURN;
    }  
  } 
  else { 
    upbnd = INFTY;
    for ( j=0, nn=n; j < nn; j++ ) { 
      lpord[j] = j;
      lpcst[j] = tcost[j] + fcost[j]/mflow[j];
    }  
    SSFCTP_ptr = lpcst;
    SSFCTP_sgn = 1; 
    qsort( lpord, nn, sizeof(int), (void*) SSFCTP_compare );   
  }  

  /* Try to reduce the problem */
  if ( preproc > 1 ) {
    first_ub = upbnd;
    for ( j=0; j < n; j++ ) first_x[j]=x[j];
    jw = ssfctp_red( n, &demand, cap, fc, tc, &lpbnd, lpcst, lpord, &q, xcur,
           &upbnd, x, &nn, &zfix, maxflow, mincflow );     
    mflow=maxflow;	
    if ( nn==0 ) {
      *objval = upbnd;
      goto RETURN;
    }	
    /* Some arc capacities may have been reduced. Recompute therefore
       the relative costs and reorder */
    for ( j=0; j < nn; j++ ) jj=lpord[j], lpcst[jj]=tc[jj]+fc[jj]/maxflow[jj];
    qsort( lpord, nn, sizeof(int), (void*) SSFCTP_compare );   
  }		

  /* Reorder the remaining arcs according to non-decreasing unit flow cost */
  for ( j=0; j < nn; j++ ) jj = lpord[j], order[j] = jj, isfree[jj] = 1;
  SSFCTP_ptr = tc;
  SSFCTP_sgn = 1;
  qsort( order, nn, sizeof(int), (void*) SSFCTP_compare ); 
   
  /* Create the "domination matrices" dom:
  
     We say that a supplier i dominates a supplier j if
    
     a) cap[i]>=cap[j] and fcost[i]+tcost[i]*cap[j] <= fcost[j]+tcost[j]*cap[j]    
     or
     b) cap[i]< cap[j] and fcost[i]+tcost[i]*cap[i] <= tcost[j]*cap[i]
    
     holds. If a) or b) holds, then there is an optimal solution x such that
     x[i]=0 => x[j] < cap[j]. If a) holds and additionally tcost[i] > tcost[j],
     then the implication x[j] > 0 => x[i] > 0 is valid.
  */     
  err    = 1001;  
  dom    = (char**) calloc( nn, sizeof(char*) ); if ( dom==NULL ) goto RETURN;
  dom[0] = (char*) calloc( nn*nn, sizeof(char) ); if ( dom[0]==NULL ) goto RETURN;
  err     = 0;
  for ( i=1; i < nn; i++ ) dom[i] = dom[i-1]+nn;
  for ( i=0; i < nn; i++ ) {
    dom[i][i] = 0;
    ii = order[i];
    for ( j=i+1; j < nn; j++ ) {
      jj = order[j];
      if ( mflow[ii] < mflow[jj] ) {
        csti = tc[ii]*mflow[ii] + fc[ii];
        cstj = tc[jj]*mflow[ii];	 
        dom[i][j] = ( csti < cstj + ZERO );
        dom[j][i] = ( cstj + fc[jj] < csti + ZERO );
      } else {
        csti = fc[ii] + tc[ii]*mflow[jj];
        cstj = fc[jj] + tc[jj]*mflow[jj];    
        dom[i][j] = ( csti < cstj + ZERO );
        dom[j][i] = ( (mflow[ii]==mflow[jj]) && (!dom[i][j]) );
      } 
    }    
  }  
  
  /* Determine "greedy solution", that is, according to the ordering above, set
     the flow on the arcs as large as possible until demand is met */
  for ( j=0; j < nn; j++ ) xcur[order[j]]=0;   
  for ( j=0,require=demand,zcur=0.0; (j < nn) && (require > 0); j++ ) {
    jj = order[j];
    xcur[jj] = MIN( mflow[jj], require );
    zcur    += fc[jj] + tc[jj]*xcur[jj];
    require -= xcur[jj];
  }   
  if ( zcur < upbnd ) 
    for ( upbnd=zcur,j=0; j < nn; j++ ) jj=order[j], x[jj]=xcur[jj];
     
  /* Apply the implicit enumeration scheme */
  do {
    bvar = nn-1;
    while ( bvar >= 0 ) {
      jj = order[bvar];
      if ( xcur[jj] > 0 ) {
        require += xcur[jj];
        zcur    -= fc[jj] + tc[jj]*xcur[jj];
        xcur[jj] = 0;
        if ( ( bvar < nn-1 ) && ( fc[jj] > 0.0 ) )
          if ( !(Dom2(bvar, nn, order, xcur, dom)) ) break;	
      }
      bvar--;
    }  
#ifdef __CHKTIM
    if ( checktim ) {
      iter++;
      iter %= 500;
      err = CPXMIP_TIME_LIM_FEAS;
      if ( ( iter==0 ) && ( getrtime() > tilim ) ) break;
      err = 0;
    }  
#endif    
    if ( bvar < 0 ) break; /* All solutions implicitely enumerated */
    
    /* Find suppliers j > bvar that dominate a supplier i < bvar with x[i] > 0.
       For such suppliers j, the setting x[i] > 0 implies that x[j] > 0. */   
    for ( j=0; j <= bvar; j++ ) isfree[order[j]]=0;
    for ( j=bvar+1, last=nforced=0, cmax=-INFTY, kmin=demand+1; j < nn; j++ ) { 
      isfree[order[j]] = 1;
      forced[last++]=j;  
    }  
    for ( i=0; i < bvar; i++ ) if ( xcur[order[i]] > 0 ) {
      last = nn - bvar -1;
      while ( nforced < last ) {
        jj = order[forced[nforced]];
        if ( dom[forced[nforced]][i] ) {
	         if ( (tc[jj] > cmax) || ( (!(tc[jj] < cmax)) && (mflow[jj] < kmin) ) )
	           cmax=tc[jj], q=nforced, kmin=mflow[jj];
          nforced++;
        }  
        else {
          itmp = forced[--last];
          forced[last] = forced[nforced];
          forced[nforced] = itmp;
	       }
      }
    }
    if ( nforced > 1 ) {/*last forced supplier is the one with largest tcost*/
      last=nforced-1;
      itmp=forced[last], forced[last]=forced[q], forced[q]=itmp;
    }  
       
    /* All suppliers from the set of forced suppliers must supply a positive
       amount. This means that all suppliers from the set forced except one
       of the ones with largest unit cost supply their full capacity.
       Note that a contradiction results, if such a supplier j is again 
       dominated by a supplier k < j with x[k]=0 */
    for ( j=0,ok=1; (j < nforced-1) && (ok); j++ ) {
      jj = forced[j];
      ii = order[jj];
      isfree[ii]= 0;
      xcur[ii]  = mflow[ii];
      zcur     += fc[ii] + tc[ii]*xcur[ii];
      require  -= xcur[ii];
      ok        = ( ( require >= 0 ) && ( zcur < upbnd ) );
      if ( ok ) ok = !( Dom1(jj, nn, bvar, order, xcur, dom) );
    }   
    
    /* Complete the partial solution by a series of free moves */
    if ( ok ) {
      tmp = zcur;
      itmp= require;
      j   = bvar;
      while ( (require > 0) && (j < nn-1) ) {
        j++, jj=order[j];
        if ( xcur[jj]==0 ) {
	         isdom[jj] = Dom1( j, nn, bvar, order, xcur, dom );
          if ( require < mflow[jj] ) {
            xcur[jj] = require;
            zcur    += fc[jj] + tc[jj]*xcur[jj];	
          }	
          else if ( ! ( isdom[jj] ) ) {
            xcur[jj] = mflow[jj];
            zcur    += fc[jj] + tc[jj]*xcur[jj];
          }
	         else
	           isfree[jj]= 0;	  	
          require -= xcur[jj];       	
        }	
      }
      /* Compute the LP bound in order to check if current branch can
         be discarded */
      if ( require > 0 ) 
        lpbnd = INFTY;
      else {	 
        if ( zcur < upbnd ) 
          for ( upbnd=zcur,j=0; j < nn; j++ ) jj=order[j], x[jj] = xcur[jj];
        for ( jj=j+1; jj < nn; jj++ ) 
	         isdom[order[jj]]=Dom1( jj,nn,bvar,order,xcur,dom );
        lpbnd = bndpart( nn, itmp, mflow, lpcst, xcur, isfree, lpord, isdom ) 
	        + tmp;
      }
      if ( LBgUB( lpbnd, upbnd ) ) { /* discard the branch and solution */
        zcur    = tmp;
	       require = itmp;
	       for ( j=bvar+1; j < nn; j++ ) { 
	         jj=order[j];
	         if ( isfree[jj] ) xcur[jj] = 0;
	       }  
      }
    }
  } while ( 1 );

RETURN:  
  *objval = upbnd + zfix;
  for ( j=nn; j < n; j++ ) i=lpord[j], x[i]=maxflow[i];
  if ( jw ) {
    fc[jw->idx] = jw->fc;
    x[jw->idx] += jw->lflow;
    free ( jw );
  }     
  if ( *objval > first_ub - ZERO ) {
    *objval = first_ub;
    for ( j=0; j < n; j++ ) x[j]=first_x[j];
  }  
  if ( zcnst < 0.0 ) {
    *objval += zcnst;
    if ( tc != tcost ) free( tc );
    if ( fc != fcost ) free( fc );
  }
  if ( dom ) {
    if ( dom[0] ) free( dom[0] );
    free ( dom );
  }
  return( err );  
	       
}

/*-----------------------------------------------------------------------------*/

int ssfctp_mt1( int n, int demand, int* cap, double* fcost, double* tcost, 
                 int preproc, int* x, double* objval ) {
/*
  Description:
  ------------
  Applies an implicit enumeration/branch-and-bound algorithm to the SSFCTP.
  The algorithm relies on ideas of the algorithm MT1 of Martello & Toth
  for the binary Knapsack problem.

  Scope: To be exported
  ------
  
  Parameters:
  -----------       
  n      : number of suppliers
  demand : the sink's demand
  cap    : array of integers of size of at least n. cap[j] is the capacity the 
           arc linking supplier j with the sink.
  fcost  : array of doubles. fcost[j] is the fixed cost of arc j
  tcost  : array of doubles. tcost[j] is the unit cost of flow on arc j
  preproc: if greater zero, a reduced cost test is applied in order to
           reduce the problem size 
  x      : array of integers. On output x[j] is the optimal flow on arc j
  objval : pointer to a double containing the cost of an optimal solution
*/  
  int    impl[n], first_x[n], lpord[n], maxflow[n], mincflow[n], tcord[n], xcur[n], 
         xlp[n];
  int    bvar, err=0, ii, itmp, j, jj, *mflow=cap, nn=n, q, require, 
         rdem, xx;
  double *fc, first_ub=INFTY, L0, L1, lobnd0, lobnd, lpcst[n], *tc, tmp, upbnd, 
         zcnst=0.0, zcur, zfix=0.0, zlp;
  Tnode  *jw=NULL;       
#ifdef __CHKTIM
  char   checktim = (tilim < TILIM);
  int    iter=0;
  if ( checktim ) setreftime();
#endif  
  
  /* Transform data such that fixed-charges and transp. cost are nonnegative */
  zcnst = ssfctp_dta ( n, demand, tcost, fcost, &tc, &fc ); 

  /* Check if problem is trivially solvable by supplying the sink's demand
     from the node that minimizes the cost */
  if ( ssfctp_trivial( n, demand, cap, fc, tc, x, objval ) ) goto RETURN;
  
  /* Solve the LP relaxation and return if the solution is integral */
  if ( ssfctp_lp( n, demand, cap, fc, tc, lpcst, lpord, &q, &zlp, xlp, 
                  x, &upbnd ) ) {
    *objval = upbnd;
    goto RETURN;
  }  

  /* Try to reduce the problem. After a call to the reduction problem,
     the remaining free variables are still sorted according to
     nondecreasing values of lpcst */
  if ( preproc ) {
    first_ub = upbnd;
    for ( j=0; j < n; j++ ) first_x[j] = x[j];
    jw = ssfctp_red( n, &demand, cap, fc, tc, &zlp, lpcst, lpord, &q, xlp,
           &upbnd, x, &nn, &zfix, maxflow, mincflow );     
    mflow=maxflow;	
    if ( nn==0 ) {
      *objval = upbnd;
      goto RETURN;
    }	
    /* Some arc capacities may have been reduced. Recompute therefore
       the relative costs and reorder */
    for ( j=0; j < nn; j++ ) jj=lpord[j], lpcst[jj]=tc[jj]+fc[jj]/maxflow[jj];
    qsort( lpord, nn, sizeof(int), (void*) SSFCTP_compare );       
  } 
     
  /* Determine an ordering of suppliers according to nondecreasing unit cost */
  SSFCTP_ptr = tc;
  SSFCTP_sgn = 1;
  for ( j=0; j < nn; j++ ) jj=lpord[j], tcord[j]=jj, impl[j]=0, xcur[jj]=0;
  qsort( tcord, nn, sizeof(int), (void*) SSFCTP_compare ); 
   
  /* Determine the improved lower bound */
  L0 = INFTY, jj = lpord[q];
  if ( q < nn-1 ) L0 = zlp + (lpcst[lpord[q+1]]-lpcst[jj])*xlp[jj];
  L1 = zlp + ( lpcst[jj] - MAX(tc[jj],lpcst[lpord[q-1]]) )*( mflow[jj] - xlp[jj] );
  lobnd0 = MIN(L0,L1);
  
   
  /* Apply the implicit enumeration scheme */
  require = 0;
  if ( ! LBgUB( lobnd0, upbnd ) ) do {
    bvar = nn-1;
    while ( bvar >= 0 ) {
      /* branching variable is the last one set to a positive value */
      jj = lpord[bvar];
      if ( xlp[jj] > 0 ) {
        require += xlp[jj];
        zlp    -= lpcst[jj]*xlp[jj];
        xlp[jj] = 0;
        if ( ( bvar < nn-1 ) && ( fc[jj] > 0.0 ) ) break;	
      }
      bvar--;
    }  
#ifdef __CHKTIM
    if ( checktim ) { /* every 500 nodes/iterations check if time limit reached */
      iter++;
      iter %= 500;
      err = CPXMIP_TIME_LIM_FEAS;
      if ( ( iter==0 ) && ( getrtime() > tilim ) ) break;
      err = 0;
    }      
#endif
    if ( bvar < 0 ) break; /* All solutions implicitely enumerated */
         
    /* Check if excluding the arc jj corresponding to bvar implies that
       other flows have to be zero or less than the capacity */
    for ( j=bvar+1, ii=jj; j < nn; j++ ) {
      if ( abs(impl[j])-1 > bvar ) impl[j]=0; /* remove old implications */
      if ( impl[j]>=0 ) {
        if ( mflow[jj=lpord[j]] <= mflow[ii] ) {
          if ( fc[ii] + tc[ii]*mflow[jj] + ZERO < fc[jj]+tc[jj]*mflow[jj] ) {
	           if ( fc[ii] + tc[ii] + ZERO < fc[jj] + tc[jj] )
	             impl[j] = -bvar-1; /* x[ii]=0 implies x[jj]=0 */
	           else
	             impl[j] = bvar+1;  /* x[ii]=0 implies x[jj] < mflow[jj] */
	         }    
        } 
        else if ( impl[j]==0 ) 
          if ( fc[ii] + tc[ii]*mflow[ii] + ZERO < tc[jj]*mflow[ii] )
	           impl[j] = bvar+1;   /* x[ii]=0 implies x[jj] < mflow[jj] */
      }
    }    
         
     /* Complete the partial LP solution by a series of free moves */
    j = bvar, tmp = zlp, itmp= require;
    while ( (require > 0) && (j < nn-1) ) {
      jj=lpord[++j];
      if ( impl[j] >= 0 ) {
        xx = mflow[jj];	if ( impl[j] > 0 ) xx--;
        xlp[jj] = MIN( require, xx );
        zlp    += xlp[jj]*lpcst[jj];
        require-= xlp[jj]; 
      }
    }
    if ( require > 0 ) 
      lobnd = INFTY; 
    else {
      L0 = INFTY;
      if ( j < nn-1 ) L0 = zlp + ( lpcst[lpord[j+1]]-lpcst[jj] )*xlp[jj];
      do { j--; } while ( (j >= 0) && (xlp[lpord[j]] == 0) );
      if ( j >= 0 )
        L1 = zlp + (lpcst[jj]-MAX(tc[jj], lpcst[lpord[j]]))*(cap[jj]-xlp[jj]);
      else
        L1 = zlp; 	
      lobnd = MIN( L0, L1 );
    } 
    if LBgUB( lobnd, upbnd ) /* discard the branch */
      for ( zlp=tmp, require=itmp, j=bvar+1; j < nn; j++ ) xlp[lpord[j]] = 0;
    else { /* Obtain the feasible solution corresponding to the arcs
              that have positive flow in the LP solution */
      for ( j=0, rdem=demand, zcur=0.0; (j < nn) && (rdem > 0); j++ ) {
        jj = tcord[j];
        if ( xlp[jj] > 0 ) {
          xcur[jj] = MIN( mflow[jj], rdem );
          rdem    -= xcur[jj];
          zcur    += fc[jj] + tc[jj]*xcur[jj];
        } else xcur[jj] = 0;
      }
      for ( jj=j; jj < nn; jj++ ) xcur[tcord[jj]] = 0;
      if ( zcur < upbnd ) {
        for ( upbnd=zcur, j=0; j < nn; j++ ) jj=lpord[j], x[jj]=xcur[jj];
        if LBgUB( lobnd0, upbnd ) break;
        if LBgUB( lobnd, upbnd )
          for ( zlp=tmp,require=itmp,j=bvar+1; j < nn; j++ ) xlp[lpord[j]] = 0;
      }	      
    }
  } while ( 1 );

RETURN:
  *objval = upbnd + zfix;
  for ( j=nn; j < n; j++ ) jj=lpord[j], x[jj]=maxflow[jj];
  if ( jw ) {
    fc[jw->idx] = jw->fc;
    x[jw->idx] += jw->lflow;
    free( jw );
  }
  if ( *objval > first_ub - ZERO ) {
    *objval = first_ub;
    for ( j=0; j < n; j++ ) x[j]=first_x[j];
  }       
  if ( zcnst < 0.0 ) {
    *objval += zcnst;
    if ( tc != tcost ) free( tc );
    if ( fc != fcost ) free( fc );
  }
  return ( err );	       
}

