/********************************************************************************

FILE      : vico.h

VERSION   : 2.9.0
CHANGE_LOG: DATE              VER.-No.  CHANGES MADE
            ---------------------------------------------------------------------
            May  12-19, 2003    1.0.0   first implementation
	    July     7, 2003    1.1.0   bug in VICuflohi removed
	    Aug      6, 2003    2.0.0   routine VIClci and VICkconf added
	    Aug      7, 2003    2.1.0   routine VICgapfn added
	    Aug     15, 2003    2.2.0   routine VICuflcov added
	    Aug     22, 2003    2.3.0   VIClci, VICkconf modified to ease handling
	    Aug     23, 2003    2.3.1   bugs in VICuflohi removed
	                                (shortest path computation now based on
					 FIFO, error in cut existence check removed)
	    Aug     26, 2003    2.3.2   minor changes, some "statics" inserted
	    Aug     27, 2003    2.3.3   minor changes: routines VICsearchcut,
	                                VICclearlst added
	    Feb     20, 2004    2.4.0	routine VICeci added
	    Feb     24, 2004    2.4.1   routine VICintsort, VICdblsort added
	    Feb     25, 2004    2.5     routine VICkpreduce added
            March   11, 2004    2.5.1   small bug in VIClci removed
            May     10, 2005    2.5.2   bug in VICkconf "(if card==n ) ... return" removed
            August   1, 2005    2.5.3   bug in kirestore() removed (all coefficients of
                                        a generated cut were modified if at least one
                                        variable was inverted!)
            May 30    , 2007    2.6.0   Started with implementing Fenchel cuts based on
                                        single-node flow structures
            June 5    , 2007    2.6.1   Different subgradient strategies for Fenchel cut
	                                generation implemented.
            June 14   , 2007    2.6.2   Different algorithms for solving the subproblem
	                                within Fenchel cut generation may optionally be
					chosen
            June 15   , 2007    2.6.3   Additional parameter "Freq" for Fenchel cut
	                                generation introduced
            June 22   , 2007    2.7.0   Additional parameter "Reduce" for Fenchel cut generation
                                        introduced. If Reduce=0 all "zero arcs" are removed and
                                        the inequality is generated for the reduced polytope
            Aug 20    , 2007    2.7.1   small bug in vicsnffenchel removed: capacities are now
                                        allowed to be zero (variable can then be ignored)
            Jan  9    , 2008    2.8.0   start to implement routine for exact knapsack separation
	    Jun 18    , 2008    2.8.1   Some bugs removed in VICkplift and VICkpsep
	    Jun 19    , 2008    2.8.2   Bug remove in VICkplift ( (t-pi) could be negative)
	    Aug 20    , 2012    2.8.3   Adjusted the uplifting in VICkplift and included possibility
	                                to exclusively fix variables of zero LP value when
					defining the reduced knapsack polytope in the exact knapsack
					separation procedure
            Dec 18    , 2012    2.9.0   Inclusion of a procedure suggested by Kaparis and Letchford (2010)
                                        to separate extended cover inequalities

LANGUAGE  : c, header file
AUTHOR    : Andreas Klose

SUBJECT   : header file for module "vico.c" (Valid Inequalities for selected
            Combinatorial Optimization problems)

NOTE      : In order to use the functions listed below from within a SUN Pascal
            program compiled with option -L using the SUN Pascal compiler,
	    compile file vico.c with option -DSUNPAS and link with libF77

	    For an example of usage, see end of this file

********************************************************************************/

#ifndef __VICO_H
#define __VICO_H

#ifndef __CPXDEFS_H
#include <ilcplex/cplex.h>
#endif

/*-------------------------------------------------------------------------------
The following data structure is used to store a cut. As an example consider the
cut x(1) + 2 x(4) + 2x(6) + 3x(8) <= 4, which is then stored as follows:

            VICcut* MyCut;

	    MyCut->sense = 'L'
	    MyCut->rhs = 4.0
	    MyCut->nzcnt = 4
	    MyCut->nzind = [ 1, 4, 6, 8 ]
	    MyCut->nzval = [ 1, 2, 2, 3 ]
	    MyCut->UsrRVal = 0.0
	    MyCut->UsrIVal = 0
	    MyCut->nextcut = NULL
-------------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

typedef struct VICcut {
  char    sense;           /* sense of inequality: 'L' means <=, 'G' is >=     */
  double  rhs;             /* right-hand side of inequality                    */
  int     nzcnt;           /* number of non-zeros in the inequality            */
  int*    nzind;           /* column/variable indices of the non-zeros         */
  double* nzval;           /* values of the non-zeros in the inequality        */
  double  UsrRVal;         /* may be used to store e.g. a dual variable        */
  int     UsrIVal;         /* may be used to store e.g. an integer flag        */
  void*   UsrDatPtr;       /* pointer to any additional user data              */
  struct  VICcut* nextcut; /* pointer to next cut                              */
} VICcut;

/*-------------------------------------------------------------------------------
The following is used for defining parameters for Fenchel cut generation
-------------------------------------------------------------------------------*/
#define SG_DEFAULT 0
#define SG_DEFLECT 1
#define SG_SMOOTH  2
#define ALG_MT1    0
#define ALG_DP     1
#define ALG_CPLEX  2

typedef struct TFENCHELopt {
  int    maxit;    /* maximum number of subgradient iterations. Default: 100 */
  int    sg_strat; /* determines the subgradient strateqy (Default: 2)
                      SG_DEFAULT: standard subgradient
		      SG_DEFLECT: subgradient deflection
		      SG_SMOOTH : exponential smoothing */
  double alpha;    /* Intial value of step size parameter (Default: 2) */
  int    H;        /* step size parameter is halfed if no improve occured in
                      lower bound after H iterations. If H=0 then the parameter
		      is continuously declined by a factor of 1.05. */
  int    CHK;      /* if > 0 then subgradient iterations are written to a file
                      named "fenchel.log". Default=0 */
  int    Algo;     /* Algorithm used for solving the single-node flow problems:
                      ALG_MT1:    ssfctp_mt1 (implicit enumeration)
		      ALG_DP:    ssfctp_dp  (dynamic programming)
		      ALG_CPLEX: cplex MIP solver
		      (Default = ALG_MT1) */
  int    Freq;     /* Frequency with which these cuts may be tried, e.g.,
                      only if no other cut found (Freq=0) or in every iteration
		      (Freq=1). Default Freg=1 */
  int    Reduce;   /* If Reduce=1, the Fenchel cut is generated for the reduced polytope
                      with all variables (x_j,y_j) with y_j=0 removed. This gives also
                      a valid inequality for the full polytope*/
} TFENCHELopt;

// Checks if a vector is almost integer
int VICisintvec( int n, double* pi );

/*-----------------------------------------------------------------------------*/
void VICsetdefaults(  );
/*------------------------------------------------------------------------------
PURPOSE:
Sets the above mentioned parameters to their above mentioned default values
------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------*/
void VICaddtolst( VICcut** first, VICcut* firstnew );
/*------------------------------------------------------------------------------
PURPOSE:
Adds a linked list of cuts to an existing linked list of cuts. The new cuts
are inserted at the top of the existing list of cuts. Let

  First -> Second -> ... -> Last -> NULL

be the list of already existing cuts and let

  Firstnew -> Secondnew -> ... -> Lastnew -> NULL

denote the linked list of new cuts. After calling the procedure, the
old list looks like

  Firstnew -> Secondnew -> ... > Lastnew -> First -> Second -> Last -> NULL

However, the memory allocated for the cuts in the new list is not copied!
Therefore, do not delete it after addition to the old list.

PARAMETERS:
- first    : pointer to the pointer to first cut in the existing list of cuts
             (if the list is empty *first must be NULL)
- firstnew : pointer to first cut in the new list of cuts (which may consists
             of just one cut)

EXAMPLE:    VICcut *MyCutsSoFar = NULL;
            VICcut *MyNewCuts = NULL;

            ProcedureForGeneratingNewCuts( &MyNewCuts );
	    if ( MyNewCuts != NULL ) VICaddtolst( &MyCutsSoFar, MyNewCuts );

-------------------------------------------------------------------------------*/

void VICfreecut ( VICcut** cut );
/*------------------------------------------------------------------------------
PURPOSE:
frees the memory allocated for a cut to which the pointer *cut points

PARAMETERS:
- cut : pointer to the pointer to the
-------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------*/
void VICfreelst ( VICcut** first );
/*------------------------------------------------------------------------------
PURPOSE:
frees the memory allocated for a linked list of cuts and empties the list

PARAMETERS:
- first : pointer to the pointer to the first cut in the list

EXAMPLE:  VICcut* MyCutList;
          many very strong and helpful cuts generated and optimum proven
	  VICfreelst( &MyCutList );
-------------------------------------------------------------------------------*/

char VICallocut( int numnz, VICcut** cut );
/*------------------------------------------------------------------------------
PURPOSE:
allocates memory required to store data of a cut

PARAMETERS:
- numnz : number of nonzero coefficients in the cut
- cut   : on completion *cut is the pointer to the cut

RETURN VALUE: 1 on success and 0 otherwise

EXAMPLE:  VICcut* MyCutPtr;
          int     numnz=1000;
	  VICallocut( numnz, &MyCutPtr );
-------------------------------------------------------------------------------*/

VICcut* VICsearchcut ( VICcut* cutlst, VICcut* cut );
/*------------------------------------------------------------------------------
PURPOSE:
searches for the cut to which the pointer "cut" is pointing to in a given
list of cuts

PARAMETERS:
- cutlst : pointer to the first cut in the linked list of cuts
- cut    : pointer to the cut which is search in the list

RETURN VALUE: NULL if the cut is not contained in the list, and the pointer
              to the cut in the list of cuts otherwise
-------------------------------------------------------------------------------*/

void VICclearlst( VICcut** first );
/*------------------------------------------------------------------------------
PURPOSE:
eliminated duplicate cuts from a linked list of cuts such that every cut
is only contained once in that list

PARAMETERS:
- first : pointer to the pointer to the first cut in the list to be cleared
-------------------------------------------------------------------------------*/

void VICsort( int n, int ascending, int doinit, int size, void* numbers, int* order );
/*------------------------------------------------------------------------------
PURPOSE: sorting of arrays of integers/doubles

PARAMETERS:
- n        : dimension of the array that has to be sorted
- ascending: sort in ascending (descending) order if ascending = 1 (0)
- doinit   : if doinit=1 it is assumed that the array "order" is not initialised
- size     : size=sizeof(int) it is assumed that numbers points to integer array
             otherwise numbers must point to array of doubles
- numbers  : pointer to array of integers/doubles of dimension of at least n
- order    : on output the sorted array is given by:
             numbers[ order[0],...,order[n-1] ]
-------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------*/
void VICcflfc( int m, int n, int* demand, int* capaci, double* x, double* y,
               VICcut** fc_cut );
/*-------------------------------------------------------------------------------
PURPOSE: tries to generate a (single) extended flow cover inequality for the
         CFLP, which is violated by the current solution (x,y). The procedure
	 may also be used to generate extended flow cover inequalities
	 for the single-node flow problem given by
	    \sum_j z(j) = d
	    z(j) <= capaci[j]*y(j)
	    0 <= z(j) <= min{d,capaci[j]}
	    y(j) \in {0,1} for all
	 In this case, call the procedure with m=1, demand = d, x=z/d

PARAMETERS:
- n      : number of potential depot sites
- m      : number of customers
- demand : pointer to an array of integers of size of at least m containing
           the customers' demands
- capaci : pointer to an array of integers of size of at least n containing
           the depot capacities
- x      : pointer to an array of doubles of size of at leat m*n containing
           the allocation part of the fractional solution, which should be
           separated by a flow cover inequality.
	   Let i=0,...,m-1 and j=0,...n-1 be the indices of customers and
	   depot sites, resp. Then x[i*n + j] is the solution value of
	   the allocation variable x(i,j), where 0 <= x(i,j) <= 1. The
	   variable x(i,j) denotes the fraction of customer i's demand met
	   from facility j.
- y      : pointer to an array of doubles of size of at least n containing
           the location part of the fractional solution, which should be
	   separated by a flow cover inequality. 0 <= y[j] <= 1
	   and x(i,j) <= y[j]
- fc_cut : pointer to a pointer to a cut. If no violated flow cover is found,
           the null pointer is returned; otherwise **fc_cut contains the cut.
	   fc_cut->nzval contains the nonzeros of the cut, and
	   fc_cut->nzind contains the column indices of the nonzeros, where
	   the index of value i*n+j corresponds to the allocation variable
	   x(i,j) and the index m*n+j corresponds to the location variable
	   y(j)

EXAMPLE:
          VICcut* MyCutList = NULL;
          VICcut* MyNewCut = NULL;
	  int     m, n, *demand=NULL, *capaci=NULL;

          Read_Problem_Data_and_Allocate_Space( m, n, demand, capaci, ... );

	  Solve_Something_like_the_LP_Relaxation_and_obtain_x_y;

	  VICcflfc( m, n, demand, capaci, x, y, &MyNewCut );
	  if ( MyNewCut != NULL ) VICaddtolst( &MyCutList, MyNewCut );
-------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------*/
void VICcflsmi( int m, int n, int* demand, int* capaci, double* x, double* y,
                VICcut** first_smi );
/*------------------------------------------------------------------------------
PURPOSE:
Generates special types of "submodular inequalities" for the CFLP using
a separation heuristic of Aardal. Several such inequalities, which cut off
the solution (x,y), may be returned.

PARAMETERS:
- n         : number of potential depot sites
- m         : number of customers
- demand    : pointer to an array of integers of size of at least m containing
              the customers' demands
- capaci    : pointer to an array of integers of size of at least n containing
              the depot capacities
- x         : pointer to an array of doubles of size of at leat m*n containing
              the allocation part of the fractional solution, which should be
              separated by a submodular inequality.
	      Let i=0,...,m-1 and j=0,...n-1 be the indices of customers and
	      depot sites, resp. Then x[i*n + j] is the solution value of
	      the allocation variable x(i,j), where 0 <= x(i,j) <= 1. The
	      variable x(i,j) denotes the fraction of customer i's demand met
	      from facility j.
- y         : pointer to an array of doubles of size of at least n containing
              the location part of the fractional solution, which should be
	      separated by a submodular inequality. x(i,j) <= y[j] <= 1
- first_smi : pointer to the first cut of a linked list of violated
              submodular inequalities found by this procedure.
	      The NULL pointer is returned if no violated submodular
	      inequality was found. Column indices of nonzeros are
	      defined in the same way as in the case of the routine
	      VICcflfc

EXAMPLE:
          VICcut* MyCutList = NULL;
          VICcut* MyNewCut = NULL;
	  int     m, n, *demand=NULL, *capaci=NULL;

          Read_Problem_Data_and_Allocate_Space( m, n, demand, capaci, ... );

	  Solve_Something_like_the_LP_Relaxation_and_obtain_x_y;

	  VICcflsmi( m, n, demand, capaci, x, y, &MyNewCut );
	  if ( MyNewCut != NULL ) VICaddtolst( &MyCutList, MyNewCut );

------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------*/
void VICuflohi( int m, int n, double* x, double* y, VICcut** first_ohi );
/*-------------------------------------------------------------------------------
PURPOSE:
Generates odd-hole inequalities violated by the solution (x,y) for the UFLP.
Several such cuts may be returned.

PARAMETERS:
- n         : number of potential depot sites
- m         : number of customers
- x         : pointer to an array of doubles of size of at leat m*n containing
              the allocation part of the fractional solution.
	      Let i=0,...,m-1 and j=0,...n-1 be the indices of customers and
	      depot sites, resp. Then x[i*n + j] is the solution value of
    	      the allocation variable x(i,j), where 0 <= x(i,j) <= 1. The
	      variable x(i,j) denotes the fraction of customer i's demand met
	      from facility j.
- y         : pointer to an array of doubles of size of at least n containing
              the location part of the fractional solution, which should be
	      separated by a submodular inequality. 0 <= y[j] <= 1, y[j]>=x[i,j]
- first_ohi : pointer to the first cut of a linked list of violated odd-hole
              inequalities found by this procedure. The NULL pointer is returned
	      if no violated inequality was found. See routine VICcflfc regarding
	      definition of column indices.

EXAMPLE:
          VICcut* MyCutList = NULL;
          VICcut* MyNewCut = NULL;

	  Solve_Something_like_the_LP_Relaxation_and_obtain_x_y;

	  VICuflohi( m, n, x, y, &MyNewCut );
	  if ( MyNewCut != NULL ) VICaddtolst( &MyCutList, MyNewCut );

-------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------*/
void VICuflcov( CPXENVptr Env, int m, int n, char SolveCov, double* x, double* y,
                VICcut** first_cut );
/*-------------------------------------------------------------------------------
PURPOSE:
Tries to find violated combinatorial inequalities for the uncapacitated
facility location problem.

Let K be a subset of the set of all customers, and let J denote a subset
of the set of potential depot sites. Define a binary matrix a(i,j) for
each i \in K and j \in J. Let b denote (a lower bound on) the minimum number
of depots j \in J required to cover each customer i \in K, that is

   b =   min \sum_{j\in J} y_j
       s.t.: \sum_{j\in J} a_{ij} y_j \ge 1 \forall i \in K
             y_j = 0,1 \forall j \in J

Then the inequality

  \sum_{i \in K} \sum_{j\in J} a_{ij} x_{ij} - \sum_{j\in J} y_j \le |K| - b

is valid for the UFLP. These inequalites generalize odd holes for the UFLP
and have been proposed by D.C. Cho et al. (1983): "On the uncapacitated
facility location problem I: Valid inequalities and facets", Mathematics
of Operations Research 8, 579-589. See also G. Cornuejols, J.-M. Thizy
(1982): "Some facets of the simple plant location polytope", Mathematical
Programming 23, 50-74.

Let (x*, y*) denote a fractional solution. In order to find a violated
inequality of the above type, the following simple heuristic is tried:

  (i)    Set J = \{ j : 0 < y*_j < 1 \}
  (ii)   Set K = \{ i : x*_{ij} > 0 for more than one j \in J \}
  (iii)  Solve the covering problem in order to find b (alternatively, find
         a lower bound on b by solving the LP relaxation of the covering
	 problem and rounding up the objective function value)
  (iv)   Check if point (x*,y*) violates the resulting inequality


PARAMETERS:
- n         : number of potential depot sites
- m         : number of customers
- x         : pointer to an array of doubles of size of at leat m*n containing
              the allocation part of the fractional solution.
	      Let i=0,...,m-1 and j=0,...n-1 be the indices of customers and
	      depot sites, resp. Then x[i*n + j] is the solution value of
    	      the allocation variable x(i,j), where 0 <= x(i,j) <= 1. The
	      variable x(i,j) denotes the fraction of customer i's demand met
	      from facility j.
- y         : pointer to an array of doubles of size of at least n containing
              the location part of the fractional solution, which should be
	      separated by a submodular inequality. 0 <= y[j] <= 1, y[j]>=x[i,j]
- first_cut : pointer to cut found by this procedure. The NULL pointer is returned
	      if no violated inequality was found. See routine VICcflfc regarding
	      definition of column indices.
-------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------*/
void VICkpreduce( int WHATRED, int TRYCLI, int n, int* cap, char sense,
                  int* weight, int* order, int* indx, VICcut** clique );
/*-------------------------------------------------------------------------------
PURPOSE: Clique generation and coefficient improvement (reduction) for a single
         knapsack inequality of the form
	    \sum_j a_j x_j <= b  or  \sum_j a_j x_j >= b where x_j=0,1 for all j
	 The procedure tries to generate a single clique inequality from the
	 knapsack inequality and to reduce the right-hand side together with
	 the coefficients of variables in the clique

PARAMETERS:
- WHATRED: Determines which coeffiecient improvement scheme is applied:
           WHAT_RED = 0 -> no coefficient improvement at all
	   WHAT_RED = 1 -> just simple coefficient improvement
	                   (increase w_j to c if x_j=1 implies all other x_k = 0)
           WHAT_RED = 2 -> only apply reduction of coefficients of
                           variables in the clique
           WHAT_RED = 3 -> do both types of improvement if possible
- TRYCLI: If TRYCLI=0 no cliques are generated, otherwise cliques
          are derived if possible. ( In order to use the coefficient reduction,
	  that is WHATRED >= 2; a clique is required and TRYCLI must equal 1)
- n     : number of (free) variables in the knapsack
- cap   : pointer to an integer containing the capacity of the
          knapsack (right-hand side of inequality). *cap is possible changed
	  (reduced)
- sense : sense of inequality ('L' for <= and 'G' for >=)
- weight: pointer to an array of integers of length of at least n containig the
          coefficients ("weights") of (free) variables in the inequality.
	  Some of the weights[j] may be changed (reduced)
- order : NULL or a pointer to an array of integers of length of at least n
          containing an ordering of the items 0,...,n-1 according to increasing
	  weights
- indx  : NULL or a pointer to an array of integers of length of at least n
          containing indices of the variables in the inequality. If indx=NULL,
	  it is assumed that the variables are numbered from 0, ..., n-1
- clique: pointer to the the generated clique cut
-------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------*/
void VIClci( int n, int cap, char sense, int* weight, int* indx, double* xlp,
             double* rco, VICcut** lci );
/*-------------------------------------------------------------------------------
PURPOSE: Heuristic procecedure of Gu, Nemhauser and Savelsbergh for generating
         lifted cover inequalties from a knapsack structure like

	      \sum_{j \in N} a_j x_j <= b
	      0 <= x_j <= 1   \forall j \in N
	      x_j \in \{0,1\} \forall j \in N

	 or
	      \sum_{j \in N} a_j x_j >= b
	      0 <= x_j <= 1   \forall j \in N
	      x_j \in \{0,1\} \forall j \in N

         See: Gu Z, Nemhauser GL, Savelsbergh MWP (1998). Lifted cover inequalities
         for 0-1 linear programs: Computation. INFORMS J. on Computing 10:427-437

PARAMETERS:
- n       : number of (free) variables appearing in the knapsack inequality
- cap     : right-hand side of the knapsack inequality
- sense   : sense (that is 'L' for <= or 'G' for >=) of the knapsack inequality
- weight  : coefficient a_j of (free) variables in the knapsack inequality
            (coefficients a_j are not restricted to be nonnegative!)
- indx    : indices of the (free) variables in the knapsack inequality. If null
            it is assumed that variables are numbered from 0 to n-1
- xlp     : LP solution of (free) variables appearing in the knapsack inequality
- rco     : absolute values of reduced costs of variables appearing in the
            knapsack inequality
- lci	  : pointer to the generated cut pointer.

EXAMPLE:
Consider the following MIP

   max 7x(0) + 3x(1) + 6x(2) + 9x(3) + 10x(4) + 6x(5) + 8x(6) + 9x(7)

   s.t.:   x(0) + 2 x(1) +   x(2)                   <= 2
         3 x(3) + 5 x(4) + 2 x(5) + 6 x(6) + 7 x(7) <= 11
	 x(j) = 0,1 for all j

Furthermore, assume that variables x(3) is fixed to one and variable
x(5) is fixed to zero. In order to derive a lifted cover inequality from
the second inequality, the procedure above has to be called in the
following way:

 n = 3
 cap = 11-3 = 8
 sense = 'L'
 weight = (5, 6, 7 )
 indx = (4, 6, 7)
 xlp = ( x(4), x(6), x(7) )
 rco = ( |redcost(4)|, |redcost(6)|, |redcost(8)| )

-------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------*/
void VICkconf( int n, int cap, char sense, int* weight, int* indx, double* xlp,
               double* rco, VICcut** kconf );
/*-------------------------------------------------------------------------------
PURPOSE:
Tries to get a (1,k)-configuration inequality using the separation heuristic of
Crowder, Johnson, Padberg in Oper. Res. 31 (1983).
For an alternative separation heuristic for (1,k)-configurations see
Carlos E. Ferreira (1997). On Combinatorial Optimization Problems arising in
Computer Systems Design. Phd Thesis, Technische Universit\E4t Berlin.

Given the Knapsack-Polytop
         X = { x : \sum_{j\in N} w[j]*x[j] <= c, x_j = 0,1 }
a (1,k)-configuration is a set NP \cup {z}, where NP \subset N, such that

(i)  \sum_{j \in NP} w[j] <= c
(ii) The set K \cup {z} is a cover with respect to N for all
     subsets K of NP with cardinality k

The corresponding (1,k)-configuration inequality is given by

        (r - k + 1)x[z] + \sum_{j\in NP} x[j] <= r, where r=|NP|

Crowder, Johnson, Padberg propose the following separation heuristic:

1. Let S \subset N be the cover, and z\in S the item with
   maximum weight. Set NP = S-{z} and k=|NP|.
2. For all j \in N-S with \sum_{l \in N-S} w[l] <= c do:
   a. Check if K \cup \{z} is a cover for any K \subseteq NP with |K|=k
   b. If this is the case build the corresponding (1,k)-configuration
      inequality and lift it.
   c. If the lifted inequality is violated, add the inequality

PARAMETERS:
- n       : number of (free) variables appearing in the knapsack inequality
- cap     : right-hand side of the knapsack inequality
- sense   : sense (that is 'L' for <= or 'G' for >=) of the knapsack inequality
- weight  : coefficient a_j of (free) variables in the knapsack inequality
            (coefficients a_j are not restricted to be nonnegative!)
- indx    : indices of the (free) variables in the knapsack inequality. If null
            it is assumed that variables are numbered from 0 to n-1
- xlp     : LP solution of (free) variables appearing in the knapsack inequality
- rco     : absolute values of reduced costs of variables appearing in the
            knapsack inequality
- kconf	  : pointer to the first cut in a linked list of generated k-configuration
            cuts

EXAMPLE: same as in case of VIClci

-------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------*/
void VICeci( int n, int cap, char sense, int* weight, int* order, int* indx,
             double* xlp, VICcut** eci );
/*-------------------------------------------------------------------------------
PURPOSE: Separation procedure of Gabrel & Minoux for finding most violated
         extended cover inequality. See:
         - V. Gabrel, M. Minoux (2002). A scheme for exact separation of
           extended cover inequalities and application to multidimensional
           knapsack problems. Oper. Res. Lett. 30:252-264

PARAMETERS:
- n       : number of (free) variables appearing in the knapsack inequality
- cap     : right-hand side of the knapsack inequality
- sense   : sense (that is 'L' for <= or 'G' for >=) of the knapsack inequality
- weight  : coefficient a_j of (free) variables in the knapsack inequality
            (coefficients a_j are not restricted to be nonnegative!)
- order   : NULL or an ordering of the items in the knapsack according to
            increasing weights
- indx    : indices of the (free) variables in the knapsack inequality. If null
            it is assumed that variables are numbered from 0 to n-1
- xlp     : LP solution of (free) variables appearing in the knapsack inequality
- eci	  : pointer to the generated cut pointer.

EXAMPLE:
Consider the following MIP

   max 7x(0) + 3x(1) + 6x(2) + 9x(3) + 10x(4) + 6x(5) + 8x(6) + 9x(7)

   s.t.:   x(0) + 2 x(1) +   x(2)                   <= 2
         3 x(3) + 5 x(4) + 2 x(5) + 6 x(6) + 7 x(7) <= 11
	 x(j) = 0,1 for all j

Furthermore, assume that variables x(3) is fixed to one and variable
x(5) is fixed to zero. In order to derive a lifted cover inequality from
the second inequality, the procedure above has to be called in the
following way:

 n = 3
 cap = 11-3 = 8
 sense = 'L'
 weight = (5, 6, 7 )
 indx = (4, 6, 7)
 xlp = ( x(4), x(6), x(7) )
 order = NULL or order = ( 0, 1, 2 )

-------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------*/
void VICecikl( int n, int cap, int do_exact, char sense, int* weight, int* indx,
               double* xlp, VICcut** eci );
/*-------------------------------------------------------------------------------
PURPOSE: Procedure for generating extended cover inequalities as suggested in
K. Kaparis, A.N. Letchford (2010). Separation algorithms for 0-1 knapsack
polytopes, Math. Prog. 124:69-91.

PARAMETERS:
- n       : number of variables appearing in the knapsack inequality
- cap     : right-hand side of the knapsack inequality
- do_exact: if equal to 1, exact separation is tried. This requires to repeatedly
            solve a 0-1 knapsack problem. If equal to 0 these knapsack problems
            are solved heuristically using a greedy method.
- sense   : sense (that is 'L' for <= or 'G' for >=) of the knapsack inequality
- weight  : coefficient a_j of variables in the knapsack inequality
- indx    : indices of the variables in the knapsack inequality. If NULL, it is
            assumed that variables are numbered from 0 to n-1
- xlp     : LP solution of (free) variables appearing in the knapsack inequality
- eci	  : pointer to the generated cut pointer (NULL if none found)

-------------------------------------------------------------------------------*/



/*-----------------------------------------------------------------------------*/
void VICgapfn( int m, int n, int IsGap, int* capaci, int** weight, double* X,
               int* indx, VICcut** fn_cut );
/*-------------------------------------------------------------------------------
PURPOSE:
Given the GAP ( or alternatively LEGAP )

       max   \sum_{i\in I} \sum_{j\in J} c_{ij} x_{ij}

       s.t.:         \sum_{i\in I} x_{ij} = 1   \forall j \in J \\

             \sum_{j\in J} w_{ij} x_{ij} <= s_i \forall i \in I \\

	             x_{ij} \in \{0,1\} \forall i,j

the procedure uses a heuristic described in Farias IR, Nemhauser GL (2001),
A family of inequalities for the generalized assignment polytope, Oper. Res.
Letters 29, for finding a violated inequality of the form

\sum_{j\in J_i} w_{ij} x_{ij} + \sum_{j\in J_i} a_j \sum_{k\ne i} x_{kj} <= s_i

where J_i \subseteq J, a_j = s_i - ( W_i - w_{ij} ), W_i = \sum_{j\in J_i} w_{ij}


PARAMETERS:
- m     : number of agents
- n     : number of jobs
- IsGap : set equal to 1 if the problem is a GAP, that is if every job has to be
          assigned to exactly one agent. Otherwise a LEGAP is assumed, that is
	  some jobs may be not assigned to any agent
- capaci: pointer to an array of integers of length of at least m containing the
          agents' capacities
- weight: pointer to an array of pointers of length of at least m such that
          weight[i][j] gives the amount of resources required by agent i
	  to perform job j
- X     : pointer to an array of doubles of length of at least m*n containing
          the fractional LP solution, X[i*n+j] must contain the LP value of
	  assignment variable x_{ij}
- indx  : if not NULL this array gives the indices of the assignment
          variables in the user's problem formulation, that is indx[i*n+j] is
	  the index of variable x_{ij}. if indx==NULL it is assumed that
	  variables are numbered from 0 to m*n-1, where i*n+j is the index
	  of variable x_{ij}
- fn_cut: pointer to the first cut in a linked list of generated cuts of this
          type
-------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------*/
//void VICsnffenchel( CPXENVptr Env, int n, int demand, int* capaci, double* x,
//                    double* y, int* indx, int* indy, VICcut** cut );
/*-------------------------------------------------------------------------------
PURPOSE:
This procedure tries to find a fenchel cutting plane that cuts of a fractional
solution (x*,z*) to the single node flow problem

        \sum_{j=1}^n x(j) = demand,
        0 <= x(j) <= z(j) for j=1,...,n,
        z(j)=0 or z(j)=capaci(j) for j=1,...,n.

Let X be the set of all feasible solutions (x,z) to the single-node flow problem
above. A fenchel cut is then given by

         \sum_{j=1}^n \pi(j) x(j) - \sum_{j=1} \lambda(j) z(j) <= \pi_0,

where (\pi, \lambda, \pi_0) is a solution to the separation problem

          \max  \pi x* - \lambda z*
          s.t.: \pi x  - \lambda y <= \pi_0 \forall (x,y)\in X
                0 <= \pi(j) <= 1.0 forall j,
                0 <= \lambda(j) <= 1.0 forall j.

The separation problem is solved by means of subgradient optimization.

PARAMETERS
- Env    : pointer to CPLEX environment as returned by CPXopenCPLEX. Env may be NULL
           if FCHELopt.Algo = ALG_MT1 or FCHELopt.Algo=ALG_DP. If Env=NULL and
	   FCHELopt.Algo = CPLEX, FCHELopt.Algo is automatically reset to
	   ALG_MT1.
- n      : number of arcs from the n sinks to the single source
- demand : the sink's demand
- capaci : pointer to an array of integers of size of at least n containing
           the arcs' capacities
- x      : pointer to an array of doubles of size of at leat n containing
           the values of the flow variables in the fractional solution
- y      : pointer to an array of doubles of size of at least n containing
           the values to the binary variables in the fractional solution,
           0 <= y[j] <= 1 (the variable z(j) is defined as z(j)=capaci(j)y(j)),
- indx   : NULL or a pointer to an integer array of size of at least n containing
           the indices of the x-variables. If indx is NULL, the indices 0,..,n-1
           are used.
- indy   : NULL or a pointer to an integer array of size of at least n containing
           the indices of the y-variables. If indy is NULL, the indices n,..., 2n-1
           are used.
- cut    : If NULL, no cut was found. Otherwise points to the generated cut.
           REMARK: Coefficients of the cut are in terms of variables (x,y) not
           (x,z)!
-------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------*/
void VICkpsep( CPXENVptr Env, int justUp, int n, int cap, char sense, int* weight,
               int* indx, double* xlp, double* rco, VICcut** cut );
/*-------------------------------------------------------------------------------
PURPOSE:
Subroutine for exact knapsack separation. Given a fractional solution x*, the most
violated valid inequality

  \pi x \le 1

is returned by optimizing over the 1-polar of the knapsack polytope, that is by
solving the separation problem

     max\{ \pi x* : \pi x^t \le 1 for each feasible solution x^t \}.

The dual of the this separation problem is solved by means of column generation.
In order to reduce the effort, the separation is performed for a reduced knapsack
polytope obtained by resp. fixing variables x_j to 0 and 1 if x*_j=0 or x^*_j=1.
Afterwards, a sequential lifting of fixed variables is applied in order to get
a valid cut.

PARAMETERS:
- Env     : pointer to the CPLEX environment as returned by CPXopenCPLEX.
            Env must be a valid pointer to the open CPLEX environment
- justUp  : if equal to one, only variables showing a value of zero in the
            LP solution are first fixed to zero. The resulting inequality
	    is then just uplifted to obtain a (strengthened) inequality
	    for the full polytope. Otherwise (justUp=0) both variables
	    showing LP-value of zero and one are fixed. The resulting
	    inequality is then first to be downlifted to obtain a valid
	    inequality, and thereafter uplifting is applied.
- n       : number of variables appearing in the knapsack inequality
- cap     : right-hand side of the knapsack inequality
- sense   : sense (that is 'L' for <= or 'G' for >=) of the knapsack inequality
- weight  : coefficient a_j of variables in the knapsack inequality
            (coefficients a_j are not restricted to be nonnegative!)
- indx    : indices of the variables in the knapsack inequality. If NULL
            it is assumed that variables are numbered from 0 to n-1
- xlp     : LP solution of variables appearing in the knapsack inequality
- rco     : absolute values of reduced costs of variables appearing in the
            knapsack inequality
- cut     : pointer to the generated cut pointer.

-------------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif



/*==============================================================================

EXAMPLE (cuts for CFLP)
-----------------------

  VICcut* All_Generated_Cuts;    // pointer to first cut in a linked list
                                    containing all cuts generated so far //

  VICcut* Cuts_of_Current_Round; // pointer to first cut in a linked list
                                    containing cuts generated in current
				    round of generating cutting planes //

  VICcut* Curr_Cut;             // pointer to first cut in a linked list
                                   of cuts generated by one of the routines
			           described above //

  int m, n;                     // number of customers and pot. depot sites //

  double **cost;                // cost[i][j] = cost of supplying all of
                                   customer i's demand from facility j //

  double *fixcost;              // fixed depot costs //

  int *demand, *capaci;         // customer demands and depot capacities //

  double *x, *y;                // pointer to fractional solution (x,y) //


  All_Generated_Cuts = NULL;
  Cuts_of_Current_Round = NULL;

  Read_CFLP_Data_and_Allocate_Mem(&m, &n, demand, capaci, cost );

  x = (double*) calloc( m*n, sizeof(double) );
  y = (double*) calloc( n, sizeof(double) );


  Solve_First_LP( Data, x, y ); // (x,y) denotes LP-solution of a CFLP //
                                // x[i*n+j] is percentage of customer i's demand
				   met from facility j in this solution, where
				   m = #customers, n =#depots //

  if ( xy_Is_Integer ) return(); // LP solution is integer //

  do {

    Cuts_of_CurrentRound = NULL;

    VICcflfc( m, n, demand, int* capaci, x, y, &Curr_Cut );
    if ( Curr_Cut ) VICaddtolst ( &Cuts_of_Current_Round, Curr_Cut );

    VICcflsmi( m, n, demand, capaci, double* x, double* y, &Curr_Cut );
    if ( Curr_Cut ) VICaddtolst ( &Cuts_of_Current_Round, Curr_Cut );

    VICuflohi( m, n, x, y, &Curr_Cut );
    if ( Curr_Cut ) VICaddtolst ( &Cuts_of_Current_Round, Curr_Cut );

    Select_Cuts_You_Think_Are_Helpful ( &Cuts_of_Current_Round );
    if ( Cuts_of_Current_Round ) {
      Add_Cuts_to_Current_LP_Relaxation( Cuts_of_Current_Round, ... );
      VICaddtolst( &All_Generated_Cuts, Cuts_of_Current_Round );
    }

  } while ( Cuts_of_Current_Round != NULL );

  // Everythink done, free memory //
  VICfreelst( &All_Generated_Cuts );

===============================================================================*/
#endif // ifdef VICO_H
