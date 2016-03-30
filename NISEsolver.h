#ifndef NISESOLVER_H_INCLUDED
#define NISESOLVER_H_INCLUDED

#include<vector>
#include<list>
#include<stdexcept>
#include<ilcplex/ilocplex.h>

#include"combo.h"
#include"vico.h"

// Used when fixing variables
enum class Fixations { ToZero , ToOne , NotFixed};
typedef std::pair<std::pair<int,int>,Fixations> FixXtrip;

typedef IloArray<IloNumArray>       IloNumMatrix;
typedef IloArray<IloNumVarArray>    IloVarMatrix;



#define NCUTS   4 //!< Number of different cuttting planes. 0=GUB, 1=LCI, 2=ECI, 3=FCP
#define GUB     0 //!< Index of the generalized upper bounds
#define LCI     1 //!< Index of the Lifted cover inequalities
#define ECI     2 //!< Index of the exted cover inequalities
#define FCP     3 //!< Index of the fenchel cutting planes



class NISEsolver{
    private:
        /**
         * @name Cplex
         * This section contains all the cplex gear needed for the algorithm to run.
         */
        ///@{
            IloEnv env;         //! Ilo-environment
            IloModel model;     //! The Ilo-model used to build the SSCFLP
            IloModel *BOCOmodel; //! Pointer to an IloModel. Used to export the cuts generated during the cutting phase of the Run method
            IloCplex cplex;     //! cplex environment to solve the model
            CPXENVptr CpxEnv;   //! Pointer to C-style CPLEX environment. Used in VICkpsep in getFCP
            int CPXerr;         //! Integer containing any CPXENV exception messages
            IloObjective OBJ;   //! Object to store the objective function
            IloNumVarArray y;   //! Vector of binary decision variables
            IloVarMatrix   x;   //! Matrix of binary decision variables
            IloVarMatrix   z;   //! Matrix of binary decision variables
            IloNumVar     f1;   //! Variable holding the value of objectve 1
            IloNumVar     f2;   //! Variable holding the value of objectve 2
        ///@}

        /**
         * @name Data
         * This section contains all the data for describing the SSCFLP.
         */
        ///@{
            int ProblemType;    //! 0=SSCFLP, 1=UFLP, 2=Knapsack, 3=Set cover, 4=FCTP
            int n;              //! The number of facilities
            int m;              //! The number of customers
            std::vector<int> d; //! Demands. d[j] = demand of customer j
            std::vector<int> s; //! Capacities. s[i] = capacity of facility i
            std::vector<std::vector<int>> c; //! Assignment costs. c[i][j] = cost of assigning customer j to facility i
            std::vector<std::vector<int>> c2;//! Cost matrix when having two objective needing matrices (fx. bi objective assingment problem)
            std::vector<int> fixed; //! Fixed opening costs. f[i] = cost of opening facility i
            std::vector<int> fixed2;//! Cost vector when having two objectives needing cost vectors (fx. bi objective knapsack problem).
            int TotalDemand;    //! Total demand. The sum of the customer demands.
        ///@}

                /**
         * @name Parameters and flags
         * This section contains a list of parameters and flags used internally in the SSCFLPsolver class.
         */
        ///@{
            bool TerminateCuttingPhase; //!< Boolean-flag used to indicate if the cutting phase should be terminated.
            bool doCuttingPhase;//!< Boolean-flag used to indicate if the cutting phase should be used.
            bool doGUB;         //!< Boolean-flag used to indicate if violated GUB constraints should be added.
            bool doLCI;         //!< Boolean-flag used to indicate if violated lifted cover inequalities should be added.
            bool doECI;         //!< Boolean-flag used to indicate if violated extend cover inequalities should be added.
            bool doFCP;         //!< Boolean-flag used to indicate if Fenchel cutting planes should be added.

            int NumCuts[NCUTS]; //!< Counts the number of cuts generated of each type
            double myZero;      //!< Below this value, and you are effectually equal to zero
            double myOne;       //!< Above this value and you are effectually equal to one
            double myTol;       //!< Tolerance for equality.
            double ObjVal;      //!< The objective function value. Used in the cutting phase.
            double OldObjVal;   //!< Objective function value of previous iteration. Used in the cutting phase.

            double* yval;       //!< Pointer to array of integers. Used to store solutions to the y-variables.
            double* xval;       //!< Pointer to array of integers. Used to store solutions to the x-variables.
        ///@}

        /**
         * @name Data structures
         * Data structures for the non-dominated frontier of BOCO and BOCO-LP
         */
            VICcut* cutlst; //! List of all cuts generated. Initialized to NULL in the constructor of the class
            std::list<std::pair<double,double> > YnLP;  //!< Set of non-dominated solutions to the BOCO-LP
            int WorstLocalNadir;
        ///@}

        /*! \brief Runs the Non-Inferior Set Estimation algorithm.
         *
         * This function runs the enhanced NISE algorithm used to strengthen the LP formulation of the BOCO problem
         * It uses a NISE algorithmic framework and calls the function CuttingPhase instead of just solving the LP
         * problems arising when changing the Lambda
         * \param doCuts boolean. Indicates if a cutting phase should be used to start with
         * \param CutModel IloModel. Model used in the BOCO solver class. Instead of parsing a list of cuts, we just append generated cuts to the CutModel
         * \param PruneTheBugger reference to a boolean. On output it is true if the current problem is feasible.
         */
        void RunNISE ( bool doCuts ,  IloModel &CutModel , bool &PruneTheBugger );

        /*! \brief Solves the model
         *
         * This funcitions solves the envoked problem using cplex and its callbacks.
         */
        void SolveModel();

        /*! \brief Generates a lifted cover inequality for the i'th capacity constraint.
         *
         * This function tries to separate a solution to the i'th capacity constraint by a lifted cover inequality
         * The procedure VIClci is used as the separation routine. It is written by Andreas Klose and can be found
         * in vico.h and vico.c
         * \param i  Integer. The index of the facility for which the capacity constraint should be separated.
         */
        VICcut* getLCI(int i);

        /*! \brief Generates an extend cover inequality for the i'th capacty constraint.
         *
         * This function tries to separate a solution to the i'th capacity constraint by a Extended Cover Inequality.
         * The procedure is implemented in VICecikl in the file vico.c.
         *  \param i Integer. The index of the facility for which the capacity constraint should be separated.
         */
        VICcut* getECI(int i);

        /*! \brief Generates a fanchel cutting plane for the i'th capacity constraint.
         *
         * This function tries to generate a Fenchel Cutting plane for the i'th capacity constraint.
         * The procedure is implemented in the function VICkpsep in vico.c
         *  \param i Integer. The index of the facility for which the capacity constraint should be separated.
         */
        VICcut* getFCP(int i);

        /*! \brief Generates a cutting plane for the total demand contraint.
         *
         * Separates a solution y* from the convex hull of integer solutions to the total demand constraint
         * that is, conv( { y\in{0,1} : sum_i s[i]*y[i]>= sum_j d[j] } )
         * Cutting planes are tried in the following order:
         *  LCI if possible,
         *  if not then a ECI if possible,
         *  else FCP if possible.
         */
        VICcut* getCutsForTD();

        /*! \brief Separates GUB constraints
         *
         * Separates Generalized Upper Bounds which are violated. Returns the violated bounds x[i][j]-y[i]<=0 as a vector of
         * index pairs (i,j).
         * \return Vector of pairs of integers. The first integer in a pair on the vector returned corresponds to the facility index while
         * the second integer corresponds to the the cutstomer index.
         */
        std::vector<std::pair<int,int> > getGUBs();

    public:
        /*!
         * Constructor of the NISEsolver class.
         */
        NISEsolver();

        /**
         * @name Model builders
         * This section contains the functions used to build the BO-ILP models we want to solve
         */
         ///@{
            /*!
             * Builds a bi-objective single source capacitated facility location problem with objectives (cx,fy). ProblemType = 0.
             * \param NumFac integer. Number of facilities.
             * \param NumCust integer. Number of customers.
             * \param demand vector of integers. d[j] is the demand of customer j.
             * \param supply vector of integers. s[i] is the supply of facility i.
             * \param FixedCost vector of integers. f[i] is the fixed cost associated with opening facility i.
             * \param VarCost vector of vectors of integers. c[i][j] is the cost of assigning customer j to facility i
             */
            void BuildSSCFLP(int NumFac, int NumCust, std::vector<int> demand, std::vector<int> supply, std::vector<int> FixedCost, std::vector<std::vector<int>> VarCost);
            /*!
             * Builds a bi-objective uncapacitated facility location problem with objective (cx,fy). ProblemType = 1.
             * \param NumFac integer. Number of facilities.
             * \param NumCust integer. Number of customers.
             * \param FixedCost vector of integers. f[i] is the fixed cost associated with opening facility i.
             * \param VarCost vector of vectors of integers. c[i][j] is the cost of assigning customer j to facility i
             */
            void BuildUFLP(int NumFac, int NumCust, std::vector<int> FixedCost, std::vector<std::vector<int>> VarCost);
            /*!
             * Builds a bi-objective knapsack problem with cost (fixed y, fixed2 y). ProlemType = 2.
             * It is assumed that all parameters are integers and that sum_i weight[i]> capacity, and that weight[i]<capacity.
             * \param NumOfItems integer. The number of items in the knapsack problem. NumOfItems>0.
             * \param Cost1 vector of integers. Cost1[i] is the cost of item i in objective 1. OBS: Coefficients need to be negative!
             * \param Cost2 vector of integers. Cost2[i] is the cost of item i in objective 2. OBS: Coefficients need to be negative!
             * \param weight vector of integers. weight[i] is the weight of item i. weight[i]>0.
             * \param capacity integer. The capacity of the knapsack problem. capacity>0.
             */
            void BuildKnapsack(int NumOfItems, std::vector<int> Cost1, std::vector<int> Cost2, std::vector<int> weight, int capacity);
            void BuildSetCover(); //! Builds a bi-objective set covering problem with cost (cx,fy) where y indicate if a column is covered. ProblemTyoe = 3.
            void BuildFCTP();   //! Builds a bi-objective fixed charge transportation problem with cost (cx, c2y
         ///@}

        /*! \brief Runs the NISE algorithm
         *
         * Runs the NISE algorithm. Its function is to be the interphase to the user.
         * \param doCuts boolean. Indicates if cuts should be added at the extremepoints
         *
         */
        std::vector<std::pair<double,double>> RUN(bool doCuts, IloModel &CutModel , bool &PruneTheBugger);

        /*!
         * Set the bounds on the variables according to the argument vector
         * \param fixed vector of pairs. The first entry of the pair is a pair of the indices of the x-variable and the second is either zero or one.
         * That is, variable x[ fixed[t].first.first ][ fixed[t].first.second ] can be fixed to fixed[t].second
         */
        void fixXVariables(std::vector< std::vector<int > > fixed);

        /*!
         * Fixes the Y variables according to the vector fixed
         * \param fixed vector of integers. fixed[i]=0 => y[i]=0, fixed[i]=1 => y[i]=1, and fixed[i]=2 0<=y[i]<=1
         */
        void fixYvariables ( std::vector<int> fixed );

        /*!
         * Sets the bounds on the objective functions
         */
        void setObjFuncBounds( double LB1, double LB2, double UB1, double UB2){
            f1.setBounds(LB1,UB1);
            f2.setBounds(LB2,UB2);
        }

                /*! \brief Runs the cutting plane phase of the cut and branch algorithm.
         *
         * This function initiate the cutting plane phase of the cut-and-branch algorithm. We continuously solve the
         * LP-relaxation of the SSCFLP and separate the LP-solution from the capacity structures using three types
         * of cutting planes
         *     Lifter Cover Inequalitites (LCI):
         *       - See: Gu Z, Nemhauser GL, Savelsbergh MWP (1998). Lifted cover inequalities
         *         for 0-1 linear programs: Computation. INFORMS J. on Computing 10:427-437
         *
         *     Extended cover inequalities (ECI):
         *       - K. Kaparis, A.N. Letchford (2010). Separation algorithms for 0-1 knapsack polytopes.
         *         Mathematical Programming. 124:69-91
         *
         *     Fenchel cutting planes (FCP):
         *       - K. Kaparis, A.N. Letchford (2010). Separation algorithms for 0-1 knapsack polytopes.
         *         Mathematical Programming, 124:69-91
         *
         * The procedure works as follows: Go through all capacity constraints and check if they are "binding enough".
         * If so, first try to separate the LP solution using an LCI. If not possible, try with a ECI. If not possible
         * either generate a FCP or prove that no violated cutting plane exists.
         */
        void CuttingPhase( IloModel & CutModel);

        /*!
         * Sets the IloModel used to export cuts
         */
        void setBOCOmodel (  IloModel *CutModel ) { BOCOmodel = CutModel; }

        /*!
         * Small function that prints out the number of cuts generated during the lifespan of the NISEsolver object
         */
        void printNcuts(){
            std::cout << "Number of GUBs : " << NumCuts[0] << std::endl;
            std::cout << "Number of LCIs : " << NumCuts[1] << std::endl;
            std::cout << "Number of ECIs : " << NumCuts[2] << std::endl;
            std::cout << "Number of FCPs : " << NumCuts[3] << std::endl;
        }

        /*!
         * Returns the VICcut list of cuts
         */
        inline
        VICcut* getCuts(){ return cutlst; }

        /*!
         * Destructor of the NISEsolver class
         */
        ~NISEsolver();
};

#endif // NISESOLVER_H_INCLUDED
