#ifndef BOCOSOLVER_H_INCLUDED
#define BOCOSOLVER_H_INCLUDED

#include<ilcplex/ilocplex.h>
#include<iostream>
#include<vector>
#include<stdexcept>
#include<algorithm>
#include<cmath>
#include<time.h>
#include "NISEsolver.h"
#include "data.h"
#include "myNodeDAta.h"
#include"NonDominatedSet.h"
#include<fstream>
#include<string>
#include<chrono>

typedef std::chrono::high_resolution_clock CPUclock;
//#include<boost/geometry.hpp>

typedef IloArray<IloNumArray>       IloNumMatrix;
typedef IloArray<IloNumVarArray>    IloVarMatrix;


class BOCOsolver{
    private:
        IloEnv env;
        IloModel model;
        IloCplex cplex;
        IloObjective OBJ;

        IloModel PrunerModel;
        IloCplex PrunerCplex;
        IloRange FirstCoord;
        IloRange SecondCoord;
        IloNumVarArray s;
        IloNumVarArray PrY;
        IloVarMatrix PrX;
        IloNumVar Prf1;
        IloNumVar Prf2;

        //! Scaling of the objectives
        int lambda1;
        int lambda2;

        VICcut* cutlst; //! List of cutting planes found by the NISEsolver

        double TheWorstLocalNadirPoint;


        std::vector<std::pair<double,double>> UBset;
        NonDominatedSet *nonDomSet;

        /*!
         * Builds the model for a single source capacitated facility location problem
         */
        void buildSSCFLP();

        /*!
         * Builds the model for an uncapacitated facility location problem.
         * I first change the data in the DATA object such that s[i]=m and d[j]=1 whereafter it calls buildSSCFLP
         * Finally, it adds the implied variable upper bounds x[i][j]-y[i]\leq 0
         */
        void buildUFLP();


    public:
        IloNumVarArray y;
        IloVarMatrix x;
        IloNumVar f1;
        IloNumVar f2;

        unsigned long LinearPruned;
        unsigned long PNPOLYpruned;
        unsigned long PNPOLYnotPruned;
        unsigned long LBupdates;
        unsigned long BOLPsolved;
        unsigned long SolutionsFound;
        unsigned long NumNodes;
        unsigned long MaxNodesInTree;
        unsigned long NumberOfParetoBranches ;
        unsigned long NumberOfOrdPB;
        unsigned long NumberOfGPB;
        unsigned long FractionalNodesCplexDecides;
        unsigned long FractioanlNodesIDecide;
        unsigned long NumCuts;
        double CuttingPhaseTime;
        int gcdOfC2;

        data *DATA;     // Pointer to data class holding the data of the instance
        /*!
         * \brief Constructor of the BCOsolver class
         *
         * Constructor of the BOCOsolver class. Creates an instance of the BOCOsolver class.
         * \param TheData data class. Copies the data from the data class TheData to the object Data contained as a member of the BOCOsolver class
         */
        BOCOsolver(const data &TheDataObject);

        /*!
         * Destructor of the BOCOsolver class. Releases memory allocated during runtime
         */
        ~BOCOsolver();

        /*!
         * Function running the bi--objective algorithm. The interface to the user.
         */
        void run();

        /*!
         * Adds a point p to UBSet if it is non--dominated by all points in UBset.
         * The set UBset is updated if it happens that the new point p dominated solutions in UBset.
         * \param p pair of ints. Contains the first and second coordinate of the point we which to inset in UBset.
         */
        void updateUBset(std::pair<double,double> p);

        /*!
         * Returns true if the outcome vector of the current node is dominated
         * \param Dominators vector of pairs. Contains on output the points which dominates the outcome vector of the current branching node (if any).
         * \param p pair of doubles. Contains the point which we test is dominated or not.
         */
        bool isItDominated(std::vector<std::pair<double,double>> &Dominators, std::pair<double,double> p);

        /*!
         * Returns true if all local nadir points lies below the linear lower bound set defined by the weighted sum scalarization
         * \param p pair of doubles. Contains the current outcome vector.
         */
        bool canItBePrunedByLinearLowerBoundSet(std::pair<double,double> p, std::vector<std::pair<double,double>> &DominatedNadirs);



        /*!
         * Returns the sum-value of the worst local nadir point (z_1+z_2) where z is the worst local nadir point
         */
        inline
        double getTheWorstLocalNadirPoint() { return TheWorstLocalNadirPoint; };

        /*!
         * Initializes the UBset with the lexicographic minimizers of the two objectives.
         * Finds the lexicographic minimizers by scaling the objectives -> seems to perform poorly
         */
        void initializeUBset();

        /*!
         * Initializes the UBSet with the lexicographix minimizers of the two objectives.
         * Finds the lexicographic minimizers by solving two problems per point.
         */
        void initializeUBset2();

        /*!
         * Returns the upper left point of the upper bound set. That is lexmin( (c¹x,c²x) : x\in X)
         */
        std::pair<double,double> getZul( ) { return *nonDomSet->UBsetBegin(); }

        /*!
         * Returns the lower left point of the upper bound set. That is lexmin( (c²x,c¹x) : x\in X)
         */
        std::pair<double,double> getZlr( ) { return *( std::prev( nonDomSet->UBsetEnd( ) ) ); }

        /*!
         * Copies the the current upper bound set to the argument vector
         * \param vector of pairs of doubles. Equals the current upper bound set on output
         */
        void getUBset ( std::vector< std::pair< double,double> > & CopyUBset );


        /*!
         *  returns the size of the UBset
         */
        inline
        size_t getUBsetSize(){ return nonDomSet->getSize(); }

        /*!
         * Returns the points on the UBset which is dominated by the linear lower bound set
         */
        std::vector<std::pair<double,double>> getDominatedPoints( double LB );

        /*!
         * Given two points z1 and z2 defining a straight line in R^2, and a point z^* and a slope of line passing through z^*, the function
         * returns the intersection point of the two defined straight lines. If the two lines are parallel the point (-1,-1) is returned.
         * \param z1 First point defining a facet of the frontier
         * \param z2 Second point defining a facet of the frontier
         * \param zstar outcome vector of the current LP solution. It is a point on the new facet.
         * \param slope The slope of the search direction of the LP
         * \return Returns an intersection point of the line defined by z1 and z2 and the line defined by zstar and slope
         */
        std::pair< double, double> getIntersectionPoint(std::pair<double,double> z1, std::pair<double,double> z2, std::pair<double,double> zstar, double slope);

        /*! \brief Method that updates the lower bound set.
         *
         * The method returns an updated lower bound set consisting of points defining the nondominated extreme points of the polygon in which solutions might lie.
         * \param OldLB vector of pairs of doubles. The lower bound which need to be updated.
         * \param doBOLP boolean. If true the bi-objective LP relaxation of the current node is solved.
         * \param LPvalue double. The optimal objective function value of scalarized LP relaxation of the current node.
         * \return An updated lower bound set.
         */
        std::vector< std::pair<double,double> > updateLBset(    const std::vector< std::pair<double,double> > &OldLB, bool doBOLP, double LPvalue , std::pair<double,double> LPoutcome,
                                                                const std::vector< std::vector<int> > &FixX, const std::vector<int> &FixY,
                                                                std::pair<double,double> Lbounds, std::pair<double,double> Ubounds , bool &PruneTheBugger );
         /**
         * @name This section contains functions and parameters for the solution of bi-objective LPs
         */
         ///@{
            NISEsolver BiObjSolver; //! Instance of the NISEsolver class. The class contains methods for solving bi-objective LPs.
            /*!
             * Builds a SSCFLP problem in the NISEsolver instance BiObjSolver
             */
            void setupANiseSSCFLP();
            /*!
             * Runs the no--inferior set estimation algorithm in order to solve a bi objective lp.
             * \param doCuts boolean. IF true, cuts are generated along the efficient frontier of the LP relaxation.
             * The cuts are storred in
             */
            std::vector<std::pair<double,double>> RunNISE(bool doCuts, bool& PruneTheBugger);

            /*!
             * Calls the two functions fixXvariables and fixYvariables for the NISEsolver
             */
            void FixNICEVariables(std::vector<int> fixY, std::vector< std::vector<int> > fixX );

            /*!
             * Sets the bounds on the objective function varaibels in the NISEsolver
             */
            void setNISEbounds( double LB1, double LB2, double UB1, double UB2);

            /*!
             * Function solving a feasibility/pruning IP.
             * It solves the problem
             * min  s1 + s2
             * st.: c1x - s1 <= z1
             *      c2x - s2 <= z2
             *      $x \in \mathcal{X}(\eta)$
             * where $\mathcal{X}(\eta)$ is the feasible set of solutions to the current branching node, and the point (z_1,z_2) is a local nadir point.
             * If the above program shows a strictly positive solution value of each nadir point, the node can be pruned by bound. Otherwise, it cannot!
             * If the parameter CheckAllNadirs is set to true, then all Nadirs points are check. Otherwise, the procedure returns as soon as a local Nadir
             * results in a solution of value zero.
             * \return true if the node can be pruned based on bound. False otherwise.
             * \param CheckAllNadirs boolean. True if all Nadirs should be checked. False is method should return as soon as a nadir point which is dominated is found
             * \param Nadir vector of pairs of doubles. Contains on input all local nadir points. If CheckAllNadirs is set to true, Nadirs contains on output all dominated local nadir points
             * \param FixX vector of vectors of ints. Passes the fixation of variables of current node. FixX[i][j] = 0 <=> x[i][j] = 0, FixX[i][j] = 1 <=> x[i][j]=1, FixX[i][j]=2 <=> 0<=x[i][j]<=1
             * \param FixY vector of ints. Passes the fixation of variables of the current node. FixY[i] = 0 <=> y[i]=0, FixY[i] =1 <=> y[i]=1, FixY[i] =2 <=> 0<=y[i]<=1
             */
            bool implicitPruner (   bool CheckAllNadirs , std::vector<std::pair<double,double>> &Nadirs, const std::vector<std::vector<int>> &FixX, const std::vector<int> &FixY ,
                                    IloNum LB1,  IloNum LB2, IloNum UB1, IloNum UB2 );
            inline
            void addToFront ( std::pair<double,double>  p ) { nonDomSet->addToFrontier( p ); }
        ///@}
};

#endif // BOCOsolver_H_INCLUDED
