#ifndef MYNODEDATA_H_INCLUDED
#define MYNODEDATA_H_INCLUDED

#include<ilcplex/ilocplex.h> //! Including the ilocplex header
#include<vector>
#include<stdexcept>

class myNodeData:public IloCplex::BranchCallbackI::NodeData{
    private:
        std::vector< std::pair< double, double > > LBset;   //! The lower bound set of the node
        int TimesSinseLastSolved;                           //! The number of branches created since we last solved the BO-LP
        bool doNISE;                                        //! Flag stating whether or not to solve the the BO-LP from scrats
    public:
        /*! \brief Constructor of the myNodeData class
         *
         * Constructor of the myNodeData class. Copies the points int the parameter LBSET to the internal LBset, sets the TimesSinseLastSolved to TimesSinse and
         * sets doNISE to DONISE.
         * \param LBSET vector of pairs of doubles. Contains the extreme points of the farthers lower bound set + additional updates
         * \param TimesSinse interger. Contains the number of branches created on the current nodes path sinse the BO-LP was last solved
         * \param DONISE boolean. Flags if the BO-LP should be solved for this node
         */
        myNodeData ( const std::vector< std::pair< double , double > > &LBSET , int TimesSinse , bool DONISE );

        /*!
         * Returns an iterator to the first element of the vector holding the lower bound set
         */
        std::vector< std::pair<double,double> >::iterator getLBsetBegin(){ return LBset.begin( ); }

        /*!
         * Returns an iterator pointing to LBset.end( )
         */
        std::vector< std::pair<double,double> >::iterator getLBsetEnd(){ return LBset.end( ); }

        /*!
         * Returns the lower bound set of the parent node
         */
        std::vector< std::pair< double , double > > getLBset(){ return LBset;}
        /*!
         * Returns the flag doNISE indicating if we should resolve the BOLP
         */
        bool iShouldSolveBOLP(){ return doNISE;}
        /*!
         * Returns the number of branches since the last time we have solved the BOLP
         */
        int getBranchesSinseLastSolve(){ return TimesSinseLastSolved; }
};
#endif // MYNODEDATA_H_INCLUDED
