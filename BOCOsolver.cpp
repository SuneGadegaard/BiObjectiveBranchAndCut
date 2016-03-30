#include "BOCOsolver.h"


bool ResolveBOLPIfAllBelow = true;
bool doUpdateOfLbSet = true;
bool PrintFrontiers = true;
bool doLPbasedPruneCheck = false;
bool PerformEPB = true;
bool CheckAllNadirs;
bool doBoundSetBasedBaranching = true;  // This means we do implicit branching instead
bool DoCutsAtRoot = true;
using namespace std::chrono;

std::ofstream myfile;

int gcd(int a, int b) {
    return b == 0 ? a : gcd(b, a % b);
}

std::string printPair(const std::pair<double,double> p){
    std::string strPair = std::to_string(p.first);
    strPair += ("\t" + std::to_string(p.second) + "\n");
    return strPair;
}


// New test for dominace based on linear programming!
// Nadir contains on output the nadir points which are dominated!
bool canIPrune ( const std::vector<std::pair<double,double>> &LB, std::vector<std::pair<double,double>> &Nadir )
{
    try
    {
        IloInt numVars = LB.size ( );
        size_t NumNadirs = Nadir.size ( );
        std::vector<std::pair<double,double>> DominatedNadirs;

        IloEnv env;
        IloModel model = IloModel ( env );
        IloCplex cplex = IloCplex ( model );
        IloNumVarArray lambda = IloNumVarArray ( env , numVars , 0 , 1 , ILOFLOAT );
        IloNumVarArray slack = IloNumVarArray ( env , 2 , 0 , IloInfinity , ILOFLOAT );
        IloRange FirstCoord, SecondCoord;
        IloExpr firstLHS(env), secondLHS(env), ConvexLHS(env);

        cplex.setOut( env.getNullStream ( ) );
        cplex.setWarning ( env.getNullStream ( ) );

        model.add ( IloMinimize( env , slack[0] + slack[1] ) );
        //build the constraints
        for ( int l=0; l<numVars; ++l )
        {
            firstLHS += (LB[l].first)*lambda[l];
            secondLHS += (LB[l].second)*lambda[l];
            ConvexLHS += lambda[l];
        }

        FirstCoord = ( firstLHS - slack[0] <= Nadir[0].first );
        SecondCoord = ( secondLHS - slack[1] <= Nadir[1].second );

        model.add ( FirstCoord );
        model.add ( SecondCoord );
        model.add ( ConvexLHS == 1 );


        for ( size_t n = 0; n<NumNadirs; ++n )
        { // Loop through all nadir points and see if the lower bound set dominates any of them
            FirstCoord.setUB( Nadir[n].first );
            SecondCoord.setUB( Nadir[n].second );
            if ( cplex.solve () )
            {
                if ( cplex.getObjValue( ) <= 0.0 )
                {
                    // The nadir is dominated
                    DominatedNadirs.push_back( Nadir[n] );
                }
            }
        }
        // Overwrite the Nadir vector with only the dominated ones
        Nadir.clear ( );
        Nadir = DominatedNadirs;

        // Cleaning up
        FirstCoord.end ( );
        SecondCoord.end ( );
        firstLHS.end ( );
        secondLHS.end ( );
        ConvexLHS.end ( );
        lambda.end ( );
        slack.end ( );
        cplex.end ( );
        model.end ( );
        env.end ( );

        return ( DominatedNadirs.empty ( ) );


    }
    catch ( std::exception &e )
    {
        std::cerr << "Error in the canIPrune : " << e.what ( ) << "\n";
        exit ( EXIT_FAILURE );
    }
}


// ---------------------------------------------------------------------------------------------//
// PIP function returning 1 if the point defined by (testx,testy) is in the polygon defined by the
// points in the two arrays vertx,verty. Same as below fuction, just using raw pointers as arguments
int pnpoly(int nvert, double *vertx, double *verty, double testx, double testy)
{
  int i, j, c = 0;
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
	 (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
       c = !c;
  }
  return c;
}

// ---------------------------------------------------------------------------------------------//
// PIP function returning 1 if the point testPoint is in the poylogon defined by the points in the poly-vector
// Same as above fuinction, just using vectors and pairs as arguments
int pnpoly(std::vector<std::pair<double,double>> &poly, std::pair<double,double> &testPoint)
{
  int nvert = poly.size();
  int i, j, c = 0;
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((poly[i].second>testPoint.second) != (poly[j].second>testPoint.second)) &&
	 (testPoint.first < (poly[j].first-poly[i].first) * (testPoint.second-poly[i].second) / (poly[j].second-poly[i].second) + poly[i].first) )
       c = !c;
  }
  return c;
}
// ---------------------------------------------------------------------------------------------//
ILOMIPCALLBACK1(Terminator, BOCOsolver&, solver ){
    try{
        IloNum BestLB = getBestObjValue();
        double Nadir = 0.0 ;//solver.getTheWorstLocalNadirPoint();
        double WorstNadir = 0.0;
        std::vector<std::pair<double,double>> UBset;
        solver.getUBset( UBset );
        for ( auto it = std::next ( UBset.begin() ); it!= UBset.end(); ++it )
        {
            Nadir = it->first + std::prev( it )->second;
            if ( Nadir > WorstNadir ) WorstNadir = Nadir;
        }
        if ( BestLB >= WorstNadir ) abort ( );
    }catch ( IloException &ie ){
        std::cerr << "IloException in Terminator : " << ie.getMessage() << std::endl;
        exit ( EXIT_FAILURE );
    }
}
// ---------------------------------------------------------------------------------------------//
ILOUSERCUTCALLBACK1(BoundCutter, BOCOsolver &, solver){
    try{
        int n = solver.DATA->getNumFaci();
        int m = solver.DATA->getNumCust();

        for ( int i=0; i<n; ++i )
        {
            for ( int j=0; j<m; ++j )
            {
                if (  getValue ( solver.x[i][j]) - getValue ( solver.y[i] )  > 0.01 ) add ( solver.x[i][j]-solver.y[i] <= 0 );
            }
        }
    }catch(IloException &ie){
        std::cerr << "IloException in BoundCutter : " << ie.getMessage() << std::endl;
    }
}

// ---------------------------------------------------------------------------------------------//
ILOINCUMBENTCALLBACK1(LexiUpdater, BOCOsolver &, solver){
    try{
        std::pair<double,double> p(0.0,0.0);
        p.first     = getValue(solver.f1);
        p.second    = getValue(solver.f2);
        ++solver.SolutionsFound;
        solver.updateUBset(p);
    }catch(std::exception &e){
        std::cerr << "Exception in the LexiUpdater incumbent callback : " << e.what() << std::endl;
    }catch(IloException &ie){
        std::cerr << "IloException in the LexiUpdater incumbent callback : " << ie.getMessage() << std::endl;
    }
}

/*
 * Implementation of the branch call back used when updating the LB set by use of the linear programming relaxation
 */
ILOBRANCHCALLBACK1(BOCObrancher, BOCOsolver &, solver){
    try{

        // Parameters and data containers for the callback
        bool            canBePruned = true,     // Flag used to indicate whether the node can be pruned or not. Assume no!
                        PruneTheBugger = false; // Flag used to indicate that the node can be pruned
        const double    myOne       = 0.9999,   // Tolerance for one. If binary variable x>=myOne then x=1
                        myZero      = 1.0-myOne;// Tolerances for zero. If binary x<=myZero then x=0

        std::pair< double, double > firstIntersection, secondIntersection; // Intersection points between the lower bound set of father node and of the linear lower bound of the node

        std::vector< std::pair<double,double>>  DominatedNadirs,    // Vector of double-pairs used to store Nadir points dominated by the linear lower bound
                                                UBset,              // Vector of double-pairs used to store the current upper bound set
                                                LBset,              // Vector of double-pairs used to store the current lower bound set
                                                poly,               // Vector of double-pairs used to store all the points defining the polygon
                                                DefiningIntersections, // Vector of double-pairs used to store the points defining the line that is intersected by the linear lower bound
                                                Dominators,         // Vector used to store outcome vectors dominating the current outcome vector. Used when doing Pareto branching
                                                Nadirs;
        std::vector< std::pair<double,double>>::iterator UpdateAfterMe;
        solver.getUBset( UBset );   // Copy the upper bound set

        /*if (  getNnodes() >= NumNodes  )
        {   // Write the current upper bound set to file
            std::string FileName = "UBset_" + std::to_string ( NumNodes );
            myfile.open( FileName );
            for ( auto it = UBset.begin ( ); it!=UBset.end ( ); ++it )
            {
                myfile << it->first << "\t" << it->second << "\n";
            }
            myfile.close( );
            NumNodes += 1000;
        }*/

        // Update max number of nodes in tree
        if ( solver.MaxNodesInTree < getNremainingNodes( ) ) solver.MaxNodesInTree = getNremainingNodes ( );

        // Retrieve the current outcome vector
        IloNum  f1 = getValue( solver.f1 ), // Value of objective 1 in the LP solution
                f2 = getValue( solver.f2 ), // Value of objective 2 in the LP solution
                LB = getObjValue ( ),       // The value of the scalarized LP solution
                UB1 = getUB ( solver.f1 ),  // The upper bound on objective 1
                UB2 = getUB ( solver.f2 ),  // The upper bound on objective 2
                LB1 = getLB ( solver.f1 ),  // The lower bound on objective 1
                LB2 = getLB ( solver.f2 );  // The lower bound on Objective 2


        // We start by checking if the node can be pruned by the linear lower bound alone
        for ( auto it=++UBset.begin(); it!=UBset.end(); ++it ){
            if ( it->first + std::prev(it)->second > LB ){ // If the Nadir point is above the linear lower bound set, add it to the potential Nadirs to check
                DominatedNadirs.push_back( std::pair<double,double>( it->first , std::prev(it)->second ) );
                // At the same time check if last added Nadir is in fact dominated by (f1,f2) as then we definitely have Nadirs in the polygon
                if ( f1<=DominatedNadirs.back().first - 0.9 && f2<=DominatedNadirs.back().second-0.9 ) canBePruned = false;
            }
            if ( it->first <= LB1 && it->second <= LB2 ){ // If we have a point on UBset which dominates this nodes ideal point, the node can be pruned!
                prune();
                return;
            }
        }
        if( DominatedNadirs.empty() ){ // If we have no dominated Nadirs we can prune the node, all solutions in this subproblem is dominated!
            solver.LinearPruned++;
            prune();
            return;
        }else if (UBset.begin()->first <= LB1 && UBset.begin()->second <= LB2 ){
            // We did not check this first point against the ideal point in the loop before, so we do it now!
            prune();
            return;
        }

        /* If it could not be pruned by the linear lower bound or ideal pointand there is still hope (doPOLY=true), we retrive the information from the NodeData object
           so we know if we should solve the bi objective LP relaxation in order to get a lower bound set */
        NodeData *ND = getNodeData();
        if ( getNnodes() > 0 ){
            if ( nullptr == ND ){
                std::cout << "Number of nodes " << getNnodes() << std::endl;
                throw std::runtime_error( "It appears that the current node has no NodeData object attached. Terminating the procedure!\n" );
            }
        }

        myNodeData *mND = (myNodeData*) ND;

        // Get the bounds of the variables and tell the NISEsolver which ones to fix
        std::vector<std::vector<int> > FixX;
        std::vector<int> FixY;
        for ( int i=0; i<solver.DATA->getNumFaci() ; ++i){
            // First we check if the upper bound is zero and the if the lower bound is one
            if ( getUB ( solver.y[i] ) <= myZero ) FixY.push_back( 0 );     // Fixed to o
            else if ( getLB ( solver.y[i] ) >= myOne ) FixY.push_back( 1 ); // Fixed to 1
            else FixY.push_back( 2 );                                       // Not fixed
            // Then we check the x-variables
            FixX.push_back( std::vector<int>() );
            for ( int j=0 ; j<solver.DATA->getNumCust() ; ++j){
                if ( getUB ( solver.x[i][j] ) <= myZero ) FixX[i].push_back( 0 );
                else if ( getLB ( solver.x[i][j] ) >= myOne ) FixX[i].push_back( 1 );
                else FixX[i].push_back ( 2 );
            }
        }
        // Call the function that updates the lower bound set
        std::pair<double,double> Outcome(f1,f2);
        std::pair<double,double> LBs(LB1,LB2);
        std::pair<double,double> UBs(UB1,UB2);

        if ( mND != nullptr ){
            if ( doUpdateOfLbSet )
            {
                LBset = solver.updateLBset( mND->getLBset() ,  mND->iShouldSolveBOLP()  , LB , Outcome , FixX, FixY , LBs , UBs , PruneTheBugger);
            }
            else
            {
                LBset = solver.updateLBset( mND->getLBset() ,  false  , LB , Outcome , FixX, FixY , LBs , UBs , PruneTheBugger);
            }
        }else{
            LBset = solver.updateLBset( LBset , true , LB , Outcome , FixX, FixY , LBs , UBs , PruneTheBugger );
        }
        if ( PruneTheBugger ){
            std::cout << "I am pruning based on the PruneTheBugger flag\n";
            prune();
            return;
        }

        // Now copy the lower bound set into the poly-vector and augment with the corner points defined by the the bounding box rect(z^ul,z^lr)
        poly = LBset;
        poly.push_back( std::pair<double,double> ( solver.getZlr().first , LBset.back().second    ) );
        poly.push_back( std::pair<double,double> ( solver.getZlr().first , solver.getZul().second ) );
        poly.push_back( std::pair<double,double> ( LBset.begin()->first  , solver.getZul().second ) );

        // Check each of the Nadir points above the linear lower bound line if they are in the polygon defined by the set point in the vector poly
        if ( !doLPbasedPruneCheck )
        {
            for ( auto it=DominatedNadirs.begin();  it!=DominatedNadirs.end(); ++it ){
                std::pair<double,double> p(it->first + 0.001 , it->second + 0.001 );
                //if ( p.first <= LBset.begin()->first || p.second <= LBset.back().second ) ; // The point p is definetely not in the polygon!
                if ( pnpoly( poly, p ) )
                {
                    canBePruned = false; // pnpoly returns 1 the point p lies in the polygon defined by poly
                    Nadirs.push_back ( std::pair<double,double> ( it->first, it->second ) );
                }
                if ( Nadirs.size ( ) >= 3 ) break;
            }

            if (!canBePruned ) solver.PNPOLYnotPruned++;
            else{
                prune();
                solver.PNPOLYpruned++;
                return;
            }
        }
        else
        {
            if ( canIPrune( LBset , DominatedNadirs ) )
            {
                prune ( );
                ++solver.PNPOLYpruned;
                return;
            }else ++solver.PNPOLYnotPruned;
        }



        /* Ok, we could not prune the node. Instead we do some branching! First we check if we can perform Pareto branching.
           Then, if that could not be done, we try to branch based on the dominated Nadir points. Finally, if that did not succeed
           either we perform branching based on a variable dichotomy */

        if ( solver.isItDominated ( Dominators , std::pair<double,double>(f1,f2) ) ){
            /* Enter this if-scope if and only if the current outcome vector of the LP-relaxation is dominated by a solution on the upper bound set.
               If that is the case, we perform Pareto-branching ala Stidsen, Andersen and Dammann "A Branch and Bound Algorithm for a Class of
               Biobjective Mixed Integer Programs" */

            // Create arrays to hold the branching information
            IloNumVarArray  vars ( getEnv () );
            IloNumArray     bounds ( getEnv () );
            IloCplex::BranchDirectionArray dirs ( getEnv( ) );
            // Fill in the variables and create node data

            vars.add ( solver.f1 );
            vars.add ( solver.f2 );




            if ( PerformEPB && ( !Nadirs.empty() ) && ( Nadirs.size () <= 2 ) )
            {   // We can perform generalized pareto branching (GPB)
                // When doing GPB we always branch down
                double f1_bound, f2_bound;
                if ( Nadirs.size ( ) == 1 )
                {
                    // We do not need to care about disjoining the branches
                    dirs.add ( IloCplex::BranchDown );
                    dirs.add ( IloCplex::BranchDown );
                    if ( Nadirs.begin ( )->first <= Dominators.begin ( )->first )
                    {
                        // If the Nadir point lies to the west of the first dominator we can use that info
                        f1_bound = std::min ( Nadirs.begin()->first , Dominators.begin()->first - 1.0  );
                        f2_bound = Nadirs.begin()->second;
                    }
                    else
                    {
                        // Then it must lie to the right of the last dominator.
                        f1_bound = Nadirs.begin ( )->first;
                        f2_bound = std::min ( Nadirs.begin ( )->second , Dominators.begin ( )->second -1.0 );
                    }
                    bounds.add ( f1_bound);
                    bounds.add ( f2_bound );
                    NodeData* mND1 = new ( getEnv( ) ) myNodeData ( LBset, 0 , true );
                    makeBranch ( vars, bounds, dirs, getObjValue ( ) , mND1 );
                    ++solver.NumberOfGPB;
                    ++solver.NumberOfParetoBranches;
                }
                else
                {
                    // We need the branches to be disjoint. So we add a lower bound on the second objective depending second Nadir
                    // SOrt the Nadirs, just to be sure
                    std::sort ( DominatedNadirs.begin(), DominatedNadirs.end ( ) );

                    // Check the location of the first dominated nadir point
                    if ( Nadirs.begin ( )->first <= Dominators.begin ( )->first )
                    {
                        // If the Nadir point lies to the west of the first dominator we can use that info
                        f1_bound = std::min ( Nadirs.begin()->first , Dominators.begin()->first - 1.0  );
                        f2_bound = Nadirs.begin()->second;
                    }
                    else
                    {
                        // Then it must lie to the right of the last dominator.
                        f1_bound = Nadirs.begin ( )->first;
                        f2_bound = std::min ( Nadirs.begin ( )->second , Dominators.begin ( )->second -1.0 );
                    }

                    // Make the first branch
                    vars.add ( solver.f2 );
                    dirs.add ( IloCplex::BranchDown );
                    dirs.add ( IloCplex::BranchDown );
                    dirs.add ( IloCplex::BranchUp );

                    bounds.add ( f1_bound );
                    bounds.add ( f2_bound );
                    bounds.add ( Nadirs.back ( ).second +1 );
                    NodeData* mND1 = new ( getEnv( ) ) myNodeData ( LBset, 0 , true );
                    makeBranch ( vars, bounds, dirs, getObjValue ( ) , mND1 );
                    ++solver.NumberOfGPB;
                    ++solver.NumberOfParetoBranches;

                    // Check the location of the second dominated nadir point
                    if ( Nadirs.back( ).first <= Dominators.back( ).first )
                    {
                        // If the Nadir point lies to the west of the first dominator we can use that info
                        f1_bound = std::min ( Nadirs.back( ).first , Dominators.back( ).first - 1.0  );
                        f2_bound = Nadirs.back( ).second;
                    }
                    else
                    {
                        // Then it must lie to the right of the last dominator.
                        f1_bound = Nadirs.back( ).first;
                        f2_bound = std::min ( Nadirs.back( ).second , Dominators.back( ).second -1.0 );
                    }
                    // Make the second branch
                    vars.clear ( ); bounds.clear ( ); dirs.clear ( );
                    vars.add ( solver.f1 );
                    vars.add ( solver.f2 );
                    bounds.add ( f1_bound );
                    bounds.add ( f2_bound );
                    dirs.add ( IloCplex::BranchDown );
                    dirs.add ( IloCplex::BranchDown );
                    NodeData* mND2 = new ( getEnv( ) ) myNodeData ( LBset, 0 , true );
                    makeBranch ( vars, bounds, dirs, getObjValue ( ) , mND2 );
                    ++solver.NumberOfGPB;
                    ++solver.NumberOfParetoBranches;
                }
            }
            else
            {
                // Do ordinary pareto branching
                // First branch to the north west of the first point on the Dominators vector
                ++solver.NumberOfOrdPB;
                ++solver.NumberOfParetoBranches;
                NodeData* mND1 = new ( getEnv( ) ) myNodeData ( LBset, 0 , true );
                bounds.add ( Dominators.begin( )->first  -1.0 );
                bounds.add ( Dominators.begin( )->second +1.0 );
                dirs.add ( IloCplex::BranchDown );
                dirs.add ( IloCplex::BranchUp );
                makeBranch ( vars, bounds, dirs, getObjValue( ) , mND1 );
                // Prepare for the second branch
                ++solver.NumberOfOrdPB;
                ++solver.NumberOfParetoBranches;;
                NodeData* mND2 = new ( getEnv( ) ) myNodeData( LBset, 0 ,true );
                bounds.clear( ); dirs.clear( ); // Clearing old information
                bounds.add ( Dominators.back( ).first  +1.0 );
                bounds.add ( Dominators.back( ).second -1.0);
                dirs.add ( IloCplex::BranchUp );
                dirs.add ( IloCplex::BranchDown );
                makeBranch ( vars, bounds, dirs, getObjValue( ), mND2 );
            }
            // Free memory of the IloArrays
            vars.end();
            bounds.end();
            dirs.end();
            return;
        }
        else
        {
            // Okay, we have no more fancy ideas, so we just let cplex decide. We retrieve Cplex proposal and creates the nodes accordingly
            IloInt NrBranches = getNbranches ( );
            if ( 2 == NrBranches )
            { // Cplex wants to create two branches, so we do as we are told
                IloNumVarArray vars ( getEnv( ) );
                IloNumArray bounds ( getEnv( ) );
                IloCplex::BranchDirectionArray dirs ( getEnv( ) );
                // Get the branching info of the first branch ( index 0 )
                getBranch ( vars, bounds , dirs , 0 );
                // Create the node data object appended to the first node
                NodeData* mND1 = new ( getEnv( ) ) myNodeData ( LBset, 0 , false );
                // make the first branch
                makeBranch ( vars , bounds , dirs , getObjValue() , mND1 );
                // Get the branching info for the second branch (index 1 )
                getBranch ( vars , bounds , dirs , 1 );
                // Create the node data object appended to the second node
                NodeData* mND2 = new ( getEnv( ) ) myNodeData ( LBset, 0 , true );
                // Make the second branch
                makeBranch ( vars , bounds , dirs , getObjValue() , mND2 );

                // Free memory
                vars.end();
                bounds.end();
                dirs.end();
                return;
            }else if ( 1 == NrBranches ){ // Cplex wants to create one branch, so we do as we are told
                IloNumVarArray vars ( getEnv( ) );
                IloNumArray bounds ( getEnv( ) );
                IloCplex::BranchDirectionArray dirs ( getEnv( ) );
                // Get the branching information of the node Cplex wants to create
                getBranch ( vars , bounds , dirs , 0 );
                // Create the node data object appended to the node
                NodeData* mND = new ( getEnv( ) ) myNodeData ( LBset, 0 , false );
                // Make the branch
                makeBranch( vars , bounds , dirs , getObjValue() , mND );

                // Free memory allocated to IloArrays
                vars.end( );
                bounds.end( );
                dirs.end( );
                return;
            }else{ // We can only make an educated guess on what Cplex would have done. We branch on the variable with the largest pseudo cost score. Priority given to indicator variables
                // Iterate through all variables and get their pseudo up/down costs
                bool foundY = false;
                IloNum  psUp,  // Pseudo up cost
                        psDown, // Pseudo down cost
                        psScore, // weighted score
                        bestYScore = -IloInfinity, // Best score for y variables
                        bestXScore = -IloInfinity, // best score for x variables
                        VarVal; // Value of the variable in the LP solution
                IloNumVar BranchVar( getEnv( ) ); // variable used to branch on
                IloExpr noGood( getEnv( ) ); // IloExpression used to store the no-good constraint we might need to use if all variables turn out to be integral

                for ( int i=0; i<solver.DATA->getNumFaci(); ++i ){
                    VarVal = getValue ( solver.y[i] );
                    if ( ( VarVal > myZero) && ( VarVal < myOne ) ){ // Look for fractional y-variables before looking for fractional x-variables
                        psDown  = getDownPseudoCost( solver.y[i] ); // Get the pseudocost of branching down
                        psUp    = getUpPseudoCost( solver.y[i] );   // Get the pseudocost of branching up
                        psScore = 0.83*IloMin ( psDown , psUp ) + 0.17*IloMax( psDown , psUp); // Get the score ( based on "Branching rules revisited" by Achterberg et al.)
                        if ( psScore > bestYScore ){ // If score improves, store the variable and update best
                            bestYScore = psScore;
                            BranchVar = solver.y[i];
                            foundY = true;
                        }
                    }
                    if ( !foundY && VarVal > myZero ){ // Only look for fractional x-variables if we have not found fractional y and if y[i]>0
                        for ( int j=0; j<solver.DATA->getNumCust(); ++j )
                        {
                            VarVal = getValue( solver.x[i][j] );
                            if ( ( VarVal > myZero ) && ( VarVal < myOne ) )
                            { // x[i][j] is fractional
                                psDown  = getDownPseudoCost( solver.x[i][j] );
                                psUp    = getUpPseudoCost( solver.x[i][j] );
                                psScore = 0.83*IloMin ( psDown , psUp ) +  0.17*IloMax ( psDown , psUp );
                                if ( psScore > bestXScore )
                                {
                                    bestXScore = psScore;
                                    BranchVar = solver.x[i][j];
                                }
                            }
                        }
                    }

                }

                if ( IloMax ( bestYScore , bestXScore ) >= myZero ){
                    // Create the node data object we need to add to the first branching node
                    NodeData* mND1 = new ( getEnv( ) ) myNodeData ( LBset, 0 , false );
                    // Make the first branch
                    makeBranch ( BranchVar , 0 , IloCplex::BranchDirection::BranchDown , getObjValue() , mND1 );
                    // Create the node data object we need to add to the second branchign node
                    NodeData* mND2 = new ( getEnv( ) ) myNodeData ( LBset, 0 , false );
                    // Make the second branch
                    makeBranch ( BranchVar , 1 , IloCplex::BranchDirection::BranchUp , getObjValue() , mND2 );
                }else{
                    // Create the node data we need to append to the single branching node with the no-good inequality
                    NodeData* mND = new ( getEnv( ) ) myNodeData ( LBset, 0 , false );
                    // Make the no-good branch
                    makeBranch( (noGood>=1) , getObjValue() , mND );
                }
                BranchVar.end();
                noGood.end();
                return;
            }
        }

    }catch(std::exception &e){
        std::cerr << "Exception in the BOCObrancher callback : " << e.what( ) << std::endl;
        exit ( EXIT_FAILURE );
    }catch(IloException &ie){
        std::cerr << "IloException in the BOCObrancher callback : " << ie.getMessage( ) << std::endl;
        exit ( EXIT_FAILURE );
    }
}


//---------------------------------------------------------------------------------------------------------------------
ILOBRANCHCALLBACK1 ( ImplicitBrancher, BOCOsolver &, solver )
{
    try
    {
        bool canBePruned = true;    // Assume we can prune the ndoe
        int     n   = solver.DATA->getNumFaci ( ),  // Retrieve the number of facilities
                m   = solver.DATA->getNumCust ( );  // Retrieve the number of customers
        IloNum  f1  = getValue( solver.f1 ), // Value of objective 1 in the LP solution
                f2  = getValue( solver.f2 ), // Value of objective 2 in the LP solution
                LB  = getObjValue ( ),       // The value of the scalarized LP solution
                UB1 = getUB ( solver.f1 ),  // The upper bound on objective 1
                UB2 = getUB ( solver.f2 ),  // The upper bound on objective 2
                LB1 = getLB ( solver.f1 ),  // The lower bound on objective 1
                LB2 = getLB ( solver.f2 );  // The lower bound on Objective 2


        std::vector<std::pair<double,double>> UBset;
        std::vector<std::pair<double,double>> Nadirs;
        std::vector<std::vector<int>> FixX ( n );
        std::vector<int> FixY ( n );


        solver.getUBset ( UBset );

        // Update max number of nodes in tree
        if ( solver.MaxNodesInTree < getNremainingNodes( ) ) solver.MaxNodesInTree = getNremainingNodes ( );

        /*if (  getNnodes() >= NumNodes  )
        {   // Write the current upper bound set to file
            std::string FileName = "UBset_" + std::to_string ( NumNodes );
            myfile.open( FileName );
            for ( auto it = UBset.begin ( ); it!=UBset.end ( ); ++it )
            {
                myfile << it->first << "\t" << it->second << "\n";
            }
            myfile.close( );
            NumNodes += 1000;
        }*/

        // Loop through the upper bound set and create the Nadirs points which are relevant. That is nadir points in {(LB1,LB2) + R_\geq^2}
        for ( auto it= std::next( UBset.begin ( ) ); it!=UBset.end ( ); ++it )
        {
            // Check if the Nadir point lies above the linear lower bound
            if ( it->first + std::prev ( it )->second > LB )
            {
                // Check if the Nadir point lies in the box containing this current subproblem
                if ( it->second >= LB2 && std::prev ( it )->first >= LB1)
                //if ( true )
                {
                    Nadirs.push_back( std::pair<double,double>( it->first , std::prev( it )->second) );
                    // At the same time check if last added Nadir is in fact dominated by (f1,f2) as then we definitely have dominated Nadirs
                    if ( ( f1 <= Nadirs.back().first ) && ( f2 <= Nadirs.back().second ) ) canBePruned = false;
                }
            }

            // If we have a point on UBset which dominates this nodes ideal point, the node can be pruned!
            if ( it->first <= LB1 && it->second <= LB2 ){
                std::cout << "Pruning based on ideal point\n";
                prune();
                return;
            }
        }

        // If we have no Nadir points lying above the lower bound defined by half plane {y : y1+y2>=LB } we can prune the ndoe
        if ( Nadirs.empty ( ) )
        {
            prune ( );
            return;
        }
        else if ( ( UBset[0].first <= LB1 ) && ( UBset[0].second <= LB2 ) )
        {
            // We did not make this test for the first point on the UBset in the loop above
            std::cout << "Pruning based on ideal point\n";
            prune ( );
            return;
        }

        //=================================================================================
        // The node could not be pruned by any easy means, so we try the implicit pruner
        // instead
        //=================================================================================

        if ( canBePruned )
        { // If there is still hope of pruning the node, we gather the info on variable fixations and hand that to the implicitPruner function
            for ( int i=0; i<n; ++i )
            {
                FixX[i] = ( std::vector<int>( m ) );
                if ( getUB( solver.y[i] ) <= 0.001 ) FixY[i] = 0;
                else if ( getLB ( solver.y[i] ) >= 0.999 ) FixY[i] = 1;
                else FixY[i] = 2;

                for ( int j=0; j<m; ++j )
                {
                    if ( getUB ( solver.x[i][j] ) <= 0.001 ) FixX[i][j] = 0;
                    else if ( getLB ( solver.x[i][j] ) >= 0.999 ) FixX[i][j] = 1;
                    else FixX[i][j] = 2;
                }

            }

            // The implicit pruner returns true if the node can be pruned.
            if ( solver.implicitPruner( CheckAllNadirs , Nadirs , FixX , FixY , LB1 , LB2 , UB1 , UB2  ) )
            {
                prune ( );
                // We want to check all Nadir points as soon as we have pruned the first node, but no earlier
                // Only check all, if we are actually performing EPB
                if ( PerformEPB && false == CheckAllNadirs ) CheckAllNadirs = true;
                return;
            }
        }
        // =================================================================================
        // Okay, it was not possible to prune the node, so we need to branch!
        // =================================================================================

        // First thing to do, is to see if we can do some Pareto branching! To that end, check if (f1,f2) is dominated by any points on UBset
        std::vector<std::pair<double,double>> Dominators;
        for ( auto it = UBset.begin ( ); it!=UBset.end ( ); ++it )
        {
            if ( ( it->first <= f1 ) && ( it->second <= f2 ) )
            {
                Dominators.push_back( *it );
            }
        }

        // Create arrays to hold the branching information
        IloNumVarArray  vars ( getEnv () );
        IloNumArray     bounds ( getEnv () );
        IloCplex::BranchDirectionArray dirs ( getEnv( ) );

        // If the vector of Dominators is non-empty, we can perform pareto branching!
        if ( !Dominators.empty ( ) )
        //if ( false )
        {
            // Fill in the variables and create node data
            vars.add ( solver.f1 );
            vars.add ( solver.f2 );
            // If we are checking all the Nadir points, and doing EPB, we are able to perform generalized Pareto branching (GPB)
            // if we have less than three dominated nadirs
            if ( PerformEPB && CheckAllNadirs && ( !Nadirs.empty ( ) ) && ( Nadirs.size ( ) <= 2 ) )
            {
                // Check if a maximum of two new nodes is enough to perform the GPB
                if ( ( !Nadirs.empty ( ) ) && ( Nadirs.size ( ) < 3 ) )
                {
                     ++solver.NumberOfGPB;

                    if ( Nadirs.size ( ) == 1 )
                    {
                        // We do not need to care about disjoining the branches
                        dirs.add ( IloCplex::BranchDown );
                        dirs.add ( IloCplex::BranchDown );
                        bounds.add ( Nadirs.begin ( )->first - 0.001 );
                        bounds.add ( Nadirs.begin ( )->second - 0.001 );
                        makeBranch ( vars, bounds, dirs, getObjValue ( ) );
                        ++solver.NumberOfGPB;
                        ++solver.NumberOfParetoBranches;
                    }
                    else
                    {
                        // We need the branches to be disjoint. So we add a lower bound on the second objective depending second Nadir
                        // SOrt the Nadirs, just to be sure
                        std::sort ( Nadirs.begin(), Nadirs.end ( ) );

                        // Make the first branch
                        vars.add ( solver.f2 );
                        dirs.add ( IloCplex::BranchDown );
                        dirs.add ( IloCplex::BranchDown );
                        dirs.add ( IloCplex::BranchUp );

                        bounds.add ( Nadirs.begin ( )->first - 0.001 );
                        bounds.add ( Nadirs.begin ( )->second - 0.001 );
                        bounds.add ( Nadirs.back ( ).second +1 );
                        makeBranch ( vars, bounds, dirs, getObjValue ( )  );
                        ++solver.NumberOfGPB;
                        ++solver.NumberOfParetoBranches;

                        // Make the second branch
                        vars.clear ( ); bounds.clear ( ); dirs.clear ( );
                        vars.add ( solver.f1 );
                        vars.add ( solver.f2 );
                        bounds.add ( Nadirs.back ( ).first - 0.001 );
                        bounds.add ( Nadirs.back ( ).second - 0.001 );
                        dirs.add ( IloCplex::BranchDown );
                        dirs.add ( IloCplex::BranchDown );
                        makeBranch ( vars, bounds, dirs, getObjValue ( ) );
                        ++solver.NumberOfGPB;
                        ++solver.NumberOfParetoBranches;
                    }
                }
            }
            else
            {
                // If we get to this points, it means we could not perform GPB. Therefore, we perform ordinary PB
                // First branch, branches to the north west of the first dominator
                bounds.add ( Dominators.begin ( )->first - 1.0 );
                bounds.add ( Dominators.begin ( )->second + 1.0 );
                dirs.add ( IloCplex::BranchDown );
                dirs.add ( IloCplex::BranchUp );
                makeBranch ( vars , bounds , dirs , LB );
                ++solver.NumberOfParetoBranches;
                ++solver.NumberOfOrdPB;
                // Second branch branches to the south east of the first dominator
                bounds.clear ( ); dirs.clear ( );
                bounds.add ( Dominators.back ( ).first + 1.0 );
                bounds.add ( Dominators.back ( ).second - 1.0 );
                dirs.add ( IloCplex::BranchUp );
                dirs.add ( IloCplex::BranchDown );
                makeBranch ( vars , bounds , dirs , LB );
                ++solver.NumberOfParetoBranches;
                ++solver.NumberOfOrdPB;
            }
            vars.end ( );
            bounds.end ( );
            dirs.end ( );
            return;
        }
        else // We cannot perform Pareto branching, so we just branch on varibles instead
        {
            // We start by finding out how many branches cplex whats to create!
            IloInt NrBranches = getNbranches ( );
            // If cplex want to create 1 and two nodes, then he decides
            if ( NrBranches == 2 )
            {
                // Get the first branch
                getBranch ( vars , bounds , dirs , 0 );
                makeBranch( vars , bounds , dirs , LB );
                solver.FractionalNodesCplexDecides ++;

                // Clear arrays
                vars.clear ( ); bounds.clear ( ); dirs.clear ( );

                // Get the second branch!
                getBranch( vars , bounds , dirs , 1 );
                makeBranch ( vars , bounds , dirs , LB );
                solver.FractionalNodesCplexDecides ++;

                // Free memory
                vars.end();
                bounds.end();
                dirs.end();
                return;

            }
            else if ( NrBranches == 1 )
            {
                // Get the only branch cplex wants to create
                getBranch( vars , bounds , dirs , 0 );
                makeBranch( vars , bounds , dirs , getObjValue() );
                solver.FractionalNodesCplexDecides ++;

                // Free memory
                vars.end ( );
                bounds.end ( );
                dirs.end ( );
                return;
            }
            else
            { // If cplex does not know what to do, or if cplex wants to create more than two nodes, then we decide instead
                bool    FractionalVariableFound = false;
                int     bestYi = -1,
                        bestXi = -1,
                        bestXj = -1;
                IloNum  BestXscore = 0.0,
                        BestYscore = 0.0,
                        Yscore = 0.0,
                        Xscore = 0.0,
                        Yvalue = 0.0,
                        Xvalue = 0.0,
                        psUp = 0.0,
                        psDown = 0.0;
                IloExpr NoGood ( getEnv ( ) );
                for ( int i=0; i<n; ++i )
                {
                    // If the y-variable is fractional, get the pseudocost score of that variable
                    Yvalue = getValue ( solver.y[i] );
                    if ( ( Yvalue > 0.001 ) && ( Yvalue < 0.999 ) )
                    {
                        FractionalVariableFound = true;
                        psDown = getDownPseudoCost( solver.y[i] );
                        psUp = getUpPseudoCost( solver.y[i] );
                        Yscore = 0.83*IloMin ( psDown , psUp ) + 0.17*IloMax( psDown , psUp);
                        if ( Yscore > BestYscore )
                        {
                            BestYscore = Yscore;
                            bestYi = i;
                        }
                    }
                    for ( int j=0; j<m; ++j )
                    {
                        // If the x variable is fractional, get the pseudocost score of that variable and compare with the best choice
                        Xvalue = getValue ( solver.x[i][j] );
                        if ( ( Xvalue >= 0.001 ) && ( Xvalue <= 0.999 ) )
                        {
                            FractionalVariableFound = true;
                            psUp = getUpPseudoCost( solver.x[i][j] );
                            psDown = getDownPseudoCost( solver.x[i][j] );
                            Xscore = 0.83*IloMin ( psDown , psUp ) + 0.17*IloMax( psDown , psUp);
                            if ( Xscore > BestXscore )
                            {
                                BestXscore = Xscore;
                                bestXi = i;
                                bestXj = j;
                            }
                        }
                    }
                }

                // If we found a fractional variable, branch on that guy!
                if ( FractionalVariableFound )
                {
                    // If the best X score is better than the best Y score, we should branch on the x-variable
                    if ( BestXscore > BestYscore )
                    {
                        // Start by making the to-zero branch
                        vars.add ( solver.x[bestXi][bestXj] );
                        bounds.add ( 0.0 );
                        dirs.add ( IloCplex::BranchDown );
                        makeBranch( vars, bounds , dirs , LB );
                        solver.FractioanlNodesIDecide ++;
                        // Note that on the to-one branch we can also change lower bound of y[bestXi] to 1
                        //vars.add ( solver.y[bestXi] );
                        bounds.clear ( ); dirs.clear ( );
                        bounds.add( 1.0 );
                        //bounds.add( 1.0 );
                        dirs.add ( IloCplex::BranchUp );
                        //dirs.add ( IloCplex::BranchUp );
                        makeBranch( vars , bounds , dirs , LB );
                        solver.FractioanlNodesIDecide ++;
                        vars.end ( );
                        bounds.end ( );
                        dirs.end ( );
                        return;
                    }
                    else
                    {
                        // Then we branch on a y-variable. Start with the "to-zero" branch
                        vars.add ( solver.y[bestYi] );
                        bounds.add ( 0.0 );
                        dirs.add ( IloCplex::BranchDown );
                        makeBranch ( vars , bounds , dirs , LB );
                        solver.FractioanlNodesIDecide ++;
                        // Then go to the "to-one" branch
                        bounds.clear ( ); dirs.clear ( );
                        bounds.add ( 1.0 );
                        dirs.add ( IloCplex::BranchUp );
                        makeBranch ( vars , bounds , dirs , LB );
                        solver.FractioanlNodesIDecide ++;
                        vars.end ( );
                        bounds.end ( );
                        dirs.end ( );
                        return;
                    }
                }
                else
                {
                    // If we found no fractional variable, it is an integer feasible node, and we make a simple pareto branch around (f1,f2)
                    vars.add ( solver.f1 );
                    vars.add ( solver.f2 );
                    // First branch goes north west of (f1,f2)
                    bounds.add ( int ( f1 - 1.0 + 0.5 ) );
                    bounds.add ( int ( f2 + 1.0 +.5 ) );
                    dirs.add ( IloCplex::BranchDown );
                    dirs.add ( IloCplex::BranchUp );
                    makeBranch ( vars , bounds , dirs , LB );
                    solver.FractioanlNodesIDecide ++;
                    // Second branch goes south east of (f1,f2)
                    bounds.clear ( ); dirs.clear ( );
                    bounds.add ( int ( f1 + 1.0 + 0.5 ) );
                    bounds.add ( int ( f2 - 1.0 + 0.5 ) );
                    dirs.add ( IloCplex::BranchUp );
                    dirs.add ( IloCplex::BranchDown );
                    makeBranch ( vars , bounds , dirs , LB );
                    solver.FractioanlNodesIDecide ++;
                    ++solver.NumberOfParetoBranches;
                }

            }
        }
    }
    catch ( const IloException &ie )
    {
        std::cerr << "IloException in ImplicitBrancher callback : " << ie.getMessage ( ) << "\n";
        exit ( EXIT_FAILURE );
    }
    catch ( const std::exception &e )
    {
        std::cerr << "Exception in the ImplicitBrancher callback : " << e.what ( ) << "\n";
        exit ( EXIT_FAILURE );
    }
    catch ( ... )
    {
        std::cerr << "Some error I do not know what is occured! Terminating! \n" ;
        exit ( EXIT_FAILURE );
    }
}




//---------------------------------------------------------------------------------------------------------------------
/*ILOBRANCHCALLBACK1(brancher, BOCOsolver &, solver){
    try{
        const double myOne = 0.9999, myZero = 1.0-myOne;
        bool doPNPOLY = true, canBePruned = true, PruneMe = false;
        IloNumArray xvals(getEnv()), yvals(getEnv());
        IloNum f1 = getValue( solver.f1 );
        IloNum f2 = getValue( solver.f2 );
        std::vector<std::pair<double,double>> Dominators; // Vector used in isItDominated
        std::vector<std::pair<double,double>> DominatedNadirs; // vector used in canItBePruned

        // Start by checking if the node can be prune based on linear lower bound set
        if(solver.canItBePrunedByLinearLowerBoundSet(std::pair<double,double>(f1,f2), DominatedNadirs) ){
            solver.LinearPruned++;
            prune();
            return;
        }
        // Ok, it could not be pruned by the weak lower bound set. Then we test if one of the dominated Nadir
        // points are dominated by the outcome vector (f1,f2) as this will indicate that there is a Nadir point
        // in the polytope and the pnpoly should not be called
        for ( auto it = DominatedNadirs.begin() ; it!=DominatedNadirs.end() ; ++it ) if ( (f1 < it->first) && (f2 < it->second) ) doPNPOLY = false;
        if( doPNPOLY ){
            // Get the bounds of the variables and tell the NISEsolver which ones to fix
            std::vector<std::vector<int> > FixX;
            std::vector<int> FixY;
            for ( int i=0; i<solver.DATA->getNumFaci() ; ++i){
                // First we check if the upper bound is zero and the if the lower bound is one
                if ( getUB ( solver.y[i] ) <= myZero ) FixY.push_back( 0 );     // Fixed to o
                else if ( getLB ( solver.y[i] ) >= myOne ) FixY.push_back( 1 ); // Fixed to 1
                else FixY.push_back( 2 );                                       // Not fixed
                // Then we check the x-variables
                FixX.push_back( std::vector<int>() );
                for ( int j=0 ; j<solver.DATA->getNumCust() ; ++j){

                    if ( getUB ( solver.x[i][j] ) <= myZero ) FixX[i].push_back( 0 );
                    else if ( getLB ( solver.x[i][j] ) >= myOne ) FixX[i].push_back( 1 );
                    else FixX[i].push_back ( 2 );
                }
            }

            // Call the function fixing variables
            solver.FixNICEVariables( FixY, FixX );
            solver.setNISEbounds( getLB(solver.f1) , getLB(solver.f2) , getUB(solver.f1) , getUB(solver.f2) );

            std::vector<std::pair<double,double>> poly = solver.RunNISE(false , PruneMe);

            for( auto it=poly.begin()+1 ; it!=poly.end(); ){
                if( std::prev(it)->first <= it->first && std::prev(it)->second <= it->second ) poly.erase(it);
                else ++it;
            }
            // Run through the vector of Nadir points, and se if any is dominated by just one of the points on poly:
            for ( auto NadirIt = DominatedNadirs.begin() ; NadirIt!= DominatedNadirs.end() && canBePruned; ++NadirIt ){
                for ( auto polyIt = poly.begin(); polyIt!=poly.end() && canBePruned; ++polyIt ){
                    if( (polyIt->first <= NadirIt->first-0.99) && (polyIt->second  <= NadirIt->second-0.99) ){
                        canBePruned = false;
                    }
                }
            }
            if( canBePruned ){
                poly.push_back( std::pair<double,double>( solver.getZlr().first , getLB( solver.f2 ) ) );
                poly.push_back( std::pair<double,double>( solver.getZlr().first , solver.getZul().second ) );
                poly.push_back( std::pair<double,double>( poly.begin()->first , solver.getZul().second ) );

                // Foreach of the dominated nadir points, check if it is inside the polygon:
                for( auto it = DominatedNadirs.begin() ; (canBePruned && ( it!=DominatedNadirs.end() ) )  ;  ++it){
                    std::pair<double,double> p( it->first -0.99 , it->second -0.99);
                    // Test if the current dominated Nadir point lies in the Pareto box
                    //if ( ( p.first <= getUB(solver.f1) ) && ( p.first>= getLB(solver.f1) ) && ( p.second<=getUB(solver.f2) ) && ( p.second >=getLB(solver.f2) ) ){
                    if ( pnpoly( poly, p ) ){
                        canBePruned = false;
                    }
                }

                if(canBePruned){

                    solver.PNPOLYpruned++;
                    //std::cout << "Pruning based on PNPOLY! Pruning time : " << (double(clock()/double(CLOCKS_PER_SEC))  -start) << std::endl;
                    if ( PrintFrontiers ){
                        myfile.open( "LPfront.txt" );
                        for ( auto it = poly.begin(); it!=poly.end(); ++it ) myfile << it->first << "\t" << it->second << std::endl;
                        myfile.close();
                        myfile.open ( "UBset.txt" );

                        size_t UBsetSize = solver.getUBsetSize();

                        for ( auto it = solver.nonDomSet->UBsetBegin(); it!=solver.nonDomSet->UBsetEnd(); ++it )
                        {
                            myfile  << it->first << "\t" << it->second << std::endl;
                        }
                        myfile.close();
                        PrintFrontiers = false;
                    }
                    prune();
                    return;
                }else{
                    //std::cout << "Pruning using PNPOLY failed!! Pruning time : " << (double(clock()/double(CLOCKS_PER_SEC))  -start) << std::endl;
                    solver.PNPOLYnotPruned++;
                }
            }
        }


        // Okay, we could not prune the node, instead we retrieve the nodedata of the current node
        NodeData * ND = getNodeData();
        if ( nullptr == ND ){
        //    throw std::runtime_error ( "There appear to be no NodeData object attached to the current node. Procedure Terminating!\n");
        }
        // Cast the NodeData to myNodeData
        myNodeData *mND = (myNodeData*) ND;
        if ( mND->iShouldSolveBOLP() ); // Call the NISEsolver and update the poly object in

        // Check if we can perform Pareto branching. Remember, if the LP solution of the current node is integral
        // it is already in the set UBset as incumbet callback is called before branch callback

        if ( solver.isItDominated(Dominators , std::pair<double,double>(f1,f2)) ){
            // There was points which dominates the current outcome vector, so we perform Pareto branching.
            if(Dominators.size()>8) std::cout << Dominators.size() <<"-Pareto braching is performed\n";
            makeBranch(solver.f1 , (Dominators.begin()->first-1.0), IloCplex::BranchDirection::BranchDown, (f1+f2) );
            makeBranch(solver.f2 , (Dominators.back().second-1.0), IloCplex::BranchDirection::BranchDown, (f1+f2) );
            return;
        }else if ( DominatedNadirs.size() <0 ){
            // If this is the case, we can perform a generalized pareto braching stating that



        }else{
            // We cannot perform Pareto branching and we have to do something else instead
            // Just let CPLEX do the dirty work

            return;
        }

    }catch(std::runtime_error &re){
        std::cerr << "Runtime error in the branch callback brancher : " << re.what() << std::endl;
        exit ( EXIT_FAILURE );
    }catch(std::exception &e){
        std::cerr << "Exception in the branch callback brancher : " << e.what() << std::endl;
        exit ( EXIT_FAILURE );
    }catch(IloException &ie){
        std::cerr << "IloException in the branch callback brancher : " << ie.getMessage() << std::endl;
        exit ( EXIT_FAILURE );
    }
}*/


//-------------------------------------------------------------------------------------------------------------
ILOINCUMBENTCALLBACK1(rejecter, BOCOsolver &, solver){
    try{
        IloNumArray yvals(getEnv());
        IloNumArray xvals(getEnv());

        IloNum f1 , f2;
        f1 = getValue(solver.f1);
        f2 = getValue(solver.f2);

        solver.addToFront( std::pair<double,double>(f1, f2 ) );
        ++solver.SolutionsFound;
        reject();
    }catch(std::exception &e){
        std::cerr << "Exception in incumbent callback rejecter : " << e.what() << std::endl;
        exit(1);
    }catch(IloException &ie){
        std::cerr <<  "IloException in incumbent callback rejecter : " << ie.getMessage() << std::endl;
        exit(1);
    }
}
/*------------------------------------------------------------------------------------------------------
    Here begins the implementation of the BOCOsolver class defined in BOCOsolver.h
    The callback implemented above is utalized by the IloCplex instance cplex.
------------------------------------------------------------------------------------------------------*/

BOCOsolver::BOCOsolver(const data &TheDataObject):cutlst(nullptr),TheWorstLocalNadirPoint(double(INT_MAX)){
    try{
        DATA = new data(TheDataObject);
        nonDomSet = new NonDominatedSet();

        LinearPruned                = 0;
        PNPOLYpruned                = 0;
        PNPOLYnotPruned             = 0;
        LBupdates                   = 0;
        BOLPsolved                  = 0;
        SolutionsFound              = 0;
        NumNodes                    = 0;
        MaxNodesInTree              = 0;
        NumberOfParetoBranches      = 0;
        NumberOfOrdPB               = 0;
        NumberOfGPB                 = 0;
        FractionalNodesCplexDecides = 0;
        FractioanlNodesIDecide      = 0;
        NumCuts                     = 0;
        CheckAllNadirs              = false;
    }catch(std::exception &e){
        std::cerr << "Error in the constructor of the BOCOsolver class : "  << e.what() << std::endl;
        exit(1);
    }
}

BOCOsolver::~BOCOsolver(){
    if ( f1.getImpl ( ) ) f1.end ( );
    if ( f2.getImpl ( ) ) f2.end ( );
    if ( y.getImpl ( ) ) y.end ( );
    if ( x.getImpl () )
    {
        if ( DATA != nullptr )
        {
            for ( int i=0; i<DATA->getNumFaci(); ++i )
            {
                if ( x[i].getImpl ( ) ) x[i].end  ( );
            }
        }
    }
    if ( s.getImpl ( ) ) s.end ( );
    if ( PrunerCplex.getImpl ( ) ) PrunerCplex.end ( );
    if ( PrunerModel.getImpl ( ) ) PrunerModel.end ( );
    if ( cplex.getImpl ( ) ) cplex.end ( );
    if ( model.getImpl ( ) ) model.end ( );
    if ( env.getImpl ( ) ) env.end ( );
}

/*------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------*/
void BOCOsolver::buildSSCFLP(){
    model = IloModel(env);
    cplex = IloCplex(model);

    PrunerModel = IloModel ( env );
    PrunerCplex = IloCplex ( PrunerModel );

    f1 = IloNumVar(env);
    f2 = IloNumVar(env);
    y = IloNumVarArray(env, DATA->getNumFaci ( ) , 0 , 1 , ILOBOOL);
    x = IloVarMatrix ( env , DATA->getNumCust ( ) );

    s = IloNumVarArray ( env , 2 , 0 , IloInfinity , ILOFLOAT ); // The slack variables used in the implicitPruner strategy
    PrX = IloVarMatrix ( env , DATA->getNumFaci ( ) );
    PrY = IloNumVarArray ( env , DATA->getNumFaci ( ) , 0 , 1 , ILOFLOAT );
    Prf1 = IloNumVar ( env , 0 , IloInfinity , ILOFLOAT);
    Prf2 = IloNumVar ( env , 0 , IloInfinity , ILOFLOAT);

    OBJ = IloMinimize(env, f1+f2);
    model.add(OBJ);
    PrunerModel.add ( IloMinimize ( env , s[0] + s[1] ) );

    IloExpr cst(env), obj1(env), obj2(env), PrCst(env), PrObj1(env), PrObj2(env);

    for( int i=0 ; i < DATA->getNumFaci() ; ++i ){
        x[i] = IloNumVarArray( env , DATA->getNumCust() , 0 , 1 , ILOBOOL);
        PrX[i] = IloNumVarArray ( env , DATA->getNumCust() , 0 , 1 , ILOFLOAT );

        for( int j=0 ; j < DATA->getNumCust() ; ++j ){
            PrCst  += DATA->getDemand(j)*PrX[i][j];
            PrObj1 += DATA->getC(i,j)*PrX[i][j];
            PrunerModel.add ( PrX[i][j] - PrY[i] <= 0 );


            cst+=DATA->getDemand(j)*x[i][j];
            obj1+=DATA->getC(i,j)*x[i][j];
            model.add(x[i][j]-y[i]<=0);

        }
        model.add( cst - DATA->getCapacity(i)*y[i] <= 0);
        obj2+=DATA->getFixedCost(i)*y[i];
        cst.clear();
        cplex.setPriority( y[i] , 1 );

        PrunerModel.add ( PrCst - DATA->getCapacity(i)*PrY[i] <= 0);
        PrObj2 += DATA->getFixedCost(i)*PrY[i];
        PrCst.clear ( );

    }

    gcdOfC2 = int ( DATA->getFixedCost(0) );
    int g = 1;
    while ( gcdOfC2 > 1 && g < DATA->getNumFaci() )
    {
        gcdOfC2 = gcd( gcdOfC2 , DATA->getFixedCost ( g++) );
    }

    model.add(f1 == obj1);
    model.add(f2 == obj2);

    PrunerModel.add ( Prf1 == PrObj1 );
    PrunerModel.add ( Prf2 == PrObj2 );

    FirstCoord = ( PrObj1 - s[0] <= IloInfinity );
    SecondCoord = ( PrObj2 - s[1] <= IloInfinity );
    PrunerModel.add ( FirstCoord );
    PrunerModel.add ( SecondCoord );

    for( int j=0 ; j < DATA->getNumCust() ; ++j ) {
        for( int i=0; i < DATA->getNumFaci() ; ++i)
        {
            cst += x[i][j];
            PrCst += PrX[i][j];
        }
        model.add(cst == 1);
        cst.clear();

        PrunerModel.add ( PrCst == 1 );
        PrCst.clear ( );
    }


    std::cout << "solving pruner Cplex:\n";
    PrunerCplex.setOut( env.getNullStream ( ) );
    PrunerCplex.setWarning( env.getNullStream ( ) );
    PrunerCplex.solve ( );

    std::cout << "After solving PrunerCplex in buildSSCFLP\n";

    cplex.setParam ( IloCplex::ClockType , 1 );
    cplex.setParam ( IloCplex::Param::TimeLimit , 3700 );

    PrCst.end ( );
    cst.end();
    PrObj1.end ( );
    PrObj2.end ( );
    obj1.end ( );
    obj2.end ( );
    setupANiseSSCFLP();
}

/*------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------*/
void BOCOsolver::updateUBset(std::pair<double,double> p){
    try{
        nonDomSet->addToFrontier( p );
    }catch(std::runtime_error &re){
        std::cerr << "Runtime error in updateUBSet in the BOCOsolver class : " << re.what() << std::endl;
        exit(1);
    }catch(std::exception &e){
        std::cerr << "Exception in updateUBset in the BOCOsolver class : " << e.what() << std::endl;
        exit(1);
    }catch(IloException &ie){
        std::cerr << "IloException in updateUBset in the BOCOsolver class : " << ie.getMessage() << std::endl;
        exit(1);
    }
}

/*------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------*/

bool BOCOsolver::isItDominated(std::vector<std::pair<double,double>> & Dominators, std::pair<double,double> p){
    try{
        Dominators.clear();
        for ( auto it = nonDomSet->UBsetBegin(); it!=nonDomSet->UBsetEnd(); ++it ){
            if ( (it->first <= p.first) && (it->second <= p.second) ){
                Dominators.push_back(*it);
            }
        }

        return (Dominators.size()>=1);
    }catch(std::runtime_error &re){
        std::cerr << "Runtime error in isItDominated in the BOCOsolver class : " << re.what() << std::endl;
        exit(1);
    }catch(std::exception &e){
        std::cerr << "Exception in isItDominated in the BOCOsolver class : " << e.what() << std::endl;
        exit(1);
    }
}

/*------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------*/
bool BOCOsolver::canItBePrunedByLinearLowerBoundSet(std::pair<double,double> p, std::vector<std::pair<double,double>> &DominatedNadirs){
    try{
        double ObjVal = p.first+p.second + 0.0001;
        DominatedNadirs.clear();
        for ( auto it=UBset.begin(); it!=UBset.end(); ++it ){
            if ( (std::next(it)->first) + (it->second) >= ObjVal ){
                DominatedNadirs.push_back(std::pair<double,double>(std::next(it)->first,it->second));
            }
        }
        return (UBset.size() == 0);
    }catch(std::runtime_error &re){
        std::cerr << "Runtime error in canItBePrunedByLinearLowerBoundSet in the BOCOsolver class : " << re.what() << std::endl;
        exit(1);
    }catch(std::exception &e){
        std::cerr << "Exception in canItBePrunedByLinearLowerBoundSet in the BOCOsolver class : " << e.what() << std::endl;
        exit(1);
    }
}
/*------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------*/
/* double BOCOsolver::getTheWorstLocalNadirPoint(){
    try{
        double nadirPoint, oldWorst;
        oldWorst = TheWorstLocalNadirPoint;
        TheWorstLocalNadirPoint = 0.0;
        for ( auto it=UBset.begin() ; std::next(it)!=UBset.end(); ++it ){
            nadirPoint = std::next(it)->first + it->second;
            if ( nadirPoint > TheWorstLocalNadirPoint ) TheWorstLocalNadirPoint = nadirPoint;
        }
        if ( TheWorstLocalNadirPoint < oldWorst ) std::cout << "New worst local nadir point : " << TheWorstLocalNadirPoint  << std::endl;
        return TheWorstLocalNadirPoint;
    }catch(std::exception &e){
        std::cerr << "Exception in getTheWorstLocalNadirPoint in the BOCOsolver class : " << e.what() << std::endl;
        exit ( EXIT_FAILURE );
    }
 }*/

/*------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------*/

void BOCOsolver::initializeUBset(){
    try{
        IloCplex initCplex = IloCplex ( model );
        initCplex.use ( LexiUpdater ( env,*this ) );
        initCplex.setParam ( IloCplex::LBHeur, 1 );
        initCplex.setOut ( env.getNullStream ( ) );
        initCplex.setWarning( env.getNullStream ( ) );
        std::pair<double,double> p(0,0);

        OBJ.setLinearCoef(f1,0.999999);
        OBJ.setLinearCoef(f2,0.000001);

        if( initCplex.solve() ){ // calculating z^{ul}
            p.first = static_cast<int>(initCplex.getValue(f1));
            p.second = static_cast<int>(initCplex.getValue(f2));
            std::cout << "Point ul : (" << p.first << "," << p.second <<")" << std::endl;
            //updateUBset(p);
            p.first = p.second = 0;

            // Calculating z^{lr}
            OBJ.setLinearCoef(f1,0.000001);
            OBJ.setLinearCoef(f2,0.999999);
            /*if(initCplex.solve()){
                p.first =  initCplex.getValue(f1);
                p.second = initCplex.getValue(f2);
                std::cout << "Point lr : (" << p.first << "," << p.second <<")" << std::endl;
                //updateUBset(p);

            }else throw std::runtime_error("Could not solve for zlr. Terminating!\n");*/
        }else throw std::runtime_error("Could not solve for zul. Terminating!\n");


        // Adding the bounding box of the non--dominated frontier
        f1.setBounds(UBset.begin()->first, UBset.back().first);
        f2.setBounds(UBset.back().second, UBset.begin()->second);
        // Deallocate memory occupied by the IloCplex instance
        initCplex.end();

        // Setting back the objective coefficients
        OBJ.setLinearCoef(f1,1);
        OBJ.setLinearCoef(f2,1);

    }catch(std::runtime_error &re){
        std::cerr << "Runtime error in initializeUBSet in the BOCOsolver class : " << re.what() << std::endl;
        exit(1);
    }catch(std::exception &e){
        std::cerr << "Exception in initializeUBSet in the BOCOsolver class : " << e.what() << std::endl;
        exit(1);
    }catch(IloException &ie){
        std::cerr << "IloException in initializeUBSet in the BOCOsolver class : " << ie.getMessage() << std::endl;
        exit(1);
    }
}
/*------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------*/
void BOCOsolver::initializeUBset2(){
    try{
        IloCplex initCplex = IloCplex ( model );
        initCplex.use ( LexiUpdater ( env,*this ) );
        initCplex.setParam ( IloCplex::LBHeur, 1 );
        initCplex.setParam ( IloCplex::EpAGap , 0.99 );
        initCplex.setParam ( IloCplex::EpGap  , 0.0 );
        initCplex.setOut ( env.getNullStream ( ) );
        initCplex.setWarning( env.getNullStream ( ) );

        OBJ.setLinearCoef ( f1 , 1.0 );
        OBJ.setLinearCoef ( f2 , 0.0 );
        if ( initCplex.solve( ) ){
            f1.setUB ( initCplex.getObjValue( ) );
            OBJ.setLinearCoef( f1 , 0.0 );
            OBJ.setLinearCoef( f2 , 1.0 );
            if ( !initCplex.solve( ) ){
                throw std::runtime_error ( "Could not solve the lex-min problem when f1 is fixed\n" );
            }
        }else throw std::runtime_error ( "Could not solve the lex-min problem when minimizing f1\n" );

        std::cout << "Halfways through the initialize2\n";

        f1.setUB( IloInfinity );

        if ( initCplex.solve( ) ){
            f2.setUB( initCplex.getObjValue( )+0.9 );
            OBJ.setLinearCoef ( f1 , 1.0 );
            OBJ.setLinearCoef ( f2 , 0.0 );
            if ( ! initCplex.solve( ) ){
                throw std::runtime_error ( "Could not solve the lex-min problem when f2 is fixed\n" );
            }
        }else throw std::runtime_error ( "Could not solve the lex-min problem when minimizing f2\n" );
        std::cout << "At the end of initialize 2\n";
        initCplex.end();
        std::cout << "After endning the initCplex\n";
        OBJ.setLinearCoef( f1 , 1.0 );
        OBJ.setLinearCoef( f2 , 1.0 );
        std::cout << "";
        f1.setBounds ( nonDomSet->UBsetBegin()->first ,  std::prev( nonDomSet->UBsetEnd() ) -> first );
        f2.setBounds ( std::prev( nonDomSet->UBsetEnd() ) -> second , nonDomSet->UBsetBegin()->second );
        std::cout << "Done with the initialization\n";
    }catch ( std::exception &e ){
        std::cerr <<"Exception in initializeUBset2 in the BOCOsolver class : " << e.what ( ) << std::endl;
        exit( EXIT_FAILURE );
    }catch ( IloException & ie){
        std::cerr << "IloException in initializeUBset2 in the BOCOsolver class : " << ie.getMessage() << std::endl;
        exit ( EXIT_FAILURE );
    }
}

/*------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------*/
std::vector< std::pair<double,double>> BOCOsolver::getDominatedPoints( double LB ){
    try{
        std::vector< std::pair<double,double> > DominatedPoints;
        for ( auto it=UBset.begin(); it!=UBset.end(); ++it ) if ( it->first + it->second > LB ) DominatedPoints.push_back(*it);
        return DominatedPoints;
    }catch ( std::exception &e ){
        std::cerr << "Exception in getDominatedPoints in the BOCOsolver class : " << e.what ( ) << std::endl;
        exit ( EXIT_FAILURE );
    }
}

/*------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------*/
std::pair<double,double> BOCOsolver::getIntersectionPoint(std::pair<double,double> z1, std::pair<double,double> z2, std::pair<double,double> zstar, double slope){
    try{
        std::pair<double,double> intersectionPoint;
        double  a = 0.0,    // the slope of the line passing through the two points z1 and z2
                c = 0.0,    // Intersection with y-axis of line going through z1 and z2
                b = 0.0,    // Slope of the line going through the point zstar with slope "slope"
                d = 0.0;    // Intersection with the y-axis of the line passing through zstar with slope "slope"
        const double tol = 1e-6;


        /*
         * First we calculate the parameters of the line going through the two points z1 an z2.
         * That line is denoted y=ax+c
         */
        // Calculate the slope, a, of the line passing through points z1 and z2:
        a = (z1.second - z2.second) / (z1.first - z2.first);
        // Then we calculate the y-axis intersection, c
        c = z1.second - a * z1.first;

        /*
         * Now calculate the parameters of the line through zstar with slope "slope"
         * That line is defined as y=bx+d
         */
        // The slope is easy as that one is given
        b = slope;
        // The intersection with the y-axis is then
        d = zstar.second - b * zstar.first;

        // If a is very close to the value of slope, then the two points z1 and z2 are very close to be on the line defined by the LP
        // Therefore, we might end up with numerical instability, so instead we just return the original point
        if ( ( a <= slope + tol )  && ( a >= slope + tol ) )
        {
            intersectionPoint.first = z1.first;
            intersectionPoint.second = z1.second;
        }
        else
        {
            // If the slope of the line between z1 and z2 is not close to slope, then we calculate the intersection points
            // we use the formulas x = (d-c) / (a-b), and y = a (d-c)/(a-b)+c
            intersectionPoint.first = (d-c)/(a-b);
            intersectionPoint.second = a * ( ( d - c ) / ( a - b ) ) + c;
        }
        return intersectionPoint;
    }catch(int IntException){
        std::cerr << "Two lines are parallel. Returns (-1,-1). The lower bound set is updated from scratch.\n";
        return std::pair<double,double>(-1.0,-1.0);
    }catch( std::exception &e){
        std::cerr << "Exception in getIntersectionPoint in the solver class : " << e.what() << std::endl;
        exit ( EXIT_FAILURE );
    }
}

/*============================================================================================================*/

 std::vector< std::pair<double,double> > BOCOsolver::updateLBset(   const std::vector< std::pair<double,double> > &OldLB, bool doBOLP, double LPvalue, std::pair<double,double> LPoutcome,
                                                                    const std::vector< std::vector<int> > &FixX, const std::vector<int> &FixY, std::pair<double,double> Lbounds,
                                                                    std::pair<double,double> Ubounds, bool &PruneTheBugger){
    try{
        bool    firstPointBelow = false,    // Flag indicating if the left most point is below the new facet. Assume no!
                lastPointBelow  = false;    // Flag indicating if the right most point is below the new facet. Assume no!
        std::vector< std::pair<double,double> > NewLB; // The vector containing the updated LB set
        if ( !doBOLP )
        { // Only do update if we should not solve BOLP from scratch

            if ( firstPointBelow || lastPointBelow )
            {
                // If either the first or the last is below, we resolve
                doBOLP = true;
            }
            else
            {
                // Both are above. Then we need to loop through the OldLB and check where the line intersects
                // We know that OldLB.begin() is above, so we add that guy to the NewLB vector
                ++LBupdates;
                NewLB.push_back( *OldLB.begin() );
                for ( auto it=std::next ( OldLB.begin() ); it!=OldLB.end(); ++it ){
                    if ( (std::prev(it)->first + std::prev(it)->second >= LPvalue) && (it->first + it->second < LPvalue ) )
                    {
                        // If this point is below the LPvalue line and the previous was above then we have found an intersection point
                        NewLB.push_back( getIntersectionPoint(*std::prev(it) , *it , LPoutcome , -1 ) );

                    }else if ( (it->first + it->second >= LPvalue) && (std::prev(it)->first +std::prev(it)->second < LPvalue) )
                    {
                        // If this point is above and previous was below, then we have found the second intersection point
                        NewLB.push_back( getIntersectionPoint(*std::prev(it) , *it , LPoutcome , -1 ) );

                    }
                    // If the point lies above the LPvalue-defined line, add it to the NewLB vector
                    if ( it->first + it->second >= LPvalue ) NewLB.push_back( *it );
                }
                return NewLB;
            }


            /*if( OldLB.begin()->first + OldLB.begin()->second < LPvalue ) firstPointBelow = true;
            if( OldLB.back().first + OldLB.back().second < LPvalue ) lastPointBelow = true;

            if ( firstPointBelow == lastPointBelow ){ // If both above or both below
                if ( firstPointBelow ){ // Now both are below. This means that all other point on OldLB is below as well. Therefore, the LBset should contain two points only
                    if ( ResolveBOLPIfAllBelow ){ // If we have decided to resolve the BOLP if all points are below the linear lower bound, then do so
                        doBOLP = true;
                    }else{
                        std::pair<double,double> p ( LPoutcome.first , -(OldLB.begin()->first - LPoutcome.first ) + LPoutcome.second );
                        NewLB.push_back( p );
                        p.first = -OldLB.back( ).second + LPoutcome.first + LPoutcome.second;
                        p.second= OldLB.back( ).second;
                        NewLB.push_back( p );
                        return NewLB;
                    }
                }else{ // Both are above. Then we need to loop through the OldLB and check where the line intersects
                    // We know that OldLB.begin() is above, so we add that guy to the NewLB vector
                    NewLB.push_back( *OldLB.begin() );
                    for ( auto it=++OldLB.begin(); it!=OldLB.end(); ++it ){
                        if ( (std::prev(it)->first + std::prev(it)->second >= LPvalue) && (it->first + it->second < LPvalue ) ){
                            // If this point is below the LPvalue line and the previous was above then we have found an intersection point
                            NewLB.push_back( getIntersectionPoint(*std::prev(it) , *it , LPoutcome , -1 ) );
                        }else if ( (it->first + it->second >= LPvalue) && (std::prev(it)->first +std::prev(it)->second < LPvalue) ){
                            // If this point is above and previous was below, then we have found the second intersection point
                            NewLB.push_back( getIntersectionPoint(*std::prev(it) , *it , LPoutcome , -1 ) );
                        }
                        // If the point lies above the LPvalue-defined line, add it to the NewLB vector
                        if ( it->first + it->second >= LPvalue ) NewLB.push_back( *it );
                    }
                    return NewLB;
                }
            }else{ // firstPointBelow OR(exclusive) secondPoitBelow
                if ( firstPointBelow ){ // First below but second above
                    // We start by inserting the intersection point between the vertival line straight above the first point and the LPvalue defined line
                    std::pair<double,double> p ( LPoutcome.first , -(OldLB.begin()->first - LPoutcome.first ) + LPoutcome.second );
                    NewLB.push_back( p );
                    // Then we loop through the list in order to find the last intersection point
                    for ( auto it=OldLB.begin(); it!=OldLB.end(); ++ it){
                        if ( (it->first + it->second >= LPvalue) && (std::prev(it)->first +std::prev(it)->second < LPvalue) ){
                            // If this point is above and previous was below, then we have found the second intersection point
                            NewLB.push_back( getIntersectionPoint(*std::prev(it) , *it , LPoutcome , -1 ) );
                        }
                        // If the point lies above the LPvalue-defined line, add it to the NewLB vector
                        if ( it->first + it->second >= LPvalue ) NewLB.push_back( *it );
                    }
                    return NewLB;
                }else{ // Second below but first above
                    for ( auto it=++OldLB.begin(); it!=OldLB.end(); ++it ){
                        if ( (std::prev(it)->first + std::prev(it)->second >= LPvalue) && (it->first + it->second < LPvalue ) ){
                            // If this point is below the LPvalue line and the previous was above then we have found an intersection point
                            NewLB.push_back( getIntersectionPoint(*std::prev(it) , *it , LPoutcome , -1 ) );
                        }
                        // If the point lies above the LPvalue-defined line, add it to the NewLB vector
                        if ( it->first + it->second >= LPvalue ) NewLB.push_back( *it );
                    }
                    // Add the intersection point between the horizontal line to thr right of OldLB.back()
                    NewLB.push_back( std::pair<double,double>( -OldLB.back().second+ LPoutcome.first + LPoutcome.second , OldLB.back().second ) );
                    return NewLB;
                }
            }*/
        }
        if ( doBOLP )
        { // Need new if-scope (instead of the obvious else-scope) as doBOLP might get modified in the if-scope above!
            ++BOLPsolved;
            FixNICEVariables( FixY, FixX );
            setNISEbounds( Lbounds.first , Lbounds.second, Ubounds.first , Ubounds.second );
            NewLB = RunNISE ( false , PruneTheBugger ); // Run the NISE solver without adding cuts!
        }
        return NewLB;
    }catch( std::exception &e){
        std::cout << "Exception in updateLBset in the BOCOsolver class : "  << e.what ( ) << std::endl;
        exit ( EXIT_FAILURE );
    }
}

/*------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------*/
void BOCOsolver::setupANiseSSCFLP(){
    try{
        // Hand the data of the problem to the NISEsolver and let it build a model.
        BiObjSolver.BuildSSCFLP ( DATA->getNumFaci() , DATA->getNumCust() ,
                    DATA->getDemands() , DATA->getCapacitites() , DATA->getFixedCosts() , DATA->getAssignmentCosts() );

    }catch ( std::exception &e ){
        std::cerr << "Exception in setupANiseSSCFLP in the BOCOsolver class : " << e.what() << std::endl;
        exit ( 1 );
    }catch ( IloException &ie ){
        std::cerr << "IloException in the setupANiseSSCFLP in the BOCOsolver class : " << ie.getMessage() << std::endl;
        exit ( 1 );
    }
}

/*------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------*/
std::vector<std::pair<double,double>> BOCOsolver::RunNISE(bool doCuts , bool &PruneTheBugger){
    return BiObjSolver.RUN(doCuts, model , PruneTheBugger);
}

/*------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------*/
void BOCOsolver::FixNICEVariables ( std::vector<int> fixY , std::vector< std::vector<int> > fixX ){
    BiObjSolver.fixXVariables( fixX );
    BiObjSolver.fixYvariables( fixY );
}

/*------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------*/
void BOCOsolver::setNISEbounds(double LB1, double LB2, double UB1, double UB2){
    BiObjSolver.setObjFuncBounds(LB1,LB2,UB1,UB2);
}
/*------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------*/
void BOCOsolver::run(){

    buildSSCFLP();

    cplex.setParam( IloCplex::Reduce,0);
   // cplex.setParam ( IloCplex::CutsFactor , 0.5 );
    cplex.setParam( IloCplex::ParallelMode, 1);
    cplex.setParam ( IloCplex::NodeSel , 2 );
    cplex.setParam( IloCplex::CutsFactor, 1.0);
    cplex.setWarning ( env.getNullStream ( ) );

    std::cout << "===================== Before solving root BOLP  =====================\n";
    bool ProblemInfeasible = false;

    auto Start = CPUclock::now ( );
    BiObjSolver.RUN ( DoCutsAtRoot, model, ProblemInfeasible );
    auto End = CPUclock::now ( );

    CuttingPhaseTime = duration_cast<duration<double>>( End - Start ).count ( ) ;
    cutlst = BiObjSolver.getCuts();
    if ( cutlst != NULL ){
        std::cout << "Cuts are handed to the BOCOsolver class!\n";
        for ( VICcut* cut = cutlst; cut!=NULL; cut=cut->nextcut ){
            IloExpr cst(env), PrCst ( env );
            int Iindex = cut->UsrIVal;
            for(int nz=0; nz < cut->nzcnt; nz++){
                cst+= double(cut->nzval[nz])*x[Iindex][cut->nzind[nz]];
                PrCst+= double( cut->nzval[nz] )*PrX[Iindex][cut->nzind[nz]];
            }
            model.add ( cst<= double(cut->rhs)*y[Iindex] );
            PrunerModel.add ( PrCst <= double ( cut->rhs ) * PrY[Iindex] );
            ++NumCuts;
        }
    }
    else std::cout << "No cuts are handed to the BOCOsolver\n";
    BiObjSolver.printNcuts();
    std::cout << "======================== Before initializing ========================\n";
    initializeUBset2();
    std::cout << "========================  initialized UBset  ========================\n ";
    for(auto it = UBset.begin(); it!=UBset.end(); ++it){
        std::cout << it->first << "\t" << it->second << std::endl;
    }
    std::cout << "======================== After initializing =========================\n";

    if ( doBoundSetBasedBaranching )
    {
        cplex.use ( BOCObrancher(env,*this) );
    }
    else
    {
        cplex.use ( ImplicitBrancher(env, *this) );
    }
    cplex.use ( rejecter(env,*this) );
    cplex.use ( BoundCutter(env,*this) );
    cplex.use ( Terminator(env,*this) );
    if( cplex.solve() )  std:: cout << "Objective function value is " << cplex.getObjValue() << std::endl;


    std::cout << "Cplex status is : " << cplex.getStatus()  << std::endl;

    for(auto it = nonDomSet->UBsetBegin(); it!=nonDomSet->UBsetEnd(); ++it){
        std::cout << it->first << "\t" << it->second << std::endl;
    }
    NumNodes = cplex.getNnodes();
    std::cout << "Number of nodes examined                      : " << cplex.getNnodes() << std::endl;
    std::cout << "Number of nodes prund by linear lower bound   : " << LinearPruned << std::endl;
    std::cout << "Number of node pruned by PNPOLY               : " << PNPOLYpruned << std::endl;
    std::cout << "Number of nodes not pruned by PNPOLY          : " << PNPOLYnotPruned << std::endl;
    std::cout << "Number of LB updates                          : " << LBupdates   << std::endl;
    std::cout << "Number of BOLP solved                         : " << BOLPsolved  << std::endl;
    std::cout << "Number of Pareto branches                     : " << NumberOfParetoBranches << std::endl;
    std::cout << "Number of GPB                                 : " << NumberOfGPB   << std::endl;
    std::cout << "Number of ordinary PB                         : " << NumberOfOrdPB << std::endl;
    std::cout << "Number of nodes Cplex decide                  : " << FractionalNodesCplexDecides << std::endl;
    std::cout << "Number of nodes I decide                      : " << FractioanlNodesIDecide << std::endl;
}

void BOCOsolver::getUBset ( std::vector< std::pair< double,double> > & CopyUBset )
{
    try
    {
        CopyUBset.clear ( );
        for ( auto it = nonDomSet->UBsetBegin(); it!= nonDomSet->UBsetEnd(); ++it )
        {
            CopyUBset.push_back ( *it );
        }
    }
    catch ( std::exception &e )
    {
        std::cerr << "ERROR in getUbset in the BOCOsolver class : " << e.what ( )  << std::endl;
        exit ( EXIT_FAILURE );
    }
    catch ( ... )
    {
        std::cerr << "Unforsee ERROR in getUBset in the BOCOsolver class!\n";
        exit ( EXIT_FAILURE );
    }
}


// ---------------------------------------------------------------------------------------------//
bool BOCOsolver::implicitPruner (   bool CheckAllNadirs , std::vector<std::pair<double,double>> &Nadirs, const std::vector<std::vector<int>> &FixX, const std::vector<int> &FixY ,
                                    IloNum LB1,  IloNum LB2, IloNum UB1, IloNum UB2)
{
    try
    {
        int n = DATA->getNumFaci ( ),
            m = DATA->getNumCust ( );
        std::vector<std::pair<double,double>> DominatedNadirs;

        for ( int i=0; i<n; ++i )
        {   // Loop through all variables and set their bounds according to the vectors FixX and FixY
            if ( 0 == FixY[i] ) PrY[i].setBounds ( 0 , 0 );
            else if ( 1 == FixY[i] ) PrY[i].setBounds ( 1 , 1 );
            else PrY[i].setBounds ( 0 , 1);
            for ( int j=0; j<m; ++j )
            {
                if ( 0 == FixX[i][j] ) PrX[i][j].setBounds( 0 , 0 );
                else if ( 1 == FixX[i][j] ) PrX[i][j].setBounds ( 1 , 1 );
                else PrX[i][j].setBounds ( 0 , 1 );
            }
        }

        Prf1.setBounds ( LB1 , UB1 );
        Prf2.setBounds ( LB2 , UB2 );
        // Start the main loop
        for ( auto it = Nadirs.begin ( ); it!= Nadirs.end ( ); ++it )
        {
            // Start by changing the bounds on the constraints corresponding to the coordinates of the nadir points
            FirstCoord.setBounds ( 0 , it->first );
            SecondCoord.setBounds ( 0 , it->second );
            // Solve the proble using the PrunerCplex object
            PrunerCplex.setOut( env.getNullStream ( ) );
            if ( PrunerCplex.solve ( ) )
            {
                IloNum s1 = PrunerCplex.getValue ( s[0] );
                IloNum s2 = PrunerCplex.getValue ( s[1] );

                if ( IloMax ( s1 , s2 ) == 0.0 )
                { // If the objective function value is equal to zero, the local Nadir point is dominated and the node cannot be pruned
                    if ( false == CheckAllNadirs )
                    {   // If we should not ckeck all Nadirs, we can now safely tell the program, that the node cannot be pruned!
                        return false;
                    }
                    else
                    {
                        // If we have to check all Nadirs, this particular Nadirs is put into the vector DominatedNadirs which is later copied to Nadirs
                        DominatedNadirs.push_back( *it );
                        if ( 3 <= DominatedNadirs.size ( ) )
                        {
                            break;
                        }
                    }
                }
            }
        }
        // Copy the dominated Nadirs (if any) to the Nadirs vector
        //std::cout << "Size of Dominated Nadirs : " << DominatedNadirs.size ( ) << "\n";
        //if (!DominatedNadirs.empty()) std::cout << "Nadir : " << DominatedNadirs.begin()->first << "\t" << DominatedNadirs.begin()->second << "\n";
        Nadirs.clear ( );
        Nadirs = DominatedNadirs;

        // If no dominated nadirs were found (Nadirs.empty()==true) we can prune the node, and we return true. Else we return false
        return ( Nadirs.empty ( ) );

    }
    catch ( const IloException &ie )
    {
        std::cerr << "IloException in the implicitPruner : " << ie.getMessage ( ) << "\n";
        exit ( EXIT_FAILURE );
    }
    catch ( const std::exception &e )
    {
        std::cerr << "std::exception in implicitPruner : " << e.what ( ) << "\n";
        exit ( EXIT_FAILURE );
    }
    catch ( ... )
    {
        std::cerr << "Unforseen error in implicitPruner in the BOCOsolver class!\n";
        exit ( EXIT_FAILURE );
    }
}
