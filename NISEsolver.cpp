#include"NISEsolver.h"

/*!
 * Simple function returning the euclidean distance between two points
 * \param p1 pair of doubles. The first point
 * \param p2 pair of doubles. The second point
 */
double dist(std::pair<double,double> p1, std::pair<double,double> p2){
    double x = p1.first - p2.first;
    double y = p1.second- p2.second;
    double dist = pow(x,2)+pow(y,2);
    dist = sqrt(dist);
    return dist;
}

/*
 * Implementation of the NISEsolver class declared in NISEsolver.h
 */


//-------------------------------------------------------------------------------------------------//
NISEsolver::NISEsolver():CPXerr(0),ProblemType(0),n(0),m(0),TotalDemand(0){
    try{
        // Set parameters
        myTol = 0.0001;
        myOne = 1.0-myTol;

        // Allocating memory and initialize cplex stuff
        model = IloModel(env);
        cplex = IloCplex(model);
        y = IloNumVarArray(env);
        x = IloVarMatrix(env);
        z = IloVarMatrix(env);
        f1 = IloNumVar( env, -IloInfinity, IloInfinity, ILOFLOAT );
        f2 = IloNumVar( env, -IloInfinity, IloInfinity, ILOFLOAT );
        CpxEnv = CPXopenCPLEX(&CPXerr);
        if(CPXerr) throw std::runtime_error("Could not open the cplex environment for the VICkpsep function.\n");

        // Initialize the cutlist to point to the null pointer
        cutlst = NULL;
        for ( int cuts =0; cuts<NCUTS; ++cuts ) NumCuts[cuts] = 0;
        NumCuts[0] = 0;
        NumCuts[1] = 0;
        NumCuts[2] = 0;
        NumCuts[3] = 0;
    }catch(std::runtime_error &re){
        std::cerr << "Runtime error in the constructor of the NISEsolver class : " << re.what() << std::endl;
        exit(1);
    }catch(std::exception &e){
        std::cerr << "Exception in NISEsolver constructor : " << e.what() << std::endl;
        exit(1);
    }catch(IloException &ie){
        std::cerr << "IloException in NISEsolver constructor : " << ie.getMessage() << std::endl;
        exit(1);
    }
}

//-------------------------------------------------------------------------------------------------//
void NISEsolver::BuildSSCFLP(int NumFac, int NumCust, std::vector<int> demand, std::vector<int> supply, std::vector<int> FixedCost, std::vector<std::vector<int>> VarCost){
    try{
        // Copy data
        n = NumFac;
        m = NumCust;
        for ( int i=0 ; i<n ; ++i ){
            s.push_back( supply[i] );
            fixed.push_back( FixedCost[i] );
            c.push_back( std::vector<int>(m) );
            for ( int j=0 ; j<m ; ++j ){
                c[i][j] = VarCost[i][j];
                if( 0 == i ){ d.push_back( demand[j] ); TotalDemand += d[j];}
            }
        }
        // Build the model
        IloExpr cst (env);
        IloExpr f1cst (env);
        IloExpr f2cst (env);
        IloExpr TDcst (env);
        OBJ = IloMinimize (env, f1+f2 );
        model.add ( OBJ );
        for ( int i=0 ; i<n ; ++i ){
            // Add variable to the y array
            y.add ( IloNumVar ( env, 0, 1, ILOFLOAT) );
            // Build the left hand side of the definition of the objective 2 constraint
            f2cst+=fixed[i]*y[i];
            // Add array of variables to the x matrix
            x.add ( IloNumVarArray ( env, m, 0, 1, ILOFLOAT ) );
            // Start building the total demand constraint
            TDcst+=s[i]*y[i];
            // Build the assignment constraint and the first objective function
            for ( int j=0 ; j<m ; ++j ){
                f1cst += c[i][j]*x[i][j];
                cst+=d[j]*x[i][j];
                model.add(x[i][j]-y[i]<=0);
            }
            // Add the assignment constraint
            model.add ( cst <= s[i]*y[i] );
            // Clear the cst expression
            cst.clear ();
        }
        // At the total demand constrain
        model.add ( TDcst >= TotalDemand );
        // Add the objective function constraints
        model.add ( f1 == f1cst );
        model.add ( f2 == f2cst );
        // Create and add the demand constraints
        for ( int j=0 ; j<m ; ++j ){
            for ( int i=0 ; i<n ; ++i ) cst += x[i][j];
            model.add ( cst == 1 );
            cst.clear ();
        }

        cst.end ( );
        f1cst.end ( );
        f2cst.end ( );
        TDcst.end ( );

        cplex.setOut(env.getNullStream());
        cplex.setWarning(env.getNullStream());
    }catch ( std::runtime_error &re ){
        std::cerr << "Runtime error in BuildSSCFLP of the NISEsolver class : " << re.what () << std::endl;
        exit ( 1 );
    }catch ( std::exception &e ){
        std::cerr << "Exception in BuildSSCFLP of the NISEsolver class: " << e.what () << std::endl;
        exit ( 1 );
    }catch ( IloException &ie ){
        std::cerr << "IloException in BuildSSCFLP of the NISEsolver class : " << ie.getMessage () << std::endl;
        exit ( 1 );
    }
}

//-------------------------------------------------------------------------------------------------//
void NISEsolver::RunNISE ( bool doCuts , IloModel &CutModel , bool& PruneTheBugger){
    try{
        YnLP.clear();
        double lambda;
        int maxiter = 10000;
        std::list<std::pair<double,double> >::iterator PlusIndex;
        std::list<std::pair<double,double> >::iterator MinusIndex;
        std::list<std::pair<double,double> >::iterator lit;
        std::pair<double,double> yplus;
        std::pair<double,double> yminus;
        std::pair<double,double> yul(0.0,0.0);
        std::pair<double,double> ylr(0.0,0.0);
        std::pair<double,double> p(0.0,0.0);

        // Setting the coefficients for the fixed cost objective, in order to calculate yul
        OBJ.setLinearCoef(f1,0.99999);
        OBJ.setLinearCoef(f2,0.00001);
        // Run the cutting phase for this parametized program
        //SolveModel();
        if(doCuts) CuttingPhase( CutModel );
        // Need to solve again, as cuts might have been added
        if(!cplex.solve()){
            std::cerr << "Cplex status is : " << cplex.getStatus() << std::endl;
            PruneTheBugger = true;
            return;
            //throw std::runtime_error("Cplex could not solve the problem after adding cuts when finding yul");
        }
        // Now we have a solution in cplex, and we need to fill that one in to yul
        yul.first = cplex.getValue( f1 );
        yul.second= cplex.getValue( f2 );

        // Put the yul on the list of non-dominated LP solutions
        YnLP.push_back(yul);

        // Now change the objective function, such that we find ylr
        OBJ.setLinearCoef(f1,0.00001);
        OBJ.setLinearCoef(f2,0.99999);
        // Run the cutting phase for this parametized program
        if(doCuts) CuttingPhase( CutModel );
        // Need to resolve, as cuts might have been added
        if(!cplex.solve()) throw std::runtime_error("Cplex could not solve the problem after adding cuts when finding ylr");
        ylr.first = cplex.getValue ( f1 );
        ylr.second= cplex.getValue ( f2 );

        // Put the ylr on the list of noon-dominated solutions for the BOCO-LP
        YnLP.push_back(ylr);

        // Set the coefficient of the objectives back to 1
        OBJ.setLinearCoef(f1,1.0);
        OBJ.setLinearCoef(f2,1.0);

        // Now we enter the main part of the algorithm
        PlusIndex   = YnLP.begin();
        MinusIndex  = std::next(PlusIndex,1);
        yplus       = *PlusIndex;
        yminus      = *MinusIndex;



        for ( int i=0 ; (i<maxiter) && (dist(ylr, *PlusIndex)>=myZero) ; ++i ){
            // Updating the scalar lambda
            if ( MinusIndex->first-PlusIndex->first == 0 ) break;
            lambda = (PlusIndex->second - MinusIndex->second)/(MinusIndex->first-PlusIndex->first);

            /*std::cout << "y- is (" << MinusIndex->first << "," << MinusIndex->second << ")\n";
            std::cout << "y+ is (" << PlusIndex->first << "," << PlusIndex->second << ")\n";
            std::cout << "Value of Lambda is " << lambda << std::endl;*/
            // Scale the objective function using lambda
            //for(int i=0; i<n; ++i)for(int j=0; j<m; ++j) OBJ.setLinearCoef(x[i][j], double(c[i][j])*lambda );
            OBJ.setLinearCoef(f1,lambda);
            // Solve the scalarized problem and add cuts
            if(doCuts) CuttingPhase( CutModel );
            // We need to resolve, as we might have added cuts.
            if(!cplex.solve()) throw std::runtime_error("Cplex could not solve the problem after adding cuts in the main loop.");

            //First we check if we found a new solution or not:
            if(cplex.getObjValue()<= (lambda*(PlusIndex->first)+PlusIndex->second)-myTol ){
                // The solution is new and we calculate the coefficients of the new point.
                p.first = cplex.getValue ( f1 );
                p.second = cplex.getValue( f2 );
                YnLP.insert(MinusIndex,p);
            }else{
                // The solution is not new and we go to next point
                PlusIndex = MinusIndex;
                yplus = *PlusIndex;
            }
            // We now move the yminus to the next point in the set YnLP
            MinusIndex = std::next(PlusIndex,1);
            yminus = *MinusIndex;
        }

        YnLP.sort();
        /*std::cout << "====================== output from NISEsolver begin ======================\n";
        for ( auto it = YnLP.begin(); it!=YnLP.end(); ++it ){
            std::cout << it->first << "\t" << it->second << std::endl;
        }
        std::cout << "======================  Output from NISEsolver ends ======================\n";*/
    }catch ( std::runtime_error &re ){
        std::cerr << "Runtime error in RunNISE of the NISEsolver class : " << re.what () << std::endl;
        exit ( 1 );
    }catch ( std::exception &e ){
        std::cerr << "Exception in RunNISE of the NISEsolver class: " << e.what () << std::endl;
        exit ( 1 );
    }catch ( IloException &ie ){
        std::cerr << "IloException in RunNISE of the NISEsolver class : " << ie.getMessage () << std::endl;
        exit ( 1 );
    }
}

//-------------------------------------------------------------------------------------------------//
std::vector<std::pair<double,double>> NISEsolver::RUN ( bool doCuts, IloModel &CutModel , bool &PruneTheBugger) {
    try{

        std::vector<std::pair<double,double>> ReturnFrontier;

        RunNISE ( doCuts , CutModel , PruneTheBugger );
        for( auto it = YnLP.begin() ; it != YnLP.end(); ++it){
            ReturnFrontier.push_back( *it );
        }
        return ReturnFrontier;
    }
    catch(std::exception &e){throw;}
    catch(IloException &ie){throw;}
    catch(...){throw;}
}
//-------------------------------------------------------------------------------------------------//
void NISEsolver::fixXVariables( std::vector< std::vector<int> > fixed){
    try{
        // Loop through all the x-variables and set the bounds according to the bounds described in the vector fixed
        for ( int i=0 ; i < n ; ++i )
            for ( int j=0 ; j<m ; ++j )
                if ( fixed[i][j] <= 1 ) x[i][j].setBounds ( fixed[i][j] , fixed[i][j] );
                else x[i][j].setBounds ( 0 , 1 );
    }catch( IloException &ie ){
        std::cerr << "IloException in fixXVariables in the NISEsolver class : " << ie.getMessage() << std::endl;
        exit ( 1 );
    }catch( std::exception &e ){
        std::cerr << "Exception in fixXVariables in the NISEsolver class : " << e.what() << std::endl;
        exit ( 1 );
    }
}
//-------------------------------------------------------------------------------------------------//
void NISEsolver::fixYvariables( std::vector<int> fixed ){
    try{
        // Loop through all the y-variables and set the bounds according to the bounds described in the vector fixed
        for ( int i=0 ; i < n ; ++i )
            if ( fixed[i] <= 1 ) y[i].setBounds ( fixed[i] , fixed[i] );
            else y[i].setBounds ( 0 , 1 );
    }catch( IloException &ie ){
        std::cerr << "IloException in fixYvariables in the NISEsolver class : " << ie.getMessage() << std::endl;
        exit ( 1 );
    }catch( std::exception &e ){
        std::cerr << "Exception in fixYvariables in the NISEsolver class : " << e.what() << std::endl;
        exit ( 1 );
    }
}

//-------------------------------------------------------------------------------------------------//
void NISEsolver::SolveModel(){
    try{
        if(cplex.solve()){// We solved the model
            OldObjVal = ObjVal;
            ObjVal = cplex.getObjValue();
            TerminateCuttingPhase = (ObjVal-OldObjVal < 0.9);
        }else{
            std::stringstream err;
            err << "Cplex did not solve the model. Staturs is " << cplex.getStatus() << std::endl;
            throw std::runtime_error(err.str().c_str());
        }
    }catch(std::exception &e){
        std::cout << "Exception in SolveModel in Algorithm class: " << e.what() << std::endl;
        throw e;
    }catch(IloException &i){
        std::cout << "IloException in SolveModel in Algorithm class: " << i.getMessage() << std::endl;
        throw i;
    }
}

//-------------------------------------------------------------------------------------------------//
void NISEsolver::CuttingPhase(IloModel &CutModel){
    try{
        bool CONTINUE = true;
        int Iindex, iterations=0, cutsgen=0;
        std::vector<std::pair<int,int> > GUBs;
        //VICcut* cutlst = NULL;
        VICcut* thiscut = NULL;
        VICcut* ThisCutList = NULL;


        int err = 0;
        CpxEnv = CPXopenCPLEX(&err);
        if(err) throw std::runtime_error("Could not open the cplex environment for the VICkpsep function.\n");
        while(CONTINUE && (iterations < 5) ){
            thiscut = NULL;
            ThisCutList  = NULL;
            SolveModel();
            for(int i=0; i<n; ++i){
                if(true) thiscut = getLCI(i); // Try to find a violated lifted cover inequality
                if(thiscut == NULL) // If no such found, try to find violated extedn cover inequality
                  if(true) thiscut = getECI(i);
                if(thiscut == NULL) // If no such found, try the fenchel cutting planes
                  if(true) thiscut = getFCP(i);
                if( thiscut ){ // Add the cut we found to the cut list cutlst
                    VICaddtolst(&ThisCutList , thiscut);
                    cutsgen++;
                }

            }
            GUBs = getGUBs();
            for(std::vector<std::pair<int,int> >::iterator it=GUBs.begin(); it!=GUBs.end(); ++it){
                model.add(x[it->first][it->second]-y[it->first]<=0);
            }

            if(ThisCutList!=NULL){
                VICaddtolst( &cutlst , ThisCutList );
                for(VICcut* cut = ThisCutList; cut!= NULL; cut=cut->nextcut){
                    IloExpr cst(env);
                    Iindex = cut->UsrIVal;
                    for(int nz=0; nz < cut->nzcnt; nz++){
                        cst+= double(cut->nzval[nz])*x[Iindex][cut->nzind[nz]];
                    }
                    model.add( cst<= ( double(cut->rhs)*y[Iindex] )  );
                }

            }else CONTINUE = false;
            CONTINUE = (iterations<=2) || (ObjVal - OldObjVal > 0.1);
            ++iterations;
            //VICfreelst(&cutlst);
        }

        err = CPXcloseCPLEX(&CpxEnv);
    }catch(std::exception &e){
        std::cout << "Exception in CuttingPhase in Algorithm class: " << e.what() << std::endl;
        throw;
    }catch(IloException &ie){
        std::cout << "IloException in CuttingPhase in Algorithm class: " << ie.getMessage() << std::endl;
        throw;
    }
}

/*-----------------------------------------------------------------------------------------*/
VICcut* NISEsolver::getLCI(int i){
    try{

        VICcut* mycut   = NULL;

        double* xlp = new double[m];
        double* rco = new double[m];
        int*    idx = new int[m];
        int*    dem = new int[m];

        for(int j=0; j<m; ++j){
            idx[j] = j;
            xlp[j] = cplex.getValue(x[i][j]);
            rco[j] = double(c[i][j])/double(d[j]);
            dem[j] = d[j];
        }

        VIClci(m, s[i], 'L', dem, idx, xlp, rco, &mycut);

        if(mycut!=NULL){
            mycut->UsrIVal = i;
            NumCuts[1]++;
        }

        delete[] rco;
        delete[] xlp;
        delete[] idx;
        delete[] dem;
        return mycut;
    }catch(std::exception &e){
        std::cout << "Exception in getLCI in Algorithm class: " << e.what() << std::endl;
        throw;
    }catch(IloException &ie){
        std::cout << "IloException in getLCI in Algorithm class: " << ie.getMessage() << std::endl;
        throw;
    }
}

/*-----------------------------------------------------------------------------------------*/
VICcut* NISEsolver::getECI(int i){
    try{
        VICcut* mycut = NULL;
        int*    idx = new int[m];
        int*    dem = new int[m];
        double* xlp = new double[m];

        for(int j=0; j<m; ++j){
            idx[j] = j;
            xlp[j] = cplex.getValue(x[i][j]);
            dem[j] = d[j];
        }

        VICecikl(m, s[i], 0, 'L', dem, idx, xlp, &mycut);


        if(mycut!=NULL){
            mycut->UsrIVal = i;
            NumCuts[2]++;
        }

        delete[] idx;
        delete[] dem;
        delete[] xlp;

        return mycut;

    }catch(std::exception &e){  std::cerr << "Exception in getECI in Algorithm class: " << e.what() << std::endl; throw;
    }catch(IloException &ie){   std::cerr << "IloException in getECI in Algorithm class: " << ie.getMessage() << std::endl; throw 1;}
}

//-----------------------------------------------------------------------------------------//
VICcut* NISEsolver::getFCP(int i){
    try{
        VICcut* mycut = NULL;
        int*    idx = new int[m];
        int*    dem = new int[m];
        double* xlp = new double[m];
        double* rco = new double[m];

        for(int j=0; j<m; ++j){
            idx[j] = j;
            xlp[j] = cplex.getValue(x[i][j]);
            rco[j] = double(c[i][j])/double(d[j]);
            dem[j] = d[j];
        }

        if(cplex.getValue(y[i])>myTol) VICkpsep(CpxEnv, 1, m, s[i], 'L', dem, idx, xlp, rco, &mycut);

        if(mycut != NULL){
            mycut->UsrIVal=i,
            NumCuts[3]++;
        }

        delete[] idx;
        delete[] xlp;
        delete[] rco;

        return mycut;
    }catch(std::exception &e){  std::cerr << "Exception in getFCP in Algorithm class: " << e.what() << std::endl; throw;
    }catch(IloException &ie){   std::cerr << "IloException in getFCP in Algorithm class: " << ie.getMessage() << std::endl; throw;
    }catch(...){ throw; }
}

//-----------------------------------------------------------------------------------------//
VICcut* NISEsolver::getCutsForTD(){
    try{
        VICcut* mycut = NULL;

        double* ylp = new double[n];
        double* rco = new double[n];
        int*    idx = new int[n];
        int*    cap = new int[n];

        for(int i=0; i<n; ++i){
            idx[i] = i;
            rco[i] = double(fixed[i])/double(s[i]);
            ylp[i] = double( cplex.getValue(y[i]) );
            cap[i] = s[i];
        }

        VIClci(n, TotalDemand, 'G', cap, idx, ylp, rco, &mycut);
        if(NULL == mycut){
            VICecikl(n, TotalDemand, 0, 'G', cap, idx, ylp, &mycut);
            if(NULL == mycut){
                VICkpsep(CpxEnv, 0, n, TotalDemand, 'G', cap, idx, ylp, rco, &mycut);
                if(mycut) NumCuts[3]++;
            }
            else NumCuts[2]++;
        }else NumCuts[1]++;


        if(mycut){
            mycut->UsrIVal=-1;
        }

        delete[] ylp;
        delete[] rco;
        delete[] idx;

        return mycut;

    }catch(std::exception &e){
        std::cerr << "Exception in getCutsForTD in Algorihm class: " << e.what() << std::endl;
        throw;
    }catch(IloException &ie){
        std::cerr << "IloException in getCutsForTD in Algorithm class: " << ie.getMessage() << std::endl;
        throw;
    }
}

//-----------------------------------------------------------------------------------------//
std::vector<std::pair<int,int> > NISEsolver::getGUBs(){
    try{
        std::vector<std::pair<int,int> > GUBs;
        double yval,xval;
        for(int i=0; i<n; ++i)
        {
            yval = cplex.getValue(y[i]);
            if(yval>myTol)
            {
                for(int j=0; j<m; ++j)
                {
                    xval = cplex.getValue(x[i][j]);
                    if(xval-yval>=0.01)
                    {
                        GUBs.push_back(std::pair<int,int>(i,j));
                        NumCuts[0]++;
                    }
                }
            }
        }
        return GUBs;
    }catch(std::exception &e){  std::cerr << "Exception in getGUBs in Algorithm class: " << e.what() << std::endl; throw;
    }catch(IloException &ie){   std::cerr << "IloException in getGUBs inAlgortihm class: " << ie.getMessage() << std::endl; throw;}
}

//-------------------------------------------------------------------------------------------------//
NISEsolver::~NISEsolver ( ) {
    // Deallocate memory occupied by CPLEX and its variables
    if ( x.getSize() > 0 ){
        IloInt xSize = x.getSize();
        for ( int i=0; i<xSize; ++i ) x[i].end();
    }
    if ( z.getSize() > 0 ){
        IloInt zSize = z.getSize();
        for ( int i=0; i<zSize; ++i ) z[i].end();
    }
    if ( x.getImpl() )      x.end();
    if ( z.getImpl() )      z.end();
    if ( y.getImpl() )      y.end();
    if ( f1.getImpl() )     f1.end();
    if ( f2.getImpl() )     f2.end();
    if ( cplex.getImpl() )  cplex.end();
    if ( model.getImpl() )  model.end();
    if ( env.getImpl() )    env.end();
    if ( cutlst )           VICfreelst( & cutlst );
    if ( (CPXerr = CPXcloseCPLEX(&CpxEnv)) ) std::cerr    << "CPXcloseCPLEX returned non-zero error messeage in destructor of NISEsolver. Error message : "
                                                        << CPXerr << std::endl;
}
