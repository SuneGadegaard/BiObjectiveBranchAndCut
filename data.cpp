#include"data.h"

int Euclid(std::pair<double,double> p1, std::pair<double,double> p2){
    double x = p1.first-p2.first;
    double y = p1.second - p2.second;
    double dist = std::pow(x,2)+std::pow(y,2);
    dist = std::sqrt(dist);

    return int(10*dist+0.5);
}


data::data(char* filename, int Format):n(0),m(0),TD(0){
    try{
        readData(filename, Format);
    }catch(std::exception &e){
        std::cerr << "Exception in constructor of the data class : " << e.what() << std::endl;
        exit(1);
    }
}

/*------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------*/

void data::readData(char* filename, int Format){
    try{
        double anyInt = 0.0;
        std::stringstream err;
        std::ifstream data_file(filename);

        if(!data_file) throw std::runtime_error(std::string("Could not open the file: ") + filename + std::string("\n\n I will terminate now!\n\n"));

        if ( 0 == Format )
        {
            data_file >> n;
            if( !data_file ) throw std::runtime_error("Could not read the number of facilities. Termintaing!\n");

            data_file >> m;
            if( !data_file ) throw std::runtime_error("Could not read the number of customers. Termintaing!\n");

            if(n<=0 || m<=0){
                err << "Non positve values for number of facilities and customers. Provided was n=" << n << ", and m=" << m << ". Terminating\n";
                throw std::runtime_error(err.str().c_str());
            }

            for(int i=0; i<n; ++i){ // Read the capacities and fixed costs of the facilities
                data_file >> anyInt;
                if( (!data_file) || (anyInt<= 0) ){
                    err << "Could not read the capacity of facility " << i << ". Terminating!\n";
                    throw std::runtime_error( err.str() );
                }else s.push_back( anyInt );

                data_file >> anyInt;
                if( (!data_file) || (anyInt < 0) ){
                    err << "Could not read the fixed cost of facility " << i << ". Terminating!\n";
                    throw std::runtime_error( err.str() );
                }else f.push_back( anyInt );
            }

            for(int j=0; j<m; ++j){ // Reading the demands of the customers
                data_file >> anyInt;
                if( (!data_file) || (anyInt<=0) ){
                    err << "Could not read the demand of customer " << j << ". Terminating!\n";
                    throw std::runtime_error( err.str() );
                }else d.push_back( anyInt );
                TD+=d[j];
            }

            for(int i=0; i<n; ++i){ // Reading the assignment costs
                c.push_back(std::vector<int>());
                for(int j=0; j<m; ++j){
                    data_file >> anyInt;
                    if ( (!data_file) || (anyInt<0) ){
                        err << "Could not read the assignment cost of indices (" << i << "," << j << "). Terminating!\n";
                        throw std::runtime_error(err.str());
                    }
                    c[i].push_back(anyInt);
                }
            }
        }
        else if ( 1 == Format )
        {
            data_file >> m;
            if(!data_file) throw std::runtime_error("Could not read the number of customers. Termintaing!\n");

            data_file >> n;
            if(!data_file) throw std::runtime_error("Could not read the number of facilities. Termintaing!\n");

            if(n<=0 || m<=0){
                err << "Non positve values for number of facilities and customers. Provided was n=" << n << ", and m=" << m << ". Terminating\n";
                throw std::runtime_error(err.str().c_str());
            }

            for(int i=0; i<n; ++i) c.push_back ( std::vector<int>(m) );

            for(int j=0; j<m; ++j){ // Reading the assingment costs
                for(int i=0; i<n; ++i){
                    data_file >> anyInt;
                    if(!data_file || anyInt < 0){
                        err << "Could not read the assignment cost of index (" << i << "," << j << "j). Terminating!\n";
                        throw std::runtime_error(err.str().c_str());
                    }
                    c[i][j] = int(anyInt);
                }
            }

            TD = 0;
            for(int j=0; j<m; ++j){ // Reading the demands
                data_file >> anyInt;
                if( ( !data_file ) || ( anyInt<=0 ) ){
                    err << "Could not read the demand of customer " << j << ". Terminating!\n";
                    throw std::runtime_error(err.str().c_str());
                }
                d.push_back( int(anyInt) );
                TD += d[j];
            }

            for( int i=0; i<n; ++i)
            {   // Reading the fixed opening costs
                data_file >> anyInt;
                if(!data_file || anyInt<0 )
                {
                    err << "Could not read the fixed opening cost of facility " << i << ". Terminating!\n";
                    throw std::runtime_error(err.str().c_str());
                }
                f.push_back( int( anyInt ) );
            }

            for(int i=0; i<n; ++i)
            {
                data_file >> anyInt;
                if(!data_file || anyInt<=0 )
                {
                    err << "Could not read the capacity of facility " << i << ". Terminating!\n";
                    throw std::runtime_error(err.str().c_str());
                }
                s.push_back( int( anyInt ) );
            }
        }
    }catch(std::exception &e){
        std::cerr << "Exception in readData of the data class : " << e.what() << std::endl;
        exit(1);
    }
}

/*************************************************************************************************************/
void data::generateRandomData ( double BoxLB, double BoxUB, int dLB, int dUB, int sLB, int sUB, int fLB, int fUB , int seed )
{
    try
    {
        int FixedScal = 0;
        double  TotalDemand   = 0.0,
                TotalCapacity = 0.0,
                TheRatio      = 0.0,
                Scale         = 0.0;
        std::vector<std::pair<double,double>> FacPos ( n );
        std::vector<std::pair<double,double>> CustPos ( m );
        std::vector<int> MinC ( n );
        std::pair<double,double> p;
        G::result_type TheSeed = seed;
        G generator ( TheSeed );
        UniReal Position    = UniReal ( BoxLB , BoxUB );
        UniReal Ratio       = UniReal ( 1.5, 4.0 );
        UniInt  Fixed       = UniInt ( 1, 5 );
        UniInt  Demand      = UniInt ( dLB , dUB );
        UniInt  Capacity    = UniInt ( sLB , sUB );

        for ( int i=0; i<n; ++i )
        {
            p.first = Position( generator );
            p.second = Position ( generator );
            FacPos[i] = p;
            MinC[i] = INT_MAX;
        }

        for ( int j=0; j<m; ++j )
        {
            p.first = Position( generator );
            p.second = Position ( generator );
            CustPos[j] = p;
        }

        c = std::vector< std::vector< int > > ( n );
        f = std::vector< int > ( n );
        s = std::vector< int > ( n );
        for ( int i=0; i<n; ++i )
        {
            c[i] = std::vector< int > ( m );
            for ( int j=0; j<m; ++j )
            {
                c[i][j] = Euclid( FacPos[i] , CustPos[j] );
                if ( c[i][j] < MinC[i] ) MinC[i] = c[i][j];
            }
            FixedScal = Fixed( generator );
            //f[i] = FixedScal * 100;
            f[i] = ( MinC[i] + 1 ) * 10;
            s[i] = Capacity ( generator );
            TotalCapacity += double ( s[i] );
        }

        d = std::vector< int > ( m );
        TD = 0;
        for ( int j = 0; j<m; ++j )
        {
            d[j] = Demand ( generator );
            TD += d[j];
            TotalDemand += double ( d[j] );
        }
        TheRatio = Ratio ( generator );
        std::cout << "The ratio is : " << TheRatio << "\n";
        Scale = TheRatio * ( TotalDemand / TotalCapacity );
        for( int i=0; i<n; ++i )
        {
            s[i] = int ( Scale * s[i] + 0.5 );
        }



    }catch(std::exception &e)
    {
        std::cerr << "Exception in the generateRandomData of the data class : " << e.what ( ) << std::endl;
        exit ( EXIT_FAILURE );
    }
}

/*************************************************************************************************************/
void data::convertToUFLP ( )
{
    try{
        std::vector<int> Cmin ( n );
        for ( int i=0; i<n; ++i )
        {
            Cmin [ i ] = INT_MAX;
            for ( int j=0; j<m; ++j )
            {
                c[i][i] = d[j] * c[i][j];
                if ( c[i][j] < Cmin[i] ) Cmin[i] = c[i][j];
            }
            f[i] = Cmin[i] * 10;
            s[i] = m;
        }
        for ( int j=0; j<m; ++j ) d[j] = 1;
    }
    catch ( std::exception & e )
    {
        std::cerr << "Exception in convertToUFLP in data class : " << e.what ( ) << std::endl;
        exit ( EXIT_FAILURE );
    }
}
