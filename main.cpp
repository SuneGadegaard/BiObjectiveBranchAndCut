#include <iostream>
#include"BOCOsolver.h"
#include"NISEsolver.h"
#include"data.h"
#include<time.h>
#include<stdlib.h>
#include<chrono>

using namespace std;
using namespace std::chrono;
typedef std::chrono::high_resolution_clock CPUclock;

int main(int argc, char** argv)
{
    try
    {

        data SSCFLPdata = data ( argv[1] , 1 );

        std::string FileName = "SummaryExpEPBNoBugs.txt";
        std::ofstream OutputFile;
        OutputFile.open ( FileName , std::ofstream::out | std::ofstream::app);
        if ( !OutputFile )
        {
            throw std::runtime_error ( "Could not open the summaryfile\n" );
        }

        OutputFile  << "n \t m \t time \t LPprune \t PNpolyPrune \t PNPOLYnPrune \t LBupdates \t BOLPsolved \t SolutionsFound \t NumNodes \t MaxNodesInTree \t NumParBranch \t"
                    << "NumOrdPB \t NumGPB \t NumCuts \t SizeFront \t CutPhaseTime\n";
        OutputFile.close ( );

        int n;
        for ( n = 5; n<=60; n=n+5 )
        {
            for ( int trials = 1; trials<=10; ++trials )
            {
                OutputFile.open ( FileName , std::ofstream::out | std::ofstream::app);
                int m = 2*n;//std::atoi ( argv[3] );
                SSCFLPdata.setNumCustomers ( m );
                SSCFLPdata.setNumFacilities ( n );
                SSCFLPdata.generateRandomData( 0.0 , 10.0 , 5 , 10 , 10 , 20 , 0 , 10 , trials );
                //SSCFLPdata.convertToUFLP ( );

                int TotCap = 0;
                for ( int i=0; i<SSCFLPdata.getNumFaci(); ++i )
                {
                    TotCap += SSCFLPdata.getCapacity( i );
                }
                cout << "Total capacity is    " << TotCap << endl;


                BOCOsolver solver = BOCOsolver(SSCFLPdata);
                cout << "Total demand after handing to solver " << solver.DATA->getTotalDemand() << endl;
                auto Start = CPUclock::now ( );
                solver.run();
                auto End = CPUclock::now ( );


//                std::cout << "Total time used : " << End - start << std::endl;

                //OutputFile  << "n \t m \t time \t LPprune \t PNpolyPrune \t PNPOLYnPrune \t LBupdates \t BOLPsolved \t SolutionsFound \t NumNodes \t MaxNodesInTree \t NumParBranch \t"
                //    << "NumOrdPB \t NumGPB \n";
                auto TotalTime = duration_cast<duration<double>>( End - Start ).count ( );
                OutputFile  << std::to_string ( n ) << "\t"
                            << std::to_string ( m ) << "\t"
                            << TotalTime << "\t"
                            << solver.LinearPruned << "\t"
                            << solver.PNPOLYpruned << "\t"
                            << solver.PNPOLYnotPruned << "\t"
                            << solver.LBupdates << "\t"
                            << solver.BOLPsolved << "\t"
                            << solver.SolutionsFound << "\t"
                            << solver.NumNodes << "\t"
                            << solver.MaxNodesInTree << "\t"
                            << solver.NumberOfParetoBranches << "\t"
                            << solver.NumberOfOrdPB << "\t"
                            << solver.NumberOfGPB << "\t"
                            << solver.NumCuts << "\t"
                            << solver.getUBsetSize() << "\t"
                            << solver.CuttingPhaseTime << std::endl;
                OutputFile.close ( );
            }
        }

        return 0;
    }
    catch ( std::exception &e )
    {
        std::cerr << "ERROR in main : " << e.what( ) << "\n";
        exit( EXIT_FAILURE );
    }
    catch ( ... )
    {
        std::cerr << "Unknown ERROR in main\n";
        exit ( EXIT_FAILURE );
    }
}
