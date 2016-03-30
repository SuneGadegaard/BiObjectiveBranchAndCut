#include"myNodeDAta.h"

myNodeData::myNodeData( const std::vector< std::pair< double , double > > &LBSET, int TimesSinse, bool DONISE):
    LBset(LBSET),TimesSinseLastSolved(TimesSinse),doNISE(DONISE){}
