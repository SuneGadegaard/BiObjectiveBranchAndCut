#include"NonDominatedSet.h"
/*
 * Implementation of the NondominatedSet class declared in NonDominatedSet.h
 */

void NonDominatedSet::addToFrontier2(std::pair<double,double> p)
{
    try
    {
        bool    ShouldInsertPoint   = true,  // Flag used to indicate if p should be inserted in the UBset
                UpdateNadir         = false; // Flag used to indicate if we should update the Nadir list
        auto Insert = NDS.end();    // Iterator pointing to the position in which p should be inserted. Asuming it should be inserted at the end.
        if( NDS.empty() ){ // If UBset is empty, just push_back p
            NDS.push_back( p );
        }else{
            for ( auto it=NDS.begin(); it!=NDS.end() && ShouldInsertPoint; ++it )
            { // Loop through the list as long as it is not known if p should be inserted
                if ( (it->first -Tol <= p.first) && (it->second-Tol <= p.second) )
                { // The new point is dominated, and we do not need to do more!
                    ShouldInsertPoint = false;
                }
                else if( (it->first +Tol >= p.first) && (it->second +Tol >= p.second) )
                { // The new point dominates *it and we overwrite *it and check the rest of the list
                    UpdateNadir = !UpdateNadir; // We have added a new point to the frontier, so the nadir points should be updated
                    it->first  = p.first;
                    it->second =p.second;
                    ShouldInsertPoint = !ShouldInsertPoint; // We should not inser p twice
                    // Check the rest of the points on the list
                    while ( (std::next(it)!=NDS.end( )) && (it->first -Tol<= std::next(it)->first) && (it->second -Tol<=std::next(it)->second) ){
                        // Loop as long as next(it) is not end of list and it dominates next(it)
                        // At the same time remove next(it)
                        NDS.erase(std::next(it));
                    }
                }else if ( (std::next(it)!=NDS.end()) && (it->first -Tol<= p.first) && (std::next(it)->first +Tol>= p.first) ){ // If p survives the chek, p should be inserted after it
                    Insert = std::next(it);
                }
            }
            if ( ShouldInsertPoint ){ // need to find the place where p should be inserted
                UpdateNadir = !UpdateNadir;
                if ( NDS.end() == Insert ) NDS.push_back( p ); // If the position in which we should insert is in the end, just push_back
                else NDS.insert(Insert, p); // Otherwise, insert at the found position.
            }
            // Upudating the list of Nadir points
            if ( UpdateNadir ){
                NadirPoints.clear();
                for ( auto it = NDS.begin(); std::next(it)!=NDS.end(); ++it ) {
                    NadirPoints.push_back( std::pair<double,double>( std::next(it)->first , it->second ) );
                }
            }
        }
    }catch(std::exception &e){
        std::cerr << "Exception in addToFrontier in the NonDominatedSet class : " << e.what( ) << std::endl;
        exit ( EXIT_FAILURE );
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
void NonDominatedSet::addToFrontier ( std::pair<double,double> p )
{
    try
    {
        // If the Non-Domintaed Set is empty, push back p and return
        if ( NDS.empty ( ) )
        {
            NDS.push_back( p );
        }
        else
        {
            auto LastPointBeforeP = NDS.begin ( );
            // Find the first element in NDS which has a strictly smaller first coordinate compared to p
            // If p turns out to be non-dominated, p should be inserted after this point. If not, p will be dominated by this point and p can be discarded
            for ( auto it = NDS.begin ( ); it!=NDS.end ( ) ; ++it )
            {
                if ( it->first <= p. first )
                {   // As the list is sorted, we can keep overwriting LastPointBeforeP
                    LastPointBeforeP = it;
                }
                else
                {   // As the list is kept sorted, we can
                    break;
                }
            } // End for

            // We have now found the point which p should be inserted after, if it is not dominated.
            // First we check if p is in fact dominated by LastPointBeforeP
            if ( ( p.first >= LastPointBeforeP->first ) && ( p.second >= LastPointBeforeP->second) )
            {
                // p is dominated, and we can simply return without doing anything further
                return;
            }
            else
            {
                // p is not dominated, and we need to insert p in the list. Note that list::insert wants to insert BEFORE the iterator it gets as argument
                // Now we need to know, if p dominates LastPointBeforeP or not. If p dominates, we overwrite LastPointBeforeP, otherwise we insert p after LastPointBeforeP
                std::cout << "* ";
                /*if ( ( p.first <= LastPointBeforeP->first ) && ( p.second <= LastPointBeforeP->second ) )
                {
                    LastPointBeforeP->first = p.first;
                    LastPointBeforeP->second = p.second;
                    std::advance ( LastPointBeforeP , 1 );
                }*/
                if ( true )
                {
                    LastPointBeforeP = std::next ( LastPointBeforeP );
                    NDS.insert ( LastPointBeforeP , p );
                }

                // Now that p is inserted we need to check if p dominates any of the points following p
                for ( auto it = LastPointBeforeP; it!=NDS.end ( ); )
                {
                    // Check if p dominates it. If so, erase it
                    if ( ( p.first <= it->first ) && ( p.second <= it->second ) )
                    {
                        it = NDS.erase ( it );
                    }
                    else
                    {
                        // p does not dominate it, it will not dominate any points after it. Therefore we can safely return
                        return;
                    }
                }
                if ( ( p.first <= NDS.back ( ).first) && ( p.second <= NDS.back ( ).second) )
                {
                    //NDS.erase ( std::prev ( NDS.end ( ) ) );
                }
            }

        }
    }
    catch ( std::exception &e )
    {
        std::cerr << "ERROR in addToFrontier2 in the NonDominatedSet class : " << e.what ( ) << "\n";
        exit ( EXIT_FAILURE );
    }
    catch ( ... )
    {
        std::cerr << "Unknown ERROR in addToFrontier2 in the NonDominatedSet class\n";
        exit ( EXIT_FAILURE );
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
void NonDominatedSet::setSubset ( unsigned int startPt, unsigned int endPt ) {
    try{
        if ( !NDS.empty () ){
            startIt = NDS.begin ( );
            std::advance ( startIt , ( startPt-1 ) );

            if ( endPt < NDS.size ( ) ){ // if endPt is small enough, advance the iterator
                endIt = NDS.begin ( );
                std::advance ( endIt, ( endPt-1 ) );
            }
            else{ // if endPt is too large, set it to end and decrement to last element
                endIt = --NDS.end();
            }
        }else{ // Size the list is empty, we just set the iterators to lsit end
            startIt = endIt = NDS.end ( );
        }
    }catch(std::exception &e){
        std::cerr << "Exception in the setSubset function in the NonDominatedSet class : " << e.what() << std::endl;
        exit ( EXIT_FAILURE );
    }
}
