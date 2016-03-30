#ifndef NONDOMINATEDSET_H_INCLUDED
#define NONDOMINATEDSET_H_INCLUDED

//---------------------------------------------------------------------------
#include<list>      //! Provides a double linked list
#include<iostream>  //! Write to stdout
#include<fstream>   //! Write to file
#include<stdexcept> //! Provides exception types
//---------------------------------------------------------------------------

/*!
 * Class implementing the functionalities of a non-dominated set.
 */


class NonDominatedSet{
    private:
        const double Tol;
        std::list< std::pair<double,double> >::iterator startIt; //! Iterator pointing to the start of a subset of NDS
        std::list< std::pair<double,double> >::iterator endIt;  //! Iterator pointing to the end of a subset of NDS
        std::list< std::pair<double,double> > NDS; //! List of non-dominated points (Non-Dominated Set)
        std::list< std::pair<double,double> > NadirPoints; //! List containing all the Nadir points. Updated whenever addToFrontier is called
    public:
        /*!
         * Constructor of the NonDominatedSet class
         */
        NonDominatedSet():Tol(1E-4){};

        /*!
         * Adds a point to the list UBset in the right place.
         * Works in O(n) time where n is the size of the current list.
         * \param p pair of doubles. p.first is the first coordinate of the point, and p.second is the second coordinate.
         */
        void addToFrontier(std::pair<double,double> p);

        /*!
         * Adds the point p to the frontier if p is not (weakly) dominated
         * Keeps the frontier sorted in increasing order of p.first
         * \param p pair of doubles. Contains the coordinates in R^2 of the point p which we want to insert.
         */
        void addToFrontier2 ( std::pair<double,double> p );

        /*!
         * Returns the size of the NDS
         * \return The size of the non-dominated set as an unsigned int
         */
        inline
        unsigned int getSize(){ return NDS.size(); }

        /*!
         * Method providing an iterator pointing to the first point on the frontier
         * \return An iterator pointing to UBset.begin()
         */
        inline
        std::list< std::pair<double,double> >::iterator UBsetBegin(){ return NDS.begin(); };

        /*!
         * Method providing an iterator pointing to the "point right after" the last point on th frontier
         * \return An iterator pointing to UBset.end()
         */
        inline
        std::list< std::pair<double,double> >::iterator UBsetEnd(){ return NDS.end(); }

        /*!
         * Returns a pointer to the list UBset
         * \return Pointer to the UBset list
         */
        inline
        std::list< std::pair<double,double> >* getUBset(){ return &NDS; }

        /*!
         * Method that sets the start and end iterators of a subset of the NDS
         * \param startPt unsigned int. Gives the position of the first element in the subset
         * \param endPt unsigned int. Gives the position of the last element in the subset
         * \note We must have startPt<=endPt. Furthermore, if endPt>=NDS.size() we set, by default, endPt = iterator to last element on NDS.
         */
        void setSubset(unsigned int startPt, unsigned int endPt);

        /*!
         * Deletes all points on the frontier and frees memory
         */
        inline
        void clear ( ) { NDS.clear ( ); }

};

#endif // NONDOMINATEDSET_H_INCLUDED
