#ifndef DATA_H_INCLUDED
#define DATA_H_INCLUDED

#include<vector>
#include<stdexcept>
#include<iostream>
#include<fstream>
#include<sstream>
#include<limits.h>
#include<random>
#include<chrono>

typedef std::uniform_int_distribution<> UniInt; //! Uniform distribution on integers
typedef std::uniform_real_distribution<> UniReal;
typedef std::mt19937_64 G; //! Random number generator based on the 64-bit Mersenne Twister by Matsumoto and Nishimura, 2000
typedef std::chrono::high_resolution_clock myclock;

class data{
    private:
        /**
         * @name Data
         * This section contains all the data for describing the SSCFLP.
         */
        ///@{
        int n; //! Number of facility sites
        int m; //! Number of customers
        int TD; //! Total demand. TD=sum_{j=1}^n d[j]
        std::vector<std::vector<int>> c; //! Assignment cost. c[i][j] is the cost of assigning customer j to facility i
        std::vector<int> f; //! Fixed opening cost. f[i] is the cost of opening facility i
        std::vector<int> d; //! Demand. d[j] is the demand at customer j
        std::vector<int> s; //! Capacity. s[i] is the capacity of facility i
        ///@}

        /*!
         * Function taking a datafile and reads the data into the member fields.
         * \param filename pointer to char array. Contains the path of valid data file
         * \param Format integer. Indicates the format of the data file.
         */
        void readData(char* filename, int Format);
    public:
        /*! \brief Constructor of the data class
         *
         * Constructor of the data class. Takes as input a pointer to a char array containing the path to the datafile.
         * The datafile must be of the form:
         * n
         * m
         * s[1] f[2]
         * s[2] f[2]
         * ...
         * s[n] f[n]
         * d[1] d[2] ... d[m]
         * c[1][1] c[1][2] ... c[1][m]
         *...
         * c[n][1] c[n][2] ... c[n][m]
         * \param filename pointer to char array. Contains the path (relative or absolute) to the data file.
         * \param Format integer. Indicating the format of the data file. 0 = yang, holmberg, gadegaard. 1 = diaz fernandez
         */
        data(char* filename, int Format);

        /*!
         * Returns the number of facilities
         */
        int getNumFaci(){return n;}
        /*!
         * Returns the number of customers
         */
        int getNumCust(){return m;}
        /*!
         * Returns the total demand
         */
        int getTotalDemand(){return TD;}
        /*!
         * Returns a vector of vectors of integers containing the assignment costs
         */
        std::vector<std::vector<int>> getAssignmentCosts(){return c;}
        /*!
         * Returns the entry (i,j) of the assignment cost matrix
         * \param i integer. Index of the facility.
         * \param j integer. Index of the customer.
         */
        int getC(int i, int j){
          if( ( c.size()>0 ) && (0 <= i && i < n) && (0 <= j && j < m) ) return c[i][j];
          else{
            std::cout << "Index i="<<i<< ". Index j=" << j<< std::endl;
            throw std::out_of_range("Asking for c[i][j] which does not exist\n");
          }
        }
        /*!
         * Returns a vector of ints containing the fized opening costs.
         */
        std::vector<int> getFixedCosts(){return f;}
        /*!
         * Returns the fixed cost of facility i
         * \param i integer. Index of the facility.
         */
        int getFixedCost(int i){
            if( (f.size()>0) && (0<=i) && (n>i) ) return f[i];
            else throw std::out_of_range("Asking for f[i] which does not exist\n");
        }
        /*!
         * Returns a vector of ints containing the demands
         */
        std::vector<int> getDemands(){return d;}
        /*!
         * Returns the demand of customer j
         * \param j integer. Index of the customer
         */
        int getDemand(int j){
          if( (d.size()>0) && (0<=j) && (m>j) ) return d[j];
          else throw std::out_of_range("Asking for d[j] which does not exist\n");
        }
        /*!
         * Return a vector of ints containing the capacities of the facilities
         */
        std::vector<int> getCapacitites(){return s;}
        /*!
         * Returns the capacity of facility i
         * \param i integer. Index of the facility
         */
        int getCapacity(int i){
            if( (s.size()>0) && (0<=i) && (n>i) ) return s[i];
            else throw std::out_of_range("Asking for s[i] which does not exist\n");
        }

        inline
        void setNumFacilities ( int NewN ) { n = NewN; }

        inline
        void setNumCustomers ( int NewM) { m = NewM; }

        void generateRandomData ( double BoxLB, double BoxUB, int dLB, int dUB, int sLB, int sUB, int fLB, int fUB , int seed );

        void convertToUFLP ();
        // Does not need an explicit destructor as everything is STL'ed
};

#endif // DATA_H_INCLUDED
