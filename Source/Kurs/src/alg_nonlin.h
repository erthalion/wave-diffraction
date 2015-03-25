#ifndef ALG_NONLIN_H_INCLUDED
#define ALG_NONLIN_H_INCLUDED

#define g 9.8
#include "alg.h"
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <unistd.h>
#include <time.h>

/**
 * @brief Class for calculation nonlinear equation of shallow water
 *
 * @detailed This class contains specific for nonlinear equations functionali
 */
class AlgNonlin : public Alg{
    /**
     * not used
     */
    double ret[3];
    public:
     /**
     * @brief Constructor
     * @param in_height Height calculation area
     * @param in_lenght Lenght calculation area
     * @param in_nx Count of Ox node
     * @param in_ny Count of Oy node
     * @param in_eps The accuracy of calculations
     * @param in_out_name Output filename for save result
     * @param count Count iteration by time
     * @param in_t Step by time
     */
    AlgNonlin(double in_height,double in_lenght,int in_nx,int in_ny,double in_eps,char *in_out_name, int in_count, double in_t):Alg(in_height,in_lenght,in_nx,in_ny,in_eps,in_out_name,in_count,in_t){};
    
    /**
     * @brief Function for recalculation of reciduals for all equations
     */
    void resolve_r();

/**
     * @bread Function for recalculation of reciduals for one equation
     * @param mode number of equation for recalculation
     */
 
    void resolve_r(int mode);
    
    /**
     * @brief Function to find method coefficient
     * @param a array of data for first equation
     * @param b array of data for second equation
     * @param c array of data for third equation
     */
    double find(double **a,double **b,double **c);
    
    /**
     * @brief Function to find method coefficient for one equation
     * @param mode number of equation
     * @param a array of data for first equation
     * @param b array of data for second equation
     * @param c array of data for third equation
     */
    double find(int mode,double **a,double **b,double **c);

/**
     * not used
     */
    double* find_all(double **a,double **b,double **c);

/**
     * @brief Function for independent verification of reciduals
     */
    void check_result();

};

#endif // ALG_NONLIN_H_INCLUDED
