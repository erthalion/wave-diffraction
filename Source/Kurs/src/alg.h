#ifndef ALG_H_INCLUDE
#define ALG_H_INCLUDE

#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include "definition.h"

/**
 * @brief Base class for calculation
 * @detailed This class contains general functional and abstract realization algorythm
 */
class Alg
{
    protected:
        /**
         * free surface
         */
        double **pov;

        /**
         * speen by Ox
         */
        double **u;

        /**
         * speed by Oy
         */
        double **v;

        /**
         * residual for first equation
         */
        double **r1;

        /**
         * residual for second equation
         */
        double **r2;

        /** residual for third equation
        */
        double **r3;

        /**
         * for using as previous value on time
         */
        double **pov_temp1;
        /**
         * for using as previous value on time
         */
        double **u_temp1;
        /**
         * for using as previous value on time
         */
        double **v_temp1;

        /**
         * for using as previuos value in time layer (for method)
         */
        double **pov_temp2;
        /**
         * for using as previuos value in time layer (for method)
         */
        double **u_temp2;
        /**
         * for using as previuos value in time layer (for method)
         */
        double **v_temp2;

        double **temp1;
        double **temp2;
        double **temp3;

        /**
         * temporary for solving (placed here only for optimize memory allocation)
         */
        double **HBu;
        /**
         * temporary for solving (placed here only for optimize memory allocation)
         */
        double **HBv;
        /**
         * temporary for solving (placed here only for optimize memory allocation)
         */
        double **uu;
        /**
         * temporary for solving (placed here only for optimize memory allocation)
         */
        double **vv;

        /**
         * temporary for solving (placed here only for optimize memory allocation)
         */
        double **HBb;
        /**
         * temporary for solving (placed here only for optimize memory allocation)
         */
        double **HBc;
        /**
         * temporary for solving (placed here only for optimize memory allocation)
         */
        double **bb;
        /**
         * temporary for solving (placed here only for optimize memory allocation)
         */
        double **cc;

        /**
         * for calculation alfa
         */
        double **z;

        /**
         * for checking(debug)
         */
        double **check1;
        /**
         * for checking(debug)
         */
        double **check2;
        /**
         * for checking(debug)
         */
        double **check3;

        /**
         * Area parameters
         */
        double height,lenght;

        /**
         * Num of nodes
         */
        int nx,ny;

        /**
         * Steps by space/time
         */
        double hx,hy,t;

        /**
         * Current iteration
         */
        int iter;

        /**
         * Presice of method of incomplete approximation of minimal reciduals
         */
        double eps;

        /**
         * Name for output file
         */
        char out_name[256];

        /**
         * Count steps by time
         */
        int count;

        /**
         * Norm of reciduals
         */
        double mod;

        /*
         * to find average some variables(debug)
         */
        double a1,a2,a3;
        /*
         * to find average some variables(debug)
         */
        int count_a1,count_a2,count_a3;

        FILE *log;

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
         * @param tau Strp by time
         */
        Alg(double in_height,double in_lenght,int in_nx,int in_ny,double in_eps,char *in_out_name,int count,double tau);

        /**
         * @brief Recalculation of residuals for all equations
         */
        virtual void resolve_r()=0;

        /**
         * @brief Recalculation of residuals for one equation
         * @param mode Number of equation for recalculation
         */
        virtual void resolve_r(int mode)=0;

        /**
         * @brief Find iteration parameter. Is a wrapper for "find"
         * @param name Name of iteration parameter (now used to find only TAU)
         * @return Iteration parameter (TAU)
         */
        double find_coeff(int name);

        /**
         * @brief Find iteration parameter.As changed algorythm of ALFA calculation, its implementation is moved to find with parameter "mode"
         * @param a First data array for operator (to find TAU -- pov)
         * @param b Second data array for operator (to find TAU -- u)
         * @param c Third data array for operator (to find TAU -- v)
         * @return Iteration parameter (TAU)
         */
        virtual double find(double **a,double **b,double **c)=0;

        /**
         * @brief Find iteration parameter for one equation. Used to find ALFA
         * @param mode Number of equations, for which coefficient searched
         * @param a First data array for operator (to find ALFA -- z)
         * @param b Second data array for operator (to find ALFA -- z)
         * @param c Third data array for operator (to find ALFA -- z)
         * @return Iteration parameter (ALFA)
         */
        virtual double find(int mode,double **a,double **b,double **c)=0;

        /**
         * @bried Deprecated method for calculation ALFA
         */
        virtual double* find_all(double **a, double **b, double **c)=0;

        /**
         * @brief Verify the calculation of residuals
         */
        virtual void check_result()=0;

        /**
         * @brief Main calculation loop by time
         * @param out_n Frequency of write result on disk (for large calculations)
         */
        void solve(int out_n);

        /**
         * @brief Recalculation pov,u,v. Is a wrapper for resolve.
         * @param coeff Iteration coefficient (now using only for TAU)
         * @param name Name of iteration coefficient (TAU or ALFA)
         */
        void resolve_param(double coeff,int name);

        /**
         * @brief Recalculation pov,u,v with iteration parameter. After changing algorythm for ALFA, using for recaulculation with TAU.
         * @param coeff Value of iteration parameter
         * @param a First data array (for TAU - pov)
         * @param b Second data array (for TAU - u)
         * @param c Third data array (for TAU - v)
         */
        void resolve(double coeff,double **a,double **b,double **c);

        /**
         * @brief Deprecated for recalculation ALFA (now realized in resolve)
         */
        void resolve_for_alfa(double coeff,double **a,double **b,double **c);

        /**
         * @brief differential analog of derivative by x
         * @param a variable which used for derivation
         * @param i node number (for check boundaries)
         * @param j node number (for check boundaries)
         */
        double diffX(double **a,int i,int j){
            //уменьшаем шаг сетки на краях
            if(i==0)
                return (a[1][j]-a[0][j])/(hx);

            if(i==(nx-1))
                return (a[nx-1][j]-a[nx-2][j])/(hx);

            //if(i<4)
            //    return (a[i+1][j]-a[i-1][j])/(2*hx/2.0);

            //if(i>(nx-5))
            //    return (a[i+1][j]-a[i-1][j])/(2*hx/2.0);

            return (a[i+1][j]-a[i-1][j])/(2*hx);
        };

        /**
         * @brief differential analog of derivative by y
         * @param a variable which used for derivation
         * @param i node number (for check boundaries)
         * @param j node number (for check boundaries)
         */
        double diffY(double **a,int i,int j){
            //уменьшаем шаг сетки на краях
            if(j==0)
                return (a[i][1]-a[i][0])/(hy);

            if(j==ny-1)
                return (a[i][ny-1]-a[i][ny-2])/(hy);

            //if(j<10)
            //    return (a[i][j+1]-a[i][j-1])/(2*hy);

            //if(j>(ny-11))
            //    return (a[i][j+1]-a[i][j-1])/(2*hy);

            return (a[i][j+1]-a[i][j-1])/(2*hy);
        };

        /**
         * @brief For access to pov from other threads
         * @return pov
         */
        double** return_pov(){return pov;};

        /**
         * @brief For access to u from other threads
         * @return u
         */
        double** return_u(){return u;};

        /**
         * @brief For access to v from other threads
         * @return v
         */
        double** return_v(){return v;};

        /**
         * @brief Describe static function of bottom in point (i,j)
         * @param i Node number
         * @param j Node number
         * @return Bottom depth
         */
        double H(int i,int j){stair(i,j); return 0.01;};

        /**
         * @brief Describe dynamic function of bottom in point (i,j)
         * @param i Node number
         * @param j Node number
         * @return Bottom depth
         */
        double B(int i,int j,int iter){return 0;};

        /**
         * @brief To set temporary (by time or method iterations) variables
         * @param number If set pov_temp1,u_temp1,v_temp1 - 0. If set pov_temp2,u_temp2_v_temp2 - 1.
         */
        void set_temp(int number);

        /**
         * @brief Find residual norm
         * @return Residual norm
         */
        double get_norm_r();

        /**
         * @brief Deprecated function for printing (debug)
         */
        void print(double **a,char *name);
        /**
         * @brief Deprecated function for printing (debug)
         */
        void print(double **a);
        /**
         * @brief Deprecated function for printing (debug)
         */
        void print_tec(double **a,char *out);

        /**
         * @brief To set arbitrary initial data
         * @param out_u Path to file with u data
         * @param out_v Path to file with v data
         * @param out_pov Path to file with pov data
         */
        void set_u_v_pov(char *out_u,char *out_v,char *out_pov);
        /**
         * @brief To set arbitrary initial data
         * @param out_u Path to file with u data
         * @param out_v Path to file with v data
         */

        void set_u_v_pov(char *out_u,char *out_v);
        /**
         * @brief To set arbitrary initial data, which described by some function
         */
        void set_u_v_pov();
        /**
         * @brief To output pov,u,v for future use in set_u_v_pov
         */
        void get_u_v_pov(char *out_u,char *out_v,char *out_pov);

        /**
         * @brief For acces to iteration number from other threads
         * @return iteration number
         */
        int getIter();
        /**
         * @brief For acces to iteration count from other threads
         * @return iteration count
         */
        int getCount();
        /**
         * @brief For acces to residual norm from other threads
         * @return residual norm
         */
        double getMod();

        /**
         * @brief Output bottom depth in file
         * @param file_name Output file name
         */
        void getH(char *file_name);

        /**
         * @brief Dectructor
         */
        ~Alg();
};


#endif
