#include "alg_nonlin.h"
#include "alg_lin.h"
#include "comparer.h"
#include <stdlib.h>

/**
 * @author Erthalion 9erthalion6@gmail.com
 * @mainpage Program for modeling motion of waves on the surface of an ideal incompressible fluid by an incomplete approximation of minimal residual. 
 * 
 * There is three general object - solver for linear equations,solver for nonlinear equations and comparer (for comparison linear and nonlinear equations). Now to change calculation parameters, is necessary to change them in main.cpp code, or create before start:
 * @code
 * if(argc != 10)
 *   {
 *       printf("insert 10 param\n");
 *       return -1;
 *   }
 * AlgLin *A=new AlgLin(atof(argv[1]),atof(argv[2]),atoi(argv[3]),atoi(argv[4]),atof(argv[5]),argv[6],atoi(argv[7]),atof(argv[8]));
 * A->solve(atoi(argv[9]));
 * @endcode
 * 
 * There is example of start nonlinear solver:
 * @code
 * AlgNonlin *A=new AlgNonlin(1,1,30,30,0.0000001,(char*)"out_nonlin.dat",100,0.01);
 * A->solve(3);
 * @endcode
 *
 * There is example of start comparer:
 * @code
 * Comparer *cmp=new Comparer(50,50,1.0/(50-1),1.0/(50-1),80);
 * cmp->delta((char*)"out_nonlin.dat",(char*)"out_lin.dat",(char*)"delta.dat");
 * @endcode
 *
 * Main calculation process accompanied by printing thread, which displays all debugging information.
 *
 * @warning this is unix thread,not crossplatform; to use on windows,should be removed all lines of code with pthread
 *
 * We use two ways to build the project:
 * <ul>
 * <li>Compile with CMake. See Build/CMakeList.txt
 * <li>Compile with own script for Profile-guided optimization (until do not understand, how to implement is in CMake). See Build/make.sh. This script has one parameter to choise mode of compilation - 1 to collect profile statistic, 2 to compile with using exist statistic (general mode).
 * </ul>
 * Usage second way:
 * @code
 * ./make.sh 1
 * then
 * ./make.sh 2
 * @endcode
 * 
 * @warning compile whith make.sh in first mode only short-time calculation because it slowly!
 *
 * Optimization flags of gcc are used in the CMakeList.txt and make.sh:
 * @code
 * -o3 -march=native -funsafe-math-optimizations -fomit-frame-pointer -fno-stack-protector -finline-small-functions -fprefetch-loop-arrays -funsafe-loop-optimizations -funroll-loops -std=c++0x -mpreferred-stack-boundary=6 -ftracer
 *  @endcode
 *
 *  In the CMakeList.txt defined two compilation type - Release and Debug. To choise:
 *
 *  @code
 *  cmake -DCMAKE_BUILD_TYPE=Release ./
 *  make
 *  @endcode
 * 
 *  @note The code has OpenMP pragmas - parallelization of program has not yer completed. It is not effective with this pragmas (not compile with -fopenmp)
 *
 *  Actual version of program is available in the repository
 *  http://code.google.com/p/erthalion-nm/
 */
int main(int argc, char **argv)
{
    AlgNonlin *A=new AlgNonlin(1,1,30,30,0.0000001,(char*)"out_nonlin.dat",100,0.01);
    A->solve(3);
    //A->getH((char*)"h.dat");

    //AlgLin *B=new AlgLin(1,1,30,30,0.0000001,(char*)"out_lin.dat",600,0.01);
    //B->solve(3);

    //Comparer *cmp=new Comparer(50,50,1.0/(50-1),1.0/(50-1),80);
    //cmp->delta((char*)"out_nonlin.dat",(char*)"out_lin.dat",(char*)"delta.dat");
	return 0;
}
