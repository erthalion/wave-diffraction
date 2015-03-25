/**
 * \file definition.h
 * This file contains the definitions of macro, which describe forms of free surface of fluid, form of bottom and some constants.
 */

/**
 * @brief Long Ox centering wave.
 * @param height height of wave
 * @param grad degree of steepness of the ascent
 */
#define centering_wave(height,grad,i,j) height*exp(-grad*(((double)i/nx-0.5)*((double)i/nx-0.5)))

/**
 * @brief Long Ox wave.
 * @param height height of wave
 * @param grad degree of steepness of the ascent
 * @param shift offset from the origin of coordinates
 */
#define wave(height,grad,i,j,shift) height*exp(-grad*(((double)i/nx-shift)*((double)i/nx-shift)))

/**
 * @brief Small centering drop.
 * @param height height of drop
 * @param grad degree of steepness of the ascent
 */
#define centering_drop(height,grad,i,j) height*exp(-grad*(((double)i/nx-0.5)*((double)i/nx-0.5)+((double)j/ny-0.5)*((double)j/ny-0.5)))

/**
 * @brief Small drop.
 * @param height height of drop
 * @param grad degree of steepness of the ascent
 * @param shiftX offset from the origin of Ox
 * @param shiftY offset from the origin of Oy
 */
#define drop(height,grad,i,j,shiftX,shiftY) height*exp(-grad*(((double)i/nx-shiftX)*((double)i/nx-shiftX)+((double)j/ny-shiftY)*((double)j/ny-shiftY)))

/**
 * @brief Bottom function.
 *
 * @code
 * stair
 * -----
 *     |
 *     |---------
 *
 * @endcode
 */
#define stair(i,j) if(i<nx/10){return 0;}

/**
 * @brief Bottom function.
 *
 * @code
 * stair
 * -----
 *     |
 *     |---------
 *
 * @endcode
 */
#define stair_shift(i,j,shift) if(i<nx/shift){return 0;}
 
/**
 * @brief Bottom function.
 *
 * @code
 * rise
 * ----
 *     -
 *      -
 *       --------
 * @endcode
 */
#define rise(i,j) if(i<nx/3&&i>nx/5){return 0.01-(nx/3-i)*0.001;}; if(i<=nx/5){return 0.01-(nx/3-nx/5)*0.001;}

/**
 * @brief Bottom function.
 *
 * @code
 * cylinder
 *       --
 *       ||
 * ------  -------
 *  @endcode
 */
#define cylinder(i,j) if(sqrt((i*hx-0.5)*(i*hx-0.5)+(j*hy-0.5)*(j*hy-0.5))<0.05){return 0.2;}

/**
 * @brief Constants for TAU mode in solve process
 */
#define TAU 0

/**
 * @brief Constants for ALFA mode in solve process
 */
#define ALFA 1

/**
 * @brief Constants for g 
 */
#define g 9.8
	
//double result=0.01;
//if((i>100)&&(i<201))//&&(j>100)&&(j<200))
//	result=0.1+0.5*sin((i-100)*hx)*sin((j-100)*hy);
//if((i>nx/4)&&(i<3*nx/4)&&(j>ny/4)&&(j<3*ny/4))
//	result=0.05+0.05*sin((i-nx/4)*hx*3.1415)*sin((j-ny/4)*hy*3.1415);
//if(i>9*nx/10)
//	result=0;

//return 0.0073+0.005*atan((i-25)*0.1);
//if(i<nx/10)
//    result=0;
//if(i>=nx/10&&i<nx/5)
//    result=(double)i/(double)nx*0.1;

//if(sqrt((i*hx-0.7)*(i*hx-0.7)+(j*hy-0.5)*(j*hy-0.5))<0.05)
//	result=0.1;
//double result=0.007*exp(-10.0*(((double)i/nx-0.7)*((double)i/nx-0.7)));
//if(i<nx/10)
//    return 0;
//if(j<3||j>ny-4)
//    return 0.01+0.001*(rand()%10+1);
//if((i>10)&&(i<20)&&(j>10)&&(j<20))
//	return 0.01;
//if(i<nx/3&&i>nx/5)
//    return 0.01-(nx/3-i)*0.001;

//if(i<=nx/5)
//    return 0.01-(nx/3-nx/5)*0.001;
