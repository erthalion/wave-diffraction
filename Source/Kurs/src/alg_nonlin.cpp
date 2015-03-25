#include "alg_nonlin.h"
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

void AlgNonlin::resolve_r()
{
    register size_t i,j;
    //создаем дополнительные массивы
    #pragma omp parallel for collapse(2)
    for(i=0;i<nx;i++)
    {
        for(j=0;j<ny;j++)
        {
            HBu[i][j]=(H(i,j)-B(i,j,iter))*u[i][j];
            HBv[i][j]=(H(i,j)-B(i,j,iter))*v[i][j];
            uu[i][j]=u[i][j]*u[i][j];
            vv[i][j]=v[i][j]*v[i][j];
        }
    }

    //внутри области
    for(i=0; i<nx; i++) {
        for(j=0; j<ny; j++) {
            r1[i][j]=(pov[i][j]-pov_temp1[i][j])/t + diffX(HBu,i,j) + diffY(HBv,i,j)-(B(i,j,iter)-B(i,j,iter-1))/t;
            r2[i][j]=(u[i][j]-u_temp1[i][j])/t+g*diffX(pov,i,j) + v[i][j]*diffY(u,i,j) + diffX(uu,i,j)/2;
            r3[i][j]=(v[i][j]-v_temp1[i][j])/t+g*diffY(pov,i,j) + u[i][j]*diffX(v,i,j) + diffY(vv,i,j)/2;
        }
    }

}

void AlgNonlin::resolve_r(int mode)
{
    register size_t i,j;
    //создаем дополнительные массивы
    #pragma omp parallel for collapse(2)
    for(i=0;i<nx;i++)
    {
        for(j=0;j<ny;j++)
        {
            HBu[i][j]=(H(i,j)-B(i,j,iter))*u[i][j];
            HBv[i][j]=(H(i,j)-B(i,j,iter))*v[i][j];
            uu[i][j]=u[i][j]*u[i][j];
            vv[i][j]=v[i][j]*v[i][j];
        }
    }

    if(mode==0){
        //внутри области
#pragma omp parallel for collapse(2)
        for(i=0; i<nx; i++) {
            for(j=0; j<ny; j++) {
                r1[i][j]=(pov[i][j]-pov_temp1[i][j])/t + diffX(HBu,i,j) + diffY(HBv,i,j)-(B(i,j,iter)-B(i,j,iter-1))/t;
            }
        }
    }

    if(mode==1){
        //внутри области
#pragma omp parallel for collapse(2)
        for(i=0; i<nx; i++) {
            for(j=0; j<ny; j++) {
                r2[i][j]=(u[i][j]-u_temp1[i][j])/t+g*diffX(pov,i,j) + v[i][j]*diffY(u,i,j) + diffX(uu,i,j)/2;
            }
        }
    }

    if(mode==2){
        //внутри области
#pragma omp parallel for collapse(2)
        for(i=0; i<nx; i++) {
            for(j=0; j<ny; j++) {
                r3[i][j]=(v[i][j]-v_temp1[i][j])/t+g*diffY(pov,i,j) + u[i][j]*diffX(v,i,j) + diffY(vv,i,j)/2;
            }
        }
    }

    //fprintf(log,"%2.20lf\n",diffY(HBv,0,0));
}

double* AlgNonlin::find_all(double **a,double **b,double **c)
{
    double temp=0;
    double arr=0;
    double arar=0;
    double p1=0,p2=0,p3=0;
    double q1=0,q2=0,q3=0;
    double coeff1,coeff2,coeff3;

    register size_t i,j;
    //создаем дополнительные массивы
    #pragma omp parallel for collapse(2)
    for(i=0;i<nx;i++)
    {
        for(j=0;j<ny;j++)
        {
            HBb[i][j]=(H(i,j)-B(i,j,iter))*b[i][j];
            HBc[i][j]=(H(i,j)-B(i,j,iter))*c[i][j];
            bb[i][j]=b[i][j]*b[i][j];
            cc[i][j]=c[i][j]*c[i][j];
        }
    }


    //внутри области
    for(i=0; i<nx; i++) {
        for(j=0; j<ny; j++) {
            temp=(a[i][j])/t + diffX(HBb,i,j) + diffY(HBc,i,j);
            arr=temp*r1[i][j];
            arar=temp*temp;
            p1+=arr;
            q1+=arar;

            temp=(b[i][j])/t+g*diffX(a,i,j) + c[i][j]*diffY(b,i,j) + diffX(bb,i,j)/2;
            arr=temp*r2[i][j];
            arar=temp*temp;
            p2+=arr;
            q2+=arar;

            temp=(c[i][j])/t+g*diffY(a,i,j)+ b[i][j]*diffX(c,i,j) + diffY(cc,i,j)/2;
            arr=temp*r3[i][j];
            arar=temp*temp;
            p3+=arr;
            q3+=arar;
        }
    }

    coeff1=p1/q1;
    coeff2=p2/q2;
    coeff3=p3/q3;

    ret[0]=coeff1;ret[1]=coeff2;ret[2]=coeff3;
    return ret;
}

double AlgNonlin::find(double **a,double **b,double **c)
{
    double temp=0;
    double arr=0;
    double arar=0;
    double p=0;
    double q=0;
    double coeff;

    register size_t i,j;
    //создаем дополнительные массивы
    #pragma omp parallel for collapse(2)
    for(i=0;i<nx;i++)
    {
        for(j=0;j<ny;j++)
        {
            HBb[i][j]=(H(i,j)-B(i,j,iter))*b[i][j];
            HBc[i][j]=(H(i,j)-B(i,j,iter))*c[i][j];
            bb[i][j]=b[i][j]*b[i][j];
            cc[i][j]=c[i][j]*c[i][j];
        }
    }


    //внутри области
    for(i=0; i<nx; i++) {
        for(j=0; j<ny; j++) {
            temp=(a[i][j])/t + diffX(HBb,i,j) + diffY(HBc,i,j);
            arr=temp*r1[i][j];
            arar=temp*temp;
            p+=arr;
            q+=arar;

            temp=(b[i][j])/t+g*diffX(a,i,j) + c[i][j]*diffY(b,i,j) + diffX(bb,i,j)/2;
            arr=temp*r2[i][j];
            arar=temp*temp;
            p+=arr;
            q+=arar;

            temp=(c[i][j])/t+g*diffY(a,i,j)+ b[i][j]*diffX(c,i,j) + diffY(cc,i,j)/2;
            arr=temp*r3[i][j];
            arar=temp*temp;
            p+=arr;
            q+=arar;
        }
    }

    coeff=p/q;

    return coeff;
}

double AlgNonlin::find(int mode,double **a,double **b,double **c)
{
    double temp=0;
    double arr=0;
    double arar=0;
    double p=0;
    double q=0;
    double coeff;

    register size_t i,j;
    //создаем дополнительные массивы
    #pragma omp parallel for collapse(2)
    for(i=0;i<nx;i++)
    {
        for(j=0;j<ny;j++)
        {
            HBb[i][j]=(H(i,j)-B(i,j,iter))*b[i][j];
            HBc[i][j]=(H(i,j)-B(i,j,iter))*c[i][j];
            bb[i][j]=b[i][j]*b[i][j];
            cc[i][j]=c[i][j]*c[i][j];
        }
    }


    if(mode==0){
        //внутри области
        for(i=0; i<nx; i++) {
            for(j=0; j<ny; j++) {
                temp=(a[i][j])/t + diffX(HBb,i,j) + diffY(HBc,i,j);
                arr=temp*r1[i][j];

                arar=temp*temp;
                p+=arr;
                q+=arar;
 
            }
        }

        coeff=p/q;

        return coeff;
    }

    if(mode==1){
        //внутри области
        for(i=0; i<nx; i++) {
            for(j=0; j<ny; j++) {
                temp=(b[i][j])/t+g*diffX(a,i,j) + c[i][j]*diffY(b,i,j) + diffX(bb,i,j)/2;
                arr=temp*r2[i][j];

                arar=temp*temp;

                p+=arr;
                q+=arar;

            }
        }

        coeff=p/q;

        return coeff;
    }

    if(mode==2){
        //внутри области
        for(i=0; i<nx; i++) {
            for(j=0; j<ny; j++) {
                temp=(c[i][j])/t+g*diffY(a,i,j)+ b[i][j]*diffX(c,i,j) + diffY(cc,i,j)/2;
                arr=temp*r3[i][j];
                arar=temp*temp;
                p+=arr;
                q+=arar;
            }
        }

        coeff=p/q;

        return coeff;
    }

    return (double)NULL;
}

void AlgNonlin::check_result(){
    register size_t i,j;
    //создаем дополнительные массивы
    #pragma omp parallel for collapse(2)
    for(i=0;i<nx;i++)
    {
        for(j=0;j<ny;j++)
        {
            HBu[i][j]=(H(i,j)-B(i,j,iter))*u[i][j];
            HBv[i][j]=(H(i,j)-B(i,j,iter))*v[i][j];
            uu[i][j]=u[i][j]*u[i][j];
            vv[i][j]=v[i][j]*v[i][j];
        }
    }

    //внутри области
    for(int i=1; i<nx-1; i++) {
        for(int j=1; j<ny-1; j++) {
            check1[i][j]=(pov[i][j]-pov_temp1[i][j])/t + (H(i+1,j) * u[i+1][j] - H(i-1,j)*u[i-1][j])/(2*hx) + (H(i,j+1)* v[i][j+1] - H(i,j-1)*v[i][j-1])/(2*hy);
            check2[i][j]=(u[i][j]-u_temp1[i][j])/t+g*(pov[i+1][j]-pov[i-1][j])/(2*hx) + v[i][j]*(u[i][j+1]-u[i][j-1])/(2*hy)+(u[i+1][j]*u[i+1][j]-u[i-1][j]*u[i-1][j])/(4*hx);
            check3[i][j]=(v[i][j]-v_temp1[i][j])/t+g*(pov[i][j+1]-pov[i][j-1])/(2*hy) + u[i][j]*(v[i+1][j]-v[i-1][j])/(2*hx)+(v[i][j+1]*v[i][j+1]-v[i][j-1]*v[i][j-1])/(4*hy);
        }
    }

    //на границе
    //Г1 : x=0 y=[0,1]
    //i=0
    for(int j=1; j<ny-1; j++) {
        check1[0][j]=(pov[0][j]-pov_temp1[0][j])/t + ( H(1,j) * u[1][j] - H(0,j)*u[0][j])/(hx) + ( H(0,j+1)* v[0][j+1] - H(0,j-1)*v[0][j-1])/(2*hy);
        check2[0][j]=(u[0][j]-u_temp1[0][j])/t+g*(pov[1][j]-pov[0][j])/(hx) + v[0][j]*(u[0][j+1]-u[0][j-1])/(2*hy)+(u[1][j]*u[1][j]-u[0][j]*u[0][j])/(2*hx);
        check3[0][j]=(v[0][j]-v_temp1[0][j])/t+g*(pov[0][j+1]-pov[0][j-1])/(2*hy) + u[0][j]*(v[1][j]-v[0][j])/(2*hx)+(v[0][j+1]*v[0][j+1]-v[0][j-1]*v[0][j-1])/(4*hy);
    }

    //Г2 : x=1 y=[0,1]
    //i=nx-1
    for(int j=1; j<ny-1; j++) {
        check1[nx-1][j]=(pov[nx-1][j]-pov_temp1[nx-1][j])/t + (H(nx-1,j) * u[nx-1][j] - H(nx-2,j)*u[nx-2][j])/(hx) + (H(nx-1,j+1) * v[nx-1][j+1] - H(nx-1,j-1)*v[nx-1][j-1])/(2*hy);
        check2[nx-1][j]=(u[nx-1][j]-u_temp1[nx-1][j])/t+g*(pov[nx-1][j]-pov[nx-2][j])/(hx); + v[nx-1][j]*(u[nx-1][j+1]-u[nx-1][j-1])/(2*hy)+(u[nx-1][j]*u[nx-1][j]-u[nx-2][j]*u[nx-2][j])/(2*hx);
        check3[nx-1][j]=(v[nx-1][j]-v_temp1[nx-1][j])/t+g*(pov[nx-1][j+1]-pov[nx-1][j-1])/(2*hy) + u[nx-1][j]*(v[nx-1][j]-v[nx-2][j])/(2*hx)+(v[nx-1][j+1]*v[nx-1][j+1]-v[nx-1][j-1]*v[nx-1][j-1])/(4*hy);
    }

    //Г3 : x=[0,1] y=0
    //j=0
    for(int i=1; i<nx-1; i++) {
        check1[i][0]=(pov[i][0]-pov_temp1[i][0])/t + (H(i+1,0) * u[i+1][0] - H(i-1,0)*u[i-1][0])/(2*hx) + ( H(i,1) * v[i][1] - H(i,0)*v[i][0])/(hy);
        check2[i][0]=(u[i][0]-u_temp1[i][0])/t+g*(pov[i+1][0]-pov[i-1][0])/(2*hx) + v[i][0]*(u[i][1]-u[i][0])/(2*hy)+(u[i+1][0]*u[i+1][0]-u[i-1][0]*u[i-1][0])/(4*hx);
        check3[i][0]=(v[i][0]-v_temp1[i][0])/t+g*(pov[i][1]-pov[i][0])/(hy) + u[i][0]*(v[i+1][0]-v[i-1][0])/(2*hx)+(v[i][1]*v[i][1]-v[i][0]*v[i][0])/(2*hy);
    }

    //Г4 : x=[0,1] y=1
    //j=ny-1
    for(int i=1; i<nx-1; i++) {
        check1[i][ny-1]=(pov[i][ny-1]-pov_temp1[i][ny-1])/t + (H(i+1,ny-1) * u[i+1][ny-1] - H(i-1,ny-1)*u[i-1][ny-1])/(2*hx) + (H(i,ny-1) * v[i][ny-1] - H(i,ny-2)*v[i][ny-2])/(hy);
        check2[i][ny-1]=(u[i][ny-1]-u_temp1[i][ny-1])/t+g*(pov[i+1][ny-1]-pov[i-1][ny-1])/(2*hx) + v[i][ny-1]*(u[i][ny-1]-u[i][ny-2])/(2*hy)+(u[i+1][ny-1]*u[i+1][ny-1]-u[i-1][ny-1]*u[i-1][ny-1])/(4*hx);
        check3[i][ny-1]=(v[i][ny-1]-v_temp1[i][ny-1])/t+g*(pov[i][ny-1]-pov[i][ny-2])/(hy) + u[i][ny-1]*(v[i+1][ny-1]-v[i-1][ny-1])/(2*hx)+(v[i][ny-1]*v[i][ny-1]-v[i][ny-2]*v[i][ny-2])/(2*hy);
    }

    //в углах
    //x=0 y=0
    check1[0][0]=(pov[0][0]-pov_temp1[0][0])/t + (H(1,0) * u[1][0] - H(0,0)*u[0][0])/(hx) + (H(0,1)* v[0][1] - H(0,0)*v[0][0])/(hy);
    check2[0][0]=(u[0][0]-u_temp1[0][0])/t+g*(pov[1][0]-pov[0][0])/(hx) + v[0][0]*(u[0][1]-u[0][0])/(2*hy)+(u[1][0]*u[1][0]-u[0][0]*u[0][0])/(2*hx);
    check3[0][0]=(v[0][0]-v_temp1[0][0])/t+g*(pov[0][1]-pov[0][0])/(hy) + u[0][0]*(v[1][0]-v[0][0])/(2*hx)+(v[0][1]*v[0][1]-v[0][0]*v[0][0])/(2*hy);

    //x=1 y=0
    check1[nx-1][0]=(pov[nx-1][0]-pov_temp1[nx-1][0])/t + (H(nx-1,0)* u[nx-1][0] - H(nx-2,0)*u[nx-2][0])/(hx) + (H(nx-1,1) * v[nx-1][1] - H(nx-1,0)*v[nx-1][0])/(hy);
    check2[nx-1][0]=(u[nx-1][0]-u_temp1[nx-1][0])/t+g*(pov[nx-1][0]-pov[nx-2][0])/(hx) + v[nx-1][0]*(u[nx-1][1]-u[nx-1][0])/(2*hy)+(u[nx-1][0]*u[nx-1][0]-u[nx-2][0]*u[nx-2][0])/(2*hx);
    check3[nx-1][0]=(v[nx-1][0]-v_temp1[nx-1][0])/t+g*(pov[nx-1][1]-pov[nx-1][0])/(hy) + u[nx-1][0]*(v[nx-1][0]-v[nx-2][0])/(2*hx)+(v[nx-1][1]*v[nx-1][1]-v[nx-1][0]*v[nx-1][0])/(2*hy);

    //x=0 y=1
    check1[0][ny-1]=(pov[0][ny-1]-pov_temp1[0][ny-1])/t + (H(1,ny-1) * u[1][ny-1] - H(0,ny-1)*u[0][ny-1])/(hx) + (H(0,ny-1)* v[0][ny-1] - H(0,ny-2)*v[0][ny-2])/(hy);
    check2[0][ny-1]=(u[0][ny-1]-u_temp1[0][ny-1])/t+g*(pov[1][ny-1]-pov[0][ny-1])/(hx) + v[0][ny-1]*(u[0][ny-1]-u[0][ny-2])/(2*hy)+(u[1][ny-1]*u[1][ny-1]-u[0][ny-1]*u[0][ny-1])/(2*hx);
    check3[0][ny-1]=(v[0][ny-1]-v_temp1[0][ny-1])/t+g*(pov[0][ny-1]-pov[0][ny-2])/(hy) + u[0][ny-1]*(v[1][ny-1]-v[0][ny-1])/(2*hx)+(v[0][ny-1]*v[0][ny-1]-v[0][ny-2]*v[0][ny-2])/(2*hy);

    //x=1 y=1
    check1[nx-1][ny-1]=(pov[nx-1][ny-1]-pov_temp1[nx-1][ny-1])/t + (H(nx-1,ny-1) * u[nx-1][ny-1] - H(nx-2,ny-1)*u[nx-2][ny-1])/(hx) + (H(nx-1,ny-1)* v[nx-1][ny-1] - H(nx-1,ny-2)*v[nx-1][ny-2])/(hy);
    check2[nx-1][ny-1]=(u[nx-1][ny-1]-u_temp1[nx-1][ny-1])/t+g*(pov[nx-1][ny-1]-pov[nx-2][ny-1])/(hx) + v[nx-1][ny-1]*(u[nx-1][ny-1]-u[nx-1][ny-2])/(2*hy)+(u[nx-1][ny-1]*u[nx-1][ny-1]-u[nx-2][ny-1]*u[nx-2][ny-1])/(2*hx);
    check3[nx-1][ny-1]=(v[nx-1][ny-1]-v_temp1[nx-1][ny-1])/t+g*(pov[nx-1][ny-1]-pov[nx-1][ny-2])/(hy) + u[nx-1][ny-1]*(v[nx-1][ny-1]-v[nx-2][ny-1])/(2*hx)+(v[nx-1][ny-1]*v[nx-1][ny-1]-v[nx-1][ny-2]*v[nx-1][ny-2])/(2*hy);

    //внутри области
    /*for(int i=0; i<nx; i++) {
        for(int j=0; j<ny; j++) {
            check1[i][j]-=(pov[i][j]-pov_temp1[i][j])/t + diffX(HBu,i,j) + diffY(HBv,i,j)-(B(i,j,iter)-B(i,j,iter-1))/t;
            check2[i][j]-=(u[i][j]-u_temp1[i][j])/t+g*diffX(pov,i,j) + v[i][j]*diffY(u,i,j) + diffX(uu,i,j)/2;
            check3[i][j]-=(v[i][j]-v_temp1[i][j])/t+g*diffY(pov,i,j) + u[i][j]*diffX(v,i,j) + diffY(vv,i,j)/2;
        }
    }*/
}
