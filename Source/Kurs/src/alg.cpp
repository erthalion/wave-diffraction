#include "alg.h"
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <omp.h>

Alg::Alg(double in_height,double in_lenght,int in_nx,int in_ny,double in_eps,char *in_out_name, int in_count, double in_t)
{
	iter=0;
	t=in_t;
	count=in_count;
	strcpy(out_name,in_out_name);
	eps=in_eps;
	nx=in_nx;
	ny=in_ny;
	height=in_height;
	lenght=in_lenght;
	hx=height/(nx-1);
	hy=lenght/(ny-1);

	pov=new double *[nx];
	u=new double *[nx];
	v=new double *[nx];

	pov_temp1=new double *[nx];
	u_temp1=new double *[nx];
	v_temp1=new double *[nx];

	pov_temp2=new double *[nx];
	u_temp2=new double *[nx];
	v_temp2=new double *[nx];

	//temp1=new double *[nx];
	//temp2=new double *[nx];
	//temp3=new double *[nx];

	r1=new double *[nx];
	r2=new double *[nx];
	r3=new double *[nx];

    HBu=new double*[nx];
    HBv=new double*[nx];
    uu=new double*[nx];
    vv=new double*[nx];

    HBb=new double*[nx];
    HBc=new double*[nx];
    bb=new double*[nx];
    cc=new double*[nx];

    z=new double*[nx];

    check1=new double*[nx];
    check2=new double*[ny];
    check3=new double*[ny];

    for(int i=0; i<nx; i++) {
		pov[i]=new double[ny];
		u[i]=new double[ny];
		v[i]=new double[ny];

		pov_temp1[i]=new double[ny];
		u_temp1[i]=new double[ny];
		v_temp1[i]=new double[ny];

		pov_temp2[i]=new double[ny];
		u_temp2[i]=new double[ny];
		v_temp2[i]=new double[ny];

		//temp1[i]=new double[ny];
		//temp2[i]=new double[ny];
		//temp3[i]=new double[ny];

		r1[i]=new double[ny];
		r2[i]=new double[ny];
		r3[i]=new double[ny];

        HBu[i]=new double[ny];
        HBv[i]=new double[ny];
        uu[i]=new double[ny];
        vv[i]=new double[ny];
        
        HBb[i]=new double[ny];
        HBc[i]=new double[ny];
        bb[i]=new double[ny];
        cc[i]=new double[ny];

        z[i]=new double[ny];

        check1[i]=new double[ny];
        check2[i]=new double[ny];
        check3[i]=new double[ny];

        /**
         * set initial state
         */
		for(int j=0; j<ny; j++) {
            pov[i][j]=centering_drop(0.007,100.0,i,j);
			u[i][j]=0;
			v[i][j]=0;

			pov_temp1[i][j]=0;
			u_temp1[i][j]=0;
			v_temp1[i][j]=0;

			pov_temp2[i][j]=0;
			u_temp2[i][j]=0;
			v_temp2[i][j]=0;

			r1[i][j]=0;
			r2[i][j]=0;
			r3[i][j]=0;

            z[i][j]=0;
		}
	}

    log=fopen("log","w");

    count_a1=0;
    count_a2=0;
    count_a3=0;
}

Alg::~Alg(){
    for(int i=0;i<nx;i++)
    {
        delete(HBu[i]);
        delete(HBv[i]);
        delete(uu[i]);
        delete(vv[i]);

        delete(HBb[i]);
        delete(HBc[i]);
        delete(bb[i]);
        delete(cc[i]);

        delete(r1[i]);
        delete(r2[i]);
        delete(r2[i]);

        delete(pov[i]);
        delete(u[i]);
        delete(v[i]);

        delete(pov_temp1[i]);
        delete(u_temp1[i]);
        delete(v_temp1[i]);

        delete(pov_temp2[i]);
        delete(u_temp2[i]);
        delete(v_temp2[i]);

        delete(z);

        delete(check1[i]);
        delete(check2[i]);
        delete(check3[i]);
    }
   
   fclose(log);
}

double Alg::find_coeff(int name)
{
	if(name==TAU) {
		return find(r1,r2,r3);
	}

	if(name==ALFA) {
		return find(z,z,z);
	}

	return 0;
}

void Alg::resolve_param(double coeff,int name)
{
    register size_t i,j;
	if(name==TAU) {
		resolve(coeff,r1,r2,r3);
	}

    if(name==ALFA) {
        for(i=0; i<nx; i++) {
            for(j=0; j<ny; j++) {
                //if(i<3||i>nx-4||j<3||j>ny-4){
                    z[i][j]=1;

                    resolve_r(0);
                    coeff = find(0,z,z,z);

                    pov[i][j]=pov[i][j]-coeff*z[i][j];

                    resolve_r(1);
                    coeff = find(1,z,z,z);

                    u[i][j]=u[i][j]-coeff*z[i][j];

                    resolve_r(2);
                    coeff = find(2,z,z,z);

                    v[i][j]=v[i][j]-coeff*z[i][j];

                    z[i][j]=0;
                //}
            }
        }
    }

}

void Alg::resolve(double coeff,double **a,double **b,double **c)
{
    register size_t i,j;
    #pragma omp parallel for collapse(2)
	for(i=0; i<nx; i++) {
		for(j=0; j<ny; j++) {
			pov[i][j]=pov[i][j]-coeff*a[i][j];
			u[i][j]=u[i][j]-coeff*b[i][j];
			v[i][j]=v[i][j]-coeff*c[i][j];
		}
	}
}

void Alg::set_temp(int number)
{
    register size_t i,j;
	if(number==1) {
        #pragma omp parallel for collapse(2)
		for(i=0; i<nx; i++) {
			for(j=0; j<ny; j++) {
                pov_temp1[i][j]=pov[i][j];
				u_temp1[i][j]=u[i][j];
				v_temp1[i][j]=v[i][j];
			}
		}
	}

    if(number==2) {
        memcpy(pov_temp2,pov,sizeof(pov));
        memcpy(u_temp2,u,sizeof(u));
        memcpy(v_temp2,v,sizeof(v));
    }

}

double Alg::get_norm_r()
{
	register double temp=0,hx_temp=hx,hy_temp=hy;
    register size_t i,j;

	#pragma omp parallel for private(i) reduction(+:temp)
    for(i=0; i<nx; i++) {
		for(j=0; j<ny; j++) {
            temp+=(r1[i][j]*r1[i][j]+r2[i][j]*r2[i][j]+r3[i][j]*r3[i][j]);
        }
	}

	temp=sqrt(temp*hx*hy);
	return temp;
}


void Alg::print(double **a,char *name)
{
	FILE *f;
	f=fopen(name,"w");

	for(int i=0; i<nx; i++) {
		for(int j=0; j<ny; j++) {
			fprintf(f,"%f ", a[i][j]);
		}
		//fprintf(f,"\n");
	}
	fclose(f);
}

void Alg::print(double **a)
{
    for(int i=0;i<nx;i++){
        for(int j=0;j<ny;j++){
            printf("%2.10lf ",a[i][j]);
        }
        printf("\n");
    }
}

int Alg::getIter()
{
    return iter;
}

double Alg::getMod()
{
    return mod;
}

int Alg::getCount()
{
    return count;
}

void *printState(void *alg)
{
    Alg *tmp=(Alg *)alg;
    int iter = tmp->getIter();
    int count = tmp->getCount();

    while(iter<count)
    {
        /**
         * Frequency of debug output
         */
        usleep(50000);
        printf("iter %d\n",tmp->getIter());
    }
}

void Alg::getH(char *file_name){
    FILE *f=fopen(file_name,"w");

    /**
     * output in file -H!
     */
    for(int i=0; i<nx; i++) {
        for(int j=0; j<ny; j++) {
            fprintf(f,"%2.8f %2.8f %2.8f\n", hx*i, hy*j, -H(i,j));
        }
        fprintf(f,"\n");
    }
    fprintf(f,"\n\n");
    fclose(f);
}

void Alg::solve(int out_n)
{
    mod=0;
    double tau=0;
    double alfa=0;
    iter=0;

    /**
     * create printing thread
     */
    //pthread_t printThread;
    //int ret=pthread_create(&printThread,NULL,printState,this);

    FILE *f,*f_u,*f_v;
    f=fopen(out_name,"w");

    int num;
    double tau_sum=0;
    clock_t t1,t2;
    t1=clock();

    //main time loop
    while(iter<=count){
        iter++;
        printf("iter %d\n",iter);

        set_temp(1);
        mod=1;
        num=0;
        while(mod>eps&&num<=200) {

            num++;
            set_temp(2);

            resolve_r();

            tau=find_coeff(TAU);
            resolve_param(tau,TAU);

            resolve_param(0,ALFA);

            mod=get_norm_r();
            printf("mod %lf\n",mod);
        }

        //output in file
        //check_result();
        if(iter % out_n == 0){
        double x,y;
        register size_t i,j;
        for(i=0, x=0.; i<nx; i++,x+=hx) {
            for(j=0, y=0.; j<ny; j++,y+=hy) {
                fprintf(f,"%2.8f %2.8f %2.8f\n", x, y,pov[i][j]);
            }
            fprintf(f,"\n");
        }
        fprintf(f,"\n\n");
        }
	}
    
    t2=clock();
    printf("time %lf\n",(double)(t2-t1)/CLOCKS_PER_SEC);

	fclose(f);
}
