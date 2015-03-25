#include "comparer.h"
#include <stdio.h>


Comparer::Comparer(int _nx,int _ny, double _hx,double _hy,int _time):nx(_nx),ny(_ny),hx(_hx),hy(_hy),time(_time)
{
}

void Comparer::delta(char* file1,char* file2,char* resultDelta){
    FILE *f1=fopen(file1,"r");
    FILE *f2=fopen(file2,"r");
    FILE *f3=fopen(resultDelta,"w");
    
    if(f1==NULL||f2==NULL||f3==NULL)
        return;

    double **pov=new double*[nx];
    for(int i=0;i<nx;i++)
    {
        pov[i]=new double[ny];
    }

    double null=0;
    double temp=0;

    for(int iter=0;iter<time;iter++)
    {
        printf("iter for delta %d\n",iter);
        for(int i=0; i<nx; i++) {
            for(int j=0; j<ny; j++) {
                fscanf(f1,"%lf %lf %lf",&null,&null, &pov[i][j]);
            }
            //fscanf(f,"\n");
        }
        //fscanf(f,"\n\n");
        //fclose(f);

        for(int i=0; i<nx; i++) {
            for(int j=0; j<ny; j++) {
                temp=0;
                fscanf(f2,"%lf %lf %lf",&null,&null, &temp);
                pov[i][j]-=temp;
            }
            //fscanf(f,"\n");
        }
        //fscanf(f,"\n\n");
        //fclose(f);

        //output in file
        for(int i=0; i<nx; i++) {
            for(int j=0; j<ny; j++) {
                fprintf(f3,"%2.8f %2.8f %2.8f\n", hx*i, hy*j, pov[i][j]);
            }
            fprintf(f3,"\n");
        }
        fprintf(f3,"\n\n");
        //fclose(f);
    }

    fclose(f1);
    fclose(f2);
    fclose(f3);
}
