#include<algorithm>
#include<cmath>
#include<iostream>
#include<SNumeric.h>

void f2(double* x, double* res)
{
    res[0]=1.0;
    res[1]=x[0];
}

void f1(double* x, double* res)
{
    res[0]=1.0;
    res[1]=x[2];
    res[2]=x[1];
}

int main()
{
    double** pptr;
    double x0[3]={.0,1.0,1};
    double x00[2]={0,.0};

    RungeKutta(f1,x0,1.0,2,30,pptr);
    for(int i=0;i<31;++i){

        for(int k=0;k<3;++k){
            std::cout<<pptr[i][k]<<"\t";
        }
        std::cout<<"\n";
    }

    for(int i=0;i<31;++i){
        delete [] pptr[i];
    }
    delete [] pptr;

    RungeKutta(f2,x00,5.0,1,30,pptr);
    for(int i=0;i<31;++i){

        for(int k=0;k<2;++k){
            std::cout<<pptr[i][k]<<"\t";
        }
        std::cout<<"\n";
    }

    for(int i=0;i<31;++i){
        delete [] pptr[i];
    }
    return 0;
}
