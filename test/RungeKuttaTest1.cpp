#include<algorithm>
#include<cmath>
#include<iostream>
#include<SNumeric.h>

double* f2(double* x)
{
    double* res=new double[2];
    res[0]=1.0;
    res[1]=x[0];
}

int main()
{
    double** pptr;
    double x0[2]={.0,.0};
    //RungeKutta(f1,x0,8.0,2,500,pptr);
    RungeKutta(f2,x0,8.0,1,20,pptr);
    for(int i=0;i<21;++i){

        for(int k=0;k<2;++k){
            std::cout<<pptr[i][k]<<"\t";
        }
        std::cout<<"\n";
    }

    for(int i=0;i<21;++i){
        delete [] pptr[i];
    }
    delete [] pptr;
    return 0;
}
