#include<algorithm>
#include<cmath>
#include<iostream>
#include<SNumeric.h>

double* f1(double* x)
{
    double* res=new double[3];
    res[0]=1.0;
    //atan2(double y, double x)
    //It's strange that here use atan2 will work wrong.
    double theta=atan((x[2]-500.0)/(x[1]-800.0*x[0]));
    res[1]=1100.0*cos(theta);
    res[2]=1100.0*sin(theta);
    return res;
}

int main()
{
    double** pptr;
    double x0[3]={.0,.0,-1000.0};
    RungeKutta(f1,x0,8.0,2,500,pptr);

    for(int i=0;i<501;++i){

        for(int k=0;k<3;++k){
            std::cout<<pptr[i][k]<<"\t";
        }
        std::cout<<"\n";
    }

    for(int i=0;i<501;++i){
        delete [] pptr[i];
    }
    delete [] pptr;
    return 0;
}
