#include<cmath>
#include"SNumeric.h"
#include"cstdlib"

/*  Calculate Gaussian distribution.
    Programmed by Shen Weihong.
    Algorithm written by Shen Weihong( original creation, by infinite series).
*/
double Gauss(double x)
{
    if(x<)
    long p=1;
    double s=0;
    double xPow=x;
    for(int i=1;i<10;i++){
        s+=xPow/(double) (p*(2*i-1));
        p*=(-1)*2*i;
        xPow*=x*x;
    }
    s*=0.39894228;
    s+=.5;
    return s;
}


/* A simple Gaussian random number generator.
 * Codes from the Internet.
 * Algorithm by Marsaglia & Bray, 1964.
 */

double GaussRand()
{
    static double V1, V2, S;
    static int phase = 0;
    double X;
    if ( phase == 0 ) {
        do {
            double U1 = (double)rand() / RAND_MAX;
            double U2 = (double)rand() / RAND_MAX;
            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while(S >= 1 || S == 0);
        X = V1 * sqrt(-2 * log(S) / S);
    }
    else {X = V2 * sqrt(-2 * log(S) / S);}
    phase = 1 - phase;
    return X;
}

double ChiSquareDistribution(int freeDegree, double x)
{
    return 1-(RombergInt([freeDegree](double v)
    {return pow(v,freeDegree/2.0-1)*exp(-.5*v)/pow(2,.5*freeDegree)/Gamma(.5*freeDegree); },
    0,x,1e-5,10));
}

double MTCL_ChiSquareDistribution(int freeDegree, double x)
{
    double sum=0,r;
    int p=0;
    for(int i=0;i<20000;i++){
        for(int i=0;i<freeDegree;++i){
            r=GaussRand();
            sum+=r*r;
        }
        if(sum>x) p++;
        sum=0;
    }
    return (double(p/20000.0));
}
