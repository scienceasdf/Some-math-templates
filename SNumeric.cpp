#include<cmath>
#include"SNumeric.h"

/* Calculate the Gamma function.
   Programmed by Shen Weihong.
   Algorithm written by Shen Weihong( original creation).
*/
double Gamma(double s)
{
    double f=1;
    while(s>=8)
    {
        f*=--s;
    }
    while(s<1) {
        f/=s;
        s++;
    }
    f*=RombergInt([s](double x){return pow(x,s-1)*exp(-x);},1e-8,100,1e-6,10);
    return f;
}

