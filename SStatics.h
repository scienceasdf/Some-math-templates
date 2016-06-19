/*  SWTL SStatics.h
 *  Brief: Provides some basic statics templates and functions.
 *  Copyright (c) SHEN WeiHong 2016.
 *  Version 1.000
 */


#ifndef STATICS_H_INCLUDED
#define STATICS_H_INCLUDED

#include<cmath>

#include"SNumeric.h"

template<class Iter1, class Iter2> class LeastSquaredMethod
{
public:
    LeastSquaredMethod(Iter1 xIterFirst, Iter1 xIterLast, Iter2 yIterFirst){
            Iter2 yIter=yIterFirst;
            averX=0;
            averY=0;
            xSD=0;
            ySD=0;
            double sumXTY=0, sumXSQ=0, sumYSQ=0;
            N=0;
            double sumDxDy=0;

            for(Iter1 xIter=xIterFirst;xIter!=xIterLast;++xIter){
                averX+=*xIter;
                averY+=*yIter;
                N++;
                sumXTY+=(*xIter)*(*yIter);
                sumXSQ+=(*xIter)*(*xIter);
                sumYSQ+=(*yIter)*(*yIter);
                ++yIter;
            }
            coeA=(N*sumXTY-averX*averY)/(-averX*averX+N*sumXSQ);

            averX/=(double)N;
            averY/=(double)N;

            yIter=yIterFirst;
            xSD=pow(sumXSQ/N-averX*averX,.5);
            ySD=pow(sumYSQ/N-averY*averY,.5);
            coeCorre=(sumXTY/N-averX*averY)/xSD/ySD;
            xSD*=pow(N/(N-1),.5);
            xSD*=pow(N/(N-1),.5);

    }

    double coeCorre;
    double averX;
    double averY;
    double coeA;
    double coeB;
    double xSD;
    double ySD;
    int N;
};

double Gauss(double x);
double GaussRand();
double ChiSquareDistribution(int freeDegree, double x);
double MTCL_ChiSquareDistribution(int freeDegree, double x);

#endif // STATICS_H_INCLUDED
