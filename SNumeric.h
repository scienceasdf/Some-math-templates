/*  SWTL SNumeric.h
 *  Brief: Provides some basic numeric templates and functions.
 *  Copyright (c) SHEN WeiHong 2016.
 *  Version 1.000
 */

#ifndef SNUMERIC_H_INCLUDED
#define SNUMERIC_H_INCLUDED


/* Using the secant iteration method to solve equations.
   Programmed by Shen Weihong.
*/
template<class F> double secantMethod(F func, double x0, double eps)
{
    double x1=1.5*x0;
    double x,y0,y1;
    do
    {
        x=x1;
        x1=x1-func(x1)/(func(x1)-func(x0))*(x1-x0);
        x0=x;
    } while(fabs(func(x1))>eps);
    return x1;
}

template<class F> double secantMethod(F func(), double x0, double x1, double eps)
{
    double x,y0,y1;
    do
    {
        x=x1;
        x1=x1-func(x1)/(func(x1)-func(x0))*(x1-x0);
        x0=x;
    } while(fabs(func(x1))>eps);
    return x1;
}

template<class F> double secantMethod(F func, double x0, double x1, double eps, int maxIter)
{
    double x,y0,y1;
    do
    {
        x=x1;
        x1=x1-func(x1)/(func(x1)-func(x0))*(x1-x0);
        x0=x;
    } while(fabs(func(x1))>eps && (--maxIter)>0);
    return x1;
}

/* Using Romberg method to calculate intergration.
   Programmed by Shen Weihong.
*/

template<class func> double RombergInt(func f, double a, double b, double eps, long n)
{
    double sum=0, h=0;
    h=(b-a)/n;
    sum=.5*(f(a)+f(b));
    for(int k=1;k<n;k++) sum+=f(a+k*h);
    double *prevT, *T;
    T=new double[1];
    prevT=new double[1];
    T[0]=sum*h;
    int j=1;
    double pow4=1;
    do{
        delete [] prevT;
        prevT=T;
        T=new double[j+1];
        h*=.5;
        for(int k=1;k<=n;++k) sum+=f(a+(2*k-1)*h);
        T[j]=h*sum;
        pow4=1;
        for(int m=1;m<=j;++m){
            pow4*=4;
            T[j-m]=(pow4*T[j-m+1]-prevT[j-m])/(pow4-1.0);
        }
        j++;
        n*=2;
    } while(fabs(T[0]-prevT[0])>eps);
    sum=T[0];
    delete [] T;
    delete [] prevT;
    return sum;

}

/* Calculate the Gamma function.
   Programmed by Shen Weihong.
   Algorithm written by Shen Weihong( original creation).
*/
double Gamma(double s);

template<class iter1, class iter2, class numArea>
double Lagrange(iter1 xIterFirst, iter1 xIterLast, iter2 yIterFirst, numArea varX)




//Solving ODE
template<class func, class numArea, class OutputIter>
numArea RungeKutta(func f, numArea lb, numArea ub, int n, OutputIter dest)

/* Calculate the convolution of two functions at some range.
 * Programmed by SHEN Weihong.
 */
template<class func1, class func2>
double  convolution(func1 f, func2 g, double lb, double ub, double t)
{
    double result;
    result=RombergInt([t,f,g](double x) {return f(x)*g(t-x);},lb,ub,1e-6,10);
    return result;
}


template<class func1, class func2, class numArea>
numArea convolution(func1 f, func2 g, numArea lb, numArea ub, numArea t)
{
    numArea result;
    result=RombergInt([t](numArea x) {return f(x)*g(t-x);},lb,ub,1e-6,10);
}


/* Calculate the convolution for two arrays.
 * Programmed by SHEN Weihong, algorithm written by SHEN weihong.
 */
template<class iter1, class iter2>
double convolution(iter1 xIterFirst, iter1 xIterLast, iter2 yIterFirst, iter2 yIterLast, int n)
{
    double result=0;
    int an=std::distance(xIterFirst,xIterLast);
    int bn=std::distance(yIterFirst,yIterLast);
    for(int i=0;i<n;++i) ++yIterFirst;
    int a=0,b=n;
    for(int i=0;i<=n;++i){
        result+=(a<an && b<bn)?(*xIterFirst)*(*yIterFirst):0;
        ++xIterFirst; --yIterFirst;
        ++a; --b;
    }
    return result;
}


template<class iter1,class iter2, class Out>
Out convolution(iter1 xIterFirst, iter1 xIterLast, iter2 yIterFirst, iter2 yIterLast, Out res)
{
    int an=std::distance(xIterFirst,xIterLast);
    int bn=std::distance(yIterFirst,yIterLast);
    int m=an+bn-1;
    for(int i=0;i<m;++i){
        (*res++)=convolution(xIterFirst, xIterLast, yIterFirst, yIterLast, i);
    }
    return res;
}

template<class iter, class outputIter>
inline void FFT(iter first, iter last, outputIter res)
{
    int n=std::distance(first,last);
    int N=n/2;
    const double PI=3.14159265358979323846;
    if(n!=2){

        complex* temp1=new complex[N];
        complex* temp2=new complex[N];
        complex* out1=new complex[N];
        complex* out2=new complex[N];
        for(int i=0;i<N;++i){
            temp1[i]=*first;
            ++first;
            temp2[i]=*first;
            ++first;
        }
        const complex J(0,1);
        complex w=exp(-2.0*PI*J/(double) n);
        complex wk=1;
        if(n>=128){
            std::thread t2([temp2,out2,&N](){FFT(temp2,temp2+N,out2);});
            FFT(temp1,temp1+N,out1);
            delete [] temp1;
            t2.join();
            delete [] temp2;
        }

        else{
            FFT(temp1,temp1+N,out1);
            FFT(temp2,temp2+N,out2);
            delete [] temp1;
            delete [] temp2;
        }

        for(int k=0;k<N;k++){
            *res=(out1[k]+wk*out2[k]);
            wk*=w;
            ++res;
        }
        wk=1;
        for(int k=0;k<N;k++){
            *res=(out1[k]-wk*out2[k]);
            wk*=w;
            ++res;
        }
        delete [] out1; delete [] out2;
    }
    else{
        complex y1=*first;
        ++first;
        complex y2=*first;
        *res=(y1+y2);
        ++res;
        *res=(y1-y2);
    }
}

template<class iter, class outputIter>
inline void IFFT(iter first, iter last, outputIter res)
{
    int n=std::distance(first,last);
    int N=n/2;
    double s=.5;
    const double PI=3.14159265358979323846;
    if(n!=2){
        complex* temp1=new complex[N];
        complex* temp2=new complex[N];
        complex* out1=new complex[N];
        complex* out2=new complex[N];
        for(int i=0;i<N;++i){
            temp1[i]=*first;
            ++first;
            temp2[i]=*first;
            ++first;
        }
        const complex J(0,1);
        complex w=exp(2.0*PI*J/(double) n);
        complex wk=1.0;
        if(n>=128){
            std::thread t2([temp2,out2,&N](){IFFT(temp2,temp2+N,out2);});
            IFFT(temp1,temp1+N,out1);
            delete [] temp1;
            t2.join();
            delete [] temp2;
        }

        else{
            IFFT(temp1,temp1+N,out1);
            IFFT(temp2,temp2+N,out2);
            delete [] temp1;
            delete [] temp2;
        }

        for(int k=0;k<N;k++){
            *res=s*(out1[k]+wk*out2[k]);
            wk*=w;
            ++res;
        }
        wk=1.0;
        for(int k=0;k<N;k++){
            *res=s*(out1[k]-wk*out2[k]);
            wk*=w;
            ++res;
        }
        delete [] out1; delete [] out2;
    }
    else{
        complex y1=*first;
        ++first;
        complex y2=*first;
        *res=s*(y1+y2);
        ++res;
        *res=s*(y1-y2);
    }
}

#endif // SNUMERIC_H_INCLUDED
