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
template<class F, class numArea> numArea secantMethod(F func, numArea x0, double eps)
{
    numArea x1=1.5*x0;
    numArea x,y0,y1;
    do
    {
        x=x1;
        x1=x1-func(x1)/(func(x1)-func(x0))*(x1-x0);
        x0=x;
    } while(abs(func(x1))>eps);
    return x1;
}

template<class F, class numArea> numArea secantMethod(F func(), numArea x0, numArea x1, double eps)
{
    numArea x,y0,y1;
    do
    {
        x=x1;
        x1=x1-func(x1)/(func(x1)-func(x0))*(x1-x0);
        x0=x;
    } while(abs(func(x1))>eps);
    return x1;
}

template<class F, class numArea> numArea secantMethod(F func, numArea x0, numArea x1, double eps, int maxIter)
{
    numArea x,y0,y1;
    do
    {
        x=x1;
        x1=x1-func(x1)/(func(x1)-func(x0))*(x1-x0);
        x0=x;
    } while(abs(func(x1))>eps && (--maxIter)>0);
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
numArea Lagrange(iter1 xIterFirst, iter1 xIterLast, iter2 yIterFirst, numArea varX)
{
    int num=xIterLast-xIterFirst-1;
    numArea l,f;
    f=0.0;
    for(int k=0;k<=num;k++)
    {
        l=1.0;
        for(int j=0;(j<=num);j++)
        {
            l*=(j==k)?1:(varX-xIterFirst[j])/(xIterFirst[k]-xIterFirst[j]);
        }
        f+=l*yIterFirst[k];
    }
    return f;
}



//Solving ODE:y'=func(x,y)
template<class BinOp, class numArea, class OutputIter>
numArea RungeKutta(BinOp func, numArea lb, numArea ub, numArea alpha, int nn, OutputIter res)
{
    numArea x0,y0,h,k1,k2,k3,k4,x1;
    x0=lb;
    y0=alpha;
    // res[0] should be y0, which will be more convenient.
    // Thus the last number which will be written is res[nn] instead of res[nn-1].
    *res=alpha;
    ++res;
    x1=lb;
    h=(ub-x0)/nn;
    int r=1;

    for(int n=1;n<(nn+1);n++)
    {
        k1=h*func(x0,y0);
        k2=h*func(x0+h/2.0,y0+k1/2.0);
        k3=h*func(x0+h/2.0,y0+k2/2.0);
        k4=h*func(x0+h,y0+k3);
        x1+=h;
        (*res)=y0+(k1+2.0*k2+2.0*k3+k4)/6.0;
        x0=x1;
        y0=*res;
        ++res;
    }
    return *res;
}

/* Solving two order ODE
 * y1=f(x,y1,y2), y2=g(x,y1,y2), given x0 ,y10, y20.
 * Specially, it can solve equations like y"=g(x,y,y'),
 * as y1=y, y2=y', f(x,y1,y2)=y2.
 */
template<class func1, class func2, class numArea, class outputIter1, class outputIter2>
void RungeKutta(func1 f, func2 g, numArea lb, numArea ub,
                 numArea init_y1, numArea init_y2, int nn, outputIter1 res_y1, outputIter2 res_y2)
{
    double m1,m2,m3,m4,k1,k2,k3,k4,y1,y2,x0,y0,z0,b,h,x1;
    x0=lb;
    h=(ub-x0)/nn;
    x1=x0;
    y1=init_y1;
    y2=init_y2;
    *res_y1=init_y1;
    *res_y2=init_y2;
    ++res_y1;
    ++res_y2;

    for(int n=1;n<(nn+1);n++)
    {
        k1=h*f(x1,y1,y2);
        m1=h*g(x1,y1,y2);
        k2=h*f(x1+.5*h,y1+k1/2,y2+m1/2);
        m2=h*g(x1+h/2,y1+k1/2,y2+m1/2);
        k3=h*f(x1+h/2,y1+k2/2,y2+m2/2);
        m3=h*g(x1+h/2,y1+k2/2,y2+m2/2);
        k4=h*f(x1+h,y1+k3,y2+m3);
        m4=h*g(x1+h,y1+k3,y2+m3);

        y1+=(k1+2*k2+2*k3+k4)/6;
        y2+=(m1+2*m2+2*m3+m4)/6;
        *res_y1=y1;
        *res_y2=y2;
        ++res_y1;
        ++res_y2;
        x1+=h;
    }
}
/* Solving ODE like y"=f(x,y,y'), given x1, y1, x2, y2.
 * Use Runge-Kutta method, shooting method and secant method.
 */
template<class func, class numArea, class outputIter>
numArea RungeKutta(func f, numArea x1, numArea x2, numArea y1, numArea y2, int nn, outputIter res)
{
    //
    secantMethod()
}

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

/* Calculate the convolution for two arrays.
 * Like polynomials' mutiply, it writes the results in a new array( or container).
 * Programmed by SHEN Weihong, algorithm written by SHEN weihong.
 */
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

//-----------------------------------------------------------------------------


/*  Calulates the fast Fourier transform based on 2.
 *  Programmed by SHEN Weihong.( original creation).
 *  Users can only use the class std::complex<double>.
 *  The result is the same with MATLAB, which is convenient for mixed programming.
 */
template<class iter, class outputIter>
inline void FFT(iter first, iter last, outputIter res)
{
    int n=std::distance(first,last);
    int N=n/2;
    const double PI=3.14159265358979323846;
    if(n!=2){

        std::complex<double>* temp1=new std::complex<double>[N];
        std::complex<double>* temp2=new std::complex<double>[N];
        std::complex<double>* out1=new std::complex<double>[N];
        std::complex<double>* out2=new std::complex<double>[N];
        for(int i=0;i<N;++i){
            temp1[i]=*first;
            ++first;
            temp2[i]=*first;
            ++first;
        }
        const std::complex<double> J(0,1);
        std::complex<double> w=exp(-2.0*PI*J/(double) n);
        std::complex<double> wk=1;
        if(n>=1024){         //If the number is too large, we can call one more thread. And the number can be changed.
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
        std::complex<double> y1=*first;
        ++first;
        std::complex<double> y2=*first;
        *res=(y1+y2);
        ++res;
        *res=(y1-y2);
    }
}

/*  Calulates the inverse fast Fourier transform based on 2.
 *  Programmed by SHEN Weihong.( original creation).
 *  The result is the same with MATLAB, which is convenient for mixed programming.
 */
template<class iter, class outputIter>
inline void IFFT(iter first, iter last, outputIter res)
{
    int n=std::distance(first,last);
    int N=n/2;
    double s=.5;
    const double PI=3.14159265358979323846;
    if(n!=2){
        std::complex<double>* temp1=new std::complex<double>[N];
        std::complex<double>* temp2=new std::complex<double>[N];
        std::complex<double>* out1=new std::complex<double>[N];
        std::complex<double>* out2=new std::complex<double>[N];
        for(int i=0;i<N;++i){
            temp1[i]=*first;
            ++first;
            temp2[i]=*first;
            ++first;
        }
        const std::complex<double> J(0,1);
        std::complex<double> w=exp(2.0*PI*J/(double) n);
        std::complex<double> wk=1.0;
        if(n>=1024){
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
        std::complex<double> y1=*first;
        ++first;
        std::complex<double> y2=*first;
        *res=s*(y1+y2);
        ++res;
        *res=s*(y1-y2);
    }
}


#endif // SNUMERIC_H_INCLUDED
