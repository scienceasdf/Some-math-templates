#ifndef SPOLYNOMIAL_H
#define SPOLYNOMIAL_H

#include<valarray>

template <class iter, class numArea>
numArea polyCalc(iter first, iter last, numArea x)
{
    --last;
    numArea result=*last;
    while(first!=last){
        --last;
        result=result*x+*last;
    }
    return result;
}

template <class inputIter, class OutputIter>
void polyTan(inputIter first, inputIter last, OutputIter res)
{
    int i=1;
    while(first!=last){
        ++first;
        *res=i*(*first);
        ++res; ++i;
    }
    --res;
    *res=0;
}

template <class inputIter, class OutputIter>
void polyInt(inputIter first, inputIter last, OutputIter res)
{
    int i=1;
    *res=0;
    ++res;
    while(first!=last){
        *res=(*first)/i;
        ++res; ++first; ++i;
    }
}

template <class OutputIter>
void getBin(int n, OutputIter res)
{
  OutputIter next=res;
  while(n!=0){
    *next=n%2;
    n/=2;
    ++next;
  }
}

template<class numArea>
class polynomial{
public:
    std::valarray<numArea> coeArray;
    polynomial()
}
template<class T>   // Written without any tests.
class polynomial{
private:
    T* array;
    int degree;
public:
    polynomial() {array=new T; T=0;}
    polynomial(T val);
    polynomial(const T val, size_t n); // Construct a (n-1) degree polynomial, all coefficients of val.
    polynomial(const T* p, size_t n);   // Construct a (n-1) degree polynomial, coefficients from the T pointer starts and on.
    polynomial(const polynomial& poly); // Copy a polynomial.
    ~polynomial() {delete [] array;}

    T operator[] (size_t n) const {return array[n];}
    T& operator[] (size_t n) {return array[n];}

    polynomial& operator*=(const T arg) {std::for_each(begin(),end(),[arg] (T& s){s*=arg;})
    polynomial& operator/=(const T arg) {std::for_each(begin(),end(),[arg] (T& s){s/=arg;})

    int getDegree() {return degree;}
    void reDegree(size_t n) { delete [] array; array=new T[n+1];}   
    // Caution! This operation will delete all the former elements!
                                     
    T* begin() {return array;}
    T* end() {return array+degree+1;}

    T calc(T val) {return polyCalc(begin(),end(),val);}

};

template<class T>
polynomial<T>::polynomial(T val)
{
    array=new T;
    *array=val;
}

template<class T>
polynomial<T>::polynomial(const T val, size_t n)
{
    array=new T[n];
    for(int i=0;i<n;++i){
        array[n]=val;
    }
}

template<class T>
polynomial<T>::polynomial(const T* p, size_t n)
{
    array=new T[n];
    for(int i=0;i<n;++i){
        array[n]=p[n];
    }
}

template<class T>
polynomial<T>::polynomial(const polynomial& poly)
{
    int n=poly.getDegree()+1;
    array=new T[n];
    for(int i=0;i<n;++i){
        array[n]=poly[n];
    }
}

template<class T>
T polynomial<T>::operator [](size_t n)
{
    return array[n];
}

template<class T>
T& polynomial<T>::operator [](size_t n)
{
    return array[n];
}

#endif // SPOLYNOMIAL_H
