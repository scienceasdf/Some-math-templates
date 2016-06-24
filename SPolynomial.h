#ifndef SPOLYNOMIAL_H
#define SPOLYNOMIAL_H

#include<valarray>
#include"SNumeric.h"
#include<algorithm>

//Calcultes the value of a polynomial on some number.
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

//Calulate the tangent of a polynomial,
// writing the coefficients into a container.
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

//Calulate the intergration of a polynomial,
// writing the coefficients into a container.
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


/* A function to write the binary representation of a numbeer into a container.
 * Programmed by SHEN Weihong.
 */
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

template<class T>   // Written without any tests.
class polynomial{
private:
    T* array;
    int degree;
    size_t size;
    void reDegree(size_t n){
        polynomial temp(this);
        delete [] array;
        array=new T[n];
        copy(temp.begin(),temp.end(),begin());
        for(int i=size;i<n;++i){
            array[i]=0;
        }
        size=n;
    }
public:
    polynomial():degree(0),size(1) {array=new T; T=0;}
    polynomial(T val);
    polynomial(size_t n,const T val); // Construct a (n-1) degree polynomial, all coefficients of val.
    polynomial(const T* p, size_t n);   // Construct a (n-1) degree polynomial, coefficients from the T pointer starts and on.
    polynomial(const polynomial& poly); // Copy a polynomial.
    ~polynomial() {delete [] array;}

    T operator[] (size_t n) const;
    T& operator[] (size_t n);

    polynomial& operator*=(const T arg) {std::for_each(begin(),end(),[arg] (T& s){s*=arg;})
    polynomial& operator/=(const T arg) {std::for_each(begin(),end(),[arg] (T& s){s/=arg;})
    polynomial& operator+=(const polynomial<T>& poly);
    polynomial& operator-=(const polynomial<T>& poly);
    polynomial& operator*=(const polynomial<T>& poly);
    
    polynomial& operator/=(const polynomial<T>& poly);
    polynomial& operator%=(const polynomial<T>& poly);
    int getDegree() {return degree;}
    size_t getSize() {return size;}

    T* begin() {return array;}
    T* end() {return array+size;}

    T calc(T val) {return polyCalc(begin(),end(),val);}

};

template<class T>
polynomial<T>::polynomial(T val)
{
    array=new T;
    *array=val;
}

template<class T>
polynomial<T>::polynomial(size_t n,const T val)
{
    array=new T[n];
    for(int i=0;i<n;++i){
        array[n]=val;
    }
    size=n;
    degree=n-1;     // Here should be changed if val==0;
}

template<class T>
polynomial<T>::polynomial(const T* p, size_t n)
{
    array=new T[n];
    for(int i=0;i<n;++i){
        array[n]=p[n];
    }
    size=n;
    degree=n-1;     // Here should be changed if T[n-1]==0 or like this.
}

template<class T>
polynomial<T>::polynomial(const polynomial& poly)
{
    int n=poly.getDegree()+1;
    array=new T[n];
    for(int i=0;i<n;++i){
        array[n]=poly[n];
    }
    size=poly.getSize();
    degree=poly.getDegree();
}

template<class T>
T polynomial<T>::operator [](size_t n)
{
    if(n<size()){
        return array[n];
    }
    else{
        reDegree(n+1);
        return 0;
    }
}

template<class T>
T& polynomial<T>::operator [](size_t n)
{
    if(n<size()){
        return array[n];
    }
    else{
        reDegree(n+1);
        return array[n];
    }
}

template<class T>
polynomial<T> operator *(polynomial<T>& lhs, polynomial<T>& rhs)
{
    int n=lhs.getDegree()+rhs.getDegree()+1;
    polynomial<T> temp(n,1.0);
    convolution(lhs.begin(),lhs.end(),rhs.begin(),rhs.end(),temp.begin());
    return temp;
}

template<class T>
polynomial<T> operator+(polynomial<T>& lhs, polynomial<T>& rhs)
{
    int n=std::max(lhs.getDegree,rhs.getDegree);
    polynomial temp(n,1.0);
    
    
}

template<class T>
polynomial<T> operator-(polynomial<T>& lhs, polynomial<T>& rhs)
{
    int n=std::max(lhs.getDegree,rhs.getDegree);
    polynomial temp(n,1.0);
    
    
}

#endif // SPOLYNOMIAL_H
