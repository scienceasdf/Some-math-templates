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

template<class T>   
class polynomial{
private:
    T* array;
    size_t degree;
    size_t size;
    bool degUpdated;
    void reDegree(size_t n){
        int m=getDegree();
        T* temp=new T[m+1];
        std::copy(begin(),end(),temp);
        delete [] array;
        array=new T[n+1];   // Degree of n, and size of (n+1).
        std::copy(temp,temp+m+1,begin());
        for(int i=size;i<n+1;++i){
            array[i]=0;
        }
        size=n+1;
        degUpdated=false;
        delete temp;
    }
public:
    polynomial():degree(0),size(1),degUpdated(true) {array=new T; *array=0;}
    polynomial(T val);
    polynomial(size_t n,const T val); // Construct a (n-1) degree polynomial, all coefficients of val.
    polynomial(const T* p, size_t n);   // Construct a (n-1) degree polynomial, coefficients from the T pointer starts and on.
    polynomial(polynomial& poly); // Copy a polynomial.
    ~polynomial() {delete [] array;}

    T operator[] (size_t n) const;
    T& operator[] (size_t n);

    polynomial<T>& operator*=(const T arg) {std::for_each(begin(),end(),[arg] (T& s){s*=arg;});}
    polynomial<T>& operator/=(const T arg) {std::for_each(begin(),end(),[arg] (T& s){s/=arg;});}
    polynomial<T>& operator+=(const polynomial<T>& poly);
    polynomial<T>& operator-=(const polynomial<T>& poly);
    polynomial<T>& operator*=(const polynomial<T>& poly);

    polynomial& operator/=(const polynomial<T>& poly);
    polynomial& operator%=(const polynomial<T>& poly);
    size_t getDegree();
    size_t getSize() {return size;}

    T* begin() {return array;}
    T* end() {return array+size;}

    T calc(T val) {return polyCalc(begin(),end(),val);}

};

template<class T>
size_t polynomial<T>::getDegree()   // Checked
{
    if(degUpdated==true){
        return degree;
    }
    else{
        T* p=end(); int i=0;
        do{
            --p;
            ++i;
        }while(*p==(T)0);
        i=size-i;
        degUpdated=true;
        return i;
    }
}

template<class T>       // Checked
polynomial<T>::polynomial(T val)
{
    array=new T;
    *array=val;
    size=1; degree=0;
    degUpdated=(val!=0);
}

template<class T>
polynomial<T>::polynomial(size_t n, const T val)    // Checked
{
    array=new T[n];
    for(int i=0;i<n;++i){
        array[i]=val;
    }
    size=n;
    degree=n-1;     // Here should be changed if val==0;
    degUpdated=(val!=(T)0);
}

template<class T>   // Checked
polynomial<T>::polynomial(const T* p, size_t n)
{
    array=new T[n];
    for(int i=0;i<n;++i){
        array[i]=p[i];
    }
    size=n;
    degree=n-1;     // Here should be changed if T[n-1]==0 or like this.
    degUpdated=false;
}


template<class T>
polynomial<T>::polynomial(polynomial& poly)     // Checked
{
    size_t n=poly.getDegree()+1;
    array=new T[n];
    for(int i=0;i<n;++i){
        array[i]=poly[i];
    }
    size=n;
    degree=n-1;
}

template<class T>
T polynomial<T>::operator [](size_t n) const    // Checked
{
    if(n<size()){
        return array[n];
    }
    else{
        reDegree(n);
        degUpdated=false;
        return 0;
    }
}

template<class T>
T& polynomial<T>::operator [](size_t n)     // Checked
{
    degUpdated=false;
    if(n<size){
        return array[n];
    }
    else{
        reDegree(n);
        return array[n];
    }
}

template<class T>
polynomial<T>& polynomial<T>::operator+=(const polynomial& rhs)     // Checked
{
    int p=degree;
    int q=rhs.degree;
    if(p>q){
        for(int i=0;i<=q;++i){
            array[i]+=rhs.array[i];
        }
    }
    if(p<q){
        reDegree(q);
        for(int i=0;i<=q;++i){
            array[i]+=rhs.array[i];
        }
    }
    if(p==q){
        for(int i=0;i<=q;++i){
            array[i]+=rhs.array[i];
        }
        degUpdated=false;   // It's possible that the degree will decrease.
    }
    return *this;
}

template<class T>
polynomial<T>& polynomial<T>::operator-=(const polynomial& rhs)     // Checked
{
    int p=degree;
    int q=rhs.degree;
    if(p>q){
        for(int i=0;i<=q;++i){
            array[i]-=rhs.array[i];
        }
    }
    if(p<q){
        reDegree(q);
        for(int i=0;i<=q;++i){
            array[i]-=rhs.array[i];
        }
    }
    if(p==q){
        for(int i=0;i<=q;++i){
            array[i]-=rhs.array[i];
        }
        degUpdated=false;   // It's possible that the degree will decrease.
    }
    return *this;
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
