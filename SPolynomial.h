#ifndef SPOLYNOMIAL_H
#define SPOLYNOMIAL_H

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

#endif // SPOLYNOMIAL_H
