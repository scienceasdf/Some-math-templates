#ifndef SLINEAR_H
#define SLINEAR_H
/*Chase Method, solving tridiagnol linear equations.
 */
double* root(double* a, double* b, double* c, double* d, int n);

template<class aIter, class bIter, class cIter, class OutputIter>
void chaseMethod(aIter aIterFirst, aIter aIterLast, bIter bIterFirst,
                 cIter cIterFirst, OutputIter dest)


#endif // SLINEAR_H
