/*Chase Method, solving tridiagnol linear equations.
 */
double* root(double* a, double* b, double* c, double* d, int n)
{
    double* alpha=new double[n-1];
    double* beta;
    beta=new double[n];
    beta[0]=b[0];
    for(int i=0;i<=n-2;++i){
        alpha[i]=a[i]/beta[i];
        beta[i+1]=b[i+1]-alpha[i]*c[i];
    }
    double* x=new double[n];
    x[0]=d[0];
    for(int i=1;i<=n-1;++i) {x[i]=d[i]-alpha[i-1]*(x[i-1]); }
    x[n-1]=beta[n-1];
    for(int i=n-1;i>=0;--i) {x[i]=(x[i]-c[i]*x[i+1])/beta[i];}
    delete [] alpha;
    delete [] beta;
    return x;
}
