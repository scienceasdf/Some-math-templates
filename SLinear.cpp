/*Chase Method, solving tridiagnol linear equations.
 */
double* root(double* a, double* b, double* c, double* d, int n)
{
    if(n==1){
        double* x=new double[1];
        x[0]=d[0]/b[0];
        return x;
    }
    else{
        double* alpha=new double[n-1];
        double* beta;
        beta=new double[n];
        beta[0]=b[0];
        for(int i=0;i<=n-2;++i){
            alpha[i]=a[i]/beta[i];
            beta[i+1]=b[i+1]-alpha[i]*(c[i]);
        }
        double* x=new double[n],*y=new double[n];
        y[0]=d[0];
        for(int i=1;i<=n-1;++i) y[i]=d[i]-alpha[i-1]*(y[i-1]);
        x[n-1]=beta[n-1];
        for(int i=n-1;i>=0;--i) x[i]=(y[i]-c[i]*x[i+1])/beta[i];
        delete [] y;
        delete [] alpha;
        delete [] beta;
        return x;
    }
}

