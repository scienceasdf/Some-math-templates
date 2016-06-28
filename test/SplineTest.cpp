#include<vector>
#include<algorithm>

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

template<class T>
class spline{
public:
    std::vector<T> vecX;
    std::vector<T> vecY;
    spline(T* xFirst, T* xLast, T* yFirst);
    spline(std::vector<T>::iterator xFirst,
           std::vector<T>::iterator xLast,
           std::vector<T>::iterator yFirst);
    ~spline(){delete [] m;}

    void addPoint(T x, T y) {vecX.insert(x); vecY.insert(y); sort(); isUpdated=false;}
    T getPos(T x);
private:
    bool isUpdated;
    T* m;
    void sort();
};

template<class T> T spline<T>::getPos(T x)
{
    if(isUpdated==false){
        int n=vecX.size()-1;
        T* alpha=new T[n+1];
        T* beta=new T[n+1];
        T* a=new T[n]; T* h=new T[n];
        T* b=new T[n+1]; T* c=new T[n];
        h[0]=vecX[1]-vecX[0];
        alpha[0]=1.0; alpha[n]=.0;
        beta[0]=3.0/h[0]*(vecY[1]-vecY[0]);
        for(int i=0; i<n; ++i){
            h[i]=vecX[i+1]-vecX[i];
            alpha[i]=h[i-1]/(h[i-1]+h[i]);
            beta[i]=3.0*((1.0-alpha[i])/h[i-1]*(vecY[i]-vecY[i-1])+alpha[i]/h[i]*(vecY[i+1]-vecY[i]);
            a[i-1]=1.0-alpha[i];
            c[i]=alpha[i];
            b[i]=2.0;
        }
        a[n-1]=1.0;
        beta[n]=3.0/h[n-1]*(vecY[n]-vecY[n-1]);
        
        m=root(a,b,c,beta,n+1);
        delete [] alpha; delete [] beta; delete [] a;
        delete [] b; delete [] c;
    }
    int i=n/2, min=0, max=n;
    while((x-vecX[i])*(x-vecX[i+1])>0){
        if((x-vecX[i])<0){
            max=i;
            i=(i+min)/2;
        }
        else{
            min=i;
            i=(i+max)/2;
        }
    }
    
    T s;
    s=(1.0+2.0*(x-x[i])/(vecX[i+1]-vecX[i]))*((x-vecX[i+1])/(vecX[i]-vecX[i+1]))*((x-vecX[i+1])/(vecX[i]-vecX[i+1]))*vecY[i];
    s+=(1.0+2.0*(x-x[i+1])/(vecX[i]-vecX[i+1]))*((x-vecX[i])/(vecX[i+1]-vecX[i]))*((x-vecX[i])/(vecX[i+1]-vecX[i]))*vecY[i+1];
    s+=(x-vecX[i])*((x-vecX[i+1])/(vecX[i]-vecX[i+1]))*((x-vecX[i+1])/(vecX[i]-vecX[i+1]))*m[i];
    s+=(x-vecX[i+1])*((x-vecX[i])/(vecX[i+1]-vecX[i]))*((x-vecX[i])/(vecX[i+1]-vecX[i]))*m[i+1];
    
    return s;
}

template<class T> void spline<T>::sort()
{
    const size_t n=vecX.size();
    
    for(int gap=n/2; 0<gap; gap/2)
        for(int i=gap; i<n; i++)
            for(int j=i-gap; 0<=j; j-=gap)
                if(vecX[j+gap]<vecX[j]){
                    std::swap(vecX[j],vecX[j+gap]);
                    std::swap(vecY[j],vecY[j+gap]);
                }
}
