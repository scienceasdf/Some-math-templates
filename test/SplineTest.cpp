#include<vector>
#include<algorithm>
#include<iostream>

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


template<class T>
class spline{
public:
    std::vector<T> vecX;
    std::vector<T> vecY;
    spline() {}
    spline(T* xFirst, T* xLast, T* yFirst);
    /*spline(std::vector<T>::iterator xFirst,
           std::vector<T>::iterator xLast,
           std::vector<T>::iterator yFirst);*/
    ~spline(){delete [] m;}

    void addPoint(T x, T y) {vecX.push_back(x); vecY.push_back(y); sort(); isUpdated=false;}
    T getPos(T x);
private:
    bool isUpdated;
    T* m;
    void sort();
};

template<class T> T spline<T>::getPos(T x)
{
    int n=vecX.size()-1;
    if(isUpdated==false){
        T* alpha=new T[n+1];
        T* beta=new T[n+1];
        T* a=new T[n]; T* h=new T[n];
        T* b=new T[n+1]; T* c=new T[n];
        h[0]=vecX[1]-vecX[0];
        alpha[0]=1.0; alpha[n]=0.0;
        beta[0]=3.0/h[0]*(vecY[1]-vecY[0]);
        c[0]=1.0;
        b[0]=2.0; b[n]=2.0;

        for(int i=1; i<n; ++i){
            h[i]=vecX[i+1]-vecX[i];
            alpha[i]=h[i-1]/(h[i-1]+h[i]);
            beta[i]=3.0*((1.0-alpha[i])/h[i-1]*(vecY[i]-vecY[i-1])+alpha[i]/h[i]*(vecY[i+1]-vecY[i]));
            a[i-1]=1.0-alpha[i];
            c[i]=alpha[i];
            b[i]=2.0;

        }
        a[n-1]=1.0;
        beta[n]=3.0/h[n-1]*(vecY[n]-vecY[n-1]);

        m=root(a,b,c,beta,n+1);
        delete [] alpha;
        delete [] beta;
        delete [] a;
        delete [] b; delete [] c; delete [] h;
        isUpdated=true;
    }

    int i=n/2, min=0, max=n-1;
    if(x<=vecX[1]) i=0;
    else{
        if(x>=vecX[n-1]) i=n-1;
        else{
            while(((x-vecX[i])*(x-vecX[i+1]))>0){
                if((x-vecX[i])<0){
                    max=i;
                    i=(i+min)/2;
                }
                else{
                    min=i;
                    i=(i+max)/2;
                }

            }
        }
    }

    T s;
    s=(1.0+2.0*(x-vecX[i])/(vecX[i+1]-vecX[i]))*((x-vecX[i+1])/(vecX[i]-vecX[i+1]))*((x-vecX[i+1])/(vecX[i]-vecX[i+1]))*vecY[i];
    s+=(1.0+2.0*(x-vecX[i+1])/(vecX[i]-vecX[i+1]))*((x-vecX[i])/(vecX[i+1]-vecX[i]))*((x-vecX[i])/(vecX[i+1]-vecX[i]))*vecY[i+1];
    s+=(x-vecX[i])*((x-vecX[i+1])/(vecX[i]-vecX[i+1]))*((x-vecX[i+1])/(vecX[i]-vecX[i+1]))*m[i];
    s+=(x-vecX[i+1])*((x-vecX[i])/(vecX[i+1]-vecX[i]))*((x-vecX[i])/(vecX[i+1]-vecX[i]))*m[i+1];

    return s;
}

template<class T> void spline<T>::sort()
{
    const size_t n=vecX.size();

    for(int gap=n/2; 0<gap; gap/=2)
        for(int i=gap; i<n; i++)
            for(int j=i-gap; 0<=j; j-=gap)
                if(vecX[j+gap]<vecX[j]){
                    std::swap(vecX[j],vecX[j+gap]);
                    std::swap(vecY[j],vecY[j+gap]);
                }
}

int main()
{
    spline<double> s1;
    s1.addPoint(1.0,2.0);
    s1.addPoint(2.0,3);
    s1.addPoint(3,1);
    s1.addPoint(4,0);
    s1.addPoint(5,2);
    s1.addPoint(0,3);
    s1.addPoint(8,1);
    s1.addPoint(7,5);
    s1.addPoint(6,4);
    for(int i=0;i<s1.vecX.size();++i){
        std::cout<<s1.vecX[i]<<"\t"<<s1.vecY[i]<<"\n";
    }

    for(int i=0;i<55;++i){
        std::cout<<(double)i/5.0<<"\t"<<s1.getPos((double)i/5.0)<<"\n";
    }
    return 0;
}

