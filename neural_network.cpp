#include <cstdio>
#include <cstring>
#include <algorithm>
#include <vector>
using namespace std;
struct rec {
    double w1,w2,w3;
    rec(){}
    rec(double w1,double w2,double w3):w1(w1),w2(w2),w3(w3){}
    rec operator+(const rec &a)
    {
        return rec(w1+a.w1,w2+a.w2,w3+a.w3);
    }
    rec operator-(const rec &a)
    {
        return rec(w1-a.w1,w2-a.w2,w3-a.w3);
    }
    rec operator*(const double x)
    {
        return rec(w1*x,w2*x,w3*x);
    }
    double operator*(const rec &a)
    {
        return w1*a.w1+w2*a.w2+w3*a.w3;
    }
    void output()
    {
        printf("%.2f %.2f %.2f\n",w1,w2,w3);
    }
};
class NeuralNetwork {
private:
    vector<rec> vec;
    double bias,c;
    int sign(double x)
    {
        if(x>=0) return 1;
        else return -1;
    }
    int f(rec &a,rec &b)
    {
        return sign(a*b);
    }
public:
    NeuralNetwork(){}
    NeuralNetwork(double w1,double w2,double w3,double b,double c):bias(b),c(c)
    {
        vec.push_back(rec(w1,w2,w3));
    }
    void study()
    {
        int n;
        double x1,x2;
        int y;
        scanf("%d",&n);
        for(int i=1;i<=n;++i)
        {
            scanf("%lf%lf%d",&x1,&x2,&y);
            rec &w=vec[vec.size()-1];
            rec x=rec(x1,x2,bias);
            if(f(w,x)==y)
            {
                //w.output();
                continue;
            }
            rec ww=w+x*(double)(y-sign(w*x))*c;
            vec.push_back(ww);
            ww.output();
        }
        rec &a=vec[vec.size()-1];
        //printf("%.2f %.2f %.2f \n",a.w1,a.w2,a.w3);
    }
};
int main()
{
    freopen("neural.in","r",stdin);
    NeuralNetwork a = NeuralNetwork(0.75,0.5,-0.6,1,0.2);
    a.study();
}
