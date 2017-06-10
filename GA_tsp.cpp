#include <cstdio>
#include <cstring>
#include <algorithm>
#include <vector>
#include <cmath>
#include <ctime>

typedef long long ll;

struct Point{
    ll x,y;
    Point(){}
    Point(ll x,ll y):x(x),y(y){}
};

struct Unit{
    std::vector<int> gene;
    Unit(){}
    Unit(std::vector<int> &a):gene(a){}
};
class GAForTsp {
private:
    std::vector<Point> P;//点集

    std::vector<std::vector<ll> > dis;//任意两点间的距离

    double mutate_probability,across_probability;

    int N;//点的数量

    ll get_fit(Unit &u);//计算

    void create_initial_populaiton(Unit population[],int n);//生成初始群体

    void across(Unit population[],int n);//基因交叉

    void select(Unit population[],int n);

    void mutate(Unit population[],int n);

public:
    GAForTsp()
    {
        srand(time(0));
    }

    GAForTsp(double p1,double p2):across_probability(p1),mutate_probability(p2)
    {
        srand(time(0));
    }

    void input();//读入点集，并求出任意两点间的距离

    void evolve(int n,int m);//进化，n表示群体数量，m表示进化代数


};
ll GAForTsp::get_fit(Unit &u)
{
    ll fit=0;
    int x,y;
    for(int i=0;i<u.gene.size();++i)
    {
        x=u.gene[i];
        if(i==u.gene.size()-1) y=u.gene[0];
        else y=u.gene[i+1];
        fit+=dis[x][y];
//        printf("%d %d %lld  ???\n",x,y,dis[x][y]);
    }
//    printf("%lld ???\n",fit);
    return fit;
}
void GAForTsp::input()
{
    scanf("%d",&N);
    ll x,y,dx,dy,tij;
    double rij;
    for(int i=0;i<N;++i)
    {
        scanf("%lld%lld",&x,&y);
        P.push_back(Point(x,y));
    }
//    while(scanf("%lld%lld",&x,&y)!=EOF)
//        P.push_back(Point(x,y));
    for(int i=0;i<N;++i)
    {
        std::vector<ll> vec;
        for(int j=0;j<N;++j)
        {
            dx=P[i].x-P[j].x;
            dy=P[i].y-P[j].y;
            rij=sqrt((dx*dx+dy*dy)/10.0);
            tij=(ll)rij;
//            printf("%d %d %lld %lld %.3f ??\n",i,j,dx,dy,rij);
            if(abs(tij-rij)<1e-8) vec.push_back(tij);
            else vec.push_back(tij+1);
        }
        dis.push_back(vec);
    }
}

void output(Unit population[],int n)
{
    for(int i=0;i<n;++i)
    {
        for(int j=0;j<population[i].gene.size();++j)
            printf("%d ",population[i].gene[j]);
        printf("\n");
    }
}
void GAForTsp::evolve(int n,int m)
{
    Unit population[n];
    create_initial_populaiton(population,n);
    ll fit,min_fit=1e18;
    int cnt=0;
    while(m--)
    {
        across(population,n);
//        output(population,n);
//        printf("........\n");
        select(population,n);
        mutate(population,n);
//        output(population,n);
//        printf("--------------\n");
        min_fit=1e18;
        for(int i=0;i<n;++i)
        {
            fit=get_fit(population[i]);
            if(fit<min_fit) min_fit=fit;
        }
        ++cnt;
        if(cnt!=100) continue;
        cnt=0;
        printf("ans %d : %lld\n",m,min_fit);
    }
    min_fit=1e18;
    for(int i=0;i<n;++i)
    {
        fit=get_fit(population[i]);
        if(fit<min_fit) min_fit=fit;
    }
    printf("ans : %lld\n",min_fit);
}

void GAForTsp::create_initial_populaiton(Unit population[],int n)
{
    std::vector<int> permulation;
    for(int i=0;i<N;++i) permulation.push_back(i);
    for(int i=0;i<n;++i)
    {
        std::random_shuffle(permulation.begin(),permulation.end());
        population[i]=Unit(permulation);
    }
//    output(population,n);
//    printf("++++++++++++++\n");
}

void output_unit(Unit &u)
{
    for(int i=0;i<u.gene.size();++i) printf("%d ",u.gene[i]);
}
void GAForTsp::across(Unit populaiton[],int n)
{
    double p;
    ll min_fit=1e18,fit;
    int min_id,j,l,r,x;
    for(int i=0;i<n;++i)
    {
        fit=get_fit(populaiton[i]);
        if(fit<min_fit)
        {
            min_fit=fit;min_id=i;
        }
    }
    int fai[N],faj[N];
    for(int i=0;i<n;++i)
    {
        p=rand()/(RAND_MAX+1.0);
        j=i+1;
        if(j==n) j=0;
        if(i==min_id||j==min_id||p>across_probability) continue;//最优个体不进行基因交叉
        l=rand()%N;r=rand()%N;
        if(l>r) std::swap(l,r);
        Unit &ui=populaiton[i],&uj=populaiton[j];
//        output_unit(ui);printf("  (i)\n");
//        output_unit(uj);printf("  (j)\n");
        memset(fai,-1,sizeof(fai));
        memset(faj,-1,sizeof(faj));
        for(int k=l;k<=r;++k)
        {
            std::swap(ui.gene[k],uj.gene[k]);
            fai[ui.gene[k]]=k;
            faj[uj.gene[k]]=k;
        }
//        output_unit(ui);printf("  (i)swap\n");
//        output_unit(uj);printf("  (j)swap\n");
        for(int k=0;k<N;++k)
        {
            if(k==l)
            {
                k=r;continue;
            }
            x=ui.gene[k];
            while(fai[x]!=-1) x=uj.gene[fai[x]];
            ui.gene[k]=x;
            x=uj.gene[k];
            while(faj[x]!=-1) x=ui.gene[faj[x]];
            uj.gene[k]=x;
        }
//        output_unit(ui);
//        printf("!!!!!\n");
//        output_unit(uj);
//        printf("  ?????\n");
    }
}


void GAForTsp::select(Unit population[], int n)
{
    ll min_fit=1e18,fit;
    int min_id,id;
    double fit_value[n];
    for(int i=0;i<n;++i)
    {
        fit=get_fit(population[i]);
        fit_value[i]=fit;
        if(fit<min_fit)
        {
            min_fit=fit;min_id=i;
        }
    }
    std::swap(population[0],population[min_id]);
    std::swap(fit_value[0],fit_value[min_id]);
    for(int i=2;i<n;++i) fit_value[i]+=fit_value[i-1];
    for(int i=1;i<n;++i) fit_value[i]/=fit_value[n-1];
    double p;
    Unit a[n];a[0]=population[0];
    for(int i=1;i<n;++i)
    {
        p=rand()/(RAND_MAX+1.0);
        id=std::lower_bound(fit_value+1,fit_value+n,p)-fit_value;
        a[i]=population[i];
    }
    for(int i=1;i<n;++i) population[i]=a[i];
}

void GAForTsp::mutate(Unit population[],int n)
{
    double p;
    int x,y;
    for(int i=1;i<n;++i)
    {
        p=rand()/(RAND_MAX+1.0);
        if(p>mutate_probability) continue;
        x=rand()%N;y=rand()%N;
        std::swap(population[i].gene[x],population[i].gene[y]);
    }
}

int main()
{
    freopen("GA_tsp.out","w",stdout);
    GAForTsp g=GAForTsp(0.65,0.0001);
    g.input();
    g.evolve(1000,100000);
}
/*
4
0 0
0 10
10 0
10 10
*/
