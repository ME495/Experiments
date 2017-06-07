#include <cstdio>
#include <cstring>
#include <algorithm>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <array>
typedef double (*p_fun)(double);
class Genetic {
private:
    double left, right;//左区间，右区间

    double eps, actual_eps;//要求的精度，实际的精度

    int code_length;//编码长度

    int n;//选取的群体中个体数量

    double across_probability, mutate_probobility;//交叉概率

    double mx_x, mx;

    p_fun f, fit;

    int get_binary_digit(int x,int k);

    void output_binary(int x);

    double decode(int code);

    void create_initial_population(int population[]);

    void across(int population[]);

    void mutate(int population[]);

    void select(int population[]);

public:
    Genetic(){}
    Genetic(double left, double right, double eps, int n, p_fun f, p_fun fit, double p1, double p2):
        left(left), right(right), eps(eps), n(n), f(f), fit(fit), across_probability(p1), mutate_probobility(p2)
    {
        double length = log((right - left) / eps) / log(2) ;
        code_length = ceil(length);
        if(std::abs(length-code_length)<1e-8) ++code_length;
        actual_eps = (right - left) / ( ( 1 << code_length ) - 1);
        mx = -1e9;
        srand(time(0));
    }

    void debug(int code);

    double evolve(int m);

};

int Genetic::get_binary_digit(int x, int k)
{
    return (x >> k) & 1;
}

void Genetic::output_binary(int x)
{
    for(int i=code_length-1;i>=0;--i)
        printf("%d",(x>>i)&1);
}

double Genetic::decode(int code)
{
    int code1 = code&1, b1, b2;
    for(int i=1; i<code_length; ++i)
    {
        b1 = get_binary_digit(code, i);
        b2 = get_binary_digit(code1, i-1);
        code1 |= (b1 ^ b2) << i;
    }
    return left + actual_eps * code1;
}

void Genetic::create_initial_population(int population[])
{
    for(int i=0;i<n;++i)
        population[i]=rand()%(1<<code_length);
}

void Genetic::across(int population[])
{
    int j, k, x, y, t;
    double p;
    for(int i=0;i<n;++i)
    {
        p = rand() / (RAND_MAX + 1.0);
        if(p>across_probability) continue;
        j = i + 1;
        if( j == n ) j = 0;
        k = rand()%code_length + 1;
        x = (1<<k) - 1;
        y = (1<<code_length) - 1;
        y ^= x;
        t = population[i]&x;
        population[i] &= y;
        population[i] |= population[j]&x;
        population[j] &= y;
        population[j] |= t;
    }
}

void Genetic::select(int population[])
{
    int a[n];
    double f_value[n], fit_value[n];
    double x, p;
    for(int i=0;i<n;++i)
    {
        x = decode(population[i]);
        f_value[i] = f(x);
        fit_value[i] = fit(x);
        if(i) fit_value[i]+=fit_value[i-1];
        if(f_value[i]>mx)
        {
            mx=f_value[i];mx_x=x;
        }
    }
    for(int i=0;i<n;++i) fit_value[i] /= fit_value[n-1];
    for(int i=0;i<n;++i)
    {
        p = rand() / (RAND_MAX + 1.0);
        int k = std::lower_bound(fit_value, fit_value+n, p) - fit_value;
        a[i] = population[k];
    }
    for(int i=0;i<n;++i) population[i] = a[i];
}

void Genetic::mutate(int population[])
{
    double p, x, f_value;
    for(int i=0;i<n;++i)
    {
        p = rand() / (RAND_MAX + 1.0);
        if(p>mutate_probobility) continue;
        int k = rand()%code_length;
        population[i] ^= 1<<k;
        x = decode(population[i]);
        f_value = f(x);
        if(f_value>mx)
        {
            mx=f_value;mx_x=x;
        }
    }
}

double Genetic::evolve(int m)
{
    int population[n];
    double f_value[n], fit_value[n];
    double x, p;
    create_initial_population(population);
    while(m--)
    {
        across(population);
        select(population);
        mutate(population);
    }
    return mx_x;
}

void Genetic::debug(int code)
{
    printf("%.2f  ??\n",f(decode(code)));
    decode(code);
}

double f(double x)
{
    return x * cos(x);
//    return x*(100000-x);
}

double fit(double x)
{
    return f(x);
}
int main()
{
    //printf("%.2f  ??\n",acos(-1)/2.0);
    printf("%.8f  !!!\n",f(acos(-1)/4.0));
    Genetic a = Genetic(0, acos(-1)/4.0, 0.01, 30, f, fit, 0.65, 0.001);
//    printf("%.8f  !!!\n",f(50000));
//    Genetic a = Genetic(0, 100000, 0.01, 100, f, fit, 0.65, 0.001);
    printf("%.8f \n", f(a.evolve(100)));
    //a.debug(1);
    return 0;
}
