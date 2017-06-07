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

    p_fun f, fit;

    int get_binary_digit(int x,int k);

    void output_binary(int x);

    double decode(int code);

    void create_initial_population(int population[]);

    void across(int population[]);

    void mutate(int population[]);

    void select(int population[]);

	void across2(int &p1, int &p2);
	
public:
    Genetic(){}
    Genetic(double left, double right, double eps, int n, p_fun f, p_fun fit, double p1, double p2):
        left(left), right(right), eps(eps), n(n), f(f), fit(fit), across_probability(p1), mutate_probobility(p2)
    {
        double length = log((right - left) / eps) / log(2) ;
        code_length = ceil(length);
        if(std::abs(length-code_length)<1e-8) ++code_length;
        actual_eps = (right - left) / ( ( 1 << code_length ) - 1);
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

void Genetic::across2(int &p1, int &p2)
{
	int k = rand()%code_length + 1;
    int x = (1<<k) - 1;
    int y = (1<<code_length) - 1;
    //output_binary(p1);printf(" p1\n");
    //output_binary(p2);printf(" p2\n");
    y ^= x;
    p1 &= y;
    p1 |= p2&x;
    //output_binary(p1);printf(" ++p1\n");
    //output_binary(p2);printf(" ++p2\n");
}

void Genetic::across(int population[])
{
    int j, k, t;
    double x, y, p;
    for(int i=0;i<n;++i)
    {
    	j = i + 1;
        if( j == n ) j = 0;
    	x=decode(population[i]);y=decode(population[j]);
        if(fit(x)>fit(y)) std::swap(population[i],population[j]);
        p = rand() / (RAND_MAX + 1.0);
        if(p>across_probability) continue;
        across2(population[i],population[j]);
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
    double p, x, f_value, mx_f=-1e12, mx_id;
    for(int i=0;i<n;++i)
    {
    	x=decode(population[i]);
    	f_value=f(x);
    	if(f_value>mx_f)
    	{
    		mx_f=f_value;
    		mx_id=i;
		}
	}
    for(int i=0;i<n;++i)
    {
    	if(i==mx_id) continue;
        p = rand() / (RAND_MAX + 1.0);
        if(p>mutate_probobility) continue;
        int k = rand()%code_length;
        population[i] ^= 1<<k;
        x = decode(population[i]);
        f_value = f(x);
    }
}

double Genetic::evolve(int m)
{
    int population[n];
    double f_value[n], fit_value[n];
    double x, p;
    create_initial_population(population);
    int cnt=0;
	while(m--)
    {
        across(population);
        //select(population);
        //mutate(population);
        for(int i=0;i<code_length;++i)
        {
        	bool flag=true;
        	for(int j=1;j<n;++j)
				if(((population[j]^population[j-1])>>i)&1)
				{
        			flag=false;break;
				}
			if(flag)
			{
				int k=rand()%n;
				int z=0;
				double f_max=fit(decode(population[z])),ff;
				for(int j=1;j<n;++j)
				{
					ff=fit(decode(population[j]));
					if(ff>f_max)
					{
						f_max=ff;z=j;
					}
				}
				z+=1+rand()%(n-1);
				if(z>n) z-=n;
				population[z]^=1<<i;
			}
		}
		/*++cnt;
        printf("(%d)\n",cnt);
        for(int i=0;i<n;++i)
        {
        	double x=decode(population[i]);
        	printf("%d ",i);
        	printf(" %.3f ",x); 
        	output_binary(population[i]);
        	printf(" %.8f",fit(x));
        	printf(" %.8f\n",f(x));
		}
		printf("--------------------\n");*/
    }
    double mx_x, mx_f=-1e12, ff;
    for(int i=0;i<n;++i)
    {
    	x=decode(population[i]);
    	ff=f(x);
    	if(ff>mx_f)
    	{
    		mx_f=ff;
    		mx_x=x;
		}
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
    //return x*x/10.0;
}

double fit(double x)
{
    return f(x);
}
int main()
{
    //printf("%.2f  ??\n",acos(-1)/2.0);
    Genetic a = Genetic(0, acos(-1)/4.0, 0.01, 20, f, fit, 0.7, 0.001);
//    printf("%.8f  !!!\n",f(50000));
    //Genetic a = Genetic(0, 100000, 0.01, 100, f, fit, 0.65, 0.001);
    printf("%.8f \n", f(a.evolve(100)));
    printf("%.8f  !!!\n",f(acos(-1)/4.0));
    //a.debug(1);
    return 0;
}
