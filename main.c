#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double nevyazka(int f_num, double h, double* next, int M);
void solve(double *next, double* right, double h, double tau, int M);
double func_u(int num, double x, double t);
double func_f(int num, double x, double t);


int main(void)
{

    int N, M;
    int f_num, n_num, m_num;


    double *curr, *next;
    double *right;
    double h, tau;


    double t=0;
    double x=0;

    int step;
    int i;


    int test_mode = 0;

    double eps = 0;

    for(f_num = 1; f_num <= 4; ++f_num)
    {
        printf("Nomer funkcii - %d:****************************************\n", f_num);
        for(N = 10; N<=1000; N*=10)
        {
            for(M = 10; M<=1000; M*=10)
            {


                h = 1./M;
                tau = 1./N;
                printf("N = %d, M = %d, ", N, M);

                curr = (double*)malloc(sizeof(double)*(M+1));
                next = (double*)malloc(sizeof(double)*(M+1));

                for(i=0; i<M+1; ++i)
                {
                    curr[i]=0;
                    next[i]=0;
                }
                x = h;
                for(i=1; i<M; ++i)
                {
                    curr[i] = func_u(f_num, x, 0);
                    x+=h;
                }

                for(step=0; step<N; ++step)
                {
                    t = step*tau;
                    right = (double*)malloc(sizeof(double)*(M-1));

                    for(i=0; i<M-1; ++i)
                    {
                        right[i]=0;
                    }


                    x=h;
                    for(i=0; i<M-1; ++i)
                    {
                        right[i] = func_f(f_num, x, t+tau) + 1./(6*tau)*(curr[i+2] + 4*curr[i+1] + curr[i]);
                        x+=h;
                    }

                    solve(next, right, h, tau, M);

                    for(i=0; i<M+1; ++i)
                    {
                        curr[i]=next[i];
                    }
                }

                eps = nevyazka(f_num, h, next, M);
                printf("eps = %0.3e\n", eps);


            }
            printf("\n");
        }
        printf("\n");
    }




    return 0;
}


void solve(double *ans, double* right, double h, double tau, int N)
{
    double *ksi,*eta,*y;
    double a,b,c;
    a = 1./(h*h) - 1./(6*tau);
    c = 2./(h*h)+ 2./(3.*tau);
    b = 1./(h*h) - 1./(6*tau);

    ksi = (double*)malloc(sizeof(double)*(N-1));
    eta = (double*)malloc(sizeof(double)*(N-1));
    y  = (double*)malloc(sizeof(double)*(N-1));

    for(int i=0; i<(N-2); ++i)
    {
        ksi[i]=0;
        eta[i]=0;
    }

    ksi[N-2]=a/c;
    eta[N-2]=right[N-2]/c;

    for( int i=N-3; i>=1; --i)
    {
        ksi[i] = a/(c-b*ksi[i+1]);
    }

    for( int i=N-3; i>=0; --i)
    {
        eta[i] = (right[i] +b*eta[i+1]) / (c-b*ksi[i+1]) ;
    }

    y[0] = eta[0];

    for(int i=0; i<N-2; i++)
    {
        y[i+1] = ksi[i+1]*y[i]+eta[i+1];
    }
       for(int i=1; i<=N-1; ++i)
    {
        ans[i] = y[i-1];
    }

}

double nevyazka(int f_num, double h, double* next, int M)
{
    double *rez;
    double x = h;
    int i;
    double summ = 0;

    rez = (double*)malloc(sizeof(double)*(M-1));
    for(i=0; i<M-1; ++i)
    {
        rez[i]=0;
    }



    for(i=1; i<=M-1; ++i)
    {
        rez[i-1] = next[i] - func_u(f_num, x, 1);
        x+=h;
    }


    for(i=0; i<M-1; ++i)
    {
        summ = summ+rez[i]*rez[i];
    }
    summ = summ*h;


    return sqrt(summ);
}

double func_u(int num, double x, double t)
{
    double res;
    double pi = M_PI;
    double e = M_E;

    if(num == 1)
    {
        res = sqrt(30)*t*x*(1-x);
    }

    if(num == 2)
    {
        res = sqrt(105)*pow(e, (-1)*t+1)*x*x*(1-x);
    }

    if(num == 3)
    {
        res = sqrt(630)*t*x*x*(1-x)*(1-x);
    }

    if(num == 4)
    {
        res = sqrt(2)*pow(e, (-1)*pi*pi*(t-1))*sin(pi*x);
    }
    return res;
}

double func_f(int num, double x, double t)
{
    double res;
    double pi = M_PI;
    double e = M_E;

    if(num == 1)
    {
        res = sqrt(30)*(x-x*x+2*t);
    }

    if(num == 2)
    {
        res = sqrt(105)*pow(e, 1-t)*(x*x*(x-1)-2+6*x);
    }

    if(num == 3)
    {
        res = (-1)*sqrt(630)*((-1)*x*x + 2*x*x*x - x*x*x*x + 2*t - 12*t*x + 12*t*x*x);
    }

    if(num == 4)
    {
        res = 0;
    }
return res;
}
