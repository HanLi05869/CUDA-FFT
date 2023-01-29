/**
 * @author HAN LI
 * @email  xxmlih05869@outlook.com
 * @time   2022.12.15
 */
#include<iostream>
#include<cmath>
#include<string>
#include<cstring>
#include <sys/time.h>
#include "support.h"
using namespace std;
const int maxn=1e6+5;
 
struct Complex //complex
{
    double x,y;
    Complex(double dx=0,double dy=0)
    {
        x=dx;
        y=dy;
    }
};

Complex operator +(Complex a,Complex b)
{
    return Complex(a.x+b.x,a.y+b.y);
}
Complex operator -(Complex a,Complex b)
{
    return Complex(a.x-b.x,a.y-b.y);
}
Complex operator *(Complex a,Complex b)
{
    return Complex(a.x*b.x-a.y*b.y,a.x*b.y+a.y*b.x);
}

const double pi=acos(-1.0); //PI
int limit=1,bit=0; //limit is the final length extended. limit = 1<<bit
int wz[maxn<<2];
int re[maxn<<2]; //save results
Complex a[maxn<<2],b[maxn<<2];
char s1[maxn],s2[maxn];//save two input numbers

void FFT(Complex *A,int inv)
{
    //base on the reverse bit sequience of the number,just swap the corresponding elements
    for(int mid=1;mid<limit;mid<<=1)
    {
        Complex wn(cos(pi/mid),inv*sin(pi/mid));
        for(int i=0;i<limit;i+=(mid<<1))
        {
            Complex w(1,0);
            for(int j=0;j<mid;j++,w=w*wn)
            {
                Complex t1=A[i+j];
                Complex t2=w*A[i+mid+j];
                A[i+j]=t1+t2;
                A[i+mid+j]=t1-t2;
            }
        }
    }
}

Timer timer;
int main(int argc,char* argv[])
{
    int returnflg = scanf("%s%s",s1,s2);
    int len1=strlen(s1),len2=strlen(s2);
    int len=len1+len2-2;//add times, denotes the highest order of the result polynominal
    len+=1;//n order polynominal needs n+1 points to represent it
    while(limit<len)//ensure that n is no less than the exponential power of 2
    {
        limit<<=1;
        bit++;
    }               
    //An n-digit decimal number can be viewed as an n-1 polynomial
    for(int i=len1-1,j=0;i>=0;i--,j++)
    {
        a[j].x=s1[i]-'0';
        a[j].y=0;
    }
    for(int i=len2-1,j=0;i>=0;i--,j++)
    {
        b[j].x=s2[i]-'0';
        b[j].y=0;
    }
    startTime(&timer);
    for(int i=0;i<limit;i++)
        wz[i]=(wz[i>>1]>>1)|((i&1)<<(bit-1));//uncanny
    stopTime(&timer);
    printf("%f\n",elapsedTime(timer)); fflush(stdout);
    
    for(int i=0;i<limit;i++)
        if(i<wz[i])
            swap(a[i],a[wz[i]]);
    startTime(&timer);
    FFT(a,1);
    stopTime(&timer);
    printf("%f\n",elapsedTime(timer)); fflush(stdout);

    for(int i=0;i<limit;i++)
        if(i<wz[i])
            swap(b[i],b[wz[i]]);
    startTime(&timer);
    FFT(b,1);
    stopTime(&timer);
    printf("%f\n",elapsedTime(timer)); fflush(stdout);
    
    for(int i=0;i<limit;i++)
        a[i]=a[i]*b[i];
    
    for(int i=0;i<limit;i++)
        if(i<wz[i])
            swap(a[i],a[wz[i]]);
    startTime(&timer);
    FFT(a,-1);
    stopTime(&timer);
    printf("%f\n",elapsedTime(timer)); fflush(stdout);
    memset(re,0,sizeof(re));
    for(int i=0;i<=limit;i++)
    {
        re[i]+=(int)(a[i].x/limit+0.5);
        if(re[i]>=10) //carry bit
        {
            re[i+1]+=re[i]/10;
            re[i]%=10;
            if(i==limit)
                ++limit;
        }
    }
    while(limit&&!re[limit])//exclude high order zero
        limit--;
    FILE* fp = fopen("./res","w");
    while(limit>=0)
        fprintf(fp,"%d",re[limit--]);
    fprintf(fp,"\n");
    fclose(fp);
    return 0;
}
