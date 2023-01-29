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
#define BLOCK_SIZE 1024
using namespace std;
const int maxn=1e6+5; 

struct Complex //complex
{
    double x,y;
    Complex()
    {
		x=0.0;
		y=0.0;
    }
    
    void set(double dx,double dy)
    {
    	x=dx;
		y=dy;
	}
	
    __device__
	Complex(double dx,double dy)
    {
        x=dx;
        y=dy;
    }
};

Complex operator +(Complex a,Complex b)
{
	Complex tmp;
	tmp.set(a.x+b.x,a.y+b.y);
    return tmp;
}
Complex operator -(Complex a,Complex b)
{
	Complex tmp;
	tmp.set(a.x-b.x,a.y-b.y);
    return tmp;
}
Complex operator *(Complex a,Complex b)
{
	Complex tmp;
	tmp.set(a.x*b.x-a.y*b.y,a.x*b.y+a.y*b.x);
    return tmp;
}

__device__
void gpuset(Complex& it,double dx=0,double dy=0)
{
	it.x=dx;
	it.y=dy;
}
__device__
Complex gpuadd(const Complex& a,const Complex& b)
{
	return Complex(a.x+b.x,a.y+b.y);
}

__device__
Complex gpusub(const Complex& a,const Complex& b)
{
	return Complex(a.x-b.x,a.y-b.y);
}

__device__
Complex gpumul(const Complex& a,const Complex& b)
{
	return Complex(a.x*b.x-a.y*b.y,a.x*b.y+a.y*b.x);
}

const double pi=acos(-1.0); //PI
int limit=1,bit=0; //limit is the final length extended. limit = 1<<bit
int wz[maxn<<2];
int re[maxn<<2]; //save results
Complex a[maxn<<2],b[maxn<<2];
char s1[maxn],s2[maxn];//save two input numbers

Timer timer;
Complex* da,*db,*dc;
int* dwz;

__global__
void FFT(Complex* __restrict__ a,int limit,int flag)
{
	for(int stride = 2; stride <= limit; stride <<= 1) 
	{
		const double _PI = acos(-1.0);
		int validx = (threadIdx.x + blockIdx.x * blockDim.x) * stride;
		for (int j = 0; j < (stride >> 1); j++) 
		{
			if (validx + j + (stride >> 1) < limit) 
			{
				Complex wn(cos(2.0*_PI*flag*j/stride), sin(2.0*_PI*flag*j/stride));
				Complex u = a[j + validx];
				Complex v = gpumul(wn, a[j + validx + stride / 2]);
				a[j + validx] = gpuadd(u, v);
				a[j + validx + (stride >> 1)] = gpusub(u, v);
			}
		}
		__syncthreads();
	}
}


__global__
void bitreverse(int* __restrict__ dwz, int bits, int limit){
	int i = (blockIdx.x * blockDim.x + threadIdx.x);
	int idx = i;
	int r = 0;
    do {
        r += i % 2 << --bits;
    } while (i /= 2);
    *(dwz+idx) = r;
}

int main(int argc,char* argv[])
{
    int returnflag = scanf("%s%s",s1,s2);
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
    //for(int i=0;i<limit;i++)
    //    wz[i]=(wz[i>>1]>>1)|((i&1)<<(bit-1));
    const unsigned int numBlocks = ceil(1.0 * limit / BLOCK_SIZE);
    
    cudaMalloc((void**)&dwz,sizeof(int) * limit);
    cudaDeviceSynchronize();
    cudaMemcpy(dwz, wz, sizeof(int) * (limit), cudaMemcpyHostToDevice);
    cudaDeviceSynchronize();
    startTime(&timer);
	bitreverse <<< numBlocks, BLOCK_SIZE >>> (dwz, bit, limit);
	cudaDeviceSynchronize();
	stopTime(&timer); 
	printf("%f\n", elapsedTime(timer));
    cudaDeviceSynchronize();
    cudaMemcpy(wz, dwz, sizeof(int) * (limit), cudaMemcpyDeviceToHost);
    cudaFree(dwz);
    
    
    for(int i=0;i<limit;i++)
        if(i<wz[i])
            swap(a[i],a[wz[i]]);
    for(int i=0;i<limit;i++)
        if(i<wz[i])
            swap(b[i],b[wz[i]]);     
	
    cudaMalloc((void**)&da, sizeof(Complex) * (limit));
    cudaMalloc((void**)&db, sizeof(Complex) * (limit));
    cudaDeviceSynchronize();
    
    cudaMemcpy(da, a, sizeof(Complex) * (limit), cudaMemcpyHostToDevice);
    cudaMemcpy(db, b, sizeof(Complex) * (limit), cudaMemcpyHostToDevice);
    cudaDeviceSynchronize();
    
    
   	
    startTime(&timer);
    FFT <<< numBlocks, BLOCK_SIZE >>> (da,limit,1);
	cudaDeviceSynchronize();
	stopTime(&timer); 
	printf("%f\n", elapsedTime(timer));
	cudaMemcpy(a, da, sizeof(Complex) * (limit), cudaMemcpyDeviceToHost);

    startTime(&timer);
	FFT <<< numBlocks, BLOCK_SIZE >>> (db,limit,1);
	cudaDeviceSynchronize();
	stopTime(&timer); 
	printf("%f\n", elapsedTime(timer));
	cudaMemcpy(b, db, sizeof(Complex) * (limit), cudaMemcpyDeviceToHost);
    
    for(int i=0;i<limit;i++)
        a[i]=a[i]*b[i];
    
    for(int i=0;i<limit;i++)
        if(i<wz[i])
            swap(a[i],a[wz[i]]);
    
    cudaMemcpy(da, a, sizeof(Complex) * (limit), cudaMemcpyHostToDevice);
    startTime(&timer);
	FFT <<< numBlocks, BLOCK_SIZE >>> (da,limit,-1);
	cudaDeviceSynchronize();
	stopTime(&timer); 
	printf("%f\n", elapsedTime(timer));
	cudaMemcpy(a, da, sizeof(Complex) * (limit), cudaMemcpyDeviceToHost);

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
    cudaFree(da);
    cudaFree(db);
    return 0;

}
