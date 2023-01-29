#include<iostream>
#include<cmath>
#include<string>
#include<cstring>
#include <sys/time.h>
#include "support.h"
const int MAX = 1e6+10;
char input_1[MAX];
char input_2[MAX];
char result[MAX<<1];

int tmp1[MAX];
int tmp2[MAX];
int tmp3[MAX<<1]; 

Timer timer;   
int main()
{
	int returnflg = scanf("%s%s",input_1,input_2);
	
	for(int i = 0;input_1[i]!='\0';i++)
		tmp1[strlen(input_1)-i-1] = input_1[i] - '0';
	for(int i = 0;input_2[i]!='\0';i++)
		tmp2[strlen(input_2)-i-1] = input_2[i] - '0';
	int size1 = strlen(input_1);
	int size2 = strlen(input_2);
	
	startTime(&timer);
	for(int i = 0;i < size1;i++)
		for(int j = 0;j < size2;j++)
			tmp3[i+j] += tmp1[i] * tmp2[j];
	
	for(int i = 0;i <= size1 + size2 -2;i++){
		tmp3[i + 1] += tmp3[i] / 10;
		tmp3[i]      =  tmp3[i] % 10; 
	}
	for(int i = 0;i <= size1 + size2 - 1;i++)
	{
		result[size1 + size2 - 1 - i] = tmp3[i] + '0';
	}
	int i, t;
	for(i = 0;;i++)
		if(result[i] != '0')
		{
			t = i;
			break;
		}
	while(i < size1 + size2)
	{
		result[i - t] = result[i];
		i += 1;
	}
	result[i - t] = '\0';
	stopTime(&timer); 
	printf("%f\n", elapsedTime(timer));
	
	return 0; 
	
}
