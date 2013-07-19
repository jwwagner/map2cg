#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <signal.h>
#include <limits.h>
#include <time.h>

void foo(int*);
void bar(int*);

int main(int argc,char *argv[])
{
	int* test = malloc(10*sizeof(int));
	int i;
	
	for(i = 0; i < 10; i++)
		{
		test[i] = i;
		printf("main: test[%d] = %d\n", i, test[i]);
		}

	printf("pass to foo\n\n");
	foo(test);
	
	for(i = 0; i < 10; i++)
		{
		printf("main: test[%d] = %d now\n", i, test[i]);
		}

	free(test);

	return 0;
}

void foo(int* test)
{
	int i;
	
	for(i = 0; i < 3; i++)
		{
		printf("test[%d] = %d\n", i, test[i]);
		bar(&test[i]);
		printf("verify: test[%d] = %d\n", i, test[i]);
		}

}

void bar(int* test)
{
	*test += 1;
	printf("bar: value of test is now %d\n", *test);
}