#include<stdio.h>
#include "../lib/f2c.h"
void printM(doublereal* Matrix, int M, int N, int m, int n) {
	int i,j, blockRow, blockColumn, mRow, mColumn;
	
	printf("\n\x1b[31mYour matrix is:\n\x1b[0m");
	
	for(i=0; i<M; i++) {
		for(j=0; j<N; j++) {
			if(j==0) printf("| ");

			blockRow = ((j*N+i)/n)%m;
			blockColumn = ((j*N+i)/M)/n;
			mRow = (j*N+i)%n;
			mColumn = ((j*N+i)/M)%n;

			if(Matrix[j*N+i] >= 0) printf(" ");
			if(blockRow == blockColumn) {
				if(Matrix[j*N+i] == 0)
					printf("%.4f ", Matrix[j*N+i]);
				else if(mRow == mColumn)
					printf("\x1b[34m%.4f\x1b[0m ", Matrix[j*N+i]);
				else
					printf("\x1b[34m%.4f\x1b[0m ", Matrix[j*N+i]);
			} else if(blockRow == blockColumn-1) {
				if(Matrix[j*N+i] == 0)
					printf("%.4f ", Matrix[j*N+i]);
				else
					printf("\x1b[31m%.4f\x1b[0m ", Matrix[j*N+i]);
			} else if(blockRow == blockColumn+1) {
				if(Matrix[j*N+i] == 0)
					printf("%.4f ", Matrix[j*N+i]);
				else
					printf("\x1b[31m%.4f\x1b[0m ", Matrix[j*N+i]);
			} else {
				printf("%.4f ", Matrix[j*N+i]);
			}

			if(j==N-1) printf("|");
		}
		printf("\n");
	}
	
}

void printSolution(doublereal* x, integer m, integer n, doublereal step, doublereal T){
	int i;
	for(i=0;i<m*n;i++) {
		integer index = m*n-1-i;
		if(i%n==0) printf("\n");
		if(x[index] < 10) printf(" ");		
		if(x[index]>=(T-step)) printf("\x1B[31m%.4lf \x1B[0m",x[index]);
		if(x[index]<(T-step) && x[index]>=(T-2*step)) printf("\x1B[33m%.4lf \x1B[0m",x[index]);
		if(x[index]<(T-2*step) && x[index]>=(T-3*step)) printf("\x1B[32m%.4lf \x1B[0m",x[index]);
		if(x[index]<(T-3*step) && x[index]>=(T-4*step)) printf("\x1B[34m%.4lf \x1B[0m",x[index]);
		if(x[index]<(T-4*step)) printf("\x1B[36m%.4lf \x1B[0m",x[index]);
	}
	printf("\n");
}

void printDraw(doublereal* x, integer m, integer n, doublereal step, integer T){
	int i;
	for(i=0;i<m*n;i++) {
		integer index = m*n-1-i;
		if(i%n==0) printf("\n");		
		if(x[index]>=(T-step)) printf("\x1B[31m* \x1B[0m");
		if(x[index]<(T-step) && x[index]>=(T-2*step)) printf("\x1B[33m* \x1B[0m");
		if(x[index]<(T-2*step) && x[index]>=(T-3*step)) printf("\x1B[32m* \x1B[0m");
		if(x[index]<(T-3*step) && x[index]>=(T-4*step)) printf("\x1B[34m* \x1B[0m");
		if(x[index]<(T-4*step)) printf("\x1B[36m* \x1B[0m");
	}
	printf("\n");
}
