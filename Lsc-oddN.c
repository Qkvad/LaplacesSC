/*=====================================================================================================================
 *
 * Author : Petra Brčić
 * College: PMF-MO (DS - Primijenjena Matematika)
 * Project: Laplace equation on semicircle - methods used: cg (to add: Jacobi, Gauss-Seidel, Gauss eliminations, ..)
 *
 *
 * BUILD: gcc Lsc-oddN.c -o Lsc-oddN_exe -lblas -llapack -lm
 * RUN  : ./Lsc-oddN_exe
 *
 ====================================================================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "../lib/f2c.h"
#include "../lib/fblaswr.h"
#include "../lib/clapack.h"
#include <math.h>
#include "cg.h"
#include "ispis.h"
#include "red.h"
#define PI (3.141592653589793)
#define RED   "\x1B[31m"
#define GRN   "\x1B[32m"
#define YEL   "\x1B[33m"
#define BLU   "\x1B[34m"
#define CYN   "\x1B[36m"
#define RESET "\x1B[0m"


void LUfSolve(doublereal* A, doublereal* b, doublereal* x, integer dim,integer m, integer n){	
	integer info;
	//printM(A,dim,dim,m,n); 
	integer* ipiv=malloc(dim*sizeof(integer));
	doublereal* work=malloc(dim*sizeof(doublereal));
	dgetrf_(&dim,&dim,A,&dim,ipiv,&info);		
	if(info==0) dgetri_(&dim,A,&dim,ipiv,work,&dim,&info);
	if(info==0){
		char trans='N';
		integer inc=1;
		doublereal alpha=1.0,beta=0.0;
		dgemv_(&trans,&dim,&dim,&alpha,A,&dim,b,&inc,&beta,x,&inc);
	}
}

void postavljanjeMatrice(integer m, integer n, doublereal* A, doublereal delta_r, doublereal delta_fi){
	integer M=m*n,i;
	doublereal r=delta_r;
	doublereal r2=delta_r*delta_r, fi2=delta_fi*delta_fi,l=delta_r*2,u=r;
	for(i=0;i<M;i++) {
		if(i%n==0 && i!=0) r+=delta_r;
		A[i+M*i]=2*(1+(r2/(r*r*fi2))); //dijagonala
		if(i%n!=0 && i!=0) A[i+M*(i-1)]=-(r2/(r*r*fi2)); //ispod dijagonale
		if(i!=0 && i%n!=0) A[i-1+M*i]=-(r2/(r*r*fi2));  //iznad dijagonale
		if(i%n==0 && i!=0) l+=delta_r;
		if(i%n==0 && i>n) u+=delta_r;
		if(i<(M-n)) A[i+M*i+n]=-(1-delta_r/(2*l));  //dijagonalni blokovi ispod dijagonale
		if(i>n-1) A[i+M*i-n]=-(1+delta_r/(2*u));  //dijagonalni blokovi iznad dijagonale
	}
}

void postavljanjeMatrice1(integer m, doublereal* A, doublereal delta_r, doublereal delta_fi){
	integer i;
	doublereal r=delta_r;
	doublereal r2=delta_r*delta_r, fi2=delta_fi*delta_fi,l=delta_r*2,u=r;
	for(i=0;i<m;i++) {
		r+=delta_r;
		A[i+m*i]=2*(1+(r2/(r*r*fi2))); //dijagonala
		if(i!=0) A[i+m*(i-1)]=-(1-delta_r/(2*l)); //ispod dijagonale
		//if(i!=(m-1)) A[i-1+m*i]=-(1+delta_r/(2*u)); //iznad dijagonale
		A[i-1+m*i]=-(1+delta_r/(2*u)); //iznad dijagonale
		if(i!=0) l+=delta_r;
		u+=delta_r;
		r+=delta_r;
	}
}

void ispisMatrice(doublereal *A,integer n){	
	int i,j;
	printf("\n");
	for(i=0;i<n;i++){
		for(j=0;j<n;j++)
			printf("%lf ",A[i+j*n]);
		printf("\n");
	}
	printf("\n");
}

void solveGe(integer m, integer n, doublereal T, doublereal R, doublereal delta_r, doublereal delta_fi){
	integer M=m*n,i;
	doublereal* A=calloc(M*M,sizeof(doublereal));
	doublereal* b=calloc(m*n,sizeof(doublereal));
	doublereal* x=calloc(m*n,sizeof(doublereal));
	for(i=(m-1)*n;i<M;i++) b[i]=T;
	postavljanjeMatrice(m,n,A,delta_r,delta_fi);
	//printM(A,M,M,m,n);
	printf("\n");
	LUfSolve(A,b,x,M,m,n);
	doublereal step=T/5;
	printSolution(x,m,n,step,T); //ispis rjesenja
	//printDraw(x,m,n,step,T);
	printf("\n");	
	free(A);free(b);free(x);
}

void solveSy(integer m,integer oldn, doublereal T, doublereal R, doublereal delta_r, doublereal delta_fi){
	integer i;	
	integer n=oldn/2;
	//prvi sustav za n=1
	doublereal delta_fi1=PI/2;
	doublereal* matrica1=calloc(m*m,sizeof(doublereal));
	//postavljanjeMatrice12(m,matrica1);
	postavljanjeMatrice1(m,matrica1,delta_r,delta_fi1);
	//ispisMatrice(matrica1,m);
	doublereal* b1=calloc(m,sizeof(doublereal));
	for(i=(m-1);i<m;i++) b1[i]=T;
	doublereal* x1=calloc(m,sizeof(doublereal));
	LUfSolve(matrica1,b1,x1,m,m,n); //u x1 je spremljeno rjesenje sustava koje je rubni uvjet za drugi sustav

	//printf("\n"); for(i=0;i<m;i++) printf("%lf\n",x1[i]);

	//drugi sustav za n=n/2
	doublereal delta_fi2=(PI/2)/(n+1);
	doublereal* matrica2=calloc(m*m*n*n,sizeof(doublereal));			
	postavljanjeMatrice(m,n,matrica2,delta_r,delta_fi2);
	//ispisMatrice(matrica2,m*n); 
	//printM(matrica2,m*n,m*n,m,n);
	doublereal* b2=calloc(m*n,sizeof(doublereal));
	integer k=0;
	doublereal dr2=delta_r*delta_r,df2=delta_fi2*delta_fi2,r=delta_r;	
	for(i=0;i<m*n;i++){
		if((i+1)%n==0){
			b2[i]=(dr2/(r*r*df2))*x1[k];
			r+=delta_r;
			k++;
		}
	}
	
	for(i=(m-1)*n;i<m*n;i++) b2[i]=b2[i]+T;
	//printf("\n\n"); for(i=0;i<m*n;i++) printf("   %.4lf\n",b2[i]); printf("\n\n");
	doublereal* x2=calloc(m*n,sizeof(doublereal));
	LUfSolve(matrica2,b2,x2,m*n,m,n);
	doublereal step=T/5;

	printf("\n");
	printSolution(x2,m,n,step,T);
	printf("\n");
	free(matrica2); free(x1); free(b1); free(b2);  free(x2);
}

void solveGeSpecialb(integer m, integer n, doublereal T, doublereal R, doublereal delta_r, doublereal delta_fi){
	integer M=m*n,i;
	doublereal* A=calloc(M*M,sizeof(doublereal));
	doublereal* b=calloc(m*n,sizeof(doublereal));
	doublereal* x=calloc(m*n,sizeof(doublereal));
	doublereal fi=delta_fi;
	for(i=(m-1)*n;i<M;i++) {
		b[i]=T*sin(fi);
		fi+=delta_fi;
	}
	postavljanjeMatrice(m,n,A,delta_r,delta_fi);
	//printM(A,M,M,m,n);
	printf("\n");
	LUfSolve(A,b,x,M,m,n);
	doublereal step=T/5;
	printf("\nRješenje bez diskontinuiteta na rubu:\n");
	printSolution(x,m,n,step,T); //ispis rjesenja
	//printDraw(x,m,n,step,T);
	printf("\n");	
	free(A);free(b);free(x);
}

void realSolution(integer m, integer n, doublereal R, doublereal T, doublereal r, doublereal fi){
	doublereal* RSM=malloc(m*n*sizeof(doublereal));
	doublereal stepfi=fi,stepr=r,sum;
	int i,j,no,k=0;
	for(j=0;j<m;j++){
		for(i=0;i<n;i++){
			no=brojClanovaReda(stepr,stepfi,R,T);
			sum=sumaReda(stepr,stepfi,R,T,no);
			RSM[k]=sum;
			stepfi+=fi;
			k++;
		}
		stepr+=r;
		stepfi=fi;
	}
	printf("\nAnalitičko rješenje:\n");
	doublereal step=T/5;
	printSolution(RSM,m,n,step,T);
	free(RSM);
}

int main(){
	integer m,n,i;
	printf("\nUnesite m: "); scanf("%ld",&m);
	printf("Unesite n: "); scanf("%ld",&n);
	doublereal R,T;
	printf("Unesite temperaturu na zakrivljenom dijelu polukruga: "); scanf("%lf",&T);
	printf("Unesite radijus polukruga: "); scanf("%lf",&R);
	doublereal delta_r,delta_fi;
	delta_r=R/(m+1); delta_fi=PI/(n+1);	
	if(n%2!=0){
		char answerSy;
		printf("S obzirom da ste unijeli neparan n, želite li problem rješavati s obzirom na simetriju [Y/n]: ");
		scanf(" %c",&answerSy);
		if(answerSy!='n') solveSy(m,n,T,R,delta_r,delta_fi);
		else {
			solveGe(m,n,T,R,delta_r,delta_fi);
			realSolution(m,n,R,T,delta_r,delta_fi);
		}
	}
	else {
		solveGe(m,n,T,R,delta_r,delta_fi);
		realSolution(m,n,R,T,delta_r,delta_fi);
		solveGeSpecialb(m,n,T,R,delta_r,delta_fi);
	}

	return 0;
}
