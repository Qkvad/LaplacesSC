#include<stdio.h>
#include<math.h>
#define PI (3.141592653589793)

int brojClanovaReda(doublereal r, doublereal fi, doublereal R, doublereal T){
	int n=1;
	doublereal tol=1e-4, elementReda, base=r/R;
	elementReda=(1/((2*(doublereal)n)-1))*pow(base,(2*n-1))*sin((2*n-1)*fi); 
	while(elementReda>tol){
		n++;
		elementReda=(1/((2*(doublereal)n)-1))*pow(base,(2*n-1))*sin((2*n-1)*fi); 
	}
	return n;
}

doublereal sumaReda(doublereal r, doublereal fi, doublereal R, doublereal T, int n){
	doublereal sum=0, base=r/R,elementReda;
	int i;
	for(i=5;i>0;i--){
		sum+=(1/((2*(doublereal)i)-1))*pow(base,(2*i-1))*sin((2*i-1)*fi); 
	}
	return ((4*T)/PI)*sum;
}
