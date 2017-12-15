/*
 * 24-703 Group Project
 * Group member: Di Wu, Jiayi Wang, Kai Ge
 * Andrew id   : dwu1, jiayiw2, kge
 *
 * Description:
 * This file decomposes a matrix into LU format and solves 
 * the matrix equation Ku=P, thus obtains the value of u, 
 * which is also known as displacement. 
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include "Integration1.h"
#include <functional>
using namespace std;


// declaration of functions/subroutines
double K(int i, int j, int n, double h);

double P(int i, int n, double h);
// declaration of output stream

ofstream data_out("output.txt"); // output data


// the main part of the program - should always be type "int"
int main() {
	
	// variable declarations
    int i,j,k; // loop counters 
    
    const int n = 16; // number of elements
	const int N = 2 * n + 1; //number of nodes
    
    double h = 1. / (2*n); // discretization length
    
    double A[N][N]; // coefficient matrix
    double L[N][N]; // lower triangular matrix
    double U[N][N]; // upper triangular matrix
    
    double x[N]; // solution vector
    double y[N]; // intermediate solution vector
    double b[N]; // right side vector
    
    double sum; // dummy variable to keep track of sums in decomposition

	// coefficient matrix for conduction problem (you could enter or read in any square matrix of size n)
	
	for(i=0;i<N;i++) {
		for(j=0;j<N;j++) A[i][j] = 0.;
	}
	
	A[0][0] = -3./2.;	//initialize the first line of A
	A[0][1] = 2.;
	A[0][2] = -1./2.;

	A[N - 1][N - 1] = 1.;	//initialize the last line of A
	for(i=1;i<N-1;i++) { //internal nodes
		if (i % 2 ==0) {	//the even lines
			A[i][i - 2] = K(3, 1, i / 2, h);
			A[i][i - 1] = K(3, 2, i / 2, h);
			A[i][i] = K(3, 3, i / 2, h) + K(1, 1, i / 2 + 1, h);
			A[i][i + 1] = K(1, 2, i / 2 + 1, h);
			A[i][i + 2] = K(1, 3, i / 2 + 1, h);
		}
		else {	//the odd lines
			A[i][i - 1] = K(2, 1, i / 2 + 1, h);
			A[i][i] = K(2, 2, i / 2 + 1, h);
			A[i][i + 1] = K(2, 3, i / 2 + 1, h);
		}
	}
 
   // right side vector
   
   // known temperature boundary conditions
   
   b[0] = 0.;
   b[N-1] = 0.;
   
   // heat generation at internal nodes
   
   for (i = 1; i < N - 1; i++) {
	   if (i % 2 == 0) {
		   b[i] = P(3, i / 2, h) + P(1, i / 2 + 1, h);
	   }
	   else {
		   b[i] = P(2, i / 2 + 1, h);
	   }
   }
  
   // initialize L and U
   
   for(i=0;i<N;i++) {
	   for(j=0;j<N;j++){
		   L[i][j] = 0.;
		   U[i][j] = 0.;
	   }
   }
    
   for(i=0;i<N;i++) U[i][i] = 1.; // for Crout reduction
  
   // Do the LU decomposition using the Crout reduction
   
   for(i=0;i<N;i++) { // loop over pairs of L columns and U rows. The three levels of loop indicate n^3 behavior
	  
     // first, column i of L
     
   		for(j=i;j<N;j++) { // i is the row
			sum=0.;
			for(k=0;k<=i-1;k++) sum = sum + L[j][k]*U[k][i];
			L[j][i] = A[j][i] - sum;
		}
		
	// second, row i of U
		   
		  for(j=i+1;j<N;j++){ // j is the column
			   sum = 0.;
			   for(k=0;k<=i-1;k++) sum = sum + L[i][k]*U[k][j];
			   U[i][j] = (A[i][j] - sum)/L[i][i];
		   }
   }
   
   // output intermediate data to screen
   
   for(i=0;i<N;i++) {
	   for(j=0;j<N;j++) cout<<A[i][j]<<'\t';
	 cout<<endl;
   }
   cout<<endl;		   
    
   for(i=0;i<N;i++) {
	   for(j=0;j<N;j++) cout<<L[i][j]<<'\t';
	   cout<<endl;
   }
   cout<<endl;	
   
   for(i=0;i<N;i++) {
	   for(j=0;j<N;j++) cout<<U[i][j]<<'\t';
	 cout<<endl;
   }
  cout<<endl;	
   
   // solve the system of equations
   
   // could loop over a series of b vectors here once the decomposition has been done
   
   // first, find the y vector
   
   y[0] = b[0]/L[0][0];
   for (i = 1; i < N; i++) {
	   sum = 0.;
	   for(j=0;j<i;j++) sum = sum + L[i][j]*y[j];
	   y[i] = (b[i]-sum)/L[i][i];
   }
   
   for(i=0;i<N;i++) cout<<y[i]<<'\t';
   cout<<endl<<endl;;
   
   // second, find the x vector
   
   x[N-1] = y[N-1];
   for(i=1;i<N;i++) {
	   j = N-i-1;
	   sum = 0;
	   for(k=j+1;k<N;k++) sum = sum + U[j][k]*x[k];
	   x[j] = y[j] - sum;
   }
   
   // output data
   
   for(i=0;i<N;i++) cout<<x[i]<<'\t';
   for(i=0;i<N;i++) data_out<<"u"<<i<<'\t'<<x[i]<<endl;
   cout<<endl;

}	
