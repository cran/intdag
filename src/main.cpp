#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include <R_ext/Rdynload.h>
#include <R_ext/Applic.h>
#include <R_ext/BLAS.h>
#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include "def.h"
#include "dmatrix.h"

using namespace std;

const double eps_abs = 0.001;
const double eps_rel = 0.01;


extern "C"{

	/* check loop in C
	 * check whether adding a directed edge from node a to node b induces cycles in graph G.
	 * ij-th element of G indicate j --> i
	 * 
	 * Arguments
	 *
	 * node: number of nodes
	 * *G: graph, a node*node matrix, with the same coding as A_NZ
	 */
	int Cycle(int node, int *G, int a, int b)
	{
		if(a == b)  return 1;  /* no edge, no cycle  */

		int i, j, lo;
		int nBot   = 0;
		int nTop   = 0;
		int SLeng  = 1;
		int bCycle = 0;

		int *color = new int[node] ();
		int *S     = new int[node] ();

		color[a] = 1;
		S[0]     = a;

		while(SLeng > 0)
		{
			i = S[nBot];
			SLeng--;
			nBot++;
			for(j=0;j<node;j++)
			{
				lo= j*node;			 
				if(G[lo+i]==1)
				{
					if(j==b)
					{
						bCycle=1;
						break;
					}
					else if(color[j]==0)
					{
						nTop++;
						S[nTop]=j;
						SLeng++;
						color[j]=1;
					}
				}
			}
			if(bCycle==1) break;
		}

		delete [] color;
		delete [] S;

		return bCycle;
	}


	// compare function, to use in quicksort when update sigma
	int cmpfunc(const void * a, const void * b){
		return (*(int*)a - *(int*)b);
	}


	// ADMM objective function
	double objective(int m, int M, int n, double mu, double mu2, double *Sig, double *X, double *A, 
		double *B, int *A_NZ)
	{
		int i, j, pos;
		double obj = 0.0;
		double *XAT = new double[n*M] ();
		// 1. likelihood part
		// XAT = X*A^T
		dmat_C_ABT(M, n, m, X, A, XAT);
	
		for(j=0; j<m; j++){
			for(i=0; i<n; i++){
				// obj += (X[i+j*M]-XAT[i+j*M]) * (X[i+j*M]-XAT[i+j*M])/2/Sig[j] + (n*log(Sig[j])/2);
				obj += (X[i+j*M]*Sig[j]-XAT[i+j*M]) * (X[i+j*M]-XAT[i+j*M])/2;
			}
			obj -= n*log(Sig[j]);
		}
		// 2. penalty part
		for(pos=0; pos<(m*m); pos++){
			if(A_NZ[pos] != 1){
				obj += mu*ABS(B[pos]);
			}
		}
		//added for the covariate part
		for(pos=(m*m); pos<(m*M); pos++){
			if(A_NZ[pos] != 1){
				obj += mu2*ABS(B[pos]);
			}
		}
		delete [] XAT;
		return(obj);
	}
	


	// update parameter Sig: variances for nodes
	void update_Sig(int m, int M, int n, double *sigma, double *Sig, double *X, double *X2, double *A){
		int i, j;
		double tmp;
		double *XAT = new double[n*M] ();
		// XAT = X*A^T
		dmat_C_ABT(M, n, m, X, A, XAT);
		// update the variances one by one
		for(j=0; j<m; j++){
			tmp = 0.0;
			// tmp = sum(X-XA^T)^2
			for(i=0; i<n; i++){
				tmp += X[i+j*n] * XAT[i + j*n];
			}
			Sig[j] = (sqrt(tmp*tmp + 4*n*X2[j]) + tmp) / 2 / X2[j];
		}
		delete [] XAT;
	}
	
	
	
	// update parameter C
	// note that all sig_X, sigma and Sig are in squared form, i.e., variances
	void update_C(int m, int M, double sig_X, double *B, double *UC, double *sigma, double *Sig, double *C){
		int i, pos;
		int S = MIN(m, M-m);
		double sig   = *sigma;
		double *Lam  = new double[m*m] ();  // intermediate parameter for variance constraint
		double *Tem1 = new double[m*(M-m)] ();
		double *Tem2 = new double[m*(M-m)] ();
		double *U    = new double[m*m] ();  // SVD matrices
		double *D    = new double[m] ();
		double *VT   = new double[m*m] ();
	
		// 1. compute Lam matrix
		for(i=0; i<m; i++){
			// Lam[i+i*m] = pow((sig-Sig[i])/sig_X, 0.5);
			Lam[i+i*m] = pow((sig*Sig[i]-1)/sig_X, 0.5);
		}
	
		// 2. compute matrix to be decomposed
		for(pos=0; pos<m*(M-m); pos++){
			Tem1[pos] = B[pos+m*m] - UC[pos];
		}

		// Tem2 = Tem1^T*Lam
		dmat_C_ATB(m, M-m, m, Tem1, Lam, Tem2);
	
		// 3. SVD
		// Tem2 = U*D*V^T
		svd_c(M-m, m, Tem2, U, D, VT);
	
		// 4. update C
		// Tem1 = V*U^T
		dmat_C_ATBT(S, m, M-m, VT, U, Tem1);
		// C = Lam*Tem1 = Lam*V*U^T
		dmat_C_AB(m, m, M-m, Lam, Tem1, C);
	
		delete [] Lam;
		delete [] Tem1;
		delete [] Tem2;
		delete [] U;
		delete [] D;
		delete [] VT;
	}


	// update A exactly
	void update_A_exact(int m, int M, double rho2, double *A, double *B, double *U, double *XTX, double *XTX_inv){
		int m_A = M-1;
		int i, j, k;
		double *V_A   = new double[M-1] ();
		double *A_row = new double[M-1] ();

		for(k=0; k<m; k++){

			//calculate V_A
			for(i=0; i<(M-1); i++){
				if(i >= k)  j=i+1;
				else        j=i;
				V_A[i] = XTX[j*M+k] + rho2*(B[j*m+k] - U[j*m+k]); // col-major matrix, k-th row, j-th col
			}

			// calcluate A
			// A_row = XTXinv * V_A
			dmat_yAx(m_A, m_A, &XTX_inv[k*(M-1)*(M-1)], V_A, A_row);

			for(i=0; i<(M-1); i++){
				if(i >= k)  j = i+1;
				else        j = i;
				A[j*m+k] = A_row[i];
			}
		}
		// print_dmatrix(A, m, M);

		delete [] V_A;
		delete [] A_row;
	}


	// update A in convex reparametrization
	void update_A_convex(int m, int M, double rho2, double *A, double *B, double *U, double *XTX, double *XTX_inv, double *Sig){
		int m_A = M-1;
		int i, j, k;
		double *V_A   = new double[M-1] ();
		double *A_row = new double[M-1] ();

		for(k=0; k<m; k++){ // update the k-th row of A

			//calculate V_A
			for(i=0; i<(M-1); i++){
				if(i >= k)  j=i+1;
				else        j=i;
				V_A[i] = XTX[j*M+k]*Sig[k] + rho2*(B[j*m+k] - U[j*m+k]); // col-major matrix, k-th row, j-th col
			}

			// calcluate A
			// A_row = XTXinv * V_A
			dmat_yAx(m_A, m_A, &XTX_inv[k*(M-1)*(M-1)], V_A, A_row);

			for(i=0; i<(M-1); i++){
				if(i >= k)  j = i+1;
				else        j = i;
				A[j*m+k] = A_row[i];
			}
		}
		// print_dmatrix(A, m, M);

		delete [] V_A;
		delete [] A_row;
	}




	/* Function to solve the unconvex interventional DAG estimation problem
	 *
	 * Variables
	 *
	 * X:       n by M data matrix
	 * A:       m by M adjacency matrix 
	 * m:       the number of nodes
	 * M:       the number of total nodes (nodes + intervention nodes)
	 * n:       sample size
	 * lambda:  penalty for A
	 * lambda2: penalty for B (intervention part)
	 * tau:     \tau in the paper, tuning parameter for TLP
	 * A_NZ:    m by M matrix (0 indicates the coresponding element of A being 0,
	 *          1 indicates nonzero. -1 indicates a suggested direction, which will be
	 *          checked outside this function. If this suggested direction does not violate
	 *          acyclicity, -1 will be converted to 1.
	 * NZ:      the number of nonzero elements in A. ( the number of 1's in A_NZ)
	 * sigma:   the parameter in the variance constraint.
	 * tol:     tolerance level
	 * obj:     objective function value
	 * XTX:     m by m matrix X^T %*% X
	 * XTX_inv: inverse matrix of XTX, (m-1)*(m-1)*m, different matrices are for calculating different A_row
	 * B:       m by m matrix (B in the paper)
	 * LLambda: m by m dual variable matrix (\lambda in the paper)
	 * Sig:     vector of length m, the error variances sigma^2_j of each node
	 * sigma_X: the error variance of intervention nodes
	 *
	 *
	 * notice that for Xi and y, they are three dimensional and row-major arrays.
	 * A, B, U, LLambda are column-major, as in R, to facilitate the use of fortran functions.
	 *
	 */
	void DAG_int_var(double *X, double *A, double *B, double *C, double *U, int *mm, int *MM, int *nn, 
		double *Lambda, double *Lambda2, double *tautau, int *A_NZ, int *NZ, double *sigma, 
		double *sigmaX, double *Sig, double *tol, double *obj, double *XTX, double *XTX_inv,
		double *rhorho, double *rhorho2, double *rhorho3, int *maxIter)  			
	{

		// take in the values from R
		int m       = (*mm);  // num of nodes
		int M       = (*MM);  // num of nodes + num of intervention modes
		int n       = (*nn);  // num of samples
		int nonzero = (*NZ);  // num of nonzero elements in A

		double lambda   = (*Lambda);
		double lambda2  = (*Lambda2);
		double tau      = (*tautau);    // tuning parameter in TLP function
		double opts_tol = (*tol);       // tolerance
		double mu       = lambda/tau;
		double mu2      = lambda2/tau;
		double rho      = (*rhorho);    // ADMM penalty for DAG constraint
		double rho2     = (*rhorho2);   // ADMM penalty for A=B
		double rho3     = (*rhorho3);   // ADMM penalty for B=C
		double sig_X    = (*sigmaX);

		// initialize the parameters
		int i          = 0;
		int j          = 0;
		int k          = 0;
		int pos        = 0;
		int iterDC     = 0;
		int iter       = 0;  // ADMM iteration number
		
		double valueDC     = 0.0;
		double valueDCp    = 0.0;   /* previous value in DC loop */
		double Bab_temp    = 0.0;
		double obj_cur     = 0.0;   /* current objective value  */
		double obj_pre     = 0.0;   /* previous objective value  */
		double res         = 0.0;   /* residuals for primal feasibility condition w.r.t. y */
		double res_temp    = 0.0;	/* intermediate variable to update y */
		double r1;

		double *F        = new double[m*M] ();      // B matrix in the previous iter.
		double *Bab      = new double[m*M] ();    	// absolute value of B
		double *UC       = new double[m*(M-m)] ();	// dual variable U_C
		double *LLambda  = new double[m*m] ();    	// lambda matrix
		double *Mat_step = new double[m*m] ();    	// intermediate matrix when update LLambda
		double *V_step   = new double[m*m] ();    	// intermediate matrix when update LLambda
		double *Xi       = new double[m*m*m] ();  	// kisi in the notes
		double *y        = new double[m*m*m] ();  	// dual variable y

		// initialize parameters
		for(pos=0; pos<m*M; pos++){
			B[pos]   = 0;
			Bab[pos] = 0;
			U[pos]   = 0;
		}

		for(pos=0; pos<m*(M-m); pos++){
			UC[pos] = 0;
		}

		for(pos=0; pos<m*m*m; pos++){
			Xi[pos] = 0;
			y[pos]  = 0;
		}

		for(i=0; i<(m*m); i++){
			LLambda[i] = 1; 
			V_step[i]  = 0;
		}

		// initiate Mat_step and V_step
		for(i = 0; i<m; i++){  // row
		    for(j =0; j<m; j++){  // column
				if(j == 0)       Mat_step[j*m + i] = 1;
				else if(i == 0)  Mat_step[j*m + i] = 0;
				else if(i == j)  Mat_step[j*m + i] = 2/((double)m);
				else             Mat_step[j*m + i] = 1/((double)m);	      
		    }	  
		}
		 

		// calculate objective value: norm((X - X*A'), 'fro')-Mu*sum(sum(abs(A)))
		obj_cur = objective(m, M, n, mu, mu2, Sig, X, A, B, A_NZ);
		obj_pre = obj_cur;


		// beginning of DC loop
		for(iterDC = 1; iterDC <= 6; iterDC ++){ //DC loop
		  

		  	// initialize all the dual variables
		  	for(i=0; i<(m*M); i++){
		  		U[i] = 0;
		  	}

		  	for(i=0; i<(m*(M-m)); i++){
		  		UC[i] = 0;
		  	}

		  	for(i=0; i<(m*m*m); i++){
		  		y[i] = 0;
		  	}
	 		

			// beginning of ADMM loop
			for(iter=1; iter < *maxIter; iter++){ // ADMM loop
	
				// step 1: A direction, update by rows

				update_A_convex(m, M, rho2, A, B, U, XTX, XTX_inv, Sig);


				// step 2: update C matrix
				update_C(m, M, sig_X, B, UC, sigma, Sig, C);


				// step 3: B direction
				for(pos=0; pos<m*M; pos++){
					F[pos] = B[pos];
				}

				// part I, adjacency matrix
				for(pos = 0; pos<(m*m); ++pos){
					j = pos/m;
					i = pos%m;
					Bab_temp = 0.0;

					// calculate the temporary sum: M_{ijk} in notes
					for(k = 0; k<m; k++){
						Bab_temp += LLambda[k*m + i] - LLambda[k*m+j] - Xi[i+j*m+k*m*m] - y[i+j*m+k*m*m];
					}

					Bab_temp = Bab_temp + (m-1)*tau;

					// use the sum to update B
					if(A_NZ[pos] == 1) {
						Bab[pos] = tau;
						B[pos]   = A[pos] + U[pos];
					}
					else {
						Bab[pos] = MAX(0, ((rho2 * ABS(A[pos]+U[pos]) + rho*Bab_temp - mu)/(rho2+m*rho)));
						B[pos]   = Bab[pos] * (A[pos]+U[pos]>=0?1:-1);
					}
					
				}

				// part II, intervention covariates
				for(i=0; i<m; i++){
					for(j=m; j<M; j++){
						if(i == j-m){
							Bab_temp = (rho2*(A[i+j*m]+U[i+j*m])+rho3*(C[i]+UC[i]))/(rho2+rho3);
							if(A_NZ[i+j*m]==0)	B[i+j*m] = Bab_temp;
							else 				B[i+j*m] = MAX(0, Bab_temp - mu2/(rho2+rho3)) - MAX(0,-Bab_temp - mu2/(rho2+rho3));
						}
						else{
							if(A_NZ[i+j*m]==0)	B[i+j*m] = A[i+j*m] + U[i+j*m];
							else 				B[i+j*m] = MAX(0,A[i+j*m]+U[i+j*m] - mu2/rho2) - MAX(0,-A[i+j*m]-U[i+j*m] - mu2/rho2);
						}
					}
				}


				// step 4: update Lambda matrix
				//initiate V_step;
				for(i = 0; i < m; ++i) {
				  for(k = 0; k < m; ++k){
				    V_step[k*m+i]=0;	   
				    for(j=0; j<m; ++j)   V_step[k*m+i] += Bab[j*m+i] + Xi[i+j*m+k*m*m] + y[i+j*m+k*m*m];		      
				  }
				}
				for(j = 0; j < m; ++j) {
				  for(k = 0; k < m; ++k){	   
				    for(i=0; i<m; ++i)   V_step[k*m+j] -= (Bab[j*m+i] + Xi[i+j*m+k*m*m] + y[i+j*m+k*m*m]);		      
				  }
				}
				for(i = 1; i < m; ++i) {
				  for (k = 0; k < m; ++k){	  	   
				    if(i == k)   V_step[k*m+i] = (-(m-1)*tau + V_step[k*m+i])/2;	
					else         V_step[k*m+i] = (tau + V_step[k*m+i])/2;		      
				  }
				}
				//V_step(1,:) = ones(1,p);
				for(k = 0; k < m; ++k)   V_step[k*m] = 1;
	
				/* LLambda = Mat_step * V_step  */
				dmat_C_AB(m, m, m, Mat_step, V_step, LLambda);
				// end of Lambda update



				// step 7: update Xi and y
				res = 0;
				r1  = 0;  // primal residual 1
				for(i = 0; i < m; ++i) {
					for(j = 0; j < m; ++j){
				    	if(i != j){
				    		for ( k = 0; k < m; ++k){	      
				    			if(k == j){
				    				/* res_temp is the increment of dual variable y in each step */
				    				/* res is the largest  */
									res_temp        = Bab[j*m+i] + Xi[i+j*m+k*m*m] - LLambda[k*m+i] + LLambda[k*m+j];
									Xi[i+j*m+k*m*m] = MAX(0, (LLambda[k*m+i] - LLambda[k*m+j] - Bab[j*m+i]- y[i+j*m+k*m*m]));
									y[i+j*m+k*m*m]  = y[i+j*m+k*m*m] + res_temp;
									r1             += res_temp * res_temp;
									// res             = MAX(res, (res_temp*res_temp));
				    		  	}
				    		  	else{
									res_temp        = Bab[j*m+i] + Xi[i+j*m+k*m*m] - LLambda[k*m+i] - tau + LLambda[k*m+j];
									Xi[i+j*m+k*m*m] = MAX(0, (LLambda[k*m+i] + tau - LLambda[k*m+j] - Bab[j*m+i]- y[i+j*m+k*m*m]));
									y[i+j*m+k*m*m]  = y[i+j*m+k*m*m] + res_temp;
									r1             += res_temp * res_temp;
									// res             = MAX(res, (res_temp*res_temp));
				    		  	}
				    		}
				    	}
				  	}
				}  // end of Xi and y update



				// step 8: U & UC direction
				for(pos = 0; pos<(m*M); pos++){
					// U[pos] = U[pos] + ((A[pos] - B[pos])/2);  // use a 1/2 step size, inexact update
					U[pos] = U[pos] + A[pos] - B[pos];
				}

				for(i=0; i<m*(M-m); i++){
					UC[i] = UC[i] + C[i] - B[i+m*m];
				} // end of U and UC update

	
				// end of all parameter updates


				// update objective value
				obj_pre = obj_cur;
				obj_cur = objective(m, M, n, mu, mu2, Sig, X, A, B, A_NZ);

				
				
				// stopping rule for ADMM
				if(ABS((obj_cur-obj_pre)) < opts_tol*obj_pre) break;
			
			}  // end of ADMM


			/* update objective values in DC loop */
			valueDCp = valueDC;//
			valueDC  = obj_cur;


			/* update A_NZ and count num of nonzero elements */
			nonzero = 0;
			for(i=1; i<m; i++){
			    for(j=0; j<i; j++){
					if((ABS(A[j*m+i])+ABS(A[i*m+j])) >= tau*0.99){   
						if(ABS(A[j*m+i]) > tau*0.5){
							if(Cycle(m, A_NZ, j, i) == 0){
								A_NZ[j*m+i]=1;
								A_NZ[i*m+j]=0;
								nonzero++;
							}
							else{      // if (Cycle(m,A_NZ, i,j) == 0) ...
								A_NZ[j*m+i]=0;
								A_NZ[i*m+j]=0;
							}
					    }
					    else{
							if(Cycle(m, A_NZ, i, j) == 0){
								
								A_NZ[j*m+i]=0;
								A_NZ[i*m+j]=1;
								nonzero++;
							}
							else {
								A_NZ[j*m+i]=0;
								A_NZ[i*m+j]=0;
							}
						}
			    
					}
				}
			}

			// changed & z = (K-nonzero)*tau; */
      		// second part
      		for(pos=(m*m); pos < (m*M); ++pos){
      		  if(ABS(A[pos])>=tau){
      		    A_NZ[pos]=1;
      		  }
      		}	
	
			// check convergence of DC loop
			if(ABS((valueDC - valueDCp)) <= opts_tol && iterDC >= 2) break;
			
		}  // end of DC loop





		*obj = valueDC;  /* update objective value in the output list to R */

		// release memory
		delete [] F;
		delete [] Bab;
		delete [] UC;
		delete [] LLambda;
		delete [] Mat_step;
		delete [] V_step;
		delete [] Xi;
		delete [] y;
		
	}  // end of DAG function
		
	



	
	/* Function to solve the unconvex interventional DAG estimation problem
	 *
	 * Variables
	 *
	 * X:       n by M data matrix
	 * A:       m by M adjacency matrix 
	 * m:       the number of nodes
	 * M:       the number of total nodes (nodes + intervention nodes)
	 * n:       sample size
	 * lambda:  penalty for A
	 * lambda2: penalty for B (intervention part)
	 * tau:     \tau in the paper, tuning parameter for TLP
	 * A_NZ:    m by M matrix (0 indicates the coresponding element of A being 0,
	 *          1 indicates nonzero. -1 indicates a suggested direction, which will be
	 *          checked outside this function. If this suggested direction does not violate
	 *          acyclicity, -1 will be converted to 1.
	 * NZ:      the number of nonzero elements in A. ( the number of 1's in A_NZ)
	 * sigma:   m by 1 variances vector (\sigma in the paper). In this version, they are fixed at 1's.
	 * tol:     tolerance level
	 * obj:     objective function value
	 * XTX:     m by m matrix X^T %*% X
	 * XTX_inv: inverse matrix of XTX, (m-1)*(m-1)*m, different matrices are for calculating different A_row
	 * B:       m by m matrix (B in the paper)
	 * LLambda: m by m dual variable matrix (\lambda in the paper)
	 * Sig:     vector of length m, the error variances sigma^2_j of each node
	 * sigma_X: the error variance of intervention nodes
	 *
	 *
	 * notice that for Xi and y, they are three dimensional and row-major arrays.
	 * A, B, U, LLambda are column-major, as in R, to facilitate the use of fortran functions.
	 *
	 */
	void DAG_int(double *X, double *A, double *B, int *mm, int *MM, int *nn, double *Lambda, double *Lambda2,
			 double *tautau, int *A_NZ, int *NZ, double *sigma, double *sigmaX, double *Sig, double *tol,
			 double *obj, double *XTX, double *XTX_inv, double *X2, double *rhorho, double *rhorho2, int *maxIter)  			
	{

		// take in the values from R
		int m       = (*mm);  // num of nodes
		int M       = (*MM);  // num of nodes + num of intervention modes
		int n       = (*nn);  // num of samples
		int nonzero = (*NZ);  // num of nonzero elements in A

		double lambda   = (*Lambda);
		double lambda2  = (*Lambda2);
		double tau      = (*tautau);    // tuning parameter in TLP function
		double opts_tol = (*tol);       // tolerance
		double mu       = lambda/tau;
		double mu2      = lambda2/tau;
		double rho      = (*rhorho);    // ADMM penalty for DAG constraint
		double rho2     = (*rhorho2);   // ADMM penalty for A=B

		// initialize the parameters
		int i          = 0;
		int j          = 0;
		int k          = 0;
		int pos        = 0;
		int iterDC     = 0;
		int iter       = 0;  // ADMM iteration number
		
		double valueDC     = 0.0;
		double valueDCp    = 0.0;   /* previous value in DC loop */
		double Bab_temp    = 0.0;
		double obj_cur     = 0.0;   /* current objective value  */
		double obj_pre     = 0.0;   /* previous objective value  */
		double res         = 0.0;   /* residuals for primal feasibility condition w.r.t. y */
		double res_temp    = 0.0;	/* intermediate variable to update y */

		double *F        = new double[m*M] ();      // B matrix in the previous iter.
		double *Bab      = new double[m*M] ();    	// absolute value of B
		double *U        = new double[m*M] ();    	// dual variable U
		double *LLambda  = new double[m*m] ();    	// lambda matrix
		double *Mat_step = new double[m*m] ();    	// intermediate matrix when update LLambda
		double *V_step   = new double[m*m] ();    	// intermediate matrix when update LLambda
		double *Xi       = new double[m*m*m] ();  	// kisi in the notes
		double *y        = new double[m*m*m] ();  	// dual variable y

		// initialize parameters
		for(pos=0; pos<m*M; pos++){
			// B[pos]   = 0;
			Bab[pos] = 0;
			U[pos]   = 0;
		}

		for(pos=0; pos<m*m*m; pos++){
			Xi[pos] = 0;
			y[pos]  = 0;
		}

		for(i=0; i<(m*m); i++){
			LLambda[i] = 1; 
			V_step[i]  = 0;
		}

		// initiate Mat_step and V_step
		for(i = 0; i<m; i++){  // row
		    for(j =0; j<m; j++){  // column
				if(j == 0)       Mat_step[j*m + i] = 1;
				else if(i == 0)  Mat_step[j*m + i] = 0;
				else if(i == j)  Mat_step[j*m + i] = 2/((double)m);
				else             Mat_step[j*m + i] = 1/((double)m);	      
		    }	  
		}
		 

		// calculate objective value: 0.5*norm((X*Sig - X*A'), 'fro')-Mu*sum(sum(abs(A)))
		obj_cur = objective(m, M, n, mu, mu2, Sig, X, A, B, A_NZ);
		obj_pre = obj_cur;

		// beginning of DC loop
		for(iterDC = 1; iterDC <= 6; iterDC ++){ //DC loop

		  	// initialize all the dual variables
		  	for(i=0; i<(m*M); i++){
		  		U[i] = 0;
		  	}

		  	for(i=0; i<(m*m*m); i++){
		  		y[i] = 0;
		  	}
	 		

			// beginning of ADMM loop
			for(iter=1; iter < *maxIter; iter++){ // ADMM loop
	
				// step 1: A direction, update by rows

				update_A_convex(m, M, rho2, A, B, U, XTX, XTX_inv, Sig);


				// step 2: update Sig, error variances for the nodes
				update_Sig(m, M, n, sigma, Sig, X, X2, A);


				// step 3: B direction
				for(pos=0; pos<m*M; pos++){
					F[pos] = B[pos];
				}

				// part I, adjacency matrix
				for(pos = 0; pos<(m*m); ++pos){
					j = pos/m;
					i = pos%m;
					Bab_temp = 0.0;

					// calculate the temporary sum: M_{ijk} in notes
					for(k = 0; k<m; k++){
						Bab_temp += LLambda[k*m + i] - LLambda[k*m+j] - Xi[i+j*m+k*m*m] - y[i+j*m+k*m*m];
					}

					Bab_temp = Bab_temp + (m-1)*tau;

					// use the sum to update B
					if(A_NZ[pos] == 1) {
						Bab[pos] = tau;
						B[pos]   = A[pos] + U[pos];
					}
					else {
						Bab[pos] = MAX(0, ((rho2 * ABS(A[pos]+U[pos]) + rho*Bab_temp - mu)/(rho2+m*rho)));
						B[pos]   = Bab[pos] * (A[pos]+U[pos]>=0?1:-1);
					}
					
				}

				// part II, intervention covariates
				for(pos=(m*m); pos<(m*M); ++pos){
					Bab_temp = A[pos]+U[pos];
					if(A_NZ[pos]==1){
						B[pos] = Bab_temp;
						// B[pos] = A[pos] + U[pos];
					}else{
						B[pos] = MAX(0, Bab_temp - mu2/rho2) - MAX(0,-Bab_temp - mu2/rho2);//softthreshold
						// B[pos] = MAX(0,A[pos]+U[pos] - mu2/rho2) - MAX(0,-A[pos]-U[pos] - mu2/rho2);
					}
				} // end of B update



				// step 4: update Lambda matrix
				//initiate V_step;
				for(i = 0; i < m; ++i) {
				  for(k = 0; k < m; ++k){
				    V_step[k*m+i]=0;	   
				    for(j=0; j<m; ++j)   V_step[k*m+i] += Bab[j*m+i] + Xi[i+j*m+k*m*m] + y[i+j*m+k*m*m];		      
				  }
				}
				for(j = 0; j < m; ++j) {
				  for(k = 0; k < m; ++k){	   
				    for(i=0; i<m; ++i)   V_step[k*m+j] -= (Bab[j*m+i] + Xi[i+j*m+k*m*m] + y[i+j*m+k*m*m]);		      
				  }
				}
				for(i = 1; i < m; ++i) {
				  for (k = 0; k < m; ++k){	  	   
				    if(i == k)   V_step[k*m+i] = (-(m-1)*tau + V_step[k*m+i])/2;	
					else         V_step[k*m+i] = (tau + V_step[k*m+i])/2;		      
				  }
				}
				//V_step(1,:) = ones(1,p);
				for(k = 0; k < m; ++k)   V_step[k*m] = 1;
	
				/* LLambda = Mat_step * V_step  */
				dmat_C_AB(m, m, m, Mat_step, V_step, LLambda);
				// end of Lambda update



				// step 5: update Xi and y
				res = 0;
				// r1  = 0;  // primal residual 1
				for(i = 0; i < m; ++i) {
					for(j = 0; j < m; ++j){
				    	if(i != j){
				    		for ( k = 0; k < m; ++k){	      
				    			if(k == j){
				    				/* res_temp is the increment of dual variable y in each step */
				    				/* res is the largest  */
									res_temp        = Bab[j*m+i] + Xi[i+j*m+k*m*m] - LLambda[k*m+i] + LLambda[k*m+j];
									Xi[i+j*m+k*m*m] = MAX(0, (LLambda[k*m+i] - LLambda[k*m+j] - Bab[j*m+i]- y[i+j*m+k*m*m]));
									y[i+j*m+k*m*m]  = y[i+j*m+k*m*m] + res_temp;
									// r1             += res_temp * res_temp;
									// res             = MAX(res, (res_temp*res_temp));
				    		  	}
				    		  	else{
									res_temp        = Bab[j*m+i] + Xi[i+j*m+k*m*m] - LLambda[k*m+i] - tau + LLambda[k*m+j];
									Xi[i+j*m+k*m*m] = MAX(0, (LLambda[k*m+i] + tau - LLambda[k*m+j] - Bab[j*m+i]- y[i+j*m+k*m*m]));
									y[i+j*m+k*m*m]  = y[i+j*m+k*m*m] + res_temp;
									// r1             += res_temp * res_temp;
									// res             = MAX(res, (res_temp*res_temp));
				    		  	}
				    		}
				    	}
				  	}
				}  // end of Xi and y update



				// step 6:  direction
				for(pos = 0; pos<(m*M); pos++){
					U[pos] = U[pos] + A[pos] - B[pos];
				}
	
				// end of all parameter updates


				// update objective value
				obj_pre = obj_cur;
				obj_cur = objective(m, M, n, mu, mu2, Sig, X, A, B, A_NZ);


				// stopping rule for ADMM
				if(ABS((obj_cur-obj_pre)) < opts_tol*obj_pre) break;
			
			}  // end of ADMM

			/* update objective values in DC loop */
			valueDCp = valueDC;//
			valueDC  = obj_cur;

			/* update A_NZ and count num of nonzero elements */
			nonzero = 0;
			for(i=1; i<m; i++){
			    for(j=0; j<i; j++){
					if((ABS(A[j*m+i])+ABS(A[i*m+j])) >= tau*0.99){   
						if(ABS(A[j*m+i]) > tau*0.5){
							if(Cycle(m, A_NZ, j, i) == 0){
								A_NZ[j*m+i]=1;
								A_NZ[i*m+j]=0;
								nonzero++;
							}
							else{      // if (Cycle(m,A_NZ, i,j) == 0) ...
								A_NZ[j*m+i]=0;
								A_NZ[i*m+j]=0;
							}
					    }
					    else{
							if(Cycle(m, A_NZ, i, j) == 0){
								
								A_NZ[j*m+i]=0;
								A_NZ[i*m+j]=1;
								nonzero++;
							}
							else {
								A_NZ[j*m+i]=0;
								A_NZ[i*m+j]=0;
							}
						}
			    
					}
				}
			}
			// changed & z = (K-nonzero)*tau; */
      		// second part
      		for(pos=(m*m); pos < (m*M);++pos){
      		  if(ABS(A[pos])>=tau){
      		    A_NZ[pos]=1;
      		  }
      		}			
			// check convergence of DC loop
			if(ABS((valueDC - valueDCp)) <= opts_tol && iterDC >= 2) break;
			
		}  // end of DC loop



		*obj = valueDC;  /* update objective value in the output list to R */


		// release memory
		delete [] F;
		delete [] Bab;
		delete [] U;
		delete [] LLambda;
		delete [] Mat_step;
		delete [] V_step;
		delete [] Xi;
		delete [] y;
		
	}  // end of DAG function






	/* Function to estimate the DAG structure based on the observational data
	 *
	 * Arguments
	 *
	 * X:       n by m data matrix
	 * A:       m by m adjacency matrix 
	 * m:       the number of variables
	 * n:       sample size
	 * lambda:  penalty (\mu in the paper) 
	 * tau:     \tau in the paper, tuning parameter for TLP
	 * A_NZ:    m by m  matrix (0 indicates the coresponding element of A being 0,
	 *          1 indicates nonzero. -1 indicates a suggested direction, which will be
	 *          checked outside this function. If this suggested direction does not violate
	 *          acyclicity, -1 will be converted to 1.
	 * NZ:      the number of nonzero elements in A. ( the number of 1's in A_NZ)
	 * sigma:   m by 1 variances vector (\sigma in the paper). In this version, they are fixed at 1's.
	 * tol:     tolerance level
	 * obj:     objective function value
	 * XTX:     m by m matrix X^T %*% X
	 * XTX_inv: inverse matrix of XTX, (m-1)*(m-1)*m, different matrices are for calculating different A_row
	 * B:       m by m matrix (B in the paper)
	 * LLambda: m by m dual variable matrix (\lambda in the paper)
	 *
	 *
	 * notice that for Xi and y, they are three dimensional and row-major arrays.
	 * A, B, U,LLambda are column-major, to facilitate the use of fortran functions.
	 *
	 */
	void DAG_obs(double *X, double *A, int *mm, int *nn, double *Lambda,
			 double *tautau, int *A_NZ, int *NZ, double *sigma, double *tol, double *obj,
			 double *XTX, double *XTX_inv, double *rhorho, int *maxIter)  			
	{

		/* take in the values from R */
		int m       = (*mm);  /* num of nodes */
		int m_A     = m-1;    /* num of nodes - 1, for updating A */
		int n       = (*nn);  /* num of samples */
		int nonzero = (*NZ);  /* num of nonzero elements in A */

		double lambda   = (*Lambda);
		double tau      = (*tautau);    /* tuning parameter in TLP function */
		double opts_tol = (*tol);       /* tolerance */
		double mu       = lambda/tau;
		double rho      = (*rhorho);    /* ADMM penalty */

		/* initialize the parameters */
		int i          = 0;
		int j          = 0;
		int k          = 0;
		int pos        = 0;
		int iterDC     = 0;
		int iter       = 0;  /* ADMM */
		
		double valueDC     = 0.0;
		double valueDCp    = 0.0;   /* previous value in DC loop */
		double Bab_temp    = 0.0;
		double obj_cur     = 0.0;   /* current objective value  */
		double obj_pre     = 0.0;   /* previous objective value  */
		double res         = 0.0;   /* residuals for primal feasibility condition w.r.t. y */
		double res_temp    = 0.0;	/* intermediate variable to update y */
		double r2, s2;
		double eps_p2, eps_d2;

	
		// define: B, U, Bab,
		double *B        = new double[m*m] ();
		double *F        = new double[m*m] ();
		double *Bab      = new double[m*m] ();    /* absolute value of B */
		double *U        = new double[m*m] ();    /* dual variable U */
		double *XAT      = new double[n*m] ();    /* X%*%A^T */
		double *LLambda  = new double[m*m] ();    /* lambda matrix */
		double *V_A      = new double[m-1] ();    /* intermediate matrix when update A */
		double *A_row    = new double[m-1] ();    /* row vectors of A matrix */
		double *Mat_step = new double[m*m] ();    /* intermediate matrix when update LLambda */
		double *V_step   = new double[m*m] ();    /* intermediate matrix when update LLambda */
		double *Xi       = new double[m*m*m] ();  /* kisi in the notes */
		double *y        = new double[m*m*m] ();  /* dual variable y */

		for(i=0;i<(m*m);i++)  LLambda[i]=1;   /* initiate LLambda */

		/* initiate Mat_step and V_step */
		for(i = 0; i<m;i++){
		    for(j =0; j<m; j++){
				if(j == 0)       Mat_step[j*m + i] = 1;
				else if(i == 0)  Mat_step[j*m + i] = 0;
				else if(i == j)  Mat_step[j*m + i] = 2/((double)m);
				else             Mat_step[j*m + i] = 1/((double)m);	      
		    }	  
		}
		
		for(i=0; i<(m*m); i++)   V_step[i]=0;

	
		/* calculate objective value: norm((X - X*A'), 'fro')-Mu*sum(sum(abs(A))) */

		/* calculate XAT = X %*% A'  */
		dmat_C_ABT(m, n, m, X, A, XAT);

		obj_cur = 0.0;

		for(pos=0; pos<(n*m); pos++){
			obj_cur += (X[pos] - XAT[pos])*(X[pos] - XAT[pos])/2;
		}

		for(pos=0; pos<(m*m); pos++){
			if(A_NZ[pos] != 1){
				obj_cur += mu*ABS(A[pos]);
			}
		}
		 
		obj_pre = obj_cur;


		/* beginning of DC loop  */
		for(iterDC = 1; iterDC <= 6; iterDC++){ //DC loop
	    
			/* beginning of ADMM loop  */
			for(iter=0; iter < *maxIter; iter++){ // ADMM loop
	
				/* step 1: A direction, update by rows  */
				for(k=0; k<m; k++){  // k is the row number, A, B, U are column-major in C++

					//calculate V_A
					for(i=0; i<(m-1); i++){
						if(i >= k)  j=i+1;
						else        j=i;
						V_A[i] = XTX[j*m+k] + rho*(B[j*m+k] - U[j*m+k]);
					}

					/* calcluate A	*/
					/* A_row = XTXinv * V_A  */
					dmat_yAx(m_A, m_A, &XTX_inv[k*(m-1)*(m-1)], V_A, A_row);

					for(i=0; i<(m-1); i++){
						if(i >= k)  j = i+1;  // j is the column number
						else        j = i;
						A[j*m+k] = A_row[i];
					}
				}  /* end of A update */
	

				/* step 2: B direction  */
				for(pos=0; pos<(m*m); pos++){
					F[pos] = B[pos];
				}
				for(pos = 0; pos<(m*m); ++pos){
					j = pos/m;
					i = pos%m;
					Bab_temp = 0.0;

					/* calculate the temporary sum: M_{ijk} in notes  */
					for(k = 0; k<m; k++){
						Bab_temp += LLambda[k*m + i] - LLambda[k*m+j] - Xi[i*m*m+j*m+k] - y[i*m*m+j*m+k];
					}

					Bab_temp = Bab_temp + (m-1)*tau;

					/* use the sum to update B */
					if(A_NZ[pos] == 1) {  // w_{ij} = 0
						Bab[pos] = tau;
						B[pos]   = A[pos] + U[pos];
					}
					else {
						Bab[pos] = MAX(0, ((rho * ABS(A[pos]+U[pos]) + rho*Bab_temp - mu)/(rho+m*rho)));
						B[pos]   = Bab[pos] * (A[pos]+U[pos]>=0?1:-1);
					}
					
				}  /* end of B update */

				
				/* step 3: update Lambda matrix  */
				//initiate V_step;
				for ( i = 0; i < m; ++i) {
				  for ( k = 0; k < m; ++k){
				    V_step[k*m+i]=0;	   
				    for(j=0; j<m; ++j)   V_step[k*m+i] += Bab[j*m+i] + Xi[i*m*m+j*m+k] + y[i*m*m+j*m+k];		      
				  }
				}
				for ( j = 0; j < m; ++j) {
				  for ( k = 0; k < m; ++k){	   
				    for(i=0; i<m; ++i)   V_step[k*m+j] -= (Bab[j*m+i] + Xi[i*m*m+j*m+k] + y[i*m*m+j*m+k]);		      
				  }
				}
				for ( i = 1; i < m; ++i) {
				  for ( k = 0; k < m; ++k){	  	   
				    if(i == k)   V_step[k*m+i] = (-(m-1)*tau + V_step[k*m+i])/2;	
					else         V_step[k*m+i] = (tau + V_step[k*m+i])/2;		      
				  }
				}
				//V_step(1,:) = ones(1,p);
				for ( k = 0; k < m; ++k)   V_step[k*m] = 1;
	
				/* LLambda = Mat_step * V_step  */
				dmat_C_AB(m, m, m, Mat_step, V_step, LLambda);
				/* end of Lambda update */


				/* step 4 & 5: update Xi and y */
				res = 0;
				for(i = 0; i < m; ++i) {
					for(j = 0; j < m; ++j){
				    	if(i != j){
				    		for ( k = 0; k < m; ++k){	      
				    			if(k == j){
				    				/* res_temp is the increment of dual variable y in each step */
				    				/* res is the largest  */
									res_temp        = Bab[j*m+i] + Xi[i*m*m+j*m+k] - LLambda[k*m+i] + LLambda[k*m+j];
									Xi[i*m*m+j*m+k] = MAX(0, (LLambda[k*m+i] - LLambda[k*m+j] - Bab[j*m+i]- y[i*m*m+j*m+k]));
									y[i*m*m+j*m+k]  = y[i*m*m+j*m+k] + res_temp;
									res             = MAX(res, (res_temp*res_temp));
				    		  	}
				    		  	else{
									res_temp        = Bab[j*m+i] + Xi[i*m*m+j*m+k] - LLambda[k*m+i] - tau + LLambda[k*m+j];
									Xi[i*m*m+j*m+k] = MAX(0, (LLambda[k*m+i] + tau - LLambda[k*m+j] - Bab[j*m+i]- y[i*m*m+j*m+k]));
									y[i*m*m+j*m+k]  = y[i*m*m+j*m+k] + res_temp;
									res             = MAX(res, (res_temp*res_temp));
				    		  	}
				    		}
				    	}
				  	}
				}  /* end of Xi and y update */


				/* step 6: U direction */
				for(pos = 0; pos<(m*m); ++pos){
					U[pos] = U[pos] + A[pos] - B[pos];
				}  /* end of U update */
	
				/* end of all parameter updates */

				/* update objective value  */
				obj_pre = obj_cur;
	
				/* calculate XAT = XAT */
				dmat_C_ABT(m, n, m, X, A, XAT);
	
				obj_cur = 0.0;
	
				for(pos=0; pos<(n*m); pos++){
					obj_cur += (X[pos] - XAT[pos])*(X[pos] - XAT[pos])/2 ;
				}
	
				for(pos=0; pos<(m*m); pos++){
					if(A_NZ[pos] != 1){
						obj_cur += mu*ABS(A[pos]);
					}
				}
				
				// compute residuals
				r2 = 0.0;
				s2 = 0.0;

				for(pos=0; pos<m*m; pos++){
					r2 += (A[pos] - B[pos]) * (A[pos] - B[pos]);
					s2 += (B[pos] - F[pos]) * (B[pos] - F[pos]);
				}

				r2 = pow(r2, 0.5);
				s2 = pow(s2, 0.5) * rho;

				// compute the stopping criterion
				eps_p2 = 0.0;
				eps_d2 = 0.0;

				for(pos=0; pos<m*m; pos++){
					eps_p2 += A[pos] * A[pos];
					eps_d2 += B[pos] * B[pos];
				}
				eps_p2 = MAX(eps_p2, eps_d2);  // max of A norm and B norm

				eps_p2 = pow(m*m, 0.5) * eps_abs + eps_rel * pow(eps_p2, 0.5);
				eps_d2 = pow(m*m, 0.5) * eps_abs + eps_rel * pow(eps_d2, 0.5);

				if(r2 <= eps_p2 && s2 <= eps_d2) break;
			
			}  /* end of ADMM  */

			/* update objective values in DC loop */
			valueDCp = valueDC;
			valueDC  = obj_cur;

			/* update A_NZ and count num of nonzero elements */
			nonzero = 0;
			for(i=1; i<m; i++){
			    for(j=0; j<i; j++){
					if((ABS(A[j*m+i])+ABS(A[i*m+j])) >= tau*0.99){   
						if(ABS(A[j*m+i]) > tau*0.5){
							if( Cycle(m, A_NZ, j, i) == 0){
								A_NZ[j*m+i]=1;
								A_NZ[i*m+j]=0;
								nonzero++;
							}
							else{      // if (Cycle(m,A_NZ, i,j) == 0) ...
								A_NZ[j*m+i]=0;
								A_NZ[i*m+j]=0;
							}
					    }
					    else{
							if(Cycle(m, A_NZ, i, j) == 0){
								
								A_NZ[j*m+i]=0;
								A_NZ[i*m+j]=1;
								nonzero++;
							}
							else {
								A_NZ[j*m+i]=0;
								A_NZ[i*m+j]=0;
							}
						}
					}
				}
			}
			// changed & z = (K-nonzero)*tau; */
			
			/* check convergence of DC loop */
			if(ABS((valueDC - valueDCp)) <= opts_tol && iterDC >= 2) break;
			
		}  /* end of DC loop */



		*obj = valueDC;  /* update objective value in the output list to R */

		/* release memory */
		delete [] B;
		delete [] Bab;
		delete [] U;
		delete [] V_A;
		delete [] A_row;
		delete [] XAT;
		delete [] LLambda;
		delete [] Mat_step;
		delete [] V_step;
		delete [] Xi;
		delete [] y;
		
	}  /* end of DAG function */

} /* extern "C"  */

