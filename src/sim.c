#include <stdlib.h>
#include <math.h>

#define D 2

#ifdef MATLAB
	#include <mex.h>
	#define malloc mxMalloc
	#define free free
	#define calloc mxCalloc
	#define realloc mxRealloc
#endif

#define float_data float
#define int_data int16_t


void diff(double *x, double endy, double endx, int n, int BorE, double *output)
{
	int i;
	
	if (BorE){
		i = 0;
		for(; i < (2*n)-2; i++){
			output[i] = x[i+2] - x[i];
		}
		
		output[i] = endy - x[i]; i++;
		output[i] = endx - x[i];
	}
	else{
		i = 0;
		output[i] = 0.0 - x[i]; i++;
		output[i] = 0.0 - x[i]; i++;
		
		for(; i < 2*n; i++){
			output[i] = x[i-2] - x[i];
		}
	}
	
}

void forces(double *r, double *v, double *f, double *dl, 
			double *dr, double n)
{

	double L, Lo, nP, k, b, be, stretch;		
	int i, j, ii, N;
	
	k = 0.25;
	b = 0.00001;
	
	N = 2;
	be = 0.0008;
	L = 50;
	nP = (double) n;
	stretch = 1.1;
	
	Lo = L/(nP*stretch);
	
	diff(r, 0.0, 0.0, n, 0, dl);
	diff(r, 0.0, (L + L/nP), n, 1, dr);
	
	
	ii = 0;
	for(i = 0, j = 1; i < D*n; i+=D)
	{
		j = i + 1;
		
		f[j] = k*(dr[j] + dl[j] - Lo*(dr[j]/sqrt(dr[i]*dr[i] + dr[j]*dr[j])+ 
		       dl[j]/sqrt(dl[i]*dl[i] + dl[j]*dl[j])));
		       
		f[i] = k*(dr[i] + dl[i] - Lo*(dr[i]/sqrt(dr[i]*dr[i] + dr[j]*dr[j])+ 
		       dl[i]/sqrt(dl[i]*dl[i] + dl[j]*dl[j]))) - b*v[i];
		ii++;
	}	
	
}	

void sim(double *r, double *v, double *aud, 
			double h, int n, int nS, int nP, double L)

{
	int i, j, k;
	double *f1, *g1, *f2, *g2, 
		   *f3, *g3, *f4, *g4, 
		   *t1, *dl, *dr;
		   
	k = n*D;	
	    
	g1 = (double *)mxMalloc(k*sizeof(double));	   
	f2 = (double *)mxMalloc(k*sizeof(double));	   
	g2 = (double *)mxMalloc(k*sizeof(double));	   
	f3 = (double *)mxMalloc(k*sizeof(double));	   
	g3 = (double *)mxMalloc(k*sizeof(double));	   
	f4 = (double *)mxMalloc(k*sizeof(double));	   
	g4 = (double *)mxMalloc(k*sizeof(double));	   
	t1 = (double *)mxMalloc(k*sizeof(double));
	
	dl = (double *)mxMalloc(k*sizeof(double));
	dr = (double *)mxMalloc(k*sizeof(double));
	
		   
	for(j = 0; j < nS; j++){
		forces(r, v, g1, dl, dr, n);
		for(i = 0; i < k; i++)
		{
			f2[i] = v[i] + 0.5*h*g1[i];	
			t1[i] = r[i] + 0.5*h*v[i];		
		}
		forces(t1, f2, g2, dl, dr, n);
		for(i = 0; i < k; i++)
		{
			f3[i] = v[i] + 0.5*h*g2[i];
			t1[i] = r[i] + 0.5*h*f2[i];			
		}		 
		forces(t1, f3, g3, dl, dr, n);
		for(i = 0; i < k; i++)
		{
			f4[i] = v[i] + h*g3[i];	
			t1[i] = r[i] + h*f3[i];	
		}		 
		forces(t1, f4, g4, dl, dr, n);
		for(i = 0; i < k; i++)
		{
			r[i] = r[i] + (1.0/6.0)*h*(v[i] + 2.0*f2[i] + 2.0*f3[i] + f4[i]);
			v[i] = v[i] + (1.0/6.0)*h*(g1[i] + 2.0*g2[i] + 2.0*g3[i] + g4[i]);	
		}	
	}
		
	mxFree(g1); 
	mxFree(f2); 
	mxFree(g2); 
	mxFree(f3); 
	mxFree(g3); 
	mxFree(f4); 
	mxFree(g4); 
	mxFree(t1); 
	mxFree(dl);
	mxFree(dr);	
}	


#ifdef MATLAB
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{ 
	double *r, *v, *aud, h, L;
	int m, n, i, nS, nP;
    if (nrhs != 6)
        mexErrMsgTxt("6 ins pls!");
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    
    r = mxGetPr(prhs[0]);
    v = mxGetPr(prhs[1]);
    h = mxGetScalar(prhs[2]);
    nS = (int) mxGetScalar(prhs[3]);
    nP = (int) mxGetScalar(prhs[4]);
    L = mxGetScalar(prhs[5]);
   
    
    aud = (double *) mxMalloc(nS * sizeof(double));
    
    
    plhs[0] = mxCreateDoubleMatrix(1, nS, mxREAL);
    
    aud = mxGetPr(plhs[0]);
    
    
    
    sim(r, v, aud, h, n, nS, nP, L);
    	
}   
#endif
