#include <mex.h>

#include <math.h>

#define D 2


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

void forces(double *r, double *v, double *f, double m, double n)
{
	double *dl, *dr;// *mdl, *mdr;
	double L, Lo, nP, k, b, be, stretch;		
	int i, j, ii, N;
	
	
    
	dl = (double *) mxMalloc(D*n*sizeof(double));
	dr = (double *) mxMalloc(D*n*sizeof(double));
	//mdl = (double *) mxMalloc(n*sizeof(double));
	//mdr = (double *) mxMalloc(n*sizeof(double));
	
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
	for(i = 0, j = 1; i < D*n; i+=D){
		j = i + 1;
		//mdl[ii] = sqrt(dl[i]*dl[i] + dl[j]*dl[j]);
		//mdr[ii] = sqrt(dr[i]*dr[i] + dr[j]*dr[j]);
		
		f[j] = k*(dr[j] + dl[j] - Lo*(dr[j]/sqrt(dr[i]*dr[i] + dr[j]*dr[j])+ 
		       dl[j]/sqrt(dl[i]*dl[i] + dl[j]*dl[j])));
		       
		f[i] = k*(dr[i] + dl[i] - Lo*(dr[i]/sqrt(dr[i]*dr[i] + dr[j]*dr[j])+ 
		       dl[i]/sqrt(dl[i]*dl[i] + dl[j]*dl[j]))) - b*v[i];
		ii++;
	}	
	
	mxFree(dl);
	mxFree(dr);
	//mxFree(mdl);
	//mxFree(mdr);
			
}	

void rkstep(double *rin, double *vin, double *rout, 
			double *vout, double h, int m, int n)
/* 	function rkstep():  
   
	inputs:
   	rin: position vector 
   	vin: velocity vector 
   	h: step size
    m: number of columns
    n: number of rows
    
   	outputs:
   	rout: next position
   	vout: next velocity
*/
{
	int i, k;
	double *f1, *g1, *f2, *g2, 
		   *f3, *g3, *f4, *g4, *t1;
		   
	k = n*D;	
	    
	g1 = (double *)mxMalloc(k*sizeof(double));	   
	f2 = (double *)mxMalloc(k*sizeof(double));	   
	g2 = (double *)mxMalloc(k*sizeof(double));	   
	f3 = (double *)mxMalloc(k*sizeof(double));	   
	g3 = (double *)mxMalloc(k*sizeof(double));	   
	f4 = (double *)mxMalloc(k*sizeof(double));	   
	g4 = (double *)mxMalloc(k*sizeof(double));	   
	t1 = (double *)mxMalloc(k*sizeof(double));	   
	
	forces(rin, vin, g1, m, n);
	
	
	for(i = 0; i < k; i++)
	{
		f2[i] = vin[i] + 0.5*h*g1[i];	
		t1[i] = rin[i] + 0.5*h*vin[i];		
	}
	forces(t1, f2, g2, m, n);
	
	for(i = 0; i < k; i++)
	{
		f3[i] = vin[i] + 0.5*h*g2[i];
		t1[i] = rin[i] + 0.5*h*f2[i];			
	}		 
	forces(t1, f3, g3, m, n);
	
	for(i = 0; i < k; i++)
	{
		f4[i] = vin[i] + h*g3[i];	
		t1[i] = rin[i] + h*f3[i];	
	}		 
	forces(t1, f4, g4, m, n);
	for(i = 0; i < k; i++)
	{
		rout[i] = rin[i] + (1.0/6.0)*h*(vin[i] + 2.0*f2[i] + 2.0*f3[i] + f4[i]);
		vout[i] = vin[i] + (1.0/6.0)*h*(g1[i] + 2.0*g2[i] + 2.0*g3[i] + g4[i]);	
	}	
	
	mxFree(g1); 
	mxFree(f2); 
	mxFree(g2); 
	mxFree(f3); 
	mxFree(g3); 
	mxFree(f4); 
	mxFree(g4); 
	mxFree(t1); 
		
}	



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{ 
	double *r, *v, *rout, *vout, h;
	int m, n, i;
    if (nrhs != 3)
        mexErrMsgTxt("3 ins pls!");
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    if(m != D) mexErrMsgTxt("2 dimensions");
    
    
    r = mxGetPr(prhs[0]);
    v = mxGetPr(prhs[1]);
    h = mxGetScalar(prhs[2]);
   
    
    rout = (double *) mxMalloc(m*n * sizeof(double));
    vout = (double *) mxMalloc(m*n * sizeof(double));
    
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(m, n, mxREAL);
    rout = mxGetPr(plhs[0]);
    vout = mxGetPr(plhs[1]);
    
    
    rkstep(r, v, rout, vout, h, m, n);
    	
}   