
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

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


void print_vec_2d(float_data * vec, int L, char * name){
	
	int k, j;
	
	printf("%s:\n", name);
	
	for (k = 0; k < L; k += 2){
		printf("[%2f, %2f], ", vec[k], vec[k + 1]);
	}
	
	printf("\n");
	
}

void print_vec(float_data * vec, int L, char * name){
	
	int k;
	
	printf("%s:\n", name);
	
	for (k = 0; k < L; k++){
		printf("%2f, ", vec[k]);
	}
	
	printf("\n");
	
}

void diff(float_data *x, float_data endy, float_data endx, int n, int BorE, float_data *output)
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

void forces(float_data *r, float_data *v, float_data *f, float_data *dl, 
			float_data *dr, float_data *md, int n, float_data L, float_data *b,
			float_data nP, float_data Lo, float_data K)
{		
	int i, j, k;
	
	diff(r, 0.0, 0.0, n, 0, dl);
	diff(r, 0.0, L+L/nP, n, 1, dr); 
	
	/*
	print_vec_2d(dl, 2*n, "dl");
	print_vec_2d(dr, 2*n, "dr");
	*/
	
	for(i = 0, j = 1,k = 0; i < D * n; i+=D, j+=D, k++){
				 
		md[i] = sqrt(dl[i]*dl[i] + dl[j]*dl[j]);
		md[j] = sqrt(dr[i]*dr[i] + dr[j]*dr[j]);
		f[j] = nP*(K*(dr[j] + dl[j] - Lo*(dr[j]/md[j] + dl[j]/md[i])) - b[k]*v[i]);	       
		f[i] = nP*(K*(dr[i] + dl[i] - Lo*(dr[i]/md[j] + dl[i]/md[i])) - b[k]*v[i]);
		/*
		f[j] += 0.1*f[i];
		f[i] += 0.1*f[j];
		*/
	}
	/*
	print_vec_2d(md, 2*n, "md");
	print_vec_2d(f, 2*n, "f");
	*/
}	

/*
 * sim - run simulation of string.  
 * r = initial position vector
 * v = initial velocity vector
 * aud = output audio
 * damp = damping vector
 * h = Runge-Kutta step-size
 * n = number of string points
 * nS = number of iterations to run for
 * nP = number of string points
 * L = the length of the string
 * K = the spring constant
 * S = the stretch factor (String tension)
 */

void sim(float_data *r, float_data *v, float_data *aud, float_data *damp,
			float_data h, int n, int nS, float_data nP, float_data L,
			float_data K, float_data S)

{
	int i, j, k;
	static float_data *f1=0, *g1=0, *f2=0, *g2=0, 
		   *f3=0, *g3=0, *f4=0, *g4=0, 
		   *t1=0, *dl=0, *dr=0, *md=0;
	
	float_data Lo = L/(nP*S);	   
		   
	k = n * 2;
		
	g1 = (float_data *)realloc(g1, k*sizeof(float_data));	   
	f2 = (float_data *)realloc(f2, k*sizeof(float_data));	   
	g2 = (float_data *)realloc(g2, k*sizeof(float_data));	   
	f3 = (float_data *)realloc(f3, k*sizeof(float_data));	   
	g3 = (float_data *)realloc(g3, k*sizeof(float_data));	   
	f4 = (float_data *)realloc(f4, k*sizeof(float_data));	   
	g4 = (float_data *)realloc(g4, k*sizeof(float_data));	   
	t1 = (float_data *)realloc(t1, k*sizeof(float_data));
	
	dl = (float_data *)realloc(dl, k*sizeof(float_data));
	dr = (float_data *)realloc(dr, k*sizeof(float_data));
	md = (float_data *)realloc(md, k*sizeof(float_data));
	
	
	/* simulate for nS steps */	   
	for(j = 0; j < nS; j++){
		
		/* write the current sample to the audio data array */
		if (aud)
			aud[j] = r[4];
		
		/* invoke the fourth-order Runge_kutta method
		   to find the next state 
		*/
		/*
		printf("-----------------------------------------\n\n");
		*/
		forces(r, v, g1, dl, dr, md, n, L, damp, nP, Lo, K);
		/*
		print_vec_2d(g1, k, "g1");
		*/
		
		for(i = 0; i < k; i++)
		{
			f2[i] = v[i] + 0.5*h*g1[i];	
			t1[i] = r[i] + 0.5*h*v[i];		
		}
		forces(t1, f2, g2, dl, dr, md, n, L, damp, nP, Lo, K);
		for(i = 0; i < k; i++)
		{
			f3[i] = v[i] + 0.5*h*g2[i];
			t1[i] = r[i] + 0.5*h*f2[i];			
		}		 
		forces(t1, f3, g3, dl, dr, md, n, L, damp, nP, Lo, K);
		for(i = 0; i < k; i++)
		{
			f4[i] = v[i] + h*g3[i];	
			t1[i] = r[i] + h*f3[i];	
		}		 
		forces(t1, f4, g4, dl, dr, md, n, L, damp, nP, Lo, K);
		for(i = 0; i < k; i++)
		{
			r[i] = r[i] + (1.0/6.0)*h*(v[i] + 2.0*f2[i] + 2.0*f3[i] + f4[i]);
			v[i] = v[i] + (1.0/6.0)*h*(g1[i] + 2.0*g2[i] + 2.0*g3[i] + g4[i]);	
		}	
	}
	
}	


#ifdef MATLAB
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{ 
	float_data *r, *v, *aud, *b, h, L, nP, K, S;
	int m, n, nS;
    if (nrhs != 8)
        mexErrMsgTxt("8 ins pls!");
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    nP = (float_data) n;
    r = mxGetPr(prhs[0]);
    v = mxGetPr(prhs[1]);
    b = mxGetPr(prhs[2]);
    h = mxGetScalar(prhs[3]);
    nS = (int) mxGetScalar(prhs[4]);
    L = mxGetScalar(prhs[5]);
    K = mxGetScalar(prhs[6]);
    S = mxGetScalar(prhs[7]);
    aud = (float_data *) mxMalloc(nS*sizeof(float_data)); 
    plhs[0] = mxCreateDoubleMatrix(1, nS, mxREAL);
    aud = mxGetPr(plhs[0]);
	
	if (r[2] > 10) mexErrMsgTxt("step size too big"); 
   
    sim(r, v, aud, b, h, n, nS, nP, L, K, S);
    
}   
#endif
