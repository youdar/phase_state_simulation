/* =============================================================
 * AggKineticLoop.c is a functions that preform most of the calculation 
 * For the Aggregation Kinetic software. The initial input, parameter
 * normalization and result output is done at the Matlab AggKinetic program
 *
 * Program Done by Youval Dar
 * 2010 UC Davis Physics
 *
 * 
 * 
 * ================================================================ */
/* $Revision: 0.0 $ */

/* The computational routine */
# include "mex.h"
# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <time.h>  

/*******    Functions  ********/

double SumArray(mxArray *a, double *aValues)
{
    /*
	* a = prhs[?];
    * baValues = mxGetPr(a);
    * 1D arrays only  !!!
    * Need to convert from 2D to 1D according the the mxArrays scheme
    */
    int i;
    double tot = 0.0;
    mwSize aLength;
	
	aLength = mxGetNumberOfElements(a);
	
    for (i = 0; i< aLength; ++i)
    {
        tot += aValues[i];
    }
    return tot;
}


double SumArrayRow(double *b, int row, int ColStart, int rowSize, int colSize)
{
    /*
    * Sum elements of one row, starting at ColStart
    * b = mxGetPr(a);
    * rowSize is the number of rows in the matrix
    * colSize is the number of colomns in th ematrix
    * Need to convert from 2D to 1D according the the mxArrays scheme
    * 1D arrays only  !!! rowSize*colSize is the number of elements in the array
    */
    int i,iStart,iEnd;
    double tot=0.0;
    
    iStart = FindPos(row,ColStart,rowSize);
    iEnd = rowSize*colSize;
    for (i = iStart; i < iEnd; i += rowSize)
    {
        tot += b[i];
    }
    return tot;
}


int FindPos(int row,int col,int rowSize)
{
    /* the mxArrays converts 2D matrices to 1D arrays by putting colomn after colomn
    * This function converts the matlab matrix position to the location in the array  
    *  
    * mxGetM (C and Fortran) Number of rows in mxArray 
    * mxGetN (C and Fortran) Number of columns in mxArray 
    */
    return col * rowSize + row;
}


double RandUn(void)
{
    /* Uniform random generator
     * on the range (0,1]
     * RAND_MAX = 32767
     * x = rand(0), {x: 0 <= x <= RAND_MAX}  in range [0,RAND_MAX]
    */
    return (rand() + 1.0)/(RAND_MAX + 1.0);
}

int RandInt(int x)
{
    /* Uniform integer random generator
     * on the range [1,x]
     * RAND_MAX = 32767
     * x = rand(0), {x: 0 <= x <= RAND_MAX}  in range [0,RAND_MAX]
    */
    return 1 + x*rand()/(RAND_MAX + 1.0);
}

void a3Update(int agn, double aa[], double a[], int aRowSize, double y[], double cbf)
{
    /*  Update a3
    *   a3 = (*cbf)*y(4:n).*(i-3)
    *   a3Update(agn,aa,a,y,*cbf);   */
    double aold;
	int ii;
		
	ii = FindPos(2,agn,aRowSize);
    aold = a[ii];
    a[ii] = cbf*y[agn]*(agn-2);
    aa[2] +=  a[ii] - aold;
}

void a4Update(int agn, double aa[], double a[], int aRowSize, double y[], double clmf)
{
    /* update a4 
    *  a4 = 2*clmf*y(2:n)
    *  a4Update(agn,aa,a,aRowSize,y,*clmf);   */
    double aold;
	int ii;
		
	ii = FindPos(3,agn,aRowSize);
    aold  = a[ii];
    a[ii] =  2*clmf*y[agn];
    aa[3] += a[ii] - aold;
}


void a7Update(int agn, double aa[], double a[], int aRowSize, double g[], double cbm)
{
    /* update a7 
    *  a7 = cbm*g(3:n).*(i-2)
    *  a7Update(agn,aa,a,aRowSize,g,*cbm);   */
    double aold;
	int ii;
		
	ii = FindPos(6,agn,aRowSize);
    aold =  a[ii];
    a[ii]= cbm*g[agn]*(agn-1);
    aa[6] += a[ii] - aold;
}


void a8Update(int agn, double aa[], double a[], int aRowSize, double g[], double clmm)
{
    /* update a8 
    *  a8 = clmm*g(2:n)
    *  a8Update(agn,aa,a,aRowSize,g,*clmm);   */
    double aold;
	int ii;
		
	ii = FindPos(7,agn,aRowSize);
    aold = a[ii];
    a[ii]= clmm*g[agn];
    aa[7] += a[ii] - aold ;
}

void a9Update(int agn, double aa[], double a[], int aRowSize, double g[], double cam)
{
    /* update a9
    *  a9 = cam*g(1:n)
    *  a9Update(agn,aa,a,aRowSize,g,*cam);   */
    double aold;
	int ii;
		
	ii = FindPos(8,agn,aRowSize);
    aold = a[ii];
    a[ii] = cam*g[agn];
    aa[8] = a[ii] - aold;
}




/************ The main mex function **************/
void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    
    /* Input variables 
     * (t,F,M,y0t,g0t,y,g,inagg,a0,a,aa,rst,re,yn,gn,clpf,clmf,clpm,clmm,cap,cam,cn,cnm,cbf,cbm,y0,p0,rp,
     *                                 crs,crp,rfp,crfm,crfp,m,n,dfy,dfg)   */

    double *t,*F,*M,*y0t,*g0t,*y,*g,*inagg,*a0,*a,*aa,*rst,*re,*yn,*gn;
    double *clpf,*clmf,*clpm,*clmm,*cap,*cam,*cn,*cnm,*cbf,*cbm,*y0,*p0,*rp;
    double *crs,*crp,*rfp,*crfm,*crfp,*m,*nn,*dfy,*dfg;
    int i,rn,n,dy1,dy2,yi,dg01,gi,g0,dg2,dfg0;
    int agn1,agn2,agn3,ii,jj,aRowSize,aColSize,iiStart,iiEnd,ConstFactor;
    double dt,r,s,dyin,dg0in,rg,sg,ry,sy;
    mxArray *a0Data;
	mwSize aLength;

    
    /* Output variables
     * When changing arrays in the MEX code, the array in Matlab is changes as well
     * so there is no need to send is as an ouput
     */
    
    /* check for proper number of arguments */
	if(nrhs!=37) 
    {
        // mexErrMsgIdAndTxt("AggKineticLoop","Error - Check input parameters.");
        mexPrintf("AggKineticLoop - Error - Check input parameters.");
    }
    
    /* setting the input values */
    t = mxGetPr(prhs[0]);
    F = mxGetPr(prhs[1]);
    M = mxGetPr(prhs[2]);
    y0t = mxGetPr(prhs[3]);
    g0t = mxGetPr(prhs[4]);
    y = mxGetPr(prhs[5]);
    g = mxGetPr(prhs[6]);
    inagg = mxGetPr(prhs[7]); 
    a0Data = prhs[8];
	a0 = mxGetPr(a0Data);
	aLength = mxGetNumberOfElements(a0Data);
    a = mxGetPr(prhs[9]);
    aRowSize = mxGetM(prhs[9]);
    aColSize = mxGetN(prhs[9]);
    aa = mxGetPr(prhs[10]);
    rst = mxGetPr(prhs[11]);
    re = mxGetPr(prhs[12]);
    yn = mxGetPr(prhs[13]);
    gn = mxGetPr(prhs[14]);
    clpf = mxGetPr(prhs[15]);
    clmf = mxGetPr(prhs[16]);
    clpm = mxGetPr(prhs[17]);
    clmm = mxGetPr(prhs[18]);
    cap = mxGetPr(prhs[19]);
    cam = mxGetPr(prhs[20]);
    cn = mxGetPr(prhs[21]);
    cnm = mxGetPr(prhs[22]);
    cbf = mxGetPr(prhs[23]);
    cbm = mxGetPr(prhs[24]);
    y0 = mxGetPr(prhs[25]);
    p0 = mxGetPr(prhs[26]);
    rp = mxGetPr(prhs[27]);
    crs = mxGetPr(prhs[28]);
    crp = mxGetPr(prhs[29]);
    rfp = mxGetPr(prhs[30]);
    crfm = mxGetPr(prhs[31]);
    crfp = mxGetPr(prhs[32]);
    m = mxGetPr(prhs[33]);
    nn = mxGetPr(prhs[34]);
    dfy = mxGetPr(prhs[35]);
    dfg = mxGetPr(prhs[36]);
    
    /* converting to integers */
    n = *nn;            // Size of output arrays
    g0 = *p0;           // initial amount of monomers 
    
    /* Seed random numbers */
    srand(time(NULL));
    
    

    /*  Start main loop */
    /* Remember that in C array starts at 0 not at 1, like in Matlab*/
    for (i = 1; i < *m; ++i)
    {
        // Add time step
        dt = log(1/RandUn())/(*a0);
        t[i] = t[i-1]+ dt;   

        // Get the reaction type number                           
        r = RandUn();           // Random number (0,1]
        rn = 0;             // the reaction number in a (0 to 9)
        s = aa[rn]/(*a0);      // accomulating the propensities
//         while (s < r )
//         {
//             rn += 1;
//             if (rn > 9)
//             {
//                  // mexErrMsgIdAndTxt("AggKineticLoop","Error - rn>10");
//                 mexPrintf("AggKineticLoop - Error - rn>10.");
//             }
//             s += aa[rn]/(*a0);  
//         }

        // find the reaction within the reaction type
        // agn1 is the possition in the y or g arrays
        agn1 = rst[rn];                             // the first j in the correct region
        ii = FindPos(rn,agn1,aRowSize);
        if (rn!=0 && rn!=9)
        {
            s += a[ii]/(*a0) - aa[rn]/(*a0);  // return value to begining of reaction
//             while (s < r)
//             {
//                 agn1 += 1;
//                 ii += aRowSize; 
//                 if (ii > 7999)
//                 {
//                      // mexErrMsgIdAndTxt("AggKineticLoop","Error - ii>7999");
//                     mexPrintf("AggKineticLoop - Error - ii>7999.");
//                 }
//                 s += a[ii]/(*a0); 
//             }
        }
    
    
        // The deterministic in-flow
        // Adjust y[0] and g0 
        dyin += dt*(*crfp);        // add accumulated fraction of monomers
        dy1 = dyin;             // moving only whole particle
        dyin -= dy1;            // keeping track of fraction of particles

        // Out-flow calculation
        ry = (*yn)*RandUn();       // using random selection according the aggregation number distribution
        sy = y[0];
        yi = 0;                 // aggregate to be flowing out
//         while (sy<ry)
//         {
//             yi += 1;
//             if (yi > 799)
//                 {
//                      // mexErrMsgIdAndTxt("AggKineticLoop","Error - yi>799");
//                     mexPrintf("AggKineticLoop - Error - yi>799.\n");
//                 }
//             sy += y[yi];
//         }

        // Out-flow calculation -  update y[yi]
        dfy[yi] += dt*y[yi]*(*crfm);
        dy2 = dfy[yi];          // the integer number of aggregates or monomers the flow out
        dfy[yi] -= dy2;         // update the reminder of the aggregates or monomers, that did not flow out

        if (yi==0)
        {
            y[0] += dy1 - dy2;            // Update free monomers number
        }
        else
        {
            y[0] += dy1;
            y[yi] -= dy2;
        }
        *yn += dy1 - dy2;                // Update number of free particles (aggregates and monomers)
        if (y[yi]<0)
        {
            *yn -= y[yi];
            y[yi] = 0;   
        }
            
    
        // Membrane protein production
        dg0in += dt*(*crp);     // add accumulated fraction of membrane proteins
        dg01 = dg0in;           // moving only whole particle
        dg0in -= dg01;          // keeping track of fraction of particles
        // Membrane sequestration
        rg = (*gn)*RandUn();       // using random selection according the aggregation number distribution
        sg = g0;                // initially
        gi = 0;                 // aggregate to be flowing out
//         while (sg<rg)
//         {
//             sg += sg + g[gi];
//             gi += 1;
//             if (gi > 799)
//             {
//                  mexErrMsgIdAndTxt("AggKineticLoop","Error - gi>799");
//                 mexPrintf("AggKineticLoop - Error - gi>799.\n");
//                 mexPrintf("i = %i \n",i);
//                 mexPrintf("rn = %i \n",rn);
//             }
//         }
        
        
        // update g0 for the in-flow
        // update g[gi]
        if (gi>0)
        {
            gi -= 1;                        // the arrays index starts at 0, not like matlab
            dfg[gi] += dt*g[gi]*(*crs);
            dg2 = dfg[gi];                  // the integer number of aggregates or monomers the flow out
            g[gi] -= dg2;                   // update g
            g0 += dg01;                 
            dfg[gi] -= dg2;                 // update the reminder of the aggregates or monomers, that did not flow out
            if (g[gi]<0)
            {
                *gn += dg01 - dg2 - g[gi]; 
                inagg[i] = inagg[i-1] + (gi+1)*(dg2+g[gi]);  
                g[gi] = 0;
            }
            else
            {
                *gn += dg01 - dg2; 
                inagg[i] = inagg[i-1] + (gi+1)*dg2;  
            } 
        }
        else
        {
            dfg0 += dt*g0*(*crs);
            dg2 = dfg0;                 // the integer number of aggregates or monomers the flow out
            g0 += dg01 - dg2;           // update g0
            dfg0 -= dg2;                // update the reminder of the aggregates or monomers, that did not flow out 
            inagg[i] = inagg[i-1];
            if (g0<0)
            {
                *gn += dg01 - dg2 - g0; 
                g0 = 0;
            }
            else
            {
                *gn += dg01 - dg2;
            }
        }
                    
        if (*gn<0)
        {
			// mexErrMsgIdAndTxt("AggKineticLoop","Error - gn < 0");
            mexPrintf("AggKineticLoop - Error - gn<0\n");
        }
		
		// Adjust M and F for the out-flow
		M[i] -= dy2*(yi+1)*(yi>0);
		F[i] -= dg2*(gi+1);      // Note that we need to use the readjusted gi
		
		// Calculate the changes acording to the reaction that took place   
		switch (rn)
        {
			case 0:
				// Calculate y and g
				if (y[0]>2)
                {
					y[0] -= 2;
					*yn -= 1;
					y[1] += 1;
					F[i] = F[i-1]+ 2;
					// reaction propensities that need update
					// a1 = cn*y[0]*(y[0]-1)/2;     Updated at the end
					// a2 = clpf*y[0]*y(2:(n-1));   Updated at the end 
					// a4 = 2*clmf*y(2:n);
                    a4Update(1,aa,a,aRowSize,y,*clmf);
					// a5 = cap*g0*y;                   Updated at the end
					// a6 = clpm*y[0]*g(1:(n-1));       Updated at the end
					// a10 = cnm*g0*y[0]*(y[0]-1)/2     Updated at the end
                }
				else
                {
					F[i] = F[i-1];
                }
				M[i] = M[i-1]; 
				break;
             case 1:
                if (y[0]>0 && y[agn1]>0)
				{
                    // agn2 is the possition in the y or g arrays
                    agn2 = agn1+1;
                    // Calculate y and g
                    y[agn1] = y[agn1]-1;
                    y[0] -= 1;
                    y[agn2] += 1;
                    F[i] = F[i-1] + 1;
                    *yn -= 1;
                    // Calculate changes
                    if (agn1>2)
                    {
                        // a3 = (*cbf)*y(4:n).*(i-3);
                        a3Update(agn1,aa,a,aRowSize,y,*cbf);
                        a3Update(agn2,aa,a,aRowSize,y,*cbf);     
                    }
                    else if (agn2>2)
                    {
                        // a3 = (*cbf)*y(4:n).*(i-3);
                        a3Update(agn2,aa,a,aRowSize,y,*cbf);
                    }
                    // reaction propensities that need update
                    // a1 = cn*y[0]*(y[0]-1)/2;     Updated at the end
                    // a2 = clpf*y[0]*y(2:(n-1));    Updated at the end 
                    // a4 = 2*clmf*y(2:n);
                    a4Update(agn1,aa,a,aRowSize,y,*clmf);
                    a4Update(agn2,aa,a,aRowSize,y,*clmf);
                    // a5 = cap*g0*y;    Updated at the end
                    // a6 = clpm*y[0]*g(1:(n-1));      Updated at the end
                    // a10 = cnm*g0*y[0]*(y[0]-1)/2      Updated at the end
				}
                else
                {
                    F[i] = F[i-1];
                }
                M[i] = M[i-1];  
                break;
            case 2:
                if (y[agn1]>0)
                {
                    // agn2,3 are the possition in the y or g arrays
                    agn2 = RandInt(agn1-2);         // The aggregate can break at any bond, randomly
                    agn3 = agn1-agn2-1;
                    // Calculate y ang g 
                    y[agn1] -= 1;               
                    y[agn2] += 1;                   // y(ii)-> y(jj)+y(11-jj)
                    y[agn3] += 1;
                    *yn += 1;
                    // Calculate changes
                    if (agn2!=agn3)
                    {
                        // a2 = clpf*y[0]*y(2:(n-1));   Updated at the end 
                        // a4 = 2*clmf*y(2:n);
                        a4Update(agn3,aa,a,aRowSize,y,*clmf);
                        // a5 = cap*g0*y;               Updated at the end
                        if (agn3>2)
                        {
                            // a3 = (*cbf)*y(4:n).*yf;
                            a3Update(agn3,aa,a,aRowSize,y,*cbf);
                        }
                    } 
                    if (agn2>2)
                    {
                        // a3 = (*cbf)*y(4:n).*yf;
                        a3Update(agn2,aa,a,aRowSize,y,*cbf);
                    }          
                    // reaction propensities that need update
                    // a2 = clpf*y[0]*y(2:(n-1));       Updated at the end 
                    // a3 = (*cbf)*y(4:n).*yf; 
                    a3Update(agn1,aa,a,aRowSize,y,*cbf);
                    // a4 = 2*clmf*y(2:n);
                    a4Update(agn1,aa,a,aRowSize,y,*clmf);
                    a4Update(agn2,aa,a,aRowSize,y,*clmf);
                    // a5 = cap*g0*y;                   Updated at the end
                }
                // F[i] does not change
                F[i] = F[i-1];
                M[i] = M[i-1];
                break;
            case 3:
                if (y[agn1]>0)
                {
                    // agn2 is the possition in the y or g arrays
                    agn2 = agn1 - 1;
                    // Calculate y and g
                    y[agn1] -= 1;   
                    y[0] += 1;
                    y[agn2] += 1; 
                    if (agn1>1)
                    {
                        F[i] = F[i-1]- 1;
                    }
                    else
                    {
                        F[i] = F[i-1]- 2;
                    }
                    *yn += 1;
                    // Calculate changes
                    if (agn1>2)
                    {
                        // a3 = (*cbf)*y(4:n).*yf;
                        a3Update(agn1,aa,a,aRowSize,y,*cbf);
                    }
                    if (agn2>2)
                    {
                        // a3 = (*cbf)*y(4:n).*yf;
                        a3Update(agn2,aa,a,aRowSize,y,*cbf);
                    }
                    if (agn2>0)
                    {
                        // a4 = 2*clmf*y(2:n);
                        a4Update(agn2,aa,a,aRowSize,y,*clmf);
                        // a5 = cap*g0*y;    Updated at the end
                    }
                    // reaction propensities that need update 
                    // a4 = 2*clmf*y(2:n);
                    a4Update(agn1,aa,a,aRowSize,y,*clmf);
                    // a1 = cn*y[0]*(y[0]-1)/2;   Updated at the end
                    // a2 = clpf*y[0]*y(2:(n-1));    Updated at the end 
                    // a5 = cap*g0*y;    Updated at the end
                    // a6 = clpm*y[0]*g(1:(n-1));      Updated at the end
                    // a10 = cnm*g0*y[0]*(y[0]-1)/2      Updated at the end
                }
                else
                {
                    F[i] = F[i-1];
                }
                M[i] = M[i-1];
                break;
            case 4:
                if (y[agn1]>0 && g0>0)
                {
                    //Adjust y and g
                    y[agn1] -= 1;
                    g0 -= 1;
                    g[agn1] += 1;              
                    F[i] = F[i-1] - agn1;  // agn1 should be always >1
                    M[i] = M[i-1] + agn1; 
                    *yn -= 1;
                    // Calculate Changes         
                    if (agn1>2)
                    {
                        // a3 = (*cbf)*y(4:n).*yf;
                        a3Update(agn1,aa,a,aRowSize,y,*cbf);
                    }
                    if (agn1>1)
                    {
                        // a7 = cbm*g(3:n).*(i-2)
                        a7Update(agn1,aa,a,aRowSize,g,*cbm);
                    }
                    // reaction propensities that need update
                    // a2 = clpf*y[0]*y(2:(n-1));    Updated at the end 
                    // a4 = 2*clmf*y(2:n);
                    a4Update(agn1,aa,a,aRowSize,y,*clmf);
                    // a5 = cap*g0*y;    Updated at the end
                    // a6 = clpm*y[0]*g(1:(n-1));      Updated at the end
                    // a8 = clmm*g(2:n);
                    a8Update(agn1,aa,a,aRowSize,g,*clmm);
                    // a9 = cam*g(1:n);
                    a9Update(agn1,aa,a,aRowSize,g,*cam);
                    // a10 = cnm*g0*y[0]*(y[0]-1)/2      Updated at the end
                }
                else
                {
                    F[i] = F[i-1];  
                    M[i] = M[i-1];
                }  
                break;
            case 5:
                if (g[agn1]>0 && y[0]>0)
                {
                    // agn2 is the possition in the y or g arrays
                    agn2 = agn1 + 1;
                    // Calculate y and g
                    g[agn1] -= 1;
                    y[0] -= 1;
                    g[agn2] += 1;
                    M[i] = M[i-1]+ 1;   
                    *yn -= 1;
                    // Calculate changes
                    if (agn1>1)
                    {
                        // a7 = cbm*g(3:n).*(i-2)
                        a7Update(agn1,aa,a,aRowSize,g,*cbm);
                        a7Update(agn2,aa,a,aRowSize,g,*cbm);
                    }
                    else if (agn2>1)
                    {
                        // a7 = cbm*g(3:n).*(i-2)
                        a7Update(agn2,aa,a,aRowSize,g,*cbm);
                    }
                    if (agn1>0)
                    {
                        // a8 = clmm*g(2:n);
                        a8Update(agn1,aa,a,aRowSize,g,*clmm);
                    }
                    // reaction propensities that need update
                    // a1 = cn*y[0]*(y[0]-1)/2;     Updated at the end
                    // a2 = clpf*y[0]*y(2:(n-1));    Updated at the end 
                    // a6 = clpm*y[0]*g(1:(n-1));      Updated at the end
                    // a8 = clmm*g(2:n);  agn2>1
                    a8Update(agn2,aa,a,aRowSize,g,*clmm);
                    // a9 = cam*g(1:n);
                    a9Update(agn1,aa,a,aRowSize,g,*cam);
                    a9Update(agn2,aa,a,aRowSize,g,*cam);
                    // a10 = cnm*g0*y[0]*(y[0]-1)/2      Updated at the end
                }
                else
                {
                    M[i] = M[i-1];
                }
                F[i] = F[i-1];
                break;
            case 6:
                if (g[agn1]>0)
                {
                    // agn2,3 are the possition in the y or g arrays
                    agn2 = RandInt(agn1-1)-1;       // The aggregate can break at any bond, randomly (g length)
                    agn3 = agn1-agn2-1;             // ang3 relate to the y length
                    // Calculate y and g
                    g[agn1] = g[agn1]- 1; 
                    g[agn2] = g[agn2]+ 1;
                    y[agn3] = y[agn3]+ 1;  
                    F[i] = F[i-1]+agn3;
                    M[i] = M[i-1]-agn3; 
                    *yn += 1;
                    // Calculate changes
                    if (agn3>2);
                    {
                        // a3 = (*cbf)*y(4:n).*yf;
                        a3Update(agn3,aa,a,aRowSize,y,*cbf);
                    }
                    if (agn3>0)
                    {
                        // a4 = 2*clmf*y(2:n);
                        a4Update(agn3,aa,a,aRowSize,y,*clmf);
                    }
                    if (agn2>0)
                    {
                        // a8 = clmm*g(2:n);
                        a8Update(agn2,aa,a,aRowSize,g,*clmm);
                    }
                    if (agn2>1)
                    {
                        // a7 = cbm*g(3:n).*(i-2)
                        a7Update(agn2,aa,a,aRowSize,g,*cbm);
                    }  
                    // reaction propensities that need update
                    // a2 = clpf*y[0]*y(2:(n-1));    Updated at the end  
                    // a5 = cap*g0*y;    Updated at the end
                    // a6 = clpm*y[0]*g(1:(n-1));      Updated at the end
                    // a7 = cbm*g(3:n).*(i-2)
                    a7Update(agn1,aa,a,aRowSize,g,*cbm);
                    // a8 = clmm*g(2:n);
                    a8Update(agn1,aa,a,aRowSize,g,*clmm);
                    // a9 = cam*g(1:n);
                    a9Update(agn1,aa,a,aRowSize,g,*cam);
                    a9Update(agn2,aa,a,aRowSize,g,*cam);
                }
                else
                {
                    F[i] = F[i-1];
                    M[i] = M[i-1];
                }
                break;
            case 7:
                if (g[agn1]>0)
                {
                    // agn2 is the possition in the y or g arrays
                    agn2 = agn1 - 1; 
                    // Calculate y and g
                    g[agn1] -= 1;   
                    g[agn2] += 1; 
                    y[0] += 1;
                    M[i] = M[i-1]- 1;   
                    *yn += 1;
                    // Calculate changes
                    if (agn1>1)
                    {
                        // a7 = cbm*g(3:n).*(i-2)
                        a7Update(agn1,aa,a,aRowSize,g,*cbm);
                    }
                    if (agn2>1)
                    {
                        // a7 = cbm*g(3:n).*(i-2)
                        a7Update(agn2,aa,a,aRowSize,g,*cbm);
                    }
                    if (agn2>0)
                    {
                        // a8 = clmm*g(2:n);
                        a8Update(agn2,aa,a,aRowSize,g,*clmm);
                    }
                    // reaction propensities that need update
                    // a1 = cn*y[0]*(y[0]-1)/2;     Updated at the end
                    // a2 = clpf*y[0]*y(2:(n-1));    Updated at the end    
                    // a6 = clpm*y[0]*g(1:(n-1));      Updated at the end  
                    // a8 = clmm*g(2:n);
                    a8Update(agn1,aa,a,aRowSize,g,*clmm);
                    // a9 = cam*g(1:n);
                    a9Update(agn1,aa,a,aRowSize,g,*cam);
                    a9Update(agn2,aa,a,aRowSize,g,*cam);
                    // a10 = cnm*g0*y[0]*(y[0]-1)/2      Updated at the end
                }
                else
                {
                    M[i] = M[i-1];
                }
                F[i] = F[i-1];
                break;
            case 8:
                if  (g[agn1]>0)
                {
                    // Calculate y and g
                    y[agn1] += 1;
                    g0 += 1;
                    g[agn1] -= 1;
                    if (agn1>0)
                    {
                        F[i] = F[i-1]+ agn1;
                    }
                    else
                    {
                        F[i] = F[i-1];
                    }
                    M[i] = M[i-1]- agn1;
                    *yn += 1;
                    // Calculate changes
                    if (agn1>0)
                    {
                        // a2 = clpf*y[0]*y(2:(n-1));    Updated at the end 
                        // a4 = 2*clmf*y(2:n);
                        a4Update(agn1,aa,a,aRowSize,y,*clmf);
                        // a6 = clpm*y[0]*g(1:(n-1));      Updated at the end
                        // a8 = clmm*g(2:n);
                        a8Update(agn1,aa,a,aRowSize,g,*clmm);
                    }   
                    if (agn1>2)
                    {
                        // a3 = (*cbf)*y(4:n).*yf;
                        a3Update(agn1,aa,a,aRowSize,y,*cbf);
                    } 
                    if (agn1>1)
                    {
                        // a7 = cbm*g(3:n).*(i-2)
                        a7Update(agn1,aa,a,aRowSize,g,*cbm);
                    }         
                    // reaction propensities that need update
                    // a1 = cn*y[0]*(y[0]-1)/2;     Updated at the end
                    // a5 = cap*g0*y;    Updated at the end
                    // a9 = cam*g(1:n);
                    a9Update(agn1,aa,a,aRowSize,g,*cam);
                    // a10 = cnm*g0*y[0]*(y[0]-1)/2      Updated at the end
                }
                else
                {
                    F[i] = F[i-1];
                    M[i] = M[i-1];
                }
                break;
            case 9:
                if (y[0]>0 && g0>0)
                {
                    //Adjust y and g
                    y[0] -= 2;
                    g0 -= 1;
                    g[1] += 1;     
                    M[i] = M[i-1] + 2;   
                    *yn -= 2;
                    // Calculate changes
                    // reaction propensities that need update
                    // a1 = cn*y[0]*(y[0]-1)/2;      Updated at the end
                    // a2 = clpf*y[0]*y(2:(n-1));    Updated at the end 
                    // a5 = cap*g0*y;    Updated at the end
                    // a6 = clpm*y[0]*g(1:(n-1));      Updated at the end
                    // a9 = cam*g(1:n);
                    a9Update(1,aa,a,aRowSize,g,*cam);
                }
                else
                {
                    M[i] = M[i-1];
                }
                F[i] = F[i-1];      // monomers are not counted in F
                break;
            default:
                // mexErrMsgIdAndTxt("AggKineticLoop","Error - Number of reaction is out of range");
                mexPrintf("AggKineticLoop - Error - rn>9.");
                break;
        }

        // Out-flow calculation - update all propensities that effected by change in y(yi)   
        // a3 = (*cbf)*y(4:n).*yf;      (4-n)
        // a4 = 2*clmf*y(2:n);          (2-n)
        if (yi>2)
        {
            // Update a3,a4
            // a3 = (*cbf)*y(4:n).*(i-3);
            a3Update(yi,aa,a,aRowSize,y,*cbf);
            // a4 = 2*clmf*y(2:n);
            a4Update(yi,aa,a,aRowSize,y,*clmf); 
        }
        else if (yi>0)
        {
            // Update a4
            // a4 = 2*clmf*y(2:n);
            a4Update(yi,aa,a,aRowSize,y,*clmf);
        }

        // for g:
        // a7 = cbm*g(3:n).*(i-2)    (3-n)
        // a8 = clmm*g(2:n);         (2-n)
        // a9 = cam*g(1:n);          (1-n)
        if (gi>1)
        {
            // Update a7,a8,a9
            // a7 = cbm*g(3:n).*(i-2)
            a7Update(gi,aa,a,aRowSize,g,*cbm);
            // a8 = clmm*g(2:n);
            a8Update(gi,aa,a,aRowSize,g,*clmm);
            // a9 = cam*g(1:n);
            a9Update(gi,aa,a,aRowSize,g,*cam);
        }
        else if (gi==1)
        {
            // Update a8,a9
            // a8 = clmm*g(2:n);
            a8Update(gi,aa,a,aRowSize,g,*clmm);
            // a9 = cam*g(1:n);
            a9Update(gi,aa,a,aRowSize,g,*cam); 
        }
        
        // Always update 
        y0t[i] = y[0];
        g0t[i] = g0;
        // reaction propensities that need update
        // a1 = cn*y[0]*(y[0]-1)/2;  
        ii = FindPos(0,0,aRowSize);
        a[ii] = (*cn)*y[0]*(y[0]-1)/2.0;
        aa[0] = a[ii];
        // a2 = clpf*y[0]*y(2:(n-1));
        iiStart = FindPos(1,rst[1],aRowSize);
        iiEnd = FindPos(1,re[1],aRowSize);
        ConstFactor = (*clpf)*y[0];
        for (ii=iiStart; ii<=iiEnd; ii += aRowSize)
        {
			jj = ii-iiStart+1;
            a[ii] = ConstFactor*y[jj];
        }
        aa[1] = SumArrayRow(a,1,0,aRowSize,aColSize);
        // a5 = cap*g0*y;
        iiStart = FindPos(4,rst[4],aRowSize);
        iiEnd = FindPos(4,re[4],aRowSize);
        ConstFactor = (*cap)*g0;
        for (ii=iiStart; ii<=iiEnd; ii += aRowSize)
        {
			jj = ii-iiStart+1;
			a[ii] = ConstFactor*y[jj];    
        }
        aa[4] = SumArrayRow(a,4,0,aRowSize,aColSize);
        // a6 = clpm*y[0]*g(1:(n-1));
        iiStart = FindPos(5,rst[5],aRowSize);
        iiEnd = FindPos(5,re[5],aRowSize);
        ConstFactor = (*clpm)*y[0];
        for (ii=iiStart; ii<=iiEnd; ii += aRowSize)
        {
			jj = ii-iiStart;
            a[ii] = ConstFactor*g[jj];
        }
        aa[5] = SumArrayRow(a,5,0,aRowSize,aColSize);
        // a10 = cnm*g0*y[0]*(y[0]-1)/2
        ii = FindPos(9,0,aRowSize);
        a[ii] = (*cnm)*g0*y[0]*(y[0]-1)/2.0;
        aa[9] = a[ii]; 
        // Sum propensities
        *a0 = SumArray(a0Data,a0);
    } 
}