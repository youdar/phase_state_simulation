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


#define N_ELEMENTS 100000

/************ The main mex function **************/
void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    
    /* Input variables 
     * (clpf,clmf,clpm,clmm,cap,cam,cn,cnm,cbf,cbm,y0,p0,rp,crs,crp,rfp,crfm,crfp)   */

    double *t = NULL;
    double *g = NULL;
    
    
 //   double t[100000],F[100000],M[100000];
//     double y0t[10000],g0t[10000];
//     double inagg[10000];
    double y[100000],g[100000];
    double a0,aa[10],a[10][500],dfy[500],dfg[500];
    int rst[10] = {0, 1, 3, 1, 1, 0, 2, 1, 0, 0};
    int re[10] = {0, 498, 499, 499, 499, 498, 499, 499, 499, 0};
    int n = 500,m = 100000 ;
    
    
    double clpf,clmf,clpm,clmm,cap,cam,cn,cnm,cbf,cbm,rp;
    double crs,crp,rfp,crfm,crfp;
    int y0,p0,g0,rn,yn,gn;
    double dfg0 = 0.0, dg0in = 0.0, dyin = 0.0;
    double dt,r,s,rg,sg,ry,sy,aold;
    
    int i,dy1,dy2,yi,gi,dg01,dg2;
    int agn1,agn2,agn3,ii,jj;

    /* Output variables
     * When changing arrays in the MEX code, the array in Matlab is changes as well
     * so there is no need to send is as an ouput
     */
    
    t = (double *)malloc(N_ELEMENTS * sizeof(double));
    if (t = NULL)
    {
        // uh oh. failed to alloc memory.
        goto exitpoint;
        
    }
    g = (double *)malloc(N_ELEMENTS * sizeof(double));
    if (g = NULL)
    {
        // uh oh. failed to alloc memory.
        goto exitpoint;
        
    }
    
    
    /* setting the input values */
    clpf = *mxGetPr(prhs[0]);
    clmf = *mxGetPr(prhs[1]);
    clpm = *mxGetPr(prhs[2]);
    clmm = *mxGetPr(prhs[3]);
    cap = *mxGetPr(prhs[4]);
    cam = *mxGetPr(prhs[5]);
    cn = *mxGetPr(prhs[6]);
    cnm = *mxGetPr(prhs[7]);
    cbf = *mxGetPr(prhs[8]);
    cbm = *mxGetPr(prhs[9]);
    y0 = *mxGetPr(prhs[10]);
    y[1] = y0;
    yn = y0;
    p0 = *mxGetPr(prhs[11]);
    g0 = p0;  
    rp = *mxGetPr(prhs[12]);
    crs = *mxGetPr(prhs[13]);
    crp = *mxGetPr(prhs[14]);
    rfp = *mxGetPr(prhs[15]);
    crfm = *mxGetPr(prhs[16]);
    crfp = *mxGetPr(prhs[17]);

    
    // None zero initial values
    a[0][0] = cn*y[0]*(y[0]-1)/2.0;
    a[9][0] = cnm*g0*y[0]*(y[0]-1)/2.0;



    // Sum all the none zero values
    aa[0] = a[0][0];
    aa[9] = a[9][0];
    a0 = 0;
    for (ii=0; ii<10; ++ii)
    {
        a0 += aa[ii];
    }
    
    /* Seed random numbers */
    srand(time(NULL));
    
    

    /*  Start main loop */
    /* Remember that in C array starts at 0 not at 1, like in Matlab*/
    for (i = 1; i < m; ++i)
    {
        // Add time step
        dt = log(1/RandUn())/a0;
        t[i] = t[i-1]+ dt;   

        // Get the reaction type number                           
        r = RandUn();           // Random number (0,1]
        rn = 0;             // the reaction number in a (0 to 9)
        s = aa[rn]/a0;      // accomulating the propensities
        while (s < r )
        {
            rn += 1;
            if (rn > 9)
            {
                 // mexErrMsgIdAndTxt("AggKineticLoop","Error - rn>10");
                mexPrintf("AggKineticLoop - Error - rn>10.");
            }
            s += aa[rn]/a0;  
        }

        // find the reaction within the reaction type
        // agn1 is the possition in the y or g arrays
        agn1 = rst[rn];                             // the first j in the correct region
        if (rn!=0 && rn!=9)
        {
            s += a[rn][agn1]/a0 - aa[rn]/a0;  // return value to begining of reaction
            while (s < r)
            {
                agn1 += 1;
                s += a[rn][agn1]/a0; 
            }
        }
    
    
        // The deterministic in-flow
        // Adjust y[0] and g0 
        dyin += dt*crfp;        // add accumulated fraction of monomers
        dy1 = dyin;             // moving only whole particle
        dyin -= dy1;            // keeping track of fraction of particles

        // Out-flow calculation
        ry = yn*RandUn();       // using random selection according the aggregation number distribution
        sy = y[0];
        yi = 0;                 // aggregate to be flowing out
        while (sy<ry)
        {
            yi += 1;
            sy += y[yi];
        }

        // Out-flow calculation -  update y[yi]
        dfy[yi] += dt*y[yi]*crfm;
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
        yn += dy1 - dy2;                // Update number of free particles (aggregates and monomers)
        if (y[yi]<0)
        {
            yn -= y[yi];
            y[yi] = 0;   
        }
            
    
        // Membrane protein production
        dg0in += dt*crp;     // add accumulated fraction of membrane proteins
        dg01 = dg0in;           // moving only whole particle
        dg0in -= dg01;          // keeping track of fraction of particles
        // Membrane sequestration
        rg = gn*RandUn();       // using random selection according the aggregation number distribution
        sg = g0;                // initially
        gi = 0;                 // aggregate to be flowing out
        while (sg<rg)
        {
            sg += g[gi];
            gi += 1;
        }
        
        
        // update g0 for the in-flow
        // update g[gi]
        if (gi>0)
        {
            gi -= 1;                        // the arrays index starts at 0, not like matlab
            dfg[gi] += dt*g[gi]*crs;
            dg2 = dfg[gi];                  // the integer number of aggregates or monomers the flow out
            g[gi] -= dg2;                   // update g
            g0 += dg01;                 
            dfg[gi] -= dg2;                 // update the reminder of the aggregates or monomers, that did not flow out
            if (g[gi]<0)
            {
                gn += dg01 - dg2 - g[gi]; 
//                 inagg[i] = inagg[i-1] + (gi+1)*(dg2+g[gi]);  
                g[gi] = 0;
            }
            else
            {
                gn += dg01 - dg2; 
//                 inagg[i] = inagg[i-1] + (gi+1)*dg2;  
            } 
        }
        else
        {
            dfg0 += dt*g0*crs;
            dg2 = dfg0;                 // the integer number of aggregates or monomers the flow out
            g0 += dg01 - dg2;           // update g0
            dfg0 -= dg2;                // update the reminder of the aggregates or monomers, that did not flow out 
//             inagg[i] = inagg[i-1];
            if (g0<0)
            {
                gn += dg01 - dg2 - g0; 
                g0 = 0;
            }
            else
            {
                gn += dg01 - dg2;
            }
        }
                    
//         if (gn<0)
//         {
// 			// mexErrMsgIdAndTxt("AggKineticLoop","Error - gn < 0");
//             mexPrintf("AggKineticLoop - Error - gn<0\n");
//         }
		
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
					yn -= 1;
					y[1] += 1;
					F[i] = F[i-1]+ 2;
					// reaction propensities that need update
					// a1 = cn*y[0]*(y[0]-1)/2;     Updated at the end
					// a2 = clpf*y[0]*y(2:(n-1));   Updated at the end 
					// a4 = 2*clmf*y(2:n);
                    aold  = a[3][1];
                    a[3][1] = 2*clmf*y[1];
                    aa[3] += a[3][1] - aold;
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
                    yn -= 1;
                    // Calculate changes
                    if (agn1>2)
                    {
                        // a3 = cbf*y(4:n).*(i-3);
                        aold = a[2][agn1] + a[2][agn2];
                        a[2][agn1] = cbf*y[agn1]*(agn1-2);
                        a[2][agn2] = cbf*y[agn2]*(agn2-2);
                        aa[2] += a[2][agn1]+a[2][agn2] - aold;
                    }
                    else if (agn2>2)
                    {
                        // a3 = cbf*y(4:n).*(i-3);
                        aold = a[2][agn2];
                        a[2][agn2] = cbf*y[agn2]*(agn2-2);
                        aa[2] += a[2][agn2] - aold;
                    }
                    // reaction propensities that need update
                    // a1 = cn*y[0]*(y[0]-1)/2;     Updated at the end
                    // a2 = clpf*y[0]*y(2:(n-1));    Updated at the end 
                    // a4 = 2*clmf*y(2:n);
                    aold  = a[3][agn1] + a[3][agn2];
                    a[3][agn1] = 2*clmf*y[agn1];
                    a[3][agn2] = 2*clmf*y[agn2];
                    aa[3] += a[3][agn1] + a[3][agn2] - aold;
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
                    yn += 1;
                    // Calculate changes
                    if (agn2!=agn3)
                    {
                        // a2 = clpf*y[0]*y(2:(n-1));   Updated at the end 
                        // a4 = 2*clmf*y(2:n);
                        aold  = a[3][agn3];
                        a[3][agn3] = 2*clmf*y[agn3];
                        aa[3] += a[3][agn3] - aold;
                        // a5 = cap*g0*y;               Updated at the end
                        if (agn3>2)
                        {
                            // a3 = cbf*y(4:n).*yf;
                            aold = a[2][agn3] + a[2][agn2];
                            a[2][agn3] = cbf*y[agn3]*(agn3-2);
                            aa[2] += a[2][agn3] - aold;
                        }
                    } 
                    if (agn2>2)
                    {
                        // a3 = cbf*y(4:n).*yf;
                        aold = a[2][agn2];
                        a[2][agn2] = cbf*y[agn2]*(agn2-2);
                        aa[2] += a[2][agn2] - aold;
                    }          
                    // reaction propensities that need update
                    // a2 = clpf*y[0]*y(2:(n-1));       Updated at the end 
                    // a3 = cbf*y(4:n).*yf; 
                    aold = a[2][agn1];
                    a[2][agn1] = cbf*y[agn1]*(agn1-2);
                    aa[2] += a[2][agn1] - aold;
                    // a4 = 2*clmf*y(2:n);
                    aold  = a[3][agn1] + a[3][agn2];
                    a[3][agn1] = 2*clmf*y[agn1];
                    a[3][agn2] = 2*clmf*y[agn2];
                    aa[3] += a[3][agn1] + a[3][agn2] - aold;
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
                    yn += 1;
                    // Calculate changes
                    if (agn1>2)
                    {
                        // a3 = cbf*y(4:n).*yf;
                        aold = a[2][agn1];
                        a[2][agn1] = cbf*y[agn1]*(agn1-2);
                        aa[2] += a[2][agn1] - aold;
                    }
                    if (agn2>2)
                    {
                        // a3 = cbf*y(4:n).*yf;
                        aold = a[2][agn2];
                        a[2][agn2] = cbf*y[agn2]*(agn2-2);
                        aa[2] += a[2][agn2] - aold;
                    }
                    if (agn2>0)
                    {
                        // a4 = 2*clmf*y(2:n);
                        aold  = a[3][agn2];
                        a[3][agn2] = 2*clmf*y[agn2];
                        aa[3] += a[3][agn2] - aold;
                        // a5 = cap*g0*y;    Updated at the end
                    }
                    // reaction propensities that need update 
                    // a4 = 2*clmf*y(2:n);
                    aold  = a[3][agn1];
                    a[3][agn1] = 2*clmf*y[agn1];
                    aa[3] += a[3][agn1] - aold;
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
                    yn -= 1;
                    // Calculate Changes         
                    if (agn1>2)
                    {
                        // a3 = cbf*y(4:n).*yf;
                        aold = a[2][agn1];
                        a[2][agn1] = cbf*y[agn1]*(agn1-2);
                        aa[2] += a[2][agn1] - aold;
                    }
                    if (agn1>1)
                    {
                        // a7 = cbm*g(3:n).*(i-2)
                        aold =  a[6][agn1];
                        a[6][agn1]= cbm*g[agn1]*(agn1-1);
                        aa[6] += a[6][agn1] - aold;
                    }
                    // reaction propensities that need update
                    // a2 = clpf*y[0]*y(2:(n-1));    Updated at the end 
                    // a4 = 2*clmf*y(2:n);
                    aold  = a[3][agn1];
                    a[3][agn1] = 2*clmf*y[agn1];
                    aa[3] += a[3][agn1] - aold;
                    // a5 = cap*g0*y;    Updated at the end
                    // a6 = clpm*y[0]*g(1:(n-1));      Updated at the end
                    // a8 = clmm*g(2:n);
                    aold = a[7][agn1];
                    a[7][agn1]= clmm*g[agn1];
                    aa[7] += a[7][agn1] - aold;
                    // a9 = cam*g(1:n);
                    aold = a[8][agn1];
                    a[8][agn1]= cam*g[agn1];
                    aa[8] += a[8][agn1] - aold;
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
                    yn -= 1;
                    // Calculate changes
                    if (agn1>1)
                    {
                        // a7 = cbm*g(3:n).*(i-2)
                        aold =  a[6][agn1] + a[6][agn2];
                        a[6][agn1]= cbm*g[agn1]*(agn1-1);
                        a[6][agn2]= cbm*g[agn2]*(agn2-1);
                        aa[6] += a[6][agn1] + a[6][agn2] - aold;
                    }
                    else if (agn2>1)
                    {
                        // a7 = cbm*g(3:n).*(i-2)
                        aold =  a[6][agn2];
                        a[6][agn2]= cbm*g[agn2]*(agn2-1);
                        aa[6] += a[6][agn2] - aold;
                    }
                    if (agn1>0)
                    {
                        // a8 = clmm*g(2:n);
                        aold = a[7][agn1];
                        a[7][agn1]= clmm*g[agn1];
                        aa[7] += a[7][agn1] - aold;
                    }
                    // reaction propensities that need update
                    // a1 = cn*y[0]*(y[0]-1)/2;     Updated at the end
                    // a2 = clpf*y[0]*y(2:(n-1));    Updated at the end 
                    // a6 = clpm*y[0]*g(1:(n-1));      Updated at the end
                    // a8 = clmm*g(2:n);  agn2>1
                    aold = a[7][agn2];
                    a[7][agn2]= clmm*g[agn2];
                    aa[7] += a[7][agn2] - aold;
                    // a9 = cam*g(1:n);
                    aold = a[8][agn1] + a[8][agn2];
                    a[8][agn1]= cam*g[agn1];
                    a[8][agn2]= cam*g[agn2];
                    aa[8] += a[8][agn1] + a[8][agn2] - aold;
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
                    yn += 1;
                    // Calculate changes
                    if (agn3>2);
                    {
                        // a3 = cbf*y(4:n).*yf;
                        aold = a[2][agn3] + a[2][agn2];
                        a[2][agn3] = cbf*y[agn3]*(agn3-2);
                        aa[2] += a[2][agn3] - aold;
                    }
                    if (agn3>0)
                    {
                        // a4 = 2*clmf*y(2:n);
                        aold  = a[3][agn3];
                        a[3][agn3] = 2*clmf*y[agn3];
                        aa[3] += a[3][agn3] - aold;
                    }
                    if (agn2>0)
                    {
                        // a8 = clmm*g(2:n);
                        aold = a[7][agn2];
                        a[7][agn2]= clmm*g[agn2];
                        aa[7] += a[7][agn2] - aold;
                    }
                    if (agn2>1)
                    {
                        // a7 = cbm*g(3:n).*(i-2)
                        aold =  a[6][agn2];
                        a[6][agn2]= cbm*g[agn2]*(agn2-1);
                        aa[6] += a[6][agn2] - aold;
                    }  
                    // reaction propensities that need update
                    // a2 = clpf*y[0]*y(2:(n-1));    Updated at the end  
                    // a5 = cap*g0*y;    Updated at the end
                    // a6 = clpm*y[0]*g(1:(n-1));      Updated at the end
                    // a7 = cbm*g(3:n).*(i-2)
                    aold =  a[6][agn1];
                    a[6][agn1]= cbm*g[agn1]*(agn1-1);
                    aa[6] += a[6][agn1] - aold;
                    // a8 = clmm*g(2:n);
                    aold = a[7][agn1];
                    a[7][agn1]= clmm*g[agn1];
                    aa[7] += a[7][agn1] - aold;
                    // a9 = cam*g(1:n);
                    aold = a[8][agn1] + a[8][agn2];
                    a[8][agn1]= cam*g[agn1];
                    a[8][agn2]= cam*g[agn2];
                    aa[8] += a[8][agn1] + a[8][agn2] - aold;
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
                    yn += 1;
                    // Calculate changes
                    if (agn1>1)
                    {
                        // a7 = cbm*g(3:n).*(i-2)
                        aold =  a[6][agn1];
                        a[6][agn1]= cbm*g[agn1]*(agn1-1);
                        aa[6] += a[6][agn1] - aold;
                    }
                    if (agn2>1)
                    {
                        // a7 = cbm*g(3:n).*(i-2)
                        aold =  a[6][agn2];
                        a[6][agn2]= cbm*g[agn2]*(agn2-1);
                        aa[6] += a[6][agn2] - aold;
                    }
                    if (agn2>0)
                    {
                        // a8 = clmm*g(2:n);
                        aold = a[7][agn2];
                        a[7][agn2]= clmm*g[agn2];
                        aa[7] += a[7][agn2] - aold;
                    }
                    // reaction propensities that need update
                    // a1 = cn*y[0]*(y[0]-1)/2;     Updated at the end
                    // a2 = clpf*y[0]*y(2:(n-1));    Updated at the end    
                    // a6 = clpm*y[0]*g(1:(n-1));      Updated at the end  
                    // a8 = clmm*g(2:n);
                    aold = a[7][agn1];
                    a[7][agn1]= clmm*g[agn1];
                    aa[7] += a[7][agn1] - aold;
                    // a9 = cam*g(1:n);
                    aold = a[8][agn1] + a[8][agn2];
                    a[8][agn1]= cam*g[agn1];
                    a[8][agn2]= cam*g[agn2];
                    aa[8] += a[8][agn1] + a[8][agn2] - aold;
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
                    yn += 1;
                    // Calculate changes
                    if (agn1>0)
                    {
                        // a2 = clpf*y[0]*y(2:(n-1));    Updated at the end 
                        // a4 = 2*clmf*y(2:n);
                        aold  = a[3][agn1];
                        a[3][agn1] = 2*clmf*y[agn1];
                        aa[3] += a[3][agn1] - aold;
                        // a6 = clpm*y[0]*g(1:(n-1));      Updated at the end
                        // a8 = clmm*g(2:n);
                        aold = a[7][agn1];
                        a[7][agn1]= clmm*g[agn1];
                        aa[7] += a[7][agn1] - aold;
                    }   
                    if (agn1>2)
                    {
                        // a3 = cbf*y(4:n).*yf;
                        aold = a[2][agn1];
                        a[2][agn1] = cbf*y[agn1]*(agn1-2);
                        aa[2] += a[2][agn1] - aold;
                    } 
                    if (agn1>1)
                    {
                        // a7 = cbm*g(3:n).*(i-2)
                        aold =  a[6][agn1];
                        a[6][agn1]= cbm*g[agn1]*(agn1-1);
                        aa[6] += a[6][agn1] - aold;
                    }         
                    // reaction propensities that need update
                    // a1 = cn*y[0]*(y[0]-1)/2;     Updated at the end
                    // a5 = cap*g0*y;    Updated at the end
                    // a9 = cam*g(1:n);
                    aold = a[8][agn1];
                    a[8][agn1]= cam*g[agn1];
                    aa[8] += a[8][agn1] - aold;
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
                    yn -= 2;
                    // Calculate changes
                    // reaction propensities that need update
                    // a1 = cn*y[0]*(y[0]-1)/2;      Updated at the end
                    // a2 = clpf*y[0]*y(2:(n-1));    Updated at the end 
                    // a5 = cap*g0*y;    Updated at the end
                    // a6 = clpm*y[0]*g(1:(n-1));      Updated at the end
                    // a9 = cam*g(1:n);
                    aold = a[8][1];
                    a[8][1]= cam*g[1];
                    aa[8] += a[8][1]  - aold;
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
        // a3 = cbf*y(4:n).*yf;      (4-n)
        // a4 = 2*clmf*y(2:n);          (2-n)
        if (yi>2)
        {
            // Update a3,a4
            // a3 = cbf*y(4:n).*(i-3);
            aold = a[2][yi];
            a[2][yi] = cbf*y[yi]*(yi-2);
            aa[2] += a[2][yi] - aold;
            // a4 = 2*clmf*y(2:n); 
            aold  = a[3][yi];
            a[3][yi] = 2*clmf*y[yi];
            aa[3] += a[3][yi] - aold;
        }
        else if (yi>0)
        {
            // Update a4
            // a4 = 2*clmf*y(2:n);
            aold  = a[3][yi];
            a[3][yi] = 2*clmf*y[yi];
            aa[3] += a[3][yi] - aold;
        }

        // for g:
        // a7 = cbm*g(3:n).*(i-2)    (3-n)
        // a8 = clmm*g(2:n);         (2-n)
        // a9 = cam*g(1:n);          (1-n)
        if (gi>1)
        {
            // Update a7,a8,a9
            // a7 = cbm*g(3:n).*(i-2)
            aold =  a[6][gi];
            a[6][gi]= cbm*g[gi]*(gi-1);
            aa[6] += a[6][gi] - aold;
            // a8 = clmm*g(2:n);
            aold = a[7][gi];
            a[7][gi]= clmm*g[gi];
            aa[7] += a[7][gi] - aold;
            // a9 = cam*g(1:n);
            aold = a[8][gi];
            a[8][gi]= cam*g[gi];
            aa[8] += a[8][gi] - aold;
        }
        else if (gi==1)
        {
            // Update a8,a9
            // a8 = clmm*g(2:n);
            aold = a[7][gi];
            a[7][gi]= clmm*g[gi];
            aa[7] += a[7][gi] - aold;
            // a9 = cam*g(1:n);
            aold = a[8][gi];
            a[8][gi]= cam*g[gi];
            aa[8] += a[8][gi] - aold;
        }
        
        // Always update 
//         y0t[i] = y[0];
//         g0t[i] = g0;
        // reaction propensities that need update
        // a1 = cn*y[0]*(y[0]-1)/2;  
        a[0][0] = cn*y[0]*(y[0]-1)/2.0;
        aa[0] = a[0][0];
        // a2 = clpf*y[0]*y(2:(n-1));
        aa[1] = 0;
        for (ii=rst[1]; ii<=re[1]; ++ii)
        {
            a[1][ii] = clpf*y[0]*y[ii];
            aa[1] += a[1][ii];
        }
        // a5 = cap*g0*y;
        aa[4] = 0;
        for (ii=rst[4]; ii<=re[4]; ++ii)
        {
			a[4][ii] = cap*g0*y[ii]; 
            aa[4] += a[4][ii];
        }
        // a6 = clpm*y[0]*g(1:(n-1));
        for (ii=rst[5]; ii<=re[5]; ++ii)
        {
            a[5][ii] = clpm*y[0]*g[ii];
            aa[5] += a[5][ii];
        }
        // a10 = cnm*g0*y[0]*(y[0]-1)/2

        a[9][0] = cnm*g0*y[0]*(y[0]-1)/2.0;
        aa[9] = a[9][0]; 
        // Sum propensities
        a0 = 0;
        for (ii=0; ii<10; ++ii)
        {
            a0 += aa[ii];
        }
    } 
    
    exitpoint:
    if (t != NULL) 
    {
        free(t);
    }
    if (g != NULL)
    {
        free(g);
    }
}
