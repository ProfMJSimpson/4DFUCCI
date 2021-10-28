/*
 * cells.c - Simulate the cell agent behaviours for the IBM
 *
 * Author: Jonah Klowss, 
 * School of Mathematical Sciences, 
 * Queensland University of Technology.
 *
 * Date: 21/10/2021
 *
 * This is a MEX file for MATLAB.
*/

#include "mex.h"
#include <math.h>
#include "matrix.h"
#include "stdio.h"

// Define constants and valuable functions
#define pi          3.14159265358979323846
#define min(a,b)    (((a) < (b)) ? (a) : (b))
#define max(a,b)    (((a) < (b)) ? (b) : (a))

// RNG with larger values
int rand2();
int rseed = 0;
inline void srand2(int x) {
    rseed = x;
}
#define RAND_MAX2 ((1U << 31) - 1)
inline int rand2() {
    return rseed = (rseed * 1103515245 + 12345) & RAND_MAX2;
}

// Uniformly distributed random number, between 0 and 1 
double unif() 
{
    
    double urdm = rand2() / ((double) RAND_MAX2);
    return urdm;
    
}

int ag_choice(double *vec, double vecS)
{
    
    int ind = 0; // Initialise index
    double agc = max(0,vec[ind]); 
    
    double u = unif() * vecS;
    
    while (u > agc) {
        ind += 1;
        agc += max(0,vec[ind]);
    }
    
    return ind;
    
}

/*  3D trilinear interpolator */
double interp3d(int I, double h, double Xp, double Yp, double Zp, double *c)
{
  // Find index location of agent
  int i = ceil(Xp/h);
  int j = ceil(Yp/h);
  int k = ceil(Zp/h);
  
  // Find surrounding x y z locations
  double x1 = h*(i-1);
  double x2 = h*i;
  double y1 = h*(j-1);
  double y2 = h*j;
  double z1 = h*(k-1);
  double z2 = h*k;
  
  // Find indices and values of surrounding nutrient concentrations
  int ind1 = pow(I,2)*(k-1) + I*(i-1) + j - 1;
  int ind2 = pow(I,2)*(k-1) + I*i + j - 1;
  int ind3 = pow(I,2)*(k-1) + I*(i-1) + j;
  int ind4 = pow(I,2)*(k-1) + I*i + j;
  int ind5 = pow(I,2)*k + I*(i-1) + j - 1;
  int ind6 = pow(I,2)*k + I*i + j - 1;
  int ind7 = pow(I,2)*k + I*(i-1) + j;
  int ind8 = pow(I,2)*k + I*i + j;
  
  double c1 = c[ind1];
  double c2 = c[ind2];
  double c3 = c[ind3];
  double c4 = c[ind4];
  double c5 = c[ind5];
  double c6 = c[ind6];
  double c7 = c[ind7];
  double c8 = c[ind8];
  
  // Calculate 3D interpolation coefficients
  double a0 = -(c1*x2*y2*z2 - c2*x1*y2*z2 - c3*x2*y1*z2 + c4*x1*y1*z2 - c5*x2*y2*z1 + c6*x1*y2*z1 + c7*x2*y1*z1 - c8*x1*y1*z1)/((x1 - x2)*(y1 - y2)*(z1 - z2));
  double a1 = (c1*y2*z2 - c2*y2*z2 - c3*y1*z2 + c4*y1*z2 - c5*y2*z1 + c6*y2*z1 + c7*y1*z1 - c8*y1*z1)/((x1 - x2)*(y1 - y2)*(z1 - z2));
  double a2 = (c1*x2*z2 - c2*x1*z2 - c3*x2*z2 + c4*x1*z2 - c5*x2*z1 + c6*x1*z1 + c7*x2*z1 - c8*x1*z1)/((x1 - x2)*(y1 - y2)*(z1 - z2));
  double a3 = (c1*x2*y2 - c2*x1*y2 - c3*x2*y1 + c4*x1*y1 - c5*x2*y2 + c6*x1*y2 + c7*x2*y1 - c8*x1*y1)/((x1 - x2)*(y1 - y2)*(z1 - z2));
  double a4 = -(c1*z2 - c2*z2 - c3*z2 + c4*z2 - c5*z1 + c6*z1 + c7*z1 - c8*z1)/((x1 - x2)*(y1 - y2)*(z1 - z2));
  double a5 = -(c1*y2 - c2*y2 - c3*y1 + c4*y1 - c5*y2 + c6*y2 + c7*y1 - c8*y1)/((x1 - x2)*(y1 - y2)*(z1 - z2));
  double a6 = -(c1*x2 - c2*x1 - c3*x2 + c4*x1 - c5*x2 + c6*x1 + c7*x2 - c8*x1)/((x1 - x2)*(y1 - y2)*(z1 - z2));
  double a7 = (c1 - c2 - c3 + c4 - c5 + c6 + c7 - c8)/((x1 - x2)*(y1 - y2)*(z1 - z2));
  
  // Interpolate to agent location
  double c_p = a0 + a1*Xp + a2*Yp + a3*Zp + a4*Xp*Yp + a5*Xp*Zp + a6*Yp*Zp + a7*Xp*Yp*Zp;
  
  return c_p;
  
}

/* The main routine */
void cells(mwSize Nmax,double *X,double *Y,double *Z,double *state,double *cellN,double *D,double *M,double *cycr,double *c_p,double *c,double *C,double *rates,double *hyp,double mu,double sigma,double h,double pdeT,double eta1,double eta2,double eta3,int rng_s,mwSize I,double *Xout,double *Yout,double *Zout,double *stateout,double *Dout,double *Mout,double *cycrout,double *Oxcout,double *Ncs,double *Cout)
{
   
  // Extract values from input arrays
  int Nr = cellN[0]; int Nyel = cellN[1]; int Ng = cellN[2]; int N = cellN[3]; int Nd = cellN[4]; // Counts of each cell type and living/dead cells
  double dmax = rates[0]; double dmin = rates[1]; double Rr = rates[2]; double Ry = rates[3]; double Rg = rates[4]; double mmax = rates[5]; double mmin = rates[6]; // Per capita rates
  double c_a = hyp[0]; double c_m = hyp[1]; double c_d = hyp[2];
  
  srand2(rng_s);
  
  int counter; int oldcounter;
  
  double t = 0; // Initialise t 
  
  while (t < pdeT && (N+Nd) < Nmax) {
      oldcounter = 0;
      counter = 0;
  
      double dr = 0;
      double mr = 0;
      double cr = 0;
      for (int q=0; q<N+Nd; q++) {
          dr += D[q];
          mr += M[q];
          cr += cycr[q];
      }
      double r = cr + mr + dr;
      // Sample timestep
      double dt = -log(unif()) / r;
      t += dt;
      
      // Break and output results without doing anything if t > pdeT
      if (t > pdeT) {
          t = pdeT;
          break;
      }
      
      // Resolve occurring events
      double u = unif()*r;
      
      // Cycling event
      if (u < cr) {
          
          int ind = ag_choice(cycr,cr);
        
          // Red to yellow event
          if (state[ind] == 1) {  
              state[ind] = 2;
              Nr -= 1;
              Nyel += 1;
              cycr[ind] = Ry; // Equation (2)
              
          // Yellow to green event
          } else if (state[ind] == 2) {
              state[ind] = 3;
              Nyel -= 1;
              Ng += 1;
              cycr[ind] = Rg; // Equation (3)
              
          // Green to two reds mitosis event    
          } else if (state[ind] == 3) {
              state[ind] = 1;
              Ng -= 1;
              Nr += 2;
              
              double xp = X[ind];
              double yp = Y[ind];
              double zp = Z[ind];
              
              // Sample a random direction to proliferate in, and place the new agents sigma apart
              double theta = 2*pi*unif(); // Azimuth angle is uniform between 0 and 2*pi
              double phi = acos(1 - 2*unif()); // Zenith angle must be sampled differently
              
              double xn1 = xp + 0.5*sigma*cos(theta)*sin(phi);
              double yn1 = yp + 0.5*sigma*sin(theta)*sin(phi);
              double zn1 = zp + 0.5*sigma*cos(phi);
              
              double xn2 = xp - 0.5*sigma*cos(theta)*sin(phi);
              double yn2 = yp - 0.5*sigma*sin(theta)*sin(phi);
              double zn2 = zp - 0.5*sigma*cos(phi);
              
              X[ind] = xn1;
              Y[ind] = yn1;
              Z[ind] = zn1;
              
              // Create the new agent
              X[N+Nd] = xn2;
              Y[N+Nd] = yn2;
              Z[N+Nd] = zn2;
              state[N+Nd] = 1;
              
              // Check nutrient concentrations and update agent rates
              c_p[ind] = interp3d(I,h,xn1,yn1,zn1,c);
              
              cycr[ind] = Rr*(pow(c_p[ind],eta1))/(pow(c_a,eta1) + pow(c_p[ind],eta1)); // Equation (1)
              M[ind] = (mmax - mmin)*(pow(c_p[ind],eta2))/(pow(c_m,eta2) + pow(c_p[ind],eta2)) + mmin; // Equation (4)
              D[ind] = (dmax-dmin)*(1 - (pow(c_p[ind],eta3))/(pow(c_d,eta3) + pow(c_p[ind],eta3))) + dmin; // Equation (5)
                            
              c_p[N+Nd] = interp3d(I,h,xn2,yn2,zn2,c);
              
              cycr[N+Nd] = Rr*(pow(c_p[N+Nd],eta1))/(pow(c_a,eta1) + pow(c_p[N+Nd],eta1)); // Equation (1)
              M[N+Nd] = (mmax - mmin)*(pow(c_p[N+Nd],eta2))/(pow(c_m,eta2) + pow(c_p[N+Nd],eta2)) + mmin; // Equation (4)
              D[N+Nd] = (dmax-dmin)*(1 - (pow(c_p[N+Nd],eta3))/(pow(c_d,eta3) + pow(c_p[N+Nd],eta3))) + dmin; // Equation (5)
              
              N += 1; // Account for increase in population
              
              // Update cell density
              int iold = round(xp/h);
              int jold = round(yp/h);
              int kold = round(zp/h);
              int oldind = pow(I,2)*kold + I*iold + jold;
              
              int inew1 = round(xn1/h);
              int jnew1 = round(yn1/h);
              int knew1 = round(zn1/h);
              int newind1 = pow(I,2)*knew1 + I*inew1 + jnew1;
              if (oldind != newind1) { // Cell moves from old volume to new one
                  C[oldind] -= 1/pow(h,3);
                  C[newind1] += 1/pow(h,3);
              }
              
              int inew2 = round(xn2/h);
              int jnew2 = round(yn2/h);
              int knew2 = round(zn2/h);
              int newind2 = pow(I,2)*knew2 + I*inew2 + jnew2;
              
              C[newind2] += 1/pow(h,3); // One new cell in the finite volume surrounding node point
          }
          
      // Movement event    
      } else if (u < cr + mr) {
          
          int ind = ag_choice(M,mr);
          
          double xp = X[ind];
          double yp = Y[ind];
          double zp = Z[ind];
          
          // Sample angles and move
          double mdist = mu;
          double theta = 2*pi*unif(); // Azimuth angle is uniform between 0 and 2*pi
          double phi = acos(1 - 2*unif()); // Zenith angle must be sampled differently
          
          double xn = xp + mdist*cos(theta)*sin(phi);
          double yn = yp + mdist*sin(theta)*sin(phi);
          double zn = zp + mdist*cos(phi);
          
          X[ind] = xn;
          Y[ind] = yn;
          Z[ind] = zn;
          
          // Update local nutrient concentration for agent, and rates
          c_p[ind] = interp3d(I,h,xn,yn,zn,c);
          
          if (state[ind] == 1) {
              cycr[ind] = Rr*(pow(c_p[ind],eta1))/(pow(c_a,eta1) + pow(c_p[ind],eta1)); // Equation (1). Only update cycling rate if agent is red
          }
          M[ind] = (mmax - mmin)*(pow(c_p[ind],eta2))/(pow(c_m,eta2) + pow(c_p[ind],eta2)) + mmin; // Equation (4)
          D[ind] = (dmax-dmin)*(1 - (pow(c_p[ind],eta3))/(pow(c_d,eta3) + pow(c_p[ind],eta3))) + dmin; // Equation (5)
          
          // Account for changes in cell density
          int iold = round(xp/h);
          int jold = round(yp/h);
          int kold = round(zp/h);
          int oldind = pow(I,2)*kold + I*iold + jold;
          
          int inew = round(xn/h);
          int jnew = round(yn/h);
          int knew = round(zn/h);
          int newind = pow(I,2)*knew + I*inew + jnew;
          
          if (oldind != newind) { // Cell moves from old volume to new one
              C[oldind] -= 1/pow(h,3);
              C[newind] += 1/pow(h,3);
          }
      
      // Death event    
      } else {
          
          int ind = ag_choice(D,dr); // Choose agent
          
          // Account for subpopulation changes
          if (state[ind] == 1) {
              Nr -= 1;
          } else if (state[ind] == 2) {
              Nyel -= 1;
          } else if (state[ind] == 3) {
              Ng -= 1; 
          }
          
          // "Kill" the agent
          state[ind] = 0;
          M[ind] = 0;
          D[ind] = 0;
          cycr[ind] = 0;
          
          // Account for population change
          N -= 1;
          Nd += 1;
          
          double xp = X[ind];
          double yp = Y[ind];
          double zp = Z[ind];
          
          // Account for change in cell density
          int iold = round(xp/h);
          int jold = round(yp/h);
          int kold = round(zp/h);
          int oldind = pow(I,2)*kold + I*iold + jold;
          
          C[oldind] -= 1/pow(h,3);
          
      }
      
  }
  
  // Prepare to export data
  Ncs[0] = Nr;
  Ncs[1] = Nyel;
  Ncs[2] = Ng;
  Ncs[3] = N;
  Ncs[4] = Nd;
  Ncs[5] = rseed;
  
  for (int q=0; q<Nmax; q++) {
      Xout[q] = X[q];
      Yout[q] = Y[q];
      Zout[q] = Z[q];
      cycrout[q] = cycr[q];
      Mout[q] = M[q];
      Dout[q] = D[q];
      stateout[q] = state[q];
      Oxcout[q] = c_p[q];
  }
  
  int q = 0;
  for (q=0; q<pow(I,3); q++) {
      Cout[q] = C[q];
  }  
  
}


/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    // Initialise gateway inputs
    mwSize Nmax;        // Sizes of X,Y,state,and other per capita storage
    mwSize countsize;   // Size of cellN and, subsequently, Ncs
    double *X;          // X positions of agents
    double *Y;          // Y positions of agents
    double *Z;          // Z positions of agents
    double *state;      // Cell cycle step of agents
    double *cellN;      // Number of cells in each type
    double *D;          // Per capita death rate storage
    double *M;          // Per capita movement rate storage
    double *cycr;       // Per capita cycling rate change (depends on cell cycle stage)
    double *c_p;        // Oxygen concentration at cell locations
    double *c;          // Oxygen concentration profile
    double *C;          // Cell count storage
    double *rates;      // Values of transition rates for each cell cycle step
    double *hyp;        // Nutrient threshold level storage
    double mu;          // Movement distance size  
    double sigma;       // Proliferation spread parameter
    double h;           // Node spacing
    double pdeT;        // Length of time to run the while loop
    double eta1;        // Index of Hill function (arrest)
    double eta2;        // Index of Hill function (movement)
    double eta3;        // Index of Hill function (death)
    const mwSize *dims; // Dimensions of c and C arrays
    mwSize ndim;        // Number of dimensions for c and C arrays (3)
    mwSize I;           // Number of x, y, z-nodes
    int rng_s;          // Random number generator seed
    
    // Initialise outputs
    double *Xout;       // X positions
    double *Yout;       // Y positions
    double *Zout;       // Z positions
    double *stateout;   // Cell cycle stage storage
    double *Dout;       // Death rate storage
    double *Mout;       // Movement rate storage
    double *cycrout;    // Cycling rate storage 
    double *Oxcout;     // Oxygen concentration at cell locations
    double *Ncs;        // Number of cells in each cycle stage storage
    double *Cout;       // Cell count storage
            
    // Grab inputs
    X = mxGetDoubles(prhs[0]);
    Nmax = mxGetM(prhs[0]);
    Y = mxGetDoubles(prhs[1]);
    Z = mxGetDoubles(prhs[2]);
    state = mxGetDoubles(prhs[3]);
    cellN = mxGetDoubles(prhs[4]);
    countsize = mxGetN(prhs[4]);
    D = mxGetDoubles(prhs[5]);
    M = mxGetDoubles(prhs[6]);
    cycr = mxGetDoubles(prhs[7]);
    c_p = mxGetDoubles(prhs[8]);
    c = mxGetDoubles(prhs[9]);
    dims = mxGetDimensions(prhs[9]);
    I = dims[0];
    C = mxGetDoubles(prhs[10]);
    rates = mxGetDoubles(prhs[11]);
    hyp = mxGetDoubles(prhs[12]);
    mu = mxGetScalar(prhs[13]);
    sigma = mxGetScalar(prhs[14]);
    h = mxGetScalar(prhs[15]);
    pdeT = mxGetScalar(prhs[16]);
    eta1 = mxGetScalar(prhs[17]);
    eta2 = mxGetScalar(prhs[18]);
    eta3 = mxGetScalar(prhs[19]);
    rng_s = mxGetScalar(prhs[20]);
    
    // Set-up outputs
    plhs[0] = mxCreateDoubleMatrix(Nmax,1,mxREAL);  // pointer to output matrix
    Xout = mxGetDoubles(plhs[0]);          // Actual output matrix
    
    plhs[1] = mxCreateDoubleMatrix(Nmax,1,mxREAL);  // pointer to output matrix
    Yout = mxGetDoubles(plhs[1]);          // Actual output matrix
    
    plhs[2] = mxCreateDoubleMatrix(Nmax,1,mxREAL);  // pointer to output matrix
    Zout = mxGetDoubles(plhs[2]);          // Actual output matrix
    
    plhs[3] = mxCreateDoubleMatrix(Nmax,1,mxREAL);  // pointer to output matrix
    stateout = mxGetDoubles(plhs[3]);          // Actual output matrix
    
    plhs[4] = mxCreateDoubleMatrix(Nmax,1,mxREAL);  // pointer to output matrix
    Dout = mxGetDoubles(plhs[4]);          // Actual output matrix
    
    plhs[5] = mxCreateDoubleMatrix(Nmax,1,mxREAL);  // pointer to output matrix
    Mout = mxGetDoubles(plhs[5]);          // Actual output matrix
    
    plhs[6] = mxCreateDoubleMatrix(Nmax,1,mxREAL);  // pointer to output matrix
    cycrout = mxGetDoubles(plhs[6]);          // Actual output matrix
    
    plhs[7] = mxCreateDoubleMatrix(Nmax,1,mxREAL);  // pointer to output matrix   
    Oxcout = mxGetDoubles(plhs[7]);            // Actual output matrix
    
    plhs[8] = mxCreateDoubleMatrix(countsize+1,1,mxREAL);  // pointer to output matrix
    Ncs = mxGetDoubles(plhs[8]);          // Actual output matrix
    
    plhs[9] = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);  // pointer to output matrix
    Cout = mxGetDoubles(plhs[9]);          // Actual output matrix

    
    // /* call the computational routine */
    cells(Nmax,X,Y,Z,state,cellN,D,M,cycr,c_p,c,C,rates,hyp,mu,sigma,h,pdeT,eta1,eta2,eta3,rng_s,I,Xout,Yout,Zout,stateout,Dout,Mout,cycrout,Oxcout,Ncs,Cout);

}