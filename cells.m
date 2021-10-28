% cells.c Simulate the agent-level behaviours in the IBM
%   cells(X,Y,Z,state,cellN,D,M,cycr,c_p,c,C,rates,hyp,mu,sigma,h,pdeT,eta1,eta2,eta3,rng_s)
%   simulates the agent-level behaviours of cell cycling, migration, and
%   death
%
%   INPUTS:
%       X       : [X1;X2;...;XNmax]              : (Nmax x 1)  vector
%       Y       : [Y1;Y2;...;YNmax]              : (Nmax x 1)  vector
%       Z       : [Z1;Z2;...;ZNmax]              : (Nmax x 1)  vector
%       state   : [state1;...;stateNmax]         : (Nmax x 1)  vector
%       cellN   : [Nr,Ny,Ng,N,Nd]                : (1 x 5)  vector
%       D       : [D1;...;DNmax]                 : (Nmax x 1)  vector
%       M       : [M1;...;MNmax]                 : (Nmax x 1)  vector
%       cycr    : [cyc1;...;cycNmax]             : (Nmax x 1)  vector
%       c_p     : [c1;...;cNmax]                 : (Nmax x 1)  vector
%       c       : [c11,...,c1I;...;cI1,...,cII]  : (Nmax x 1)  matrix
%       C       : [C11,...,C1I;...;CI1,...,CII]  : (Nmax x 1)  matrix
%       rates   : [dmax,dmin,Rr,Ry,Rg,mmax,mmin] : (1 x 7)  vector
%       hyp     : [c_a,c_m,c_d]                  : (1 x 3)  vector
%       mu      : (mu)                           : scalar migration dist.
%       sigma   : (sigma)                        : scalar prolif. dist.
%       h       : (h)                            : scalar node spacing
%       pdeT    : (pdeT)                         : scalar cells.c end time
%       eta1    : (eta1)                         : scalar G1-eS Hill index
%       eta2    : (eta2)                         : scalar migr. Hill index
%       eta3    : (eta3)                         : scalar death Hill index
%       rngs    : (int)                          : rng seed
% 
%   OUTPUTS:
%       [Xout,Yout,Zout,stateout,Dout,Mout,cycrout,Oxcout,Ncs,Cout] = cells( ... )
%
%   where:
%       Xout    : X positions after pdeT hours
%       Yout    : Y positions after pdeT hours
%       Zout    : Z positions after pdeT hours
%       stateout: Agent cycle/death status after pdeT hours
%       Dout    : Agent death rates after pdeT hours 
%       Mout    : Agent migration rates after pdeT hours 
%       cycrout : Agent cycling rates after pdeT hours 
%       Oxcout  : Local agent nutrient concentrations after pdeT hours 
%       Ncs     : Number of agents per subpopulation after pdeT hours
%       Cout    : [x,y,z] agent density matrix after pdeT hours
%
%   Copyright 2021 Jonah J. Klowss
%   MEX Function