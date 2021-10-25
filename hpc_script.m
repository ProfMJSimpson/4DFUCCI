function hpc_script(sim_id)

%-% SET THE PARAMETERS FOR THE MODEL

% For testing purposes, I = 101 is recommended. The results we 
% show in the paper were calculated with high-performance computing with 
% I = 201.

% Experimental variability
% exp_N0 = [25700 35900 31000 32300 27800 28900 31000 26600 32400 29900];
% exp_rad = [232.75 260.13 247.76 251.23 238.97 242.19 247.93 235.47 251.48 244.89];

% All same initial conditions
exp_N0 = [ 30000 30000 30000 30000 30000 30000 30000 30000 30000 30000 ];
exp_rad = [245 245 245 245 245 245 245 245 245 245];

N0 = exp_N0(sim_id); % Initial population
Nmax = 300000; % Maximum number of cells for simulation
T = 240; % Final time of simulation (hours)
L = 4000; % Length of domain
I = 101; % Number of x, y, or z nodes
xmesh = linspace(0,L,I); % Generate x mesh
ymesh = linspace(0,L,I); % Generate y mesh
zmesh = linspace(0,L,I); % Generate z mesh

sigma = 12; % Proliferation spread parameter

alpha = 0.15; % Consumption-diffusion ratio
h = xmesh(2) - xmesh(1); % Node spacing parameter
c_b = 1; % Boundary nutrient concentration

r = exp_rad(sim_id); % Initial radius

mu = 12; % Movement distance
c_a = 0.4; % Critical arrest concentration
c_m = 0.5; % Critical movement concentration
c_d = 0.1; % Critical death concentration
dtsave = 2; % When to save information in simulation (MUST BE >= deltat)
pdeT = 1; % Time between PDE solutions 
eta1 = 5; % Hill index for arrest
eta2 = 5; % Hill index for movement
eta3 = 15; % Hill index for death

parms = [ alpha h c_b r mu c_a c_m c_d I pdeT eta1 eta2 eta3]; % List parameters

dmax = 2; % Maximum death rate
dmin = 0.0005; % Minimum death rate
Rr = 0.0470; % G1-eS rate (unimpeded)
Ry = 0.4898; % eS-S/G2/M rate
Rg = 0.062; % Mitosis rate
mmax = 0.12; % Maximum movement rate
mmin = mmax/2; % Minimum movement rate

rates = [ dmax dmin Rr Ry Rg mmax mmin ]; % List rate parameters

% Initialise storage data
numt = T / dtsave + 1;
tg = zeros(numt,1); % Time point storage
for i = 1:numt
    tg(i) = (i-1)*dtsave;
end

% Run script (call ibm3d.m function)
[Nvec,dvec,rvec,yvec,gvec,arrvec,radii,radarr,radnec,X,Y,Z,state,c,c_p,Xsnap,Ysnap,Zsnap,arrsnap,statesnap,csnap,Nsnap] = ibm3d(N0,Nmax,T,L,sigma,parms,rates,tg,xmesh,ymesh,zmesh,sim_id);

% Save workspace
mkdir Workspaces
filename = ['Workspaces/workspace_out_' num2str(sim_id) '.mat'];
save(filename,'Nvec','dvec','rvec','yvec','gvec','arrvec','radii','radarr','radnec','X','Y','Z','state','c','c_p','Xsnap','Ysnap','Zsnap','arrsnap','statesnap','csnap','Nsnap')

end