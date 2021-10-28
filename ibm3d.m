function [Nvec,dvec,rvec,yvec,gvec,arrvec,radii,radarr,radnec,X,Y,Z,state,c,c_p,Xsnap,Ysnap,Zsnap,arrsnap,statesnap,csnap,Nsnap] = ibm3d(N0,Nmax,T,L,sigma,parms,rates,tg,xmesh,ymesh,zmesh,sim_id)
%--------- IBM3D ---------%
% ibm3d.m initialises and performs the calculations for the IBM. Nutrient
% concentration updates are performed at the top level, and agent-level
% behaviours are performed in cells.c (see cells.m).
%
% Inputs:
%     N0: Initial population
%     Nmax: Maximum allowed population
%     T: Termination time for simulation
%     L: Domain length
%     sigma: Proliferation distance
%     parms: [alpha h c_b r mu c_a c_m c_d I pdeT eta1 eta2 eta3]
%           WHERE
%               alpha: Consumption-diffusion ratio
%               h: Node spacing
%               c_b: Nutrient concentration on the boundary
%               r: Initial radius of spheroid
%               mu: Migration distance
%               c_a: G1 arrest threshold
%               c_m: Migration threshold
%               c_d: Death threshold
%               I: Number of x, y, or z nodes
%               pdeT: Time between steady-state solutions
%               eta1: Arrest Hill index
%               eta2: Migration Hill index
%               eta3: Death Hill index
%     rates: [dmax dmin Rr Ry Rg mmax mmin]
%           WHERE
%               dmax: Maximum death rate
%               dmin: Minimum death rate
%               Rr: Maximum G1-eS rate
%               Ry: Constant eS-S/G2/M rate
%               Rg: Constant S/G2/M-mitosis rate
%               mmax: Maximum migration rate
%               mmin: Minimum migration rate
%     tg: Times to save population and radius data
%     xmesh: x in [0,L] interval divided up I times with spacing h
%     ymesh: y in [0,L] interval divided up I times with spacing h
%     zmesh: z in [0,L] interval divided up I times with spacing h
%     sim_id: Name of simulation and random number generator seed
%
% Outputs: 
%     Nvec: Time series living population data
%     dvec: Time series deaad population data
%     rvec: Time series red (all) population data
%     yvec: Time series yellow population data
%     gvec: Time series green population data
%     arrvec: Time series arrested red population data
%     radii: Time series outer radius data
%     radarr: Time series arrested radius data
%     radnec: Time series necrotic radius data
%     X: Final x positions of agents
%     Y: Final y positions of agents
%     Z: Final z positions of agents
%     state: Final states of agents
%     c: Final nutrient concentration profile
%     c_p: Final local nutrient concentration for all agents
%     Xsnap: Agent x positions every 24 hours
%     Ysnap: Agent y positions every 24 hours
%     Zsnap: Agent z positions every 24 hours
%     arrsnap: Arrested red agent [x y z] positions every 24 hours
%     statesnap: Agent cycling/dead status every 24 hours
%     csnap: Nutrient concentration profile every 24 hours
%     Nsnap: Total population every 24 hours

% Initialise seed
rng(sim_id)

% Extract parameters from storage vectors

alpha = parms(1);   % Consumption-diffusion ratio
h = parms(2);       % Node spacing
c_b = parms(3);     % Maximum nutrient (saturation) level
r = parms(4);       % Radius of initial spheroid
mu = parms(5);      % Movement distance
c_a = parms(6);     % G1 arrest threshold
c_m = parms(7);     % Critical nutrient concentration (movement)
c_d = parms(8);     % Cell death threshold
I = parms(9);       % Number of x-, y-, or z-nodes
pdeT = parms(10);   % Frequency of updating nutrient PDE profile
eta1 = parms(11);   % Hill function (arrest) index
eta2 = parms(12);   % Hill function (movement) index
eta3 = parms(13);   % Hill function (death) index

dmax = rates(1);    % Maxmimum death rate when nutrient below death threshold
dmin = rates(2);    % Intrinsic background death rate
Rr = rates(3);      % Red to yellow transition rate
Ry = rates(4);      % Yellow to green transition rate
Rg = rates(5);      % Green to two reds mitosis rate
mmax = rates(6);    % Maximum movement rate
mmin = rates(7);    % Minimum movement rate

tgsize = length(tg); % Size for vectors storing population and dead count

%% Initialise 
t = 0; % Initialise time

Nn = I^3; % Calculate the number of nodes

c_p = zeros(Nmax,1); % Oxygen concentrations at cell positions

% Initialise output data vectors
Nvec = zeros(tgsize,1);
dvec = zeros(tgsize,1);
rvec = zeros(tgsize,1);
yvec = zeros(tgsize,1);
gvec = zeros(tgsize,1);
radii = zeros(tgsize,1);
radarr = zeros(tgsize,1);
radnec = zeros(tgsize,1);
arrvec = zeros(tgsize,1);

% Initialise agent placement uniformly in sphere
rvals = 2*rand(N0,1) - 1;
elev = asin(rvals); % Elevation angle sampling

azim = 2*pi*rand(N0,1); % Azimuthal angle sampling

rads = r*(rand(N0,1)).^(1/3); % Radius sampling

[Xi,Yi,Zi] = sph2cart(azim,elev,rads); % Convert from spherical polar to cartesian coordinates
% Move to the domain centre
Xi = Xi + L/2;
Yi = Yi + L/2;
Zi = Zi + L/2;

% Put initial agent locations into the vectors for all locations
X = [Xi ; zeros(Nmax - N0,1)];
Y = [Yi ; zeros(Nmax - N0,1)];
Z = [Zi ; zeros(Nmax - N0,1)];

N = N0; % Number of living agents
Nd = 0; % Number of dead agents

Nvec(1) = N; % Living population at t = 0 is N
radii(1) = r; % Initial radius is r

C = celldens(I,xmesh,ymesh,zmesh,X,Y,Z,N,h); % Calculate agent density

% Snapshot storage initialisation (data for each day)
snapcount = T/24+1;
Xsnap = zeros(Nmax,snapcount);
Ysnap = zeros(Nmax,snapcount);
Zsnap = zeros(Nmax,snapcount);
arrsnapX = zeros(Nmax,snapcount);
arrsnapY = zeros(Nmax,snapcount);
arrsnapZ = zeros(Nmax,snapcount);
statesnap = zeros(Nmax,snapcount);
csnap = zeros(Nn,snapcount);
Nsnap = zeros(1,snapcount);

%% Initialise nutrient concentration

% Create the n^3xn^3 coefficient matrix for nutrient concentration calculation
iter = 0;
rhs = zeros(I^3,1); % Right-hand side vector
adj = zeros(I^3 - 2*I^2,1); % Adjacent diagonals
mdiag = zeros(I^3 - 2*I^2,1); % Main diagonal

for i = 1:I^2
    rhs(i) = c_b;
    
    ind = I^3 - (i-1);
    rhs(ind) = c_b;
    
end

% Construct coefficient matrix as according to Equation (9)
for j = 1:length(adj)
    iter = iter + 1;
    ind = I^2 + j; % Offset index to avoid first z-layer
    if iter < I + 1
        mdiag(j) = 1;
        rhs(ind) = c_b;
    elseif mod(iter,I) == 0
        mdiag(j) = 1;
        rhs(ind) = c_b;
        if iter == I^2
            iter = 0;
        end
    elseif mod(iter,I) == 1
        mdiag(j) = 1;
        rhs(ind) = c_b;
    elseif iter > I^2 - I + 1
        mdiag(j) = 1;
        rhs(ind) = c_b;
    else 
        adj(j) = 1/h^2;
        mdiag(j) = -6/h^2-alpha*C(ind);
    end
end

e1 = [ adj ; zeros(2*I^2,1)];
e2 = [ zeros(I^2-I,1) ; adj ; zeros(I^2+I,1) ];
e3 = [ zeros(I^2-1,1) ; adj ; zeros(I^2+1,1) ];
e4 = [ ones(I^2,1) ; mdiag ; ones(I^2,1) ];
e5 = [ zeros(I^2+1,1) ; adj ; zeros(I^2-1,1) ];
e6 = [ zeros(I^2+I,1) ; adj ; zeros(I*I-I,1) ];
e7 = [ zeros(2*I^2,1) ; adj ];

% Construct the coefficient matrix
As = spdiags([e1 e2 e3 e4 e5 e6 e7],[-I^2 -I -1:1 I I^2],Nn,Nn);

c_old = ones(Nn,1);
[c_v,flag] = gmres(As,rhs,[],1e-8,300,[],[],c_old); % Heightened restrictions to initialise better. Call flag to prevent printing of GMRES status
c_old = c_v;

c = reshape(c_v,[I I I]); 

next_c_solve = pdeT; % Set next time point at which the PDE should be solved

%% Complete initialisations

% Rate vectors (death, cycling, movement)
D = zeros(Nmax,1);
M = zeros(Nmax,1);
cycr = zeros(Nmax,1);


% Initialise cell type (red, yellow, green) storage - no initial dead cells
Nr = 0; Nyel = 0; Ng = 0; % Initialise counts of cells of each type

% Find proportion of red/yellow/green in freely cycling conditions
full_cycle = 1/Rr + 1/Ry + 1/Rg; % Time of cell cycle
redprop = 1/Rr / full_cycle; % Proportion of time spent in red
yelprop = 1/Ry / full_cycle; % Proportion of time spent in red
greprop = 1/Rg / full_cycle; % Proportion of time spent in red

% Initialise storage for indicators of cell cycle or death status
state = zeros(Nmax,1);

for q = 1:N0
    c_p(q) = interp3d(h,X(q),Y(q),Z(q),c);
    
    if c_p(q) > c_a % Identify cycling status of cell (region-wise)
        u = rand;
        if u < redprop
            state(q) = 1;
            Nr = Nr + 1;
        elseif u < redprop + yelprop
            state(q) = 2;
            Nyel = Nyel + 1;
        elseif u < redprop + yelprop + greprop
            state(q) = 3;
            Ng = Ng + 1;
        end
    else % Tuned to get an initial arrest radius matching 
        if c_p(q) > c_a - 0.137 %Rr(c) < 0.01*Rr -- consider this as an initial arrest region
            u = rand;
            if u < 0.16 % reduced chance of greens and yellows in arrest zone. Arrest calculated when green radial density drops below 20%
                u2 = rand;
                if u2 < (1/Ry)/(1/Ry + 1/Rg) % Time in yellow / time in yellow and green
                    state(q) = 2;
                    Nyel = Nyel + 1;
                else % (1/Rg)/(1/Ry + 1/Rg) condition: time in green / time in yellow and green
                    state(q) = 3;
                    Ng = Ng + 1;
                end
            else
                state(q) = 1;
                Nr = Nr + 1;
            end
        else % All red when local concentration under c_a - 0.137
            state(q) = 1;
            Nr = Nr + 1;
        end
    end
    
    % Set rates of cycling, movement, and death
    if state(q) == 1
        cycr(q) = Rr*(c_p(q)^eta1)./(c_a^(eta1) + c_p(q)^eta1); % Equation (1)
    elseif state(q) == 2
        cycr(q) = Ry; % Equation (2)
    elseif state(q) == 3
        cycr(q) = Rg; % Equation (3)
    end
    M(q) = (mmax - mmin)*(c_p(q)^eta2)./(c_m^(eta2) + c_p(q)^eta2) + mmin; % Equation (4)
    D(q) = (dmax-dmin)*(1 - c_p(q)^(eta3)./(c_d^(eta3) + c_p(q)^(eta3))) + dmin; % Equation (5)
end

% Initialise agent type counts., no initial dead cells
rvec(1) = Nr;
yvec(1) = Nyel;
gvec(1) = Ng;

% Save initial state as first snapshot
Xsnap(:,1) = X;
Ysnap(:,1) = Y;
Zsnap(:,1) = Z;
statesnap(:,1) = state;
csnap(:,1) = c_old;
Nsnap(1) = N0;
arrsnaptempX = []; arrsnaptempY = []; arrsnaptempZ = [];
centroid_i = [mean(X(1:N+Nd)) mean(Y(1:N+Nd)) mean(Z(1:N+Nd))];
radq = 0; % Initialise estimate for arrest radius at t = 0
for q = 1:N+Nd
    if state(q) == 1 && c_p(q) < c_a - 0.137 % When Rr(c) < 0.01Rr
        arrsnaptempX = [arrsnaptempX ; X(q)];
        arrsnaptempY = [arrsnaptempY ; Y(q)];
        arrsnaptempZ = [arrsnaptempZ ; Z(q)];
        
        radqn = sqrt(((X(q)-centroid_i(1)).^2 + (Y(q)-centroid_i(2)).^2 + (Z(q)-centroid_i(3)).^2));
        if radqn > radq
            radq = radqn; % Estimate initial arrest radius from placement of agents in arrested region
        end
    end
end

% Radius and count of arrested region at t = 0
radarr(1) = radq;
arrvec(1) = length(arrsnaptempX);

% Arrested agent locations at t = 0
arrsnapX(1:length(arrsnaptempX),1) = arrsnaptempX;
arrsnapY(1:length(arrsnaptempY),1) = arrsnaptempY;
arrsnapZ(1:length(arrsnaptempZ),1) = arrsnaptempZ;

%% While loop

idx = 2; % Iteration counter for saving system state (population data)
snapidx = 2; % Iteration counter for saving the snapshots of spheroid growth

rng_s = sim_id; % Set initial rng_s (C code seed) value as sim_id;

sz = 6; % Sizes of cell representation in image processing for radius calculations
thr = 0.3; % Threshold for minimum connected area size for outer radius calculations

while t < T && N+Nd < Nmax
    
    t = t + pdeT; % Take a time step forward    
    
    if t >= tg(idx) % If at next save step, store the current system state
        Nvec(idx) = N;
        dvec(idx) = Nd;
        rvec(idx) = Nr;
        yvec(idx) = Nyel;
        gvec(idx) = Ng;
        
        % Find the dead agent locations
        deadX = []; deadY = []; deadZ = [];
        for q = 1:N+Nd
            if state(q) == 0
                deadX = [deadX ; X(q)];
                deadY = [deadY ; Y(q)];
                deadZ = [deadZ ; Z(q)];
            end
        end
        
        % Radius calculations
        [r_o,r_n] = radcalc(X,Y,Z,state,L,sz,thr);
        radii(idx) = r_o;
        radnec(idx) = r_n;
        
        % Arrest radius and agent classifications
        % Calculate the spheroid centroid (mean position of all agents)
        centroid = [mean(X(1:N+Nd)) mean(Y(1:N+Nd)) mean(Z(1:N+Nd))];
        
        r_a = arrestradius(X,Y,Z,state,centroid,radii(idx));
        radarr(idx) = r_a;
        
        arrsnaptempX = []; arrsnaptempY = []; arrsnaptempZ = [];
        for q = 1:N+Nd
            agent_pos_rad = sqrt(((X(q)-centroid(1)).^2 + (Y(q)-centroid(2)).^2 + (Z(q)-centroid(3)).^2));
            if state(q) == 1 && agent_pos_rad < r_a
                arrsnaptempX = [arrsnaptempX ; X(q)];
                arrsnaptempY = [arrsnaptempY ; Y(q)];
                arrsnaptempZ = [arrsnaptempZ ; Z(q)];
            end
        end
        
        arrvec(idx) = length(arrsnaptempX);
        if arrvec(idx) == 0
            radarr(idx) = 0;
        end
        
        idx = idx + 1; % Update iteration count to next step
    end
    
    % If surpassed the next day time point
    if t >= (snapidx-1)*24
        Xsnap(:,snapidx) = X;
        Ysnap(:,snapidx) = Y;
        Zsnap(:,snapidx) = Z;
        statesnap(:,snapidx) = state;
        csnap(:,snapidx) = c_old;
        Nsnap(snapidx) = N+Nd;
        
        arrsnapX(1:length(arrsnaptempX),snapidx) = arrsnaptempX;
        arrsnapY(1:length(arrsnaptempY),snapidx) = arrsnaptempY;
        arrsnapZ(1:length(arrsnaptempZ),snapidx) = arrsnaptempZ;
        
        snapidx = snapidx + 1;
    end
    
    % After saving previous system state, break the loop if *OVER* endtime or there are no agents.
    if t > T || N == 0
        break;
    end
    
    % Update nutrient PDE and nutrient at cell locations
    [c_v,flag] = gmres(As,rhs,[],1e-6,200,[],[],c_old);
    
    c_old = c_v; % Copy the vectorised nutrient concentration across to use as guess for next solution
    c = reshape(c_v,[I,I,I]); % Reshape into I*I*I profile (3D)
    
    next_c_solve = next_c_solve + pdeT; % Update to the next solution time point
    
    % Update local nutrient concentrations at cell locations
    for q = 1:N+Nd
        c_p(q) = interp3d(h,X(q),Y(q),Z(q),c); % Interpolate to find local nutrient concentration
        
        % Update nutrient dependent rates
        if state(q) ~= 0 % If not dead
            M(q) = (mmax - mmin)*(c_p(q)^eta2)./(c_m^(eta2) + c_p(q)^eta2) + mmin; % Equation (4)
            D(q) = (dmax-dmin)*(1 - c_p(q)^(eta3)./(c_d^(eta3) + c_p(q)^(eta3))) + dmin; % Equation (5)
        end
        if state(q) == 1
            cycr(q) = Rr*(c_p(q)^eta1)./(c_a^(eta1) + c_p(q)^eta1); % Equation (1)
        end
    end
    
    % Resolve the cell behaviour events up to the next PDE solution time
    cellN = [Nr Nyel Ng N Nd]; % Cell numbers
    hyp = [c_a c_m c_d]; % Hypoxia effect levels
    [Xout,Yout,Zout,stateout,Dout,Mout,cycrout,Oxcout,Ncs,Cout] = cells(X,Y,Z,state,cellN,D,M,cycr,c_p,c,C,rates,hyp,mu,sigma,h,pdeT,eta1,eta2,eta3,rng_s);
    
    % Recover outputs to main variables
    X = Xout;
    Y = Yout;
    Z = Zout;
    cycr = cycrout;
    state = stateout;
    M = Mout;
    D = Dout;
    C = Cout;
    c_p = Oxcout;
    Nr = Ncs(1); Nyel = Ncs(2); Ng = Ncs(3); N = Ncs(4); Nd = Ncs(5); rng_s = Ncs(6);
    
    % Update coefficient matrix -- only the main diagonal will change
    % (changes to agent density) according to Equation (9)
    initer = 0;
    for j = 1:length(adj)
        initer = initer + 1;
        ind = I^2 + j; % Offset index to avoid first z-layer
        if initer < I + 1
            mdiag(j) = 1;
        elseif mod(initer,I) == 0
            mdiag(j) = 1;
            if initer == I^2
                initer = 0;
            end
        elseif mod(initer,I) == 1
            mdiag(j) = 1;
        elseif initer > I^2 - I + 1
            mdiag(j) = 1;
        else
            mdiag(j) = -6/h^2-alpha*C(ind);
        end
    end
    
    e4 = [ ones(I^2,1) ; mdiag ; ones(I^2,1) ];
    As = spdiags([e1 e2 e3 e4 e5 e6 e7],[-I^2 -I -1:1 I I^2],Nn,Nn);
    
end

arrsnap = [arrsnapX , arrsnapY , arrsnapZ];

%% Subroutines

function C = celldens(I,xvec,yvec,zvec,X,Y,Z,N,h)
% CELLDENS Calculates the cell density in the control volume around the
% nodes of the FDM mesh

C = zeros(I,I,I);

for k = 1:I
    zp = zvec(k);
    for i = 1:I
        xp = xvec(i);
        for j = 1:I
            yp = yvec(j);
            ind = I^2*(k-1) + I*(i-1) + j;
            Ccount = 0;
            
            for q = 1:N
                xag = X(q);
                yag = Y(q);
                zag = Z(q);
                
                xsp = abs(xag - xp);
                ysp = abs(yag - yp);
                zsp = abs(zag - zp);
                
                if ((xsp < h/2) && (ysp < h/2) && (zsp < h/2))
                    Ccount = Ccount + 1;
                end
            end
            
            C(ind) = Ccount/h^3;
            
        end
    end
end

end

function c_p = interp3d(h,Xp,Yp,Zp,c)
% INTERP3D Three dimensional linear interpolation of the nutrient
% concentration at the location of the cell

i = ceil(Xp/h);
j = ceil(Yp/h);
k = ceil(Zp/h);

x1 = h*(i-1);
x2 = h*i;
y1 = h*(j-1);
y2 = h*j;
z1 = h*(k-1);
z2 = h*k;

c1 = c(j,i,k);
c2 = c(j,i+1,k);
c3 = c(j+1,i,k);
c4 = c(j+1,i+1,k);
c5 = c(j,i,k+1);
c6 = c(j,i+1,k+1);
c7 = c(j+1,i,k+1);
c8 = c(j+1,i+1,k+1);

a0 = -(c1*x2*y2*z2 - c2*x1*y2*z2 - c3*x2*y1*z2 + c4*x1*y1*z2 - c5*x2*y2*z1 + c6*x1*y2*z1 + c7*x2*y1*z1 - c8*x1*y1*z1)/((x1 - x2)*(y1 - y2)*(z1 - z2));
a1 = (c1*y2*z2 - c2*y2*z2 - c3*y1*z2 + c4*y1*z2 - c5*y2*z1 + c6*y2*z1 + c7*y1*z1 - c8*y1*z1)/((x1 - x2)*(y1 - y2)*(z1 - z2));
a2 = (c1*x2*z2 - c2*x1*z2 - c3*x2*z2 + c4*x1*z2 - c5*x2*z1 + c6*x1*z1 + c7*x2*z1 - c8*x1*z1)/((x1 - x2)*(y1 - y2)*(z1 - z2));
a3 = (c1*x2*y2 - c2*x1*y2 - c3*x2*y1 + c4*x1*y1 - c5*x2*y2 + c6*x1*y2 + c7*x2*y1 - c8*x1*y1)/((x1 - x2)*(y1 - y2)*(z1 - z2));
a4 = -(c1*z2 - c2*z2 - c3*z2 + c4*z2 - c5*z1 + c6*z1 + c7*z1 - c8*z1)/((x1 - x2)*(y1 - y2)*(z1 - z2));
a5 = -(c1*y2 - c2*y2 - c3*y1 + c4*y1 - c5*y2 + c6*y2 + c7*y1 - c8*y1)/((x1 - x2)*(y1 - y2)*(z1 - z2));
a6 = -(c1*x2 - c2*x1 - c3*x2 + c4*x1 - c5*x2 + c6*x1 + c7*x2 - c8*x1)/((x1 - x2)*(y1 - y2)*(z1 - z2));
a7 = (c1 - c2 - c3 + c4 - c5 + c6 + c7 - c8)/((x1 - x2)*(y1 - y2)*(z1 - z2));

c_p = a0 + a1*Xp + a2*Yp + a3*Zp + a4*Xp*Yp + a5*Xp*Zp + a6*Yp*Zp + a7*Xp*Yp*Zp;

end

function [rad_im,necrad] = radcalc(X,Y,Z,state,L,sz,thr)

imgpx = zeros(L+1,L+1);

cent = find((abs(Z - L/2)) < 18);
liv = find(state ~= 0);
Zinds = intersect(liv,cent);
Xrs = X(Zinds);
Yrs = Y(Zinds);
Xpxs = round(Xrs);
Ypxs = round(Yrs);

non_zero = find(Xpxs,1,'last');

im_idx = sub2ind([L+1 L+1],Ypxs(1:non_zero),Xpxs(1:non_zero));
imgpx(im_idx) = 1;

nears = bwdist(imgpx);
img = nears < sz;
img = flipdim(img,1); % Flip so that low Y values correspond with bottom of image

sph_bwareaopen = bwareaopen(img,ceil(thr*sum(img,"all")));

sph_imfill = imfill(sph_bwareaopen,'holes');

% Generate mask of known cell locations
img_d = zeros(L+1,L+1);

Zdead = intersect(find(state == 0),find(abs(Z - L/2) < 18));
Xd = X(Zdead);
Yd = Y(Zdead);

Xpxs_d = round(Xd);
Ypxs_d = round(Yd);
idx_d = sub2ind([L+1 L+1],Ypxs_d,Xpxs_d);
img_d(idx_d) = 1;

nears = bwdist(img_d);
dead_mask = nears < sz;
dead_mask = flipdim(dead_mask,1);

nec_msk = ~sph_bwareaopen & sph_imfill; % Look for region present in filled image but not in unfilled image (necrotic area)

nec_clr = bwareafilt(nec_msk,1); % Remove any unnecessary peripheral segments

nec_region = regionprops(nec_clr); % Establish region of necrotic cells
if ~isempty(nec_region)
    nec_Area = nec_region.Area; % Extract area
else
    nec_Area = 0;
end

% Code to check whether necrotic area is only unfilled space or is where
% dead cells are
overlap = dead_mask & nec_msk; % Create a mask of overlap between dead cell locations and 

overlap_clr = bwareafilt(overlap,1);

overlap_region = regionprops(overlap_clr);
if ~isempty(overlap_region)
    overlapA = overlap_region.Area;
else
    overlapA = 0;
end
ov_rad = sqrt(overlapA/pi);

if ov_rad < 12
    necrad = 0;
else
    necrad = sqrt(nec_Area/pi);
end

sph_regionprops = regionprops(sph_imfill);
if ~isempty(sph_regionprops)
    Area = sph_regionprops.Area;
else
    Area = 0;
end

rad_im = sqrt(Area/pi);    

if necrad < 18 % No sense to a necrotic radius less than a cell size, also accounting for empty space
    necrad = 0;
end

end

function arr_rad = arrestradius(X,Y,Z,state,centroid,radius)

thresh = 0.2; % Threshold at which a region is considered arrested

grerads = [];
for q = 1:N+Nd
    if state(q) == 3
        rad = sqrt((X(q)-centroid(1)).^2 + (Y(q)-centroid(2)).^2 + (Z(q)-centroid(3)).^2);
        grerads = [grerads ; rad];
    end
end

rbw = 20; % (R)adial (b)in (w)idth for density calculations 

bins = 0:rbw:10*ceil(radius/10); % Round to the next 10 for the radius
data = histcounts(grerads,bins);
centres = zeros(length(bins)-1,1);
for i = 1:length(bins)-1
    centres(i) = (bins(i+1)+bins(i))/2;
end

dens1 = data./(4.*pi.*centres'.^2*(rbw)); % This approximates the radial shell volume extremely well
scaled_dens = dens1' ./ max(dens1);

% Remove the outer tail of the cell concentration distribution, as this
% leads to a better fit of the important data at the boundary of the
% arrested region
[~,index] = max(scaled_dens(5:length(scaled_dens))); % Only consider the maximum outside of the direct interior (useful for lower time)
index = index + 4; % Re-offset the index to bring it in terms of the index of the full scaled_dens array
fullcent = centres;
%Remove most central bin
centres = centres(3:index);
scaled_dens = scaled_dens(3:index);

% Maximum distance
Dmax = (fullcent(2) - fullcent(1))/2 + fullcent(end);

% Function to fit
fun = @(r,gamma) gamma(1) * exp(-exp(gamma(2) * (gamma(3) - r)));

% Residual  function
res = @(gamma) sum((scaled_dens - fun(centres,gamma)).^2);

% Start points to try
Gamma = [0.5,0.1,0.2*Dmax; 0.5,0.1,0.5*Dmax;  0.5,0.1,0.8*Dmax];

% Fit using fminsearch (for each starting point)
resids = zeros(size(Gamma,1),1);
for i = 1:size(Gamma,1)
    [gamma,resid] = fminsearch(res,Gamma(i,:),optimset('Display','none'));
    Gamma(i,:)    = gamma;
    resids(i) = resid;
end

% Find best fit
[~,idox] = min(resids);
gamma = Gamma(idox,:);
func = @(r) fun(r,gamma);

% Gompertz function
fun2  = @(r) gamma(1) * exp(-exp(gamma(2) * (gamma(3)-r)));

% Gompertz inverse function
finv = @(i) gamma(3) - log(-log(i / gamma(1))) / gamma(2);

% Check if green everywhere (i.e., phase 1)
if fun2(Dmax) < fun2(0) || fun2(0) > thresh
    arr_rad = 0;
% Otherwise, identify boundary (when fit function reaches 0.2)
else
    arr_rad = finv(thresh);
end

if arr_rad < 12 % No meaning to distinguishing an arrested region in a radius less than one cell diameter
    arr_rad = 0;
end

end

end