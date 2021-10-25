% Clean up before running
close all;
clc;

%% Establish parameters as they are in IBM (to use in figure generation)

Nmax = 300000;

T = 240;
L = 4000;
I = 101;
xmesh = linspace(0,L,I);
ymesh = linspace(0,W,I);
zmesh = linspace(0,H,I);

sigp = 12; sigm = 0;
sigs = [ sigp sigm ];

Diff = 7.2e6;
h = xmesh(2) - xmesh(1);
kappa = 1.08e6;
Omax = 1;
ri = 245;
mu = 12;
c_a = 0.4;
c_m = 0.5;
c_d = 0.1;
dtsave = 2;
minDT = 1;
eta1 = 5;
eta2 = 5;
eta3 = 15;

parms = [ Diff h kappa Omax ri mu c_a c_m c_d I minDT eta1 eta2 eta3];

dmax = 2;
dmin = 0.0005;
Rr = 0.0470;
Ry = 0.4898;
Rg = 0.0619;
mmax = 0.12;
mmin = mmax/2;
rates = [ dmax dmin Rr Ry Rg mmax mmin ];

numt = T / dtsave + 1;
tg = zeros(numt,1);
for i = 1:numt
    tg(i) = (i-1)*dtsave;
end
TgDN = 24/dtsave; % Number of time points in a day (e.g. if dtsave = 2 hours, TgDN = 12)

runcount = 10; % Change this to the number of runs performed

Ndays = T/24; % Number of days represented by T (hours)

%% Import variables

Nall = zeros(length(tg),runcount);
dall = zeros(length(tg),runcount);
redall = zeros(length(tg),runcount);
yelall = zeros(length(tg),runcount);
greall = zeros(length(tg),runcount);
arrall = zeros(length(tg),runcount);
radall = zeros(length(tg),runcount);
radarrall = zeros(length(tg),runcount);
radnecall = zeros(length(tg),runcount);
Xsnapall = zeros(Nmax,runcount*(Ndays+1));
Ysnapall = zeros(Nmax,runcount*(Ndays+1));
Zsnapall = zeros(Nmax,runcount*(Ndays+1));
statesnapall = zeros(Nmax,runcount*(Ndays+1));
csnapall = zeros(I^3,runcount*(Ndays+1));
arrsnapall = zeros(Nmax,runcount*(Ndays+1)*3);

for q = 1:runcount
    workspace_name = ['Workspaces\workspace_out_' num2str(q) '.mat'];
    load(workspace_name)
    
    % Pass variables across into storage, not including data for
    % representative figures (X,Y,Z,c, etc)
    Nall(:,q) = Nvec;
    dall(:,q) = dvec;
    redall(:,q) = rvec;
    yelall(:,q) = yvec;
    greall(:,q) = gvec;
    arrall(:,q) = arrvec;
    radall(:,q) = radii;
    radarrall(:,q) = radarr;
    radnecall(:,q) = radnec;
    Xsnapall(:,(q-1)*(Ndays+1)+1:q*(Ndays+1)) = Xsnap;
    Ysnapall(:,(q-1)*(Ndays+1)+1:q*(Ndays+1)) = Ysnap;
    Zsnapall(:,(q-1)*(Ndays+1)+1:q*(Ndays+1)) = Zsnap;
    statesnapall(:,(q-1)*(Ndays+1)+1:q*(Ndays+1)) = statesnap;
    csnapall(:,(q-1)*(Ndays+1)+1:q*(Ndays+1)) = csnap;
    arrsnapall(:,3*(q-1)*(Ndays+1)+1:3*q*(Ndays+1)) = arrsnap;
    
    % Importing the workspace overrides the previously imported variables. 
    % The last workspace is used as the representative simulation.
    disp(radii(1))
end

%% Recreate mesh for visualisations and recentre data at origin
[Xm,Ym,Zm] = meshgrid(xmesh,ymesh,zmesh);

% -- Recentre at origin -- %
% Data is processed in the IBM with 0 < (x,y,z) < L as this allows for
% easier algorithms. We visualise the data for -L/2 < (x,y,z) < L/2.

% Data in representative simulation
% Agent locations
X = X - L/2;
Y = Y - L/2;
Z = Z - L/2;

% Mesh grid
Xm = Xm - L/2;
Ym = Ym - L/2;
Zm = Zm - L/2;

% Snapshot agent locations
Xsnap = Xsnap - L/2;
Ysnap = Ysnap - L/2;
Zsnap = Zsnap - L/2;
arrsnap = arrsnap - L/2;

% Do same for data in all simulations
Xsnapall = Xsnapall - L/2;
Ysnapall = Ysnapall - L/2;
Zsnapall = Zsnapall - L/2;
arrsnapall = arrsnapall - L/2;

% Indices of X,Y,Z vectors without cells. Convert these to 0 for easier
% analysis later.
non_inds = find(Xsnapall == -2000);
Xsnapall(non_inds) = 0;
Ysnapall(non_inds) = 0;
Zsnapall(non_inds) = 0;
arr_non_inds = find(arrsnapall == -2000);
arrsnapall(arr_non_inds) = 0;

%% Find average data -- living, dead, red, yellow, green, arrested red, and radii (total, arrested, necrotic)
avgN = mean(Nall,2); avgd = mean(dall,2); 
avgred = mean(redall,2); avgyel = mean(yelall,2); avggre = mean(greall,2);
avgarr = mean(arrall,2); r_o_avg = mean(radall,2); r_a_avg = mean(radarrall,2); r_n_avg = mean(radnecall,2);

cycdata = redall - arrall;
avgcyc = mean(cycdata,2);

% Standard deviation in population calculations
stdN = std(Nall,0,2); stdG1 = std(redall,0,2); stdeS = std(yelall,0,2); stdG2 = std(greall,0,2);
stdarr = std(arrall,0,2); stdd = std(dall,0,2); stdcyc = std(cycdata,0,2);

% Standard deviation in radius calculations
stdrad = std(radall,0,2); 
stdradarr = std(radarrall,0,2);
stdradnec = std(radnecall,0,2); 

%% Prepare for representative spheroid visualisations

% Assort agents by state
reds = [];
yels = [];
gres = [];
deads = [];
for i = 1:Nvec(end) + dvec(end)
    if state(i) == 1
        reds = [ reds; X(i) Y(i) Z(i) ];
    elseif state(i) == 2
        yels = [ yels; X(i) Y(i) Z(i) ];
    elseif state(i) == 3
        gres = [ gres; X(i) Y(i) Z(i) ];
    elseif state(i) == 0
        deads = [ deads; X(i) Y(i) Z(i) ];
    end
end

% Set radius of cell representation
rad = 6;

%% Prepare calculations for visualising averaged time-series data

% Calculate error regions
Tg2 = [tg; flipud(tg)]'; % Used for plotting the error region with "fill" 
avgNtot = avgN + avgd; % Average total population
sdNtot = sqrt((std(Nall,0,2)).^2 + (std(dall,0,2)).^2); % Stdev of total population
Nbound1 = avgN + sdNtot; % Upper boundary of fill for N
Nbound2 = avgN - sdNtot; % Lower boundary of fill for N
dbound1 = avgd + std(dall,0,2); % Upper boundary of fill for d
dbound2 = avgd - std(dall,0,2); % Lower boundary of fill for d
arrbound1 = avgarr + std(arrall,0,2); % Upper boundary of fill for arrested
arrbound2 = avgarr - std(arrall,0,2); % Lower boundary of fill for arrested 
rbound1 = avgred + stdG1; % Upper boundary of fill for red
rbound2 = avgred - stdG1; % Lower boundary of fill for red
ybound1 = avgyel + stdeS; % Upper boundary of fill for yellow 
ybound2 = avgyel - stdeS; % Lower boundary of fill for yellow 
gbound1 = avggre + stdG2; % Upper boundary of fill for green 
gbound2 = avggre - stdG2; % Lower boundary of fill for green 
cycbound1 = avgcyc + stdcyc; % Upper boundary of fill for cycling red 
cycbound2 = avgcyc - stdcyc; % Lower boundary of fill for cycling red

% Get geometry of error regions prepared
inbetweenredpop = [rbound1 ; flipud(rbound2)]';
inbetweenyelpop = [ybound1 ; flipud(ybound2)]';
inbetweengrepop = [gbound1 ; flipud(gbound2)]';
inbetweenN = [Nbound1 ; flipud(Nbound2)]';
inbetweend = [dbound1 ; flipud(dbound2)]';
inbetweenarr = [arrbound1 ; flipud(arrbound2)]';
inbetweencyc = [cycbound1 ; flipud(cycbound2)]';

% Radius error marking (area within plus/minus one standard deviation),
% same process as above
bound1 = r_o_avg + stdrad;
bound2 = r_o_avg - stdrad;
bound3 = r_a_avg + stdradarr;
bound4 = r_a_avg - stdradarr;
bound5 = r_n_avg + stdradnec;
bound6 = r_n_avg - stdradnec;
inbetween1 = [bound1; flipud(bound2) ]';
inbetween2 = [bound3; flipud(bound4) ]';
inbetween3 = [bound5; flipud(bound6) ]';

%% Prepare data for day-by-day plotting (separate snapshot outputs into days)

% Initialise storage
snapreds = zeros(Nmax,3*size(statesnap,2));
snapyels = zeros(Nmax,3*size(statesnap,2));
snapgres = zeros(Nmax,3*size(statesnap,2));
snapdeads = zeros(Nmax,3*size(statesnap,2));

% Initialise red/yellow/green/dead cell count
szreds = zeros(1,size(statesnap,2));
szyels = zeros(1,size(statesnap,2));
szgres = zeros(1,size(statesnap,2));
szdeads = zeros(1,size(statesnap,2));

% Extract data and import into storage 
for j = 1:size(statesnap,2)
    dreds = []; % Daily reds
    dyels = []; % Daily yellows
    dgres = []; % Daily greens
    ddeads = []; % Daily deads
    for i = 1:Nsnap(j) % For the population size at this day
        if statesnap(i,j) == 1 % If red
            dreds = [ dreds ; Xsnap(i,j) Ysnap(i,j) Zsnap(i,j) ]; % Add to red storage
        elseif statesnap(i,j) == 2 % If yellow
            dyels = [ dyels ; Xsnap(i,j) Ysnap(i,j) Zsnap(i,j) ]; % Add to yellow storage
        elseif statesnap(i,j) == 3 % If green
            dgres = [ dgres ; Xsnap(i,j) Ysnap(i,j) Zsnap(i,j) ]; % Add to green storage
        elseif statesnap(i,j) == 0 % If dead
            ddeads = [ ddeads ; Xsnap(i,j) Ysnap(i,j) Zsnap(i,j) ]; % Add to dead storage
        end
    end
    % Import into the snapreds storage for the Figure 3 and 4 plot
    % generators
    snapreds(1:size(dreds,1),3*(j-1)+1:3*(j-1)+3) = dreds;
    snapyels(1:size(dyels,1),3*(j-1)+1:3*(j-1)+3) = dyels;
    snapgres(1:size(dgres,1),3*(j-1)+1:3*(j-1)+3) = dgres;
    snapdeads(1:size(ddeads,1),3*(j-1)+1:3*(j-1)+3) = ddeads;
    
    % Keep track of the population count for each type on each day
    szreds(j) = size(dreds,1);
    szyels(j) = size(dyels,1);
    szgres(j) = size(dgres,1);
    szdeads(j) = size(ddeads,1);
end

%% Save daily oxygen profiles as .mat files for Figure 5

% Make a directory in Figure 5 folder to save in
mkdir 'Figure 5'/Prof3D

% Save
for i = 1:Ndays+1
    c_temp = reshape(csnap(:,i),[I I I]);
    filename = ['Figure 5/Prof3D/c_' num2str(i-1) '.mat'];
    save(filename,'c_temp')
end

%% Prepare for density profile calculations

% Please note: this process is relatively computationally expensive

rmax = max(radall(1:(numt-1)/Ndays:end,:),[],2); % Get the largest radius result for all simulations on each day

rbw = 10; % Bin width for distributions

% Initialise with 100 rows to guarantee enough space for all bins

% Radial data averages
radial_data_red = zeros(100,Ndays+1);
radial_data_yel = zeros(100,Ndays+1);
radial_data_gre = zeros(100,Ndays+1);
radial_data_arr = zeros(100,Ndays+1);
radial_data_dead = zeros(100,Ndays+1);
radial_data_cyc = zeros(100,Ndays+1);

% Radial data standard deviations
radial_sd_red = zeros(100,Ndays+1);
radial_sd_yel = zeros(100,Ndays+1);
radial_sd_gre = zeros(100,Ndays+1);
radial_sd_arr = zeros(100,Ndays+1);
radial_sd_dead = zeros(100,Ndays+1);
radial_sd_cyc = zeros(100,Ndays+1);

rho_r = zeros(100,runcount);
rho_y = zeros(100,runcount);
rho_a = zeros(100,runcount);
rho_g = zeros(100,runcount);
scaled_dead = zeros(100,runcount);
rho_cyc = zeros(100,runcount);

for q = 1:Ndays+1
    % Bin EDGES for distance to the periphery (edges are passed into
    % histcounts). Upper bound rounded down to nearest 10 microns to avoid
    % collecting less dense region.
    pbin = 10*floor(rmax(q)/10) - (0:rbw:10*floor(rmax(q)/10));   
    % Set storage for agent data at day q
    dayqX = zeros(Nmax,runcount);
    dayqY = zeros(Nmax,runcount);
    dayqZ = zeros(Nmax,runcount);
    dayqstate = zeros(Nmax,runcount);
    dayqN = zeros(1,runcount);
    
    dayqreds = [];
    dayqyels = [];
    dayqgres = [];
    dayqdeads = [];
    
    % Initialise with 100 rows for enough space for all bins
    radial_day_red = zeros(100,runcount);
    radial_day_yel = zeros(100,runcount);
    radial_day_arr = zeros(100,runcount);
    radial_day_gre = zeros(100,runcount);
    radial_day_dead = zeros(100,runcount);
    radial_day_cyc = zeros(100,runcount);
    
       
    for i = 1:runcount
        % Fill with the data of the i-th run
        dayqX(:,i) = Xsnapall(:,(i-1)*(Ndays+1)+q);
        dayqY(:,i) = Ysnapall(:,(i-1)*(Ndays+1)+q);
        dayqZ(:,i) = Zsnapall(:,(i-1)*(Ndays+1)+q);
        dayqstate(:,i) = statesnapall(:,(i-1)*(Ndays+1)+q);
        dayqN(i) = find(dayqX(:,i),1,'last');
        
        dayq_runi_reds = [];
        dayq_runi_yels = [];
        dayq_runi_gres = [];
        dayq_runi_deads = [];
        
        for j = 1:dayqN(i)
            if dayqstate(j,i) == 1
                dayq_runi_reds = [ dayq_runi_reds ; dayqX(j,i) dayqY(j,i) dayqZ(j,i) ];
            elseif dayqstate(j,i) == 2
                dayq_runi_yels = [ dayq_runi_yels ; dayqX(j,i) dayqY(j,i) dayqZ(j,i) ];
            elseif dayqstate(j,i) == 3
                dayq_runi_gres = [ dayq_runi_gres ; dayqX(j,i) dayqY(j,i) dayqZ(j,i) ];
            elseif dayqstate(j,i) == 0
                dayq_runi_deads = [ dayq_runi_deads ; dayqX(j,i) dayqY(j,i) dayqZ(j,i) ];
            end
        end
        
        % Calculate the distribution 
        dayqruniNarr = find(arrsnapall(:,((i-1)*3*(Ndays+1))+q),1,'last');
        dayqruniradarr = radarrall(TgDN*(q-1)+1,i);
        [redc,yelc,grec,arrc,deadc,cycc] = radcalcs(dayq_runi_reds,dayq_runi_yels,dayq_runi_gres,arrsnapall(1:dayqruniNarr,((i-1)*3*(Ndays+1))+[q q+Ndays+1 q+2*(Ndays+1)]),dayq_runi_deads,pbin,dayqruniradarr);
        dayqrunired_distr = redc;
        dayqruniyel_distr = yelc;
        dayqrunigre_distr = grec;
        dayqruniarr_distr = arrc;
        dayqrunidead_distr = deadc;
        dayqrunicyc_distr = cycc;
        
        radial_day_red(1:length(dayqrunired_distr),i) = dayqrunired_distr;
        radial_day_yel(1:length(dayqruniyel_distr),i) = dayqruniyel_distr;
        radial_day_gre(1:length(dayqrunigre_distr),i) = dayqrunigre_distr;
        radial_day_arr(1:length(dayqruniarr_distr),i) = dayqruniarr_distr;
        radial_day_dead(1:length(dayqrunidead_distr),i) = dayqrunidead_distr;
        radial_day_cyc(1:length(dayqrunicyc_distr),i) = dayqrunicyc_distr;
    end
    
    % Averages. Don't use dead cell data.
    radial_data_red(:,q) = mean(radial_day_red,2);
    radial_data_yel(:,q) = mean(radial_day_yel,2);
    radial_data_gre(:,q) = mean(radial_day_gre,2);
    radial_data_arr(:,q) = mean(radial_day_arr,2);
    radial_data_cyc(:,q) = mean(radial_day_cyc,2);
    
    % Scaled by occupied volume
    rho_r(1:length(pbin)-1,q) = (radial_data_red(1:length(pbin)-1,q))./(4.*pi.*flip(pbin(1:end-1))'.^2*(rbw));
    rho_y(1:length(pbin)-1,q) = (radial_data_yel(1:length(pbin)-1,q))./(4.*pi.*flip(pbin(1:end-1))'.^2*(rbw));
    rho_g(1:length(pbin)-1,q) = (radial_data_gre(1:length(pbin)-1,q))./(4.*pi.*flip(pbin(1:end-1))'.^2*(rbw));
    rho_a(1:length(pbin)-1,q) = (radial_data_arr(1:length(pbin)-1,q))./(4.*pi.*flip(pbin(1:end-1))'.^2*(rbw));
    rho_cyc(1:length(pbin)-1,q) = (radial_data_cyc(1:length(pbin)-1,q))./(4.*pi.*flip(pbin(1:end-1))'.^2*(rbw));
    
    % Remove nans
    rho_r(isnan(rho_r)) = 0;
    rho_y(isnan(rho_y)) = 0;
    rho_g(isnan(rho_g)) = 0;
    rho_a(isnan(rho_a)) = 0;
    rho_cyc(isnan(rho_cyc)) = 0;
    
    % Standard deviations
    radial_sd_red(:,q) = std(radial_day_red,0,2);
    radial_sd_yel(:,q) = std(radial_day_yel,0,2);
    radial_sd_gre(:,q) = std(radial_day_gre,0,2);
    radial_sd_arr(:,q) = std(radial_day_arr,0,2);
    radial_sd_cyc(:,q) = std(radial_day_cyc,0,2);
    
end