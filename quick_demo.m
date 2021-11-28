%----------------------- EXAMPLE RUN -----------------------%
% Run this file to do a shorter, single simulation run, for the purposes of
% demonstrating the functionality of the code, and example outputs.

% Runtime for quick_demo.m is approximately 378 seconds (6.3 minutes)

%-% SET THE PARAMETERS FOR THE MODEL

% Simulation choice: 1 for short, 6 day simulation. 2 for longer, 10 day
% simulation.
choice = 1;

% All same initial conditions
exp_N0 = 30000;
exp_rad = 245;

% Acts as a seed for this demonstration and nothing more. In the main IBM,
% sim_id also establishes the initial population and radius.
sim_id = 1; 

N0 = exp_N0; % Initial population
Nmax = 300000; % Maximum number of cells for simulation
endtimes = [ 144 240 ]; % Final times of simulation (hours)
T = endtimes(choice); % Simulation end point (h)
L = 4000; % Length of domain
I = 51; % Number of x, y, or z nodes
xmesh = linspace(0,L,I); % Generate x mesh
ymesh = linspace(0,L,I); % Generate y mesh
zmesh = linspace(0,L,I); % Generate z mesh

sigma = 12; % Proliferation spread parameter

alpha = 0.15; % Consumption-diffusion ratio
h = xmesh(2) - xmesh(1); % Node spacing parameter
c_b = 1; % Boundary nutrient concentration

r = exp_rad; % Initial radius

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
TgDN = 24/dtsave; % Number of time points in a day (e.g. if dtsave = 2 hours, TgDN = 12)

runcount = 1; % Keep this at 1, 

Ndays = T/24; % Number of days represented by T (hours)

% Run script (call ibm3d.m function)
[Nvec,dvec,rvec,yvec,gvec,arrvec,radii,radarr,radnec,X,Y,Z,state,c,c_p,Xsnap,Ysnap,Zsnap,arrsnap,statesnap,csnap,Nsnap] = ibm3d(N0,Nmax,T,L,sigma,parms,rates,tg,xmesh,ymesh,zmesh,sim_id);

% Send data to "all" storage, even though only one run. 
Nall = Nvec;
dall = dvec;
redall = rvec;
yelall = yvec;
greall = gvec;
arrall = arrvec;
radall = radii;
radarrall = radarr;
radnecall = radnec;
Xsnapall = Xsnap;
Ysnapall = Ysnap;
Zsnapall = Zsnap;
statesnapall = statesnap;
csnapall = csnap;
arrsnapall = arrsnap;

%% Data analysis and preparation for figure generation
%-% Recreate mesh for visualisations and recentre data at origin
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

%-% Find average data -- living, dead, red, yellow, green, arrested red, and radii (total, arrested, necrotic)
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

%-% Prepare data for plotting daily information (separate snapshot outputs into days)
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

%-% Prepare to plot concentration profiles
rmax = max(radall(1:(numt-1)/Ndays:end,:),[],2);
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
    pbin = 10*floor(rmax(q)/10) - (0:rbw:10*floor(rmax(q)/10)); % Equation (S13)
    % Set storage for agent data at day q
    dayqX = zeros(Nmax,runcount);
    dayqY = zeros(Nmax,runcount);
    dayqZ = zeros(Nmax,runcount);
    dayqstate = zeros(Nmax,runcount);
    dayqN = zeros(1,runcount);
    
    dayqreds = [];
    dayqyels = [];
    dayqgres = [];
    
    % Initialise with 100 rows for enough space for all bins
    radial_day_red = zeros(100,runcount);
    radial_day_yel = zeros(100,runcount);
    radial_day_arr = zeros(100,runcount);
    radial_day_gre = zeros(100,runcount);
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
        
        for j = 1:dayqN(i)
            if dayqstate(j,i) == 1
                dayq_runi_reds = [ dayq_runi_reds ; dayqX(j,i) dayqY(j,i) dayqZ(j,i) ];
            elseif dayqstate(j,i) == 2
                dayq_runi_yels = [ dayq_runi_yels ; dayqX(j,i) dayqY(j,i) dayqZ(j,i) ];
            elseif dayqstate(j,i) == 3
                dayq_runi_gres = [ dayq_runi_gres ; dayqX(j,i) dayqY(j,i) dayqZ(j,i) ];
            end
        end
        
        % Calculate the distribution 
        dayqruniNarr = find(arrsnapall(:,((i-1)*3*(Ndays+1))+q),1,'last'); % Number of arrested agents at t
        dayqruniradarr = radarrall(TgDN*(q-1)+1,i); % Arrested radius at t
        [redc,yelc,grec,arrc,cycc] = radcalcs(dayq_runi_reds,dayq_runi_yels,dayq_runi_gres,arrsnapall(1:dayqruniNarr,((i-1)*3*(Ndays+1))+[q q+Ndays+1 q+2*(Ndays+1)]),pbin,dayqruniradarr);
        dayqrunired_distr = redc;
        dayqruniyel_distr = yelc;
        dayqrunigre_distr = grec;
        dayqruniarr_distr = arrc;
        dayqrunicyc_distr = cycc;
        
        radial_day_red(1:length(dayqrunired_distr),i) = dayqrunired_distr;
        radial_day_yel(1:length(dayqruniyel_distr),i) = dayqruniyel_distr;
        radial_day_gre(1:length(dayqrunigre_distr),i) = dayqrunigre_distr;
        radial_day_arr(1:length(dayqruniarr_distr),i) = dayqruniarr_distr;
        radial_day_cyc(1:length(dayqrunicyc_distr),i) = dayqrunicyc_distr;
    end
    
    % Averages. Don't use dead cell data.
    radial_data_red(:,q) = mean(radial_day_red,2);
    radial_data_yel(:,q) = mean(radial_day_yel,2);
    radial_data_gre(:,q) = mean(radial_day_gre,2);
    radial_data_arr(:,q) = mean(radial_day_arr,2);
    radial_data_cyc(:,q) = mean(radial_day_cyc,2);
    
    % Scaled by occupied volume -- Equation (S16)
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

%% Plot figures

Diam = 10; % Diameter slightly smaller than cell size, as in experimental images, we only see nucleus

%-% Equator, upper, and lower tropics without dead agents (Fig 3)
figure
set(gcf,'Position',[20 20 1800 900])
for k = 1:Ndays+1
    Narr = find(arrsnap(:,k)+L/2,1,'last'); % Number of arrested cells
    
    % Reds, yellows, and greens at Day k-1
    redk = snapreds(1:szreds(k),3*(k-1)+1:3*(k-1)+3);
    yelk = snapyels(1:szyels(k),3*(k-1)+1:3*(k-1)+3);
    grek = snapgres(1:szgres(k),3*(k-1)+1:3*(k-1)+3);
        
    % Total and necrotic radius at k-1 days
    radP = radii(TgDN*(k-1)+1);
    radN = radnec(TgDN*(k-1)+1);
        
    % Indices of agents at equator
    rinds = find(abs(redk(:,3) - 0) < 6);
    yinds = find(abs(yelk(:,3) - 0) < 6);
    ginds = find(abs(grek(:,3) - 0) < 6);
    
    % Plot 2D equator slice
    subplot(3,Ndays+1,k)
    axis equal
    axis([-0.1*L 0.1*L -0.1*L 0.1*L -0.1*L 0.1*L])
    % Set marker size to agent size (12 micron diameter)
    ax = gca;
    AR = get(gca, 'dataaspectratio');
    oldunits = get(ax,'Units');
    set(ax,'Units','points');
    pos = get(ax,'Position');
    set(ax,'Units',oldunits');
    XL = xlim(ax);
    points_per_unit = Diam*pos(3)/(XL(2) - XL(1));
    marker_size = points_per_unit.^2*pi/4;
    hold on
    if ~isempty(redk)
        scatter3(redk(rinds,1),redk(rinds,2),redk(rinds,3),marker_size,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0.8 0 0])
    end
    if ~isempty(yelk)
        scatter3(yelk(yinds,1),yelk(yinds,2),yelk(yinds,3),marker_size,'MarkerFaceColor',[1 1 0],'MarkerEdgeColor',[0.8 0.8 0])
    end
    if ~isempty(grek)
        scatter3(grek(ginds,1),grek(ginds,2),grek(ginds,3),marker_size,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 0.8 0])
    end
    view(2)
    grid on
    title(['Day ', num2str(k-1),' Equator'])
    set(gca,'TickLabelInterpreter', 'LaTeX','Fontsize' , 12,'Color',[0 0 0],'GridColor',[0.5 0.5 0.5])
    % Remove/rename tick labels
    ax.XTick = [-400 -200 0 200 400];
    ax.XTickLabel = {"","","","",""};
    ax.YTick = [-400 -200 0 200 400];
    ax.YTickLabel = {"","","","",""};
    
    % Indices for relevant upper/lower tropics
    if radN < 50 % If necrotic core is not significant in size 
        rinds_ut = find(abs(redk(:,3) + radP/2) < 6);
        yinds_ut = find(abs(yelk(:,3) + radP/2) < 6);
        ginds_ut = find(abs(grek(:,3) + radP/2) < 6);
        rinds_lt = find(abs(redk(:,3) - radP/2) < 6);
        yinds_lt = find(abs(yelk(:,3) - radP/2) < 6);
        ginds_lt = find(abs(grek(:,3) - radP/2) < 6);
    else
        rinds_ut = find(abs(redk(:,3) + radN) < 6);
        yinds_ut = find(abs(yelk(:,3) + radN) < 6);
        ginds_ut = find(abs(grek(:,3) + radN) < 6);
        rinds_lt = find(abs(redk(:,3) - radN) < 6);
        yinds_lt = find(abs(yelk(:,3) - radN) < 6);
        ginds_lt = find(abs(grek(:,3) - radN) < 6);
    end
    
    subplot(3,Ndays+1,Ndays+1+k)
    axis equal
    axis([-0.1*L 0.1*L -0.1*L 0.1*L -0.1*L 0.1*L])
    hold on
    if ~isempty(redk)
        scatter3(redk(rinds_lt,1),redk(rinds_lt,2),redk(rinds_lt,3),marker_size,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0.8 0 0])
    end
    if ~isempty(yelk)
        scatter3(yelk(yinds_lt,1),yelk(yinds_lt,2),yelk(yinds_lt,3),marker_size,'MarkerFaceColor',[1 1 0],'MarkerEdgeColor',[0.8 0.8 0])
    end
    if ~isempty(grek)
        scatter3(grek(ginds_lt,1),grek(ginds_lt,2),grek(ginds_lt,3),marker_size,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 0.8 0])
    end
    view(2)
    grid on
    title(['Day ', num2str(k-1),' Lower Tropic'])
    set(gca,'TickLabelInterpreter', 'LaTeX','Fontsize' , 12,'Color',[0 0 0],'GridColor',[0.5 0.5 0.5])
    % Remove/rename tick labels
    ax = gca;
    ax.XTick = [-400 -200 0 200 400];
    ax.XTickLabel = {"","","","",""};
    ax.YTick = [-400 -200 0 200 400];
    ax.YTickLabel = {"","","","",""};
    
    subplot(3,Ndays+1,2*(Ndays+1)+k)
    axis equal
    axis([-0.1*L 0.1*L -0.1*L 0.1*L -0.1*L 0.1*L])
    hold on
    if ~isempty(redk)
        scatter3(redk(rinds_ut,1),redk(rinds_ut,2),redk(rinds_ut,3),marker_size,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0.8 0 0])
    end
    if ~isempty(yelk)
        scatter3(yelk(yinds_ut,1),yelk(yinds_ut,2),yelk(yinds_ut,3),marker_size,'MarkerFaceColor',[1 1 0],'MarkerEdgeColor',[0.8 0.8 0])
    end
    if ~isempty(grek)
        scatter3(grek(ginds_ut,1),grek(ginds_ut,2),grek(ginds_ut,3),marker_size,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 0.8 0])
    end
    view(2)
    grid on
    title(['Day ', num2str(k-1),' Upper Tropic'])
    set(gca,'TickLabelInterpreter', 'LaTeX','Fontsize' , 12,'Color',[0 0 0],'GridColor',[0.5 0.5 0.5])
    % Remove/rename tick labels
    ax = gca;
    ax.XTick = [-400 -200 0 200 400];
    ax.XTickLabel = {"","","","",""};
    ax.YTick = [-400 -200 0 200 400];
    ax.YTickLabel = {"","","","",""};
    
end

%-% 3D and 2D spheroid images with dead agents
figure
set(gcf,'Position',[20 20 1800 900])
for k = 1:Ndays+1
    Narr = find(arrsnap(:,k)+L/2,1,'last'); % Number of arrested cells
    
    % Reds, yellows, and greens at Day k-1
    redk = snapreds(1:szreds(k),3*(k-1)+1:3*(k-1)+3);
    yelk = snapyels(1:szyels(k),3*(k-1)+1:3*(k-1)+3);
    grek = snapgres(1:szgres(k),3*(k-1)+1:3*(k-1)+3);
    deadk = snapdeads(1:szdeads(k),3*(k-1)+1:3*(k-1)+3);
    
    % Find the indices of the octant to remove
    rXinds = find(redk(:,1) < 0);
    rYinds = find(redk(:,2) < 0);
    rZinds = find(redk(:,3) > 0);
    rinds_oct = intersect(intersect(rXinds,rYinds),rZinds);

    yXinds = find(yelk(:,1) < 0);
    yYinds = find(yelk(:,2) < 0);
    yZinds = find(yelk(:,3) > 0);
    yinds_oct = intersect(intersect(yXinds,yYinds),yZinds);

    gXinds = find(grek(:,1) < 0);
    gYinds = find(grek(:,2) < 0);
    gZinds = find(grek(:,3) > 0);
    ginds_oct = intersect(intersect(gXinds,gYinds),gZinds);

    dXinds = find(deadk(:,1) < 0);
    dYinds = find(deadk(:,2) < 0);
    dZinds = find(deadk(:,3) > 0);
    dinds_oct = intersect(intersect(dXinds,dYinds),dZinds);
        
    % Indices of agents at equator
    rinds_eq = find(abs(redk(:,3) - 0) < 6);
    yinds_eq = find(abs(yelk(:,3) - 0) < 6);
    ginds_eq = find(abs(grek(:,3) - 0) < 6);
    dinds_eq = find(abs(deadk(:,3) - 0) < 6);
        
    redsSP = redk;
    redsSP(rinds_oct,:) = [];
    yelsSP = yelk;
    yelsSP(yinds_oct,:) = [];
    gresSP = grek;
    gresSP(ginds_oct,:) = [];
    deadsSP = deadk;
    deadsSP(dinds_oct,:) = [];
    
    subplot(3,Ndays+1,k)
    % Set an appropriate agent size for 3D
    view(3)
    axis equal
    axis([-0.1*L 0.1*L -0.1*L 0.1*L -0.1*L 0.1*L])
    ax = gca;
    AR = get(gca, 'dataaspectratio');
    oldunits = get(ax,'Units');
    set(ax,'Units','points');
    pos = get(ax,'Position');
    set(ax,'Units',oldunits');
    XL = xlim(ax);
    points_per_unit = Diam*pos(3)/(XL(2) - XL(1));
    marker_size = points_per_unit.^2*pi/4;
    hold on
    if ~isempty(redsSP)
        scatter3(redsSP(:,1),redsSP(:,2),redsSP(:,3),marker_size,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0.8 0 0])
    end
    if ~isempty(yelsSP)
        scatter3(yelsSP(:,1),yelsSP(:,2),yelsSP(:,3),marker_size,'MarkerFaceColor',[1 1 0],'MarkerEdgeColor',[0.8 0.8 0])
    end
    if ~isempty(gresSP)
        scatter3(gresSP(:,1),gresSP(:,2),gresSP(:,3),marker_size,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 0.8 0])
    end
    if ~isempty(deadsSP)
        scatter3(deadsSP(:,1),deadsSP(:,2),deadsSP(:,3),marker_size,'MarkerFaceColor',[0 1 1],'MarkerEdgeColor',[0 0.8 0.8])
    end
    grid on
    set(gca,'TickLabelInterpreter', 'LaTeX','Fontsize' , 12,'Color',[0 0 0],'GridColor',[0.7 0.7 0.7])
    % Plot black lines to clearly identify the removed octant
    plot3([-10,-0.92*max(abs(redk(rinds_eq,1)))],[-10 -10],[10 10],'k','Linewidth',2)
    plot3([-10 -10],[-10,-0.92*max(abs(redk(rinds_eq,1)))],[10 10],'k','Linewidth',2)
    plot3([-10 -10],[-10 -10],[10,0.92*max(abs(redk(rinds_eq,1)))],'k','Linewidth',2)
    th = pi/2:0.01:pi;
    xy = (0.92*max(abs(redk(rinds_eq,1)))-10) .* cos(th) - 10;
    zh = (0.92*max(abs(redk(rinds_eq,1)))-10) .* sin(th) + 10;
    plot3(-10*ones(length(xy)),xy,zh,'k','Linewidth',2)
    plot3(xy,-10*ones(length(xy)),zh,'k','Linewidth',2)
    plot3(xy,fliplr(xy),10*ones(length(xy)),'k','Linewidth',2)
    title(['Day ', num2str(k-1),' 3D Spheroid'])
    ax.GridAlpha = 0.25;
    % Remove/edit axis ticks
    ax.XTick = [-400 -200 0 200 400];
    ax.XTickLabel = {"","","","",""};
    ax.YTick = [-400 -200 0 200 400];
    ax.YTickLabel = {"","","","",""};
    ax.ZTick = [-400 -200 0 200 400];
    ax.ZTickLabel = {"","","","",""};
    
    subplot(3,Ndays+1,Ndays+1+k)
    axis equal
    axis([-0.1*L 0.1*L -0.1*L 0.1*L -0.1*L 0.1*L])
    % Set agent size as appropriate for 2D figure
    ax = gca;
    AR = get(gca, 'dataaspectratio');
    oldunits = get(ax,'Units');
    set(ax,'Units','points');
    pos = get(ax,'Position');
    set(ax,'Units',oldunits');
    XL = xlim(ax);
    points_per_unit = Diam*pos(3)/(XL(2) - XL(1));
    marker_size = points_per_unit.^2*pi/4;
    hold on
    if ~isempty(redk)
        scatter3(redk(rinds_eq,1),redk(rinds_eq,2),redk(rinds_eq,3),marker_size,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0.8 0 0])
    end
    if ~isempty(yelk)
        scatter3(yelk(yinds_eq,1),yelk(yinds_eq,2),yelk(yinds_eq,3),marker_size,'MarkerFaceColor',[1 1 0],'MarkerEdgeColor',[0.8 0.8 0])
    end
    if ~isempty(grek)
        scatter3(grek(ginds_eq,1),grek(ginds_eq,2),grek(ginds_eq,3),marker_size,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 0.8 0])
    end
    if ~isempty(deadk)
        scatter3(deadk(dinds_eq,1),deadk(dinds_eq,2),deadk(dinds_eq,3),marker_size,'MarkerFaceColor',[0 1 1],'MarkerEdgeColor',[0 0.8 0.8])
    end
    view(2)
    title(['Day ', num2str(k-1),' 2D Slice'])
    grid on
    set(gca,'TickLabelInterpreter', 'LaTeX','Fontsize' , 12,'Color',[0 0 0],'GridColor',[0.7 0.7 0.7])
    % Remove/edit axis ticks
    ax.XTick = [-400 -200 0 200 400];
    ax.XTickLabel = {"","","","",""};
    ax.YTick = [-400 -200 0 200 400];
    ax.YTickLabel = {"","","","",""};
    
end

%-% Agent concentration profiles
for q = 1:Ndays+1
    % Bin EDGES for distance to the periphery (edges are passed into
    % histcounts). Minus one in first term to
    pbin = 10*floor(rmax(q)/10) - (0:rbw:10*floor(rmax(q)/10)); % Equation (S13)
    
    % Scaling for the normalisation
    rho_max = max(rho_cyc,[],'all');
    
    % Calculate 1D nutrient profile 
    cdq = reshape(csnap(:,q),[I I I]);
    c1D = cdq((I-1)/2+1,(I-1)/2+1:(I-1)/2+12,(I-1)/2+1); % +12 should be enough
    xvs = Xm((I-1)/2+1,(I-1)/2+1:(I-1)/2+12,(I-1)/2+1);
    cshow = interp1(xvs,c1D,fliplr(pbin),'spline');
    
    % Cleaner without error bars
    subplot(3,Ndays+1,2*(Ndays+1)+q)
    hold on
    % Shift pbin to account for difference between bin counts (one fewer
    % than bin edges) -- Scaling with Equation (S17)
    hd1 = plot(pbin(2:end),(rho_cyc(1:length(pbin)-1,q))./(rho_max),'r-','LineWidth',2);
    hd2 = plot(pbin(2:end),(rho_y(1:length(pbin)-1,q))./(rho_max),'-','LineWidth',2,'Color',[1 0.8 0]);
    hd3 = plot(pbin(2:end),(rho_g(1:length(pbin)-1,q))./(rho_max),'g-','LineWidth',2);
    hd4 = plot(pbin(4:end),(rho_a(3:length(pbin)-1,q))./(rho_max),'r--','LineWidth',2);
    title(['Day ', num2str(q-1),' Density Profile'])
    axis([0 1.1*max(rmax) 0 1]);
    plot(pbin,cshow,'k-','LineWidth',2)
    hold off
    
end

%-% 2D and 1D nutrient profiles
load('cMap.mat')
figure
set(gcf,'Position',[20 20 1800 600])
for k = 1:Ndays+1
    
    % Find outer, arrest, and necrotic radii at Day k-1
    r_o_k = radii((k-1)*TgDN+1);
    r_a_k = radarr((k-1)*TgDN + 1);
    r_n_k = radnec((k-1)*TgDN + 1);
    
    % x and y coordinates for plotting the circle
    theta = 0:0.01:2*pi;
    xcircA = r_a_k*cos(theta);
    ycircA = r_a_k*sin(theta);
    xcircN = r_n_k*cos(theta);
    ycircN = r_n_k*sin(theta);
    
    % Reshape c on Day k-1 to the right dimensions
    c_k = reshape(csnap(:,k),[I I I]);
    
    % Proportion of domain to include in Figure (+- x,y values)
    fig_domain = 0.1;
    
    % 2D equator profile and 1D profile of c(x,0,0,t)
    c_eq = c_k(:,:,(I-1)/2+1);
    c_1D = c_k((I-1)/2+1,:,(I-1)/2+1);
    
    subplot(2,Ndays+1,k)
    set(gcf,'Renderer','Painters')
    surf(Xm(:,:,(I-1)/2+1),Ym(:,:,(I-1)/2+1),c_eq,'EdgeColor','none')
    view(2)
    shading interp
    lighting phong
    hold on
    % Plot the contour when c = c_a
    [c2,hc2] = contour3(Xm(:,:,(I-1)/2+1),Ym(:,:,(I-1)/2+1),c_eq,[c_a c_a],'LineWidth',2,'Color','r');
    plot3(c2(1,2:end),c2(2,2:end),0.05*L*ones(1,length(c2(1,:))-1),'LineWidth',2,'Color','r')
    % Plot the outine of the necrotic region
    if r_n_k > 0
        plot3(xcircN,ycircN,0.08*L*ones(1,length(xcircN)),'LineWidth',2,'Color',[1 1 1])
    end
    title(['Day ', num2str(k-1),' Nutrient 2D'])
    axis equal
    axis([-fig_domain*L fig_domain*L -fig_domain*L fig_domain*L])
    ax = gca;
    ax.XTick = [0];
    ax.XTickLabel = {"0"};
    ax.YTick = [0];
    ax.YTickLabel = {"0"};
    h = colorbar; 
    colormap(ax,cMap)
    caxis([0 1])
    set(h, 'Position', [.915 .625 .01 .260])
    
    xvec1D = linspace(-L/2,L/2,201);
    % Perform spline interpolation if necessary to make profile smoother
    c_1Dfine = interp1(Xm((I-1)/2+1,:,(I-1)/2+1),c_1D,xvec1D,'spline');
    
    subplot(2,Ndays+1,Ndays+1+k)
    plot(xvec1D,c_1Dfine,'k','LineWidth',2)
    yline(c_a,'r--','LineWidth',2);
    yline(c_d,'c--','LineWidth',2);
    patch([-r_o_k r_o_k r_o_k -r_o_k],[-0.1 -0.1 1.1 1.1],[0.8 0.8 0.4],'EdgeColor','none','FaceAlpha',0.3)
    title(['Day ', num2str(k-1),' Nutrient 1D'])
    axis([-fig_domain*L fig_domain*L 0 1])
    
end

%-% Population plots
figure
set(gcf,'Position',[20 20 1800 500])
% FIGURE 7 LEFT PANEL
% Plot average living and dead populations with shading of one standard deviation
subplot(1,3,1)
plot(tg,avgN,'k-','LineWidth',2)
hold on 
hd5 = plot(tg,avgd,'c--','LineWidth',2);
axis([0 T 0 max(avgN)*1.1]) % Allow for buffer zone at top of figure (1.1*max)
xlabel('Time (hrs)')
ylabel('Population')
title('Living and dead populations')

% FIGURE 7 MIDDLE PANEL
% Average red-type cell subpopulations (all, cycling, arrested)
subplot(1,3,2)
hd1 = plot(tg,avgred,'r:','LineWidth',2);
hold on 
hd2 = plot(tg,avgarr,'r--','LineWidth',2);
hd6 = plot(tg,avgcyc,'r-','LineWidth',2);
xlabel('Time (hrs)')
ylabel('Population')
axis([0 T 0 max(avgN)*1.1]) % Allow for buffer zone at top of figure (1.1*max)
title('Red subpopulation classes')

% FIGURE 7 RIGHT PANEL
% Average green and yellow cycling populations
subplot(1,3,3)
hold on 
plot(tg,avgyel,'-','Color',[1 0.8 0],'LineWidth',2);
plot(tg,avggre,'g-','LineWidth',2);
xlabel('Time (hrs)')
ylabel('Population')
axis([0 T 0 max(avggre)*1.1]) % Allow for buffer zone at top of figure (1.1*max)
title('Green and yellow subpopulations')

%-% Radius comparison plot
% Import experimental data and prepare for use
Days = [0 1 2 3 4 5 6 8 10];

outers = [];
inhibs = [];
necros = [];

for i = 1:length(Days)
    dy = Days(i);
    
    filepath = ['RadiusData/day' num2str(dy) '.csv'];
    Tab = readtable(filepath);
    time = Tab.Day;
    time = 24*(time - 4); % Account for offset in experimental image labels (formation time). Time is in days here.
    outer = Tab.OuterRadius;
    inhib = Tab.InhibitedRadius;
    necr = Tab.NecroticRadius;
    
    outers = [outers ; time outer];
    inhibs = [inhibs ; time inhib];
    necros = [necros ; time necr];
end

% Read IncuCyte data
incu = readtable('RadiusData/IncucyteData.csv');

% Extract the data for the 793b cell line with 10 000 cell initial density
C793b10k = incu{strcmp(incu{:, 'CellLine'}, '793b') & incu{:, 'InitialCondition'} == 10000, [4,11]};

% Find initial radius estimate
day4s = C793b10k(C793b10k(:,1)==4,:);
% Initial average radius and stddev
avgrad0 = mean(day4s(:,2));
sdrad = std(day4s(:,2));

incudata = [];

for i = 1:Ndays/0.25 + 1
    tv = 3.75 + 0.25*i;
    hr = C793b10k(C793b10k(:,1)==tv,:);
    hr(:,1) = 24*(hr(:,1) - 4);
    
    incudata = [incudata ; hr];
end

fullouter = [ outers ; incudata ]; % All outer radius information

% Calculate average outer radius data for all unique time points.
unique_outer = unique(fullouter(:,1));
outer_rad_avg = zeros(length(unique_outer),1);
outer_rad_sd = zeros(length(unique_outer),1);
for j = 1:length(unique_outer)
    inds = find(fullouter(:,1) == unique_outer(j));
    data = fullouter(inds,:);
    radii_avg = mean(data(:,2));
    radii_sd = std(data(:,2));
    outer_rad_avg(j) = radii_avg;
    outer_rad_sd(j) = radii_sd;
end

% Plot with errorbars for outer radius measurements (larger sample size)
% and scatter for inner radii (arrested + necrotic)
figure
% Experimental results
errorbar(unique_outer,outer_rad_avg,outer_rad_sd,'s','MarkerSize',5,'MarkerFaceColor','Black','Color','Black')
axis([0 T 0 1.1*max(outer_rad_avg)])
hold on
in1 = scatter(inhibs(:,1),inhibs(:,2),25,'r','filled');
ne1 = scatter(necros(:,1),necros(:,2),25,[0 0.8 0.8],'filled');
in1.MarkerFaceAlpha = 0.4;
ne1.MarkerFaceAlpha = 0.4;
% Model results
plot(tg,r_o_avg,'k-','LineWidth',2)
plot(tg,r_a_avg,'r--','LineWidth',2)
plot(tg,r_n_avg,'c--','LineWidth',2)
xlabel('Time (days)')
ylabel('Radius')
ax = gca;
ax.XTick = 0:24:T;
ax.XTickLabels = {"0","1","2","3","4","5","6","7","8","9","10"};
title('Radius comparison')