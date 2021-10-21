% Clear figures and command window
close all
clc

Diam = 10;

% Make directories to store figures
mkdir Images2D
mkdir Images3D
mkdir ConcentrationProfiles

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
    
    figure
    % Set an appropriate agent size for 3D
    view(3)
    axis equal
    axis([-0.1*L 0.1*L -0.1*W 0.1*W -0.1*H 0.1*H])
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
    ax.GridAlpha = 0.25;
    % Remove/edit axis ticks
    ax.XTick = [-400 -200 0 200 400];
    ax.XTickLabel = {"","","","",""};
    ax.YTick = [-400 -200 0 200 400];
    ax.YTickLabel = {"","","","",""};
    ax.ZTick = [-400 -200 0 200 400];
    ax.ZTickLabel = {"","","","",""};
    
    cfig = gcf;
    set(gcf,'Color',[1 1 1])
    set(gcf,'InvertHardCopy','off');
    figname = ['Images3D\day' num2str(k-1) '.pdf'];
    saveas(cfig,figname)
    
    figure
    axis equal
    axis([-0.1*L 0.1*L -0.1*W 0.1*W -0.1*H 0.1*H])
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
    grid on
    set(gca,'TickLabelInterpreter', 'LaTeX','Fontsize' , 12,'Color',[0 0 0],'GridColor',[0.7 0.7 0.7])
    % Remove/edit axis ticks
    ax.XTick = [-400 -200 0 200 400];
    ax.XTickLabel = {"","","","",""};
    ax.YTick = [-400 -200 0 200 400];
    ax.YTickLabel = {"","","","",""};
    
    cfig = gcf;
    set(gcf,'Color',[1 1 1])
    set(gcf,'InvertHardCopy','off');
    set(gcf,'renderer','Painters')
    figname2 = ['Images2D\day' num2str(k-1) '.pdf'];
    saveas(cfig,figname2)
    
end

% Concentration profiles
%% Nutrient and cell concentration profiles over time (distance from periphery) 

for q = 1:Ndays+1
    % Bin EDGES for distance to the periphery (edges are passed into
    % histcounts). Minus one in first term to
    pbin = 10*floor(rmax(q)/10) - (0:rbw:10*floor(rmax(q)/10)); 
    
    % Error bounds
    boundred1 = radial_data_red(1:length(pbin)-1,q) + radial_sd_red(1:length(pbin)-1,q);
    boundred2 = radial_data_red(1:length(pbin)-1,q) - radial_sd_red(1:length(pbin)-1,q);
    boundyel1 = radial_data_yel(1:length(pbin)-1,q) + radial_sd_yel(1:length(pbin)-1,q);
    boundyel2 = radial_data_yel(1:length(pbin)-1,q) - radial_sd_yel(1:length(pbin)-1,q);
    boundgre1 = radial_data_gre(1:length(pbin)-1,q) + radial_sd_gre(1:length(pbin)-1,q);
    boundgre2 = radial_data_gre(1:length(pbin)-1,q) - radial_sd_gre(1:length(pbin)-1,q);
    boundarr1 = radial_data_arr(3:length(pbin)-1,q) + radial_sd_arr(3:length(pbin)-1,q);
    boundarr2 = radial_data_arr(3:length(pbin)-1,q) - radial_sd_arr(3:length(pbin)-1,q);
    boundcycred1 = radial_data_cyc(1:length(pbin)-1,q) + radial_sd_cyc(1:length(pbin)-1,q);
    boundcycred2 = radial_data_cyc(1:length(pbin)-1,q) - radial_sd_cyc(1:length(pbin)-1,q);
    
    % Geometry for the "fill" function 
    perip = [pbin(2:end) fliplr(pbin(2:end))];
    
    % Scaling for the normalisation
    rho_max = max(rho_cyc,[],'all');
    
    % Scale error bounds by radial volume and scaling factor
    % Shift pbin to account for difference between bin counts (one fewer
    % than bin edges)
    scaledred1 = boundred1./(4.*pi.*flip(pbin(1:end-1))'.^2*(rbw)*rho_max);
    scaledred2 = boundred2./(4.*pi.*flip(pbin(1:end-1))'.^2*(rbw)*rho_max);
    scaledyel1 = boundyel1./(4.*pi.*flip(pbin(1:end-1))'.^2*(rbw)*rho_max);
    scaledyel2 = boundyel2./(4.*pi.*flip(pbin(1:end-1))'.^2*(rbw)*rho_max);
    scaledgre1 = boundgre1./(4.*pi.*flip(pbin(1:end-1))'.^2*(rbw)*rho_max);
    scaledgre2 = boundgre2./(4.*pi.*flip(pbin(1:end-1))'.^2*(rbw)*rho_max);
    scaledarr1 = boundarr1./(4.*pi.*flip(pbin(1:end-3))'.^2*(rbw)*rho_max);
    scaledarr2 = boundarr2./(4.*pi.*flip(pbin(1:end-3))'.^2*(rbw)*rho_max);  
    scaleddiff1 = boundcycred1./(4.*pi.*flip(pbin(1:end-1))'.^2*(rbw)*rho_max);
    scaleddiff2 = boundcycred2./(4.*pi.*flip(pbin(1:end-1))'.^2*(rbw)*rho_max);
    
    inbetweenred = [scaledred1; flipud(scaledred2) ]';
    inbetweenred(isnan(inbetweenred)) = 0;
    inbetweenyel = [scaledyel1; flipud(scaledyel2) ]';
    inbetweenyel(isnan(inbetweenyel)) = 0;
    inbetweengre = [scaledgre1; flipud(scaledgre2) ]';
    inbetweengre(isnan(inbetweengre)) = 0;
    inbetweenarrdens = [scaledarr1; flipud(scaledarr2) ]';
    inbetweenarrdens(isnan(inbetweenarrdens)) = 0;
    inbetweencycred = [scaleddiff1; flipud(scaleddiff2) ]';
    inbetweencycred(isnan(inbetweencycred)) = 0;
    
    % Calculate 1D nutrient profile 
    cdq = reshape(csnap(:,q),[I I I]);
    c1D = cdq((I-1)/2+1,(I-1)/2+1:(I-1)/2+12,(I-1)/2+1); % +12 should be enough
    xvs = Xm((I-1)/2+1,(I-1)/2+1:(I-1)/2+12,(I-1)/2+1);
    cshow = interp1(xvs,c1D,fliplr(pbin),'spline');
    
    % Cleaner without error bars
    figure
    hold on
    % Shift pbin to account for difference between bin counts (one fewer
    % than bin edges)
    hd1 = plot(pbin(2:end),(rho_cyc(1:length(pbin)-1,q))./(rho_max),'r-','LineWidth',2);
    shadecyc = fill(perip,inbetweencycred,'r');
    set(shadecyc,'FaceAlpha',0.2,'EdgeAlpha',0.2)
    hd2 = plot(pbin(2:end),(rho_y(1:length(pbin)-1,q))./(rho_max),'-','LineWidth',2,'Color',[1 0.8 0]);
    shadeyel = fill(perip,inbetweenyel,[1 0.8 0]);
    set(shadeyel,'FaceAlpha',0.2,'EdgeAlpha',0.2)
    hd3 = plot(pbin(2:end),(rho_g(1:length(pbin)-1,q))./(rho_max),'g-','LineWidth',2);
    shadegre = fill(perip,inbetweengre,'g');
    set(shadegre,'FaceAlpha',0.2,'EdgeAlpha',0.2)
    hd4 = plot(pbin(4:end),(rho_a(3:length(pbin)-1,q))./(rho_max),'r--','LineWidth',2);
    shadearr = fill(perip(3:end-2),inbetweenarrdens,'r-');
    set(shadearr,'FaceAlpha',0.2,'EdgeAlpha',0.2)
    axis([0 1.1*max(rmax) 0 1]);
    plot(pbin,cshow,'k-','LineWidth',2)
    hold off
    
    cfig = gcf;
    figurename = ['ConcentrationProfiles\day' num2str(q-1) '.pdf'];
    saveas(cfig,figurename)
    
end