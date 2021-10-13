% Clear figures and command window
close all
clc

Diam = 10;

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
    peripD = 10*ceil(maxsphrads(q)/10 - 1) - (0:rbw:10*ceil(maxsphrads(q)/10)); % Distance to the periphery
    
    % Error bounds (calculated just in case wanted, but clutters the plot)
    boundred1 = radial_data_red(1:length(peripD)-1,q) + radial_sd_red(1:length(peripD)-1,q);
    boundred2 = radial_data_red(1:length(peripD)-1,q) - radial_sd_red(1:length(peripD)-1,q);
    boundyel1 = radial_data_yel(1:length(peripD)-1,q) + radial_sd_yel(1:length(peripD)-1,q);
    boundyel2 = radial_data_yel(1:length(peripD)-1,q) - radial_sd_yel(1:length(peripD)-1,q);
    boundgre1 = radial_data_gre(1:length(peripD)-1,q) + radial_sd_gre(1:length(peripD)-1,q);
    boundgre2 = radial_data_gre(1:length(peripD)-1,q) - radial_sd_gre(1:length(peripD)-1,q);
    boundarr1 = radial_data_arr(3:length(peripD)-1,q) + radial_sd_arr(3:length(peripD)-1,q);
    boundarr2 = radial_data_arr(3:length(peripD)-1,q) - radial_sd_arr(3:length(peripD)-1,q);
    boundcycred1 = radial_data_cyc(1:length(peripD)-1,q) + radial_sd_cyc(1:length(peripD)-1,q);
    boundcycred2 = radial_data_cyc(1:length(peripD)-1,q) - radial_sd_cyc(1:length(peripD)-1,q);
    
    % Geometry for the "fill" function 
    perip2 = [peripD(1:end-1) fliplr(peripD(1:end-1))];
    
    % Scaling for the normalisation
    max_scale = max(scaled_cyc,[],'all');
    
    scaledred1 = boundred1./(4.*pi.*flip(peripD(1:end-1))'.^2*(rbw)*max_scale);
    scaledred2 = boundred2./(4.*pi.*flip(peripD(1:end-1))'.^2*(rbw)*max_scale);
    scaledyel1 = boundyel1./(4.*pi.*flip(peripD(1:end-1))'.^2*(rbw)*max_scale);
    scaledyel2 = boundyel2./(4.*pi.*flip(peripD(1:end-1))'.^2*(rbw)*max_scale);
    scaledgre1 = boundgre1./(4.*pi.*flip(peripD(1:end-1))'.^2*(rbw)*max_scale);
    scaledgre2 = boundgre2./(4.*pi.*flip(peripD(1:end-1))'.^2*(rbw)*max_scale);
    scaledarr1 = boundarr1./(4.*pi.*flip(peripD(1:end-3))'.^2*(rbw)*max_scale);
    scaledarr2 = boundarr2./(4.*pi.*flip(peripD(1:end-3))'.^2*(rbw)*max_scale);  
    scaleddiff1 = boundcycred1./(4.*pi.*flip(peripD(1:end-1))'.^2*(rbw)*max_scale);
    scaleddiff2 = boundcycred2./(4.*pi.*flip(peripD(1:end-1))'.^2*(rbw)*max_scale);
    
    inbetweenred = [scaledred1; flipud(scaledred2) ]';
    inbetweenred(isnan(inbetweenred)) = 0;
    inbetweenyel = [scaledyel1; flipud(scaledyel2) ]';
    inbetweenyel(isnan(inbetweenyel)) = 0;
    inbetweengre = [scaledgre1; flipud(scaledgre2) ]';
    inbetweengre(isnan(inbetweengre)) = 0;
    inbetweenarr = [scaledarr1; flipud(scaledarr2) ]';
    inbetweenarr(isnan(inbetweenarr)) = 0;
    inbetweencycred = [scaleddiff1; flipud(scaleddiff2) ]';
    inbetweencycred(isnan(inbetweencycred)) = 0;
    
    % Calculate 1D nutrient profile 
    cdq = reshape(csnap(:,q),[Ny Nx Nz]);
    c1D = cdq((Ny-1)/2+1,(Nx-1)/2+1:(Nx-1)/2+12,(Nz-1)/2+1); % +12 should be enough
    xvs = Xm((Ny-1)/2+1,(Nx-1)/2+1:(Nx-1)/2+12,(Nz-1)/2+1);
    cshow = interp1(xvs,c1D,fliplr(peripD),'spline');
    
    % Cleaner without error bars
    figure
    hold on
    hd1 = plot(peripD(1:end-1),(scaled_cyc(1:length(peripD)-1,q))./(max_scale),'r-','LineWidth',2);
    shadecyc = fill(perip2,inbetweencycred,'r');
    set(shadecyc,'FaceAlpha',0.2,'EdgeAlpha',0.2)
    hd2 = plot(peripD(1:end-1),(scaled_yel(1:length(peripD)-1,q))./(max_scale),'-','LineWidth',2,'Color',[1 0.8 0]);
    shadeyel = fill(perip2,inbetweenyel,[1 0.8 0]);
    set(shadeyel,'FaceAlpha',0.2,'EdgeAlpha',0.2)
    hd3 = plot(peripD(1:end-1),(scaled_gre(1:length(peripD)-1,q))./(max_scale),'g-','LineWidth',2);
    shadegre = fill(perip2,inbetweengre,'g');
    set(shadegre,'FaceAlpha',0.2,'EdgeAlpha',0.2)
    hd4 = plot(peripD(3:end-1),(scaled_arr(3:length(peripD)-1,q))./(max_scale),'r--','LineWidth',2);
    shadearr = fill(perip2(3:end-2),inbetweenarr,'r-');
    set(shadearr,'FaceAlpha',0.2,'EdgeAlpha',0.2)
    axis([0 1.1*max(maxsphrads) 0 1]);
    plot(peripD,cshow,'k-','LineWidth',2)
    hold off
    
    cfig = gcf;
    figurename = ['ConcentrationProfiles\day' num2str(q-1) '.pdf'];
    saveas(cfig,figurename)
    
end