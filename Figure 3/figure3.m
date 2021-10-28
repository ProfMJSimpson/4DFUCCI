% Clear figures and command window
close all
clc

Diam = 10; % Slightly smaller than cell size -- staining only occurs for the nucleus of the cell

% Create file directories for images
mkdir Equators
mkdir LowerTropics
mkdir UpperTropics

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
    figure
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
    set(gca,'TickLabelInterpreter', 'LaTeX','Fontsize' , 12,'Color',[0 0 0],'GridColor',[0.5 0.5 0.5])
    % Remove/rename tick labels
    ax.XTick = [-400 -200 0 200 400];
    ax.XTickLabel = {"","","","",""};
    ax.YTick = [-400 -200 0 200 400];
    ax.YTickLabel = {"","","","",""};
    
    cfig = gcf;
    set(gcf,'Color',[1 1 1])
    set(gcf,'InvertHardCopy','off');
    set(gcf,'renderer','Painters')
    figname1 = ['Equators\day' num2str(k-1) '.pdf'];
    saveas(cfig,figname1)
    
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
    
    figure
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
    set(gca,'TickLabelInterpreter', 'LaTeX','Fontsize' , 12,'Color',[0 0 0],'GridColor',[0.5 0.5 0.5])
    % Remove/rename tick labels
    ax.XTick = [-400 -200 0 200 400];
    ax.XTickLabel = {"","","","",""};
    ax.YTick = [-400 -200 0 200 400];
    ax.YTickLabel = {"","","","",""};
    
    cfig = gcf;
    set(gcf,'Color',[1 1 1])
    set(gcf,'InvertHardCopy','off');
    set(gcf,'renderer','Painters')
    figname2 = ['LowerTropics\day' num2str(k-1) '.pdf'];
    saveas(cfig,figname2)
    
    figure
    axis equal
    axis([-0.1*L 0.1*L -0.1*L 0.1*L -0.1*L 0.1*L])
    hold on
    if ~isempty(redk)
        scatter3(redk(rinds_ut,1),redk(rinds_ut,2),redk(rinds_ut,3),marker_size,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0.8 0 0])
    end
    if ~isempty(yels)
        scatter3(yelk(yinds_ut,1),yelk(yinds_ut,2),yelk(yinds_ut,3),marker_size,'MarkerFaceColor',[1 1 0],'MarkerEdgeColor',[0.8 0.8 0])
    end
    if ~isempty(gres)
        scatter3(grek(ginds_ut,1),grek(ginds_ut,2),grek(ginds_ut,3),marker_size,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 0.8 0])
    end
    view(2)
    grid on
    set(gca,'TickLabelInterpreter', 'LaTeX','Fontsize' , 12,'Color',[0 0 0],'GridColor',[0.5 0.5 0.5])
    % Remove/rename tick labels
    ax.XTick = [-400 -200 0 200 400];
    ax.XTickLabel = {"","","","",""};
    ax.YTick = [-400 -200 0 200 400];
    ax.YTickLabel = {"","","","",""};
    
    cfig = gcf;
    set(gcf,'Color',[1 1 1])
    set(gcf,'InvertHardCopy','off');
    set(gcf,'renderer','Painters')
    figname3 = ['UpperTropics\day' num2str(k-1) '.pdf'];
    saveas(cfig,figname3)
    
end