% Clear figures and command window
close all
clc

load('cMap.mat')

% Make directories in which to store figures
mkdir Nutrient1D
mkdir Nutrient2D
mkdir Nutrient3D

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
    
    figure
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
    axis equal
    axis([-fig_domain*L fig_domain*L -fig_domain*L fig_domain*L])
    ax = gca;
    ax.XTick = [0];
    ax.XTickLabel = {"0"};
    ax.YTick = [0];
    ax.YTickLabel = {"0"};
    colorbar 
    colormap(ax,cMap)
    caxis([0 1])
    
    cfig = gcf;
    figurename_2D = ['Nutrient2D\day' num2str(k-1) '_2d.pdf'];
    saveas(cfig,figurename_2D)
    
    xvec1D = linspace(-L/2,L/2,201);
    % Perform spline interpolation if necessary to make profile smoother
    c_1Dfine = interp1(Xm((I-1)/2+1,:,(I-1)/2+1),c_1D,xvec1D,'spline');
    
    figure
    plot(xvec1D,c_1Dfine,'k','LineWidth',2)
    yline(c_a,'r--','LineWidth',2);
    yline(c_d,'c--','LineWidth',2);
    patch([-r_o_k r_o_k r_o_k -r_o_k],[-0.1 -0.1 1.1 1.1],[0.8 0.8 0.4],'EdgeColor','none','FaceAlpha',0.3)
    axis([-fig_domain*L fig_domain*L 0 1])
    
    cfig2 = gcf;
    figurename_1D = ['Nutrient1D\day' num2str(k-1) '_1d.pdf'];
    saveas(cfig2,figurename_1D)
    
end