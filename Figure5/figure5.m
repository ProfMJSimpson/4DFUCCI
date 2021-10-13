% Clear figures and command window
close all
clc

load('cMap.mat')

for k = 1:Ndays+1
    
    % Find outer, arrest, and necrotic radii at Day k-1
    sphr_rad_k = radii((k-1)*TgDN+1);
    arrest_rad_k = radarr((k-1)*TgDN + 1);
    necro_rad_k = radnec((k-1)*TgDN + 1);
    
    % x and y coordinates for plotting the circle
    theta = 0:0.01:2*pi;
    xcircA = arrest_rad_k*cos(theta);
    ycircA = arrest_rad_k*sin(theta);
    xcircN = necro_rad_k*cos(theta);
    ycircN = necro_rad_k*sin(theta);
    
    % Reshape c on Day k-1 to the right dimensions
    c_k = reshape(csnap(:,k),[Ny Nx Nz]);
    
    % Proportion of domain to include in Figure (+- x,y values)
    fig_domain = 0.1;
    
    % 2D equator profile and 1D profile of c(x,0,0,t)
    c_eq = c_k(:,:,(Nz-1)/2+1);
    c_1D = c_k((Ny-1)/2+1,:,(Nz-1)/2+1);
    
    figure
    set(gcf,'Renderer','Painters')
    surf(Xm(:,:,(Nz-1)/2+1),Ym(:,:,(Nz-1)/2+1),c_eq,'EdgeColor','none')
    view(2)
    shading interp
    lighting phong
    hold on
    % Plot the contour when c = c_a
    [c2,hc2] = contour3(Xm(:,:,(Nz-1)/2+1),Ym(:,:,(Nz-1)/2+1),c_eq,[c_a c_a],'LineWidth',2,'Color','r');
    plot3(c2(1,2:end),c2(2,2:end),0.05*H*ones(1,length(c2(1,:))-1),'LineWidth',2,'Color','r')
    % Plot the outine of the necrotic region
    if necro_rad_k > 0
        plot3(xcircN,ycircN,0.08*H*ones(1,length(xcircN)),'LineWidth',2,'Color',[1 1 1])
    end
    axis equal
    axis([-fig_domain*L fig_domain*L -fig_domain*W fig_domain*W])
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
    c_1Dfine = interp1(Xm((Ny-1)/2+1,:,(Nz-1)/2+1),c_1D,xvec1D,'spline');
    
    figure
    plot(xvec1D,c_1Dfine,'k','LineWidth',2)
    yline(c_a,'r--','LineWidth',2);
    yline(c_d,'c--','LineWidth',2);
    patch([-sphr_rad_k sphr_rad_k sphr_rad_k -sphr_rad_k],[-0.1 -0.1 1.1 1.1],[0.8 0.8 0.4],'EdgeColor','none','FaceAlpha',0.3)
    axis([-fig_domain*L fig_domain*L 0 1])
    
    cfig2 = gcf;
    figurename_1D = ['Nutrient1D\day' num2str(k-1) '_1d.pdf'];
    saveas(cfig2,figurename_1D)
    
end