% Clear figures and command window
close all
clc

% Make figure storage directory
mkdir FigureOutput

% FIGURE 7 LEFT PANEL
% Plot average living and dead populations with shading of one standard deviation
figure
plot(tg,avgN,'k-','LineWidth',2)
hold on 
hd5 = plot(tg,avgd,'c--','LineWidth',2);
shadeN = fill(Tg2,inbetweenN,'k');
set(shadeN,'FaceAlpha',0.2,'EdgeAlpha',0.2)
shadedeadp = fill(Tg2,inbetweend,'c');
set(shadedeadp,'FaceAlpha',0.2,'EdgeAlpha',0.2)
axis([0 T 0 max(avgN)*1.1]) % Allow for buffer zone at top of figure (1.1*max)
xlabel('Time (hrs)')
ylabel('Population')

cfig_outer = gcf;
figurename_outer = ['FigureOutput\live-dead.pdf'];
saveas(cfig_outer,figurename_outer)

% FIGURE 7 MIDDLE PANEL
% Average red-type cell subpopulations (all, cycling, arrested)
figure
hd1 = plot(tg,avgred,'r:','LineWidth',2);
hold on 
shadered = fill(Tg2,inbetweenredpop,'r');
set(shadered,'FaceAlpha',0.2,'EdgeAlpha',0.2)
hd2 = plot(tg,avgarr,'r--','LineWidth',2);
shadearr = fill(Tg2,inbetweenarr,'r');
set(shadearr,'FaceAlpha',0.2,'EdgeAlpha',0.2)
hd6 = plot(tg,avgcyc,'r-','LineWidth',2);
shadecyc = fill(Tg2,inbetweencyc,'r');
set(shadecyc,'FaceAlpha',0.2,'EdgeAlpha',0.2)
xlabel('Time (hrs)')
ylabel('Population')
axis([0 T 0 max(avgN)*1.1]) % Allow for buffer zone at top of figure (1.1*max)

cfig_reds = gcf;
figurename_reds = ['FigureOutput\reds.pdf'];
saveas(cfig_reds,figurename_reds)

% FIGURE 7 RIGHT PANEL
% Average green and yellow cycling populations
figure
hold on 
plot(tg,avgyel,'-','Color',[1 0.8 0],'LineWidth',2);
shadeyel = fill(Tg2,inbetweenyelpop,'y');
set(shadeyel,'FaceAlpha',0.2,'EdgeAlpha',0.2)
plot(tg,avggre,'g-','LineWidth',2);
shadegre = fill(Tg2,inbetweengrepop,'g');
set(shadegre,'FaceAlpha',0.2,'EdgeAlpha',0.2)
xlabel('Time (hrs)')
ylabel('Population')
axis([0 T 0 max(avggre)*1.1]) % Allow for buffer zone at top of figure (1.1*max)

cfig_yelgre = gcf;
figurename_yelgre = ['FigureOutput\yellow-green.pdf'];
saveas(cfig_yelgre,figurename_yelgre)