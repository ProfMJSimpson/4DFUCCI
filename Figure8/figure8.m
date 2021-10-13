% Clear figures and command window
close all
clc

% Import experimental data and prepare for use
Days = [0 1 2 3 4 5 6 8 10];

outers = [];
inhibs = [];
necros = [];

for i = 1:length(Days)
    dy = Days(i);
    
    filepath = ['..\RadiusData\day' num2str(dy) '.csv'];
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
incu = readtable('..\RadiusData\IncucyteData.csv');

% Extract the data for the 793b cell line with 10 000 cell initial density
C793b10k = incu{strcmp(incu{:, 'CellLine'}, '793b') & incu{:, 'InitialCondition'} == 10000, [4,11]};

% Find initial radius estimate
day4s = C793b10k(C793b10k(:,1)==4,:);
% Initial average radius and stddev
avgrad0 = mean(day4s(:,2));
sdrad = std(day4s(:,2));

Ndays = 10;

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
errorbar(unique_outer,outer_rad_avg,outer_rad_sd,'s','MarkerSize',5,'MarkerFaceColor','Black','Color','Black')
axis([0 endtime 0 1.1*max(outer_rad_avg)])
hold on
in1 = scatter(inhibs(:,1),inhibs(:,2),25,'r','filled');
ne1 = scatter(necros(:,1),necros(:,2),25,[0 0.8 0.8],'filled');
in1.MarkerFaceAlpha = 0.4;
ne1.MarkerFaceAlpha = 0.4;
plot(Tg,avgrad,'k-','LineWidth',2)
shadered = fill(Tg2,inbetween1,'k');
set(shadered,'FaceAlpha',0.2,'EdgeAlpha',0.2)
plot(Tg,avgradarr,'r--','LineWidth',2)
shade2 = fill(Tg2,inbetween2,'r');
set(shade2,'FaceAlpha',0.2,'EdgeAlpha',0.2)
plot(Tg,avgradnec,'c--','LineWidth',2)
shade3 = fill(Tg2,inbetween3,[0 1 1]);
set(shade3,'FaceAlpha',0.2,'EdgeAlpha',0.2)
xlabel('Time (days)')
ylabel('Radius')
ax = gca;
ax.XTick = 0:24:endtime;
ax.XTickLabels = {"0","1","2","3","4","5","6","7","8","9","10"};

cfig_rad = gcf;
figurename_rad = ['RadiusPlot\radii.pdf'];
saveas(cfig_rad,figurename_rad)