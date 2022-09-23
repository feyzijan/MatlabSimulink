clear all;
close all;
clc;

% MATLAB EXERCISES I

% OLD FAITHFUL DATA
% 1. Import and boxplot
% Initialize variables.
T = readtable('faithful.csv');

% Allocate imported array to column variable names
waiting = T{:, 1};
duration = T{:, 2};
day = T{:, 3};

% boxplot by day 
figure()
boxplot(waiting, day)
xlabel('Day')
ylabel('Waiting time between successive eruptions (mins)')

% 2. histograms with different number of bins, with kernel density smooth.

figure()
nbins = [ 10 20 50 ]

for bi = 1 : 3
bins = linspace(min(waiting), max(waiting), nbins(bi));

freq = hist(waiting, bins); 
class = bins(2) - bins(1);
relfreq = freq/(sum(freq)*class);

ksestimate = ksdensity(waiting, bins);
% ksestimate = ksestimate/(sum(ksestimate)*class);

subplot(1, 3, bi);
	hold on;
	bar(bins, relfreq, 1, 'FaceColor', 'w');
	
    plot(bins, ksestimate, 'b-');
	axis square;
	set(gca, 'Box', 'On', 'FontSize', 8);
xlabel('Waiting time (mins)')
ylabel('Relative frequency');
end
hold off

% using the "histogram" command
figure()
nbins = [ 10 20 50 ]
for bi = 1 : 3
bins = linspace(min(waiting), max(waiting), nbins(bi));

ksestimate = ksdensity(waiting, bins);

subplot(1, 3, bi);
	hold on;
	histogram(bins, waiting, 1, 'FaceColor', 'w', 'Normalization', 'pdf');
	
    plot(bins, ksestimate, 'b-');
	axis square;
	set(gca, 'Box', 'On', 'FontSize', 8);
xlabel('Waiting time (mins)')
ylabel('Relative frequency');
end
hold off

% Kernel density smooths with different bandwidths
% kernel density evaluated at 100 points across the range of the data.
[f,xi,bw] = ksdensity(waiting)

figure()
plot(xi, f)
xlabel('Waiting time (mins)')
ylabel('density')
hold on

[f,xi] = ksdensity(waiting,'width',10);
plot(xi,f,'--r','LineWidth',1.5);
[f,xi] = ksdensity(waiting,'width',2);
plot(xi,f,'-.k','LineWidth',1.5);
legend('bw = default','bw = 10','bw = 2');
hold off

%3. Plots of successive waiting times by day: subplots and axis limits

% plot daily data 
% clf clear current figure
figure(3); clf;
% need day as numeric

for mi = 1 : 15
    
    loc = find( day == mi );

    subplot(3, 5, mi);
    hold on; % use hold on to allow overlaying of plots
    plot(waiting(loc), 'k-'); % line plot
    set(gca, 'Box', 'On', 'FontSize', 8);
    xlim([1 15]);
    ylim([30 120]);
    if mod(mi, 5) == 1
        ylabel('Waiting time (mins)'); % label y-axis
    end
    if mi > 2*5
        xlabel('Eruption num'); % label x-axis
    end
    title(['Day: ' num2str(mi)]);
    
end

% 4. scatter plot with linear fit
figure(4); clf;
% get number of data points
n = length(duration);
% do linear regression
lagduration = lagmatrix(duration, 1)
Y = waiting;
X = [ ones(n, 1) lagduration ];
B = regress(Y, X);

% fit linear regression
waitingest = B(1) + B(2)*lagduration;

% plot data and line of best fit
figure(4); clf;
hold on; % use hold on to allow overlaying of plots
plot(lagduration, waiting, 'ko'); % plot data
plot(lagduration, waitingest, 'r-'); % plot line of best fit
set(gca, 'Box', 'On', 'FontSize', 8);
xlabel('Previous duration (minutes)'); % label x-axis
ylabel('Waiting time (minutes)'); % label y-axis


%5. k means clustering lagged duration or lagged waiting time and waiting time to next eruption
figure(5); clf;
%lagwaiting = lagmatrix(waiting, 1)

X = [ lagduration waiting ];
K = 2;
C = kmeans(X, K);
col{1} = 'r';
col{2} = 'g';
col{3} = 'b';
col{4} = 'm';

hold on
plot(lagduration, waiting, 'k.')

for c=1:K
loc = find( C == c );
plot(lagduration(loc), waiting(loc), 'color', col{c}, 'marker', 'o', 'markersize', 4, 'LineStyle', 'None');
end
axis square;
set(gca, 'Box', 'On', 'FontSize', 8);
xlabel('last duration (minutes)')
ylabel('waiting time (minutes)')