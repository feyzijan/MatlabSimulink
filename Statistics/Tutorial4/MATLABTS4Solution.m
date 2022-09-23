clear all;
close all;
clc;

% MATLAB EXERCISES IV

%% MVN

d = 3; % number of dimensions
n = 1000; % number of random vectors to simulate

mu = [ -1 0 2 ]';
Sigma = [ 1 0.5 0; 0.5 2 0; 0 0 0.5 ];

z = randn(d, n); % generate standard normal random variates
x = repmat(mu, [1 n]) + chol(Sigma)*z; % convert to random vectors from the desired multivariate normal

figure(1);

for i = 1 : d
for j = 1 : d
   
plotno = (i-1)*d + j;
        
if i == j
subplot(d, d, plotno);
histogram( x(i, :), [-10:0.1:10] );
xlabel(['x_' num2str(i)]);
xlim([-10 10]);
            
else
            
subplot(d, d, plotno);
plot( x(i, :), x(j, :), 'k.' );
xlabel(['x_' num2str(i)]);
ylabel(['x_' num2str(j)]);
xlim([-10 10]);
ylim([-10 10]);
            
end      
end 
end

%% Bootstrapping

mu = 1;
sigma = sqrt(2);
n = 200;

sample = normrnd(mu,sigma, n, 1);

bootsamp1 = randsample(sample, n, true);
bootsamp2 = randsample(sample, n, true);
bootsamp3 = randsample(sample, n, true);

figure(2)
subplot(2,3,1);
histogram(sample);
subplot(2,3,4);
histogram(bootsamp1);
subplot(2,3,5);
histogram(bootsamp2);
subplot(2,3,6);
histogram(bootsamp3);

k = 1000;
% try different values of k
bmeans = zeros(1, k);
for i = 1:k
    bmeans(i) = mean(randsample(sample, n, true));
end

figure(3); clf;

bins = linspace(min(bmeans), max(bmeans), 20);
freq = hist(bmeans, bins); 
class = bins(2) - bins(1);
relfreq = freq/(sum(freq)*class);
normalpdf = normpdf(bins, mu, sigma/sqrt(n));
	hold on;
	bar(bins, relfreq, 1, 'FaceColor', 'w');
	plot(bins, normalpdf, 'b-');
	axis square;
	set(gca, 'Box', 'On', 'FontSize', 8);
xlabel('Mean of resamples')
ylabel('Relative frequency');

% using the bootstrp command instead of explicitly drawing resamples 
% and comparing theoretical and bootstrapped distributions of the mean by density plots :

bmeansnew = bootstrp(k, @mean, sample);

[fi, xi] = ksdensity(bmeansnew);

normalpdf1 = normpdf(xi, mu, sigma/sqrt(n));

figure(4);
hold on;
plot(xi, fi);
plot(xi, normalpdf1,'r-');

%% Bootstrapping Old Faithful waiting times

% Old Faithful dataset
% Initialize variables
T = readtable('faithful.csv');

% Allocate imported array to column variable names
waiting = T{:, 1};
duration = T{:, 2};
day = T{:, 3};

nf = length(waiting);

bootsamp1 = randsample(waiting, nf, true);
bootsamp2 = randsample(waiting, nf, true);
bootsamp3 = randsample(waiting, nf, true);

figure(5)
subplot(2,3,1);
histogram(waiting);
xlim([0 200]);
ylim ([0 100]);
subplot(2,3,4);
histogram(bootsamp1);
xlim([0 200]);
ylim ([0 100]);
subplot(2,3,5);
histogram(bootsamp2);
xlim([0 200]);
ylim ([0 100]);
subplot(2,3,6);
histogram(bootsamp3);
xlim([0 200]);
ylim ([0 100]);

%% number of bootstrap samples
k = 300;
% try different values of k

%% mean

bmeansfaith = zeros(1, k);
for i = 1:k
    bmeansfaith(i) = mean(randsample(waiting, nf, true));
end

% compare with theory (CLT)
% don't know what the true mean and variance are - use sample values

meanf = mean(waiting);
sdf = sqrt(var(waiting));

[fi, xi] = ksdensity(bmeansfaith);

normalpdf2 = normpdf(xi, meanf, sdf/sqrt(nf));
figure(6);
hold on;
plot(xi, fi);
plot(xi, normalpdf2,'r-');
xlabel('Bootstrapping distribution of mean waiting compared with theoretical distribution');

%% median

nf = length(waiting);
k = 200;

% try different values of k

bmediansfaith = zeros(1, k);

for i = 1:k
    bmediansfaith(i) = median(randsample(waiting, nf, true));
end

[fi2, xi2] = ksdensity(bmediansfaith);

figure(7);

plot(xi2, fi2);
xlabel('Bootstrapping distribution of median waiting');

%% 5% trimmed mean 

nf = length(waiting);
k = 200;

% try different values of k

btmeansfaith = zeros(1, k);

for i = 1:k
    btmeansfaith(i) = trimmean(randsample(waiting, nf, true), 10);
end

[fi3, xi3] = ksdensity(btmeansfaith);

figure(8);

plot(xi3, fi3);

xlabel('Bootstrapping distribution of 5% trimmed mean waiting');

% all in one figure

figure(9)
hold on;

plot(xi, fi, 'b-');
plot(xi2, fi2, 'r-');
plot(xi3, fi3, 'g-');

%% Bivariate

n = 1000;

x = zeros(1, n);
y = zeros(1, n);

for i = 1 : n
    
u1 = rand;
yroots = roots( [2 0 3 -5*u1] );
y(i) = yroots(3);
    
u2 = rand;
xroots = roots([ 1 2*y(i)^2 -u2*(2*y(i)^2+1) ]);
x(i) = xroots(2);
    
end

xx = linspace(0, 1, n);
yy = linspace(0, 1, n);
[XX, YY] = meshgrid(xx, yy);
f = (6/5)*( XX + YY.^2 );

figure(10); clf;
hold on;
contour(XX, YY, f, 20);
plot(x, y, 'k.');

xlim([0 1]);
ylim([0 1]);
xlabel('x');
ylabel('y');