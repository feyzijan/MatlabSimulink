clear all;
close all;
clc;

% MATLAB EXERCISES III

%% Inverse Transform Sampling (Exponential)

n = 100; % number of samples to take - try different values n
lambda = 0.5;

u = rand(1, n); % vector of n uniformly distributed random numbers
x = zeros(1, n); % vector to hold random numbers from inverse transform
y = zeros(1, n); % vector to hold random numbers from inverse transform using built in function

for i = 1 : n
    x(i) = log((1-u(i))^(-1/lambda));
end

for i = 1 : n
    y(i) = expinv(u(i),1/lambda);
end

x_true = makedist('Exponential', 'mu', 1/lambda); %make probability distribution

figure(1);

qqplot(x, x_true);
axis square;
xlabel('Theoretical');
ylabel('Sample');
title('Exponential');

figure(2);

qqplot(y, x_true);
axis square;
xlabel('Theoretical');
ylabel('Sample');
title('Exponential');

%% Inverse Transform Sampling (Poisson)

n = 100; % number of samples to take - try different values n
lambda = 500;

u = rand(1, n); % vector of n uniformly distributed random numbers
x = zeros(1, n); % vector to hold  random numbers from inverse transform

for i = 1 : n
    x(i) = poissinv(u(i),lambda);
end

x_true = makedist('Poisson', 'lambda', lambda);

figure(3);

qqplot(x, x_true);
axis square;
xlabel('Theoretical');
ylabel('Sample');
title('Poisson');

%% Box Muller Sampling
% try different values of n for all these

n = 1000;
u1 = rand(1, n);
u2 = rand(1, n);

x = sqrt(-2*log(u1)).*cos(2*pi*u2);

figure(4);
qqplot(x);
axis square;
title('Normal');