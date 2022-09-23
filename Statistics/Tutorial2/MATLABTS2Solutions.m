%MATLAB EXERCISES II

%Simulation

n = 100;
type = randi(10, [1 n]);
test = randi(10, [1 n]);
athtype = type > 1; %0 = cheat, 1 = clean
testcor = test > 2; %0 = incorrect, 1 = correct
cheatfail = (athtype == 0) & (testcor == 1); %1 = the athlete cheated and the test is correct, 0 = otherwise
cleanfail = (athtype == 1) & (testcor == 0); %1 = the athlete is clean and the test is incorrect, 0 = otherwise
ncheatfail = sum(cheatfail);
ncleanfail = sum(cleanfail);
nfail = ncheatfail + ncleanfail;
condprob = ncheatfail / nfail;

%Repeat this simulation many times (for example, 1000 times)

n = 100;
for count = 1:1000
type = randi(10, [1 n]);
test = randi(10, [1 n]);
athtype = type > 1; %0 = cheat, 1 = clean
testcor = test > 2; %0 = incorrect, 1 = correct
cheatfail = (athtype == 0) & (testcor == 1); %1 = the athlete cheated and the test is correct, 0 = otherwise
cleanfail = (athtype == 1) & (testcor == 0); %1 = the athlete is clean and the test is incorrect, 0 = otherwise
ncheatfail = sum(cheatfail);
ncleanfail = sum(cleanfail);
nfail = ncheatfail + ncleanfail;
condprob = ncheatfail / nfail;
condp(count) = condprob;
end
hist(condp)

%Poisson distribution

lambda = 4;
x = 0:20;
p = poisspdf(x, lambda);
F = poisscdf(x, lambda);

figure(1); clf;

subplot(1, 2, 1);

bar(x, p, 0.05);
axis square;
xlabel('x');
ylabel('p(x)');

subplot(1, 2, 2);

stairs(x, F);
axis square;
xlabel('x');
ylabel('F(x)');

% Generate random numbers from Poisson(4)

n = 100;
rn = poissrnd(lambda, [1 n]);

samplemean = mean(rn);
samplevar= var(rn);
 
% pmfs
nx = hist(rn, x);
nx = nx./(sum(nx));

figure(2); clf;
hold on;
bar(x, p, 0.05);

plot(x, nx, 'r*');

% cdfs
ecdf = cumsum(nx);

figure(1); clf;
hold on;
stairs(x, F);

plot(x, ecdf, 'r*');
axis([0 15 0 1]);


% Normal distribution

mu = 0; sigma2 = 1; nx = 100;
x = linspace(-10, 10, nx);
f = normpdf(x, mu, sqrt(sigma2));
F = normcdf(x, mu, sqrt(sigma2));

% Plotting pdf and cdf
figure(1); clf;
subplot(1, 2, 1);
plot(x, f, 'k-'); axis square; xlabel('x');
ylabel('f(x)');
subplot(1, 2, 2);
plot(x, F, 'k-'); axis square; xlabel('x');
ylabel('F(x)');

% Generate random numbers from N(0,1)
n = 100;
rn = normrnd(mu, sqrt(sigma2), [1 n]);

samplemean = mean(rn);
samplevariance = var(rn);

% Plotting histogram (and overlay densities)

nx = hist(rn, x);
dx = x(2)-x(1);
nx = nx./(sum(nx)*dx);
hold on;

bar(x, nx, 1, 'FaceColor', 'w');
plot(x, f, 'r-', 'LineWidth', 2);