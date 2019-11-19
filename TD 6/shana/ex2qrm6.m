clear all
clc;
close all;

%import data
data = readtable('dataQRM72.xlsx');
prices = data.PriceOrBid_AskAverage;

s = find(strcmp(data.NamesDate,"08/04/2013")); % Find index of start date

x0n = [0.3 0.6 0.5];
x0t = [0.3 0.6 0.5 5];
lbn = [10^-10,10^-10,10^-10];
lbt = [10^-10,10^-10,10^-10,100];

%log return 
R = log(prices(2:end)./prices(1:end-1));
n = length(prices);

sigt = zeros(n-s,1);
mlen = zeros(n-s,3);
mlet = zeros(n-s,4);

% Alphas for VaR(alpha)
a1 = 0.95;
a2 = 0.99;

q1n = norminv(a1);
q2n = norminv(a2);

varn = zeros(n-s,2);
vart = zeros(n-s,2);

% For each day between 04/2013 and 04/2016
for t = s:n-1
    t
    
    Rt = R(t-s+1:t);
    Xt = Rt - mean(Rt);
    sigt(t-s+1) = std(Xt); % Empirical stdv for each t
    mu = mean(Rt);
    
    % (1)
    mlen(t-s+1,:) = fmincon(@(args) L_Normal(0, Xt, sigt(t-s+1), args), x0n,[],[],[],[],lbn,[]);
    
    sigmat1n = sigmat1(mlen(t-s+1,:),0,Xt,sigt(t-s+1));
    % (2) Value at risk normal innovations
    varn(t-s+1,1) = sigmat1n*q1n - mu;
    varn(t-s+1,2) = sigmat1n*q2n - mu;
    
    % (4)
    mlet(t-s+1,:) = fmincon(@(args) L_Student(0, Xt, sigt(t-s+1), args), x0t,[],[],[],[],lbt,[]);
    
    sigmat1t = sigmat1(mlet(t-s+1,1:3),0,Xt,sigt(t-s+1));
    nu = mlet(t-s+1,4);
    c = sqrt((nu - 2)/nu);
    % (4) Value at risk student t innovations
    vart(t-s+1,1) = sigmat1t*c*tinv(a1,nu) - mu;
    vart(t-s+1,2) = sigmat1t*c*tinv(a2,nu) - mu;
end
%% (3) BREACHES
% Negative returns
neg_rets = -R(s+1:n-1);

% PLOTS of vars
figure;
plot(neg_rets)
hold on;
plot(varn(1:end-1,1))
plot(varn(1:end-1,2))
legend('Negative log returns','Var(0.95)', 'Var(0.99)')

figure;
plot(neg_rets)
hold on;
plot(vart(1:end-1,1))
plot(vart(1:end-1,2))
legend('Negative log returns','Var(0.95)', 'Var(0.99)')

% NUMBER of breaches
% Normal
length(find(neg_rets > varn(1:end-1,1)))
length(find(neg_rets > varn(1:end-1,2)))

% Student
length(find(neg_rets > vart(1:end-1,1)))
length(find(neg_rets > vart(1:end-1,2)))
%% (5) VARIANCES
var_mlen = zeros(n-s,1);
var_mlet = zeros(n-s,1);

for i = 1:n-s
    % MLE Zt NORMAL
    if 1 - mlen(i,2) - mlen(i,3) > 0 % stationarity condition
        var_mlen(i) = mlen(i,1) / (1 - mlen(i,2) - mlen(i,3));
    end
    % MLE Zt STUDENT t
    if 1 - mlet(i,2) - mlet(i,3) > 0 % stationarity condition
        var_mlet(i) = mlet(i,1) / (1 - mlet(i,2) - mlet(i,3));
    end
end

% PLOTS MLE variance
figure;
plot(var_mlen,'r') % MLE normal var
hold on;
plot(var_mlet,'g') % MLE student var
plot(sigt.^2) % empirical var
legend('Normal Zt','Student Zt', 'Empirical', 'Location','southeast')
