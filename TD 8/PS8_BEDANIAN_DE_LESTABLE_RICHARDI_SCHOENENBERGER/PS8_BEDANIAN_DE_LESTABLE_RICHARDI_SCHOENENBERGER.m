%-------------------- QRM Assignment 8 --------------------
%========== Maxime Richiardi, Ivan Schönenberger ==========
%========== Antoine Bedanian, Nicolas de Lestable==========
%----------------------------------------------------------

clc;
clear all;
close all;

% Load Data
data = csvread('Global_3_Factors_Daily2.csv');
x = data(:,2);
y = data(:,3);

% Get ECDF
[qx , vx] = ecdf(x);
[qy , vy] = ecdf(y);

x_quant = zeros(length(data),1);
y_quant = zeros(length(data),1);
%Get quantile of data
for i=1:length(data)
    x_quant(i) = max(qx(vx <= x(i))) - eps;
    y_quant(i) = max(qy(vy <= y(i))) - eps;
end

% Proceed to QMLE on Gumbel,Clayton and Frank distribution
data_quant = [x_quant , y_quant];
gumbel = @(alpha) - sum(log(copulapdf('Gumbel',data_quant,alpha)));
clayton = @(alpha) - sum(log(copulapdf('Clayton',data_quant,alpha)));
frank = @(alpha) - sum(log(copulapdf('Frank',data_quant,alpha)));

x0 = 2; A = []; b = []; Aeq =[]; beq = [];
lb = [1+eps,0+eps,-Inf];
ub = Inf;
[gumbel_par , gumbel_ML] = fmincon(gumbel,x0,A,b,Aeq,beq,lb(1),ub);
[clayton_par , clayton_ML] = fmincon(clayton,x0,A,b,Aeq,beq,lb(2),ub);
[frank_par , frank_ML] = fmincon(frank,x0,A,b,Aeq,beq,lb(3),ub);

fig = figure(1);
scatter(x_quant,y_quant,'.','r');
title('Scatter plot of pseudo sample empirical distribution');
xlabel('X');
ylabel('Y');
saveas(fig,'QRM_PS8.png');

fprintf('|-------------------------------------------|\n')
fprintf('| Distribution | Parameter value | ML value |\n')
fprintf('|    Clayton   |     %2.2s    | %2.2s |\n',clayton_par,-clayton_ML)
fprintf('|    Gumbel    |     %2.2s    | %2.2s|\n',gumbel_par,-gumbel_ML)
fprintf('|     Frank    |     %2.2s   | %2.2s |\n',frank_par,-frank_ML)
fprintf('|-------------------------------------------|\n')



