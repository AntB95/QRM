%-------------------- QRM Assignment 10 -------------------
%========== Maxime Richiardi, Ivan Schönenberger ==========
%========== Antoine Bedanian, Nicolas de Lestable==========
%----------------------------------------------------------

clc;
clear all;
close all;


%Set parameters
N = 1000;
S = 10000;
Pi = 0.01;
Rho = 0.01;

Corr = (1 - Rho) * eye(N) + Rho * ones(N);


X_norm = copularnd('Gaussian',Corr,S);
X_t3 = copularnd('t',Corr,3,S);
X_t10 = copularnd('t',Corr,10,S);

threshold = @(x) length(find(x<=Pi));
L_norm = zeros(S,1);
L_t3 = zeros(S,1);
L_t10 = zeros(S,1);
Corr_norm = zeros(N,1);
Corr_t3 = zeros(N-1,1);
Corr_t10 = zeros(N-1,1);
for i = 1:S
    L_norm(i) = threshold(X_norm(i,:));
    L_t3(i) = threshold(X_t3(i,:));
    L_t10(i) = threshold(X_t10(i,:));
    if i < N
        Bin = [X_norm(:,i:i+1) <= Pi  X_t3(:,i:i+1) <= Pi X_t10(:,i:i+1) <= Pi ];
        Corr_Bin = corr(Bin);
        Corr_norm(i) = Corr_Bin(2,1) ;
        Corr_t3(i) = Corr_Bin(4,3) ;
        Corr_t10(i) = Corr_Bin(6,5);
        if i == 1
            fprintf('The correlation of default occurence between the 2 first columns are :\n');
            fprintf('Gaussian copula: %2.4f\n',Corr_norm(1));
            fprintf('t3 copula: %2.4f\n',Corr_t3(1));
            fprintf('t10 copula: %2.4f\n',Corr_t10(1));
        end
    end
end

x = 0:1:50;
subplot(1,3,1)
histogram(L_norm,'Normalization','probability');
axis([0 50 0 1])
hold on
    line(x,binopdf(0:1:50,N,Pi),'Color','red','MarkerSize',15);
hold off
title('Default repartition for Gaussian copula');
ylabel('Probability');
xlabel('Number of default');
subplot(1,3,2)
histogram(L_t3,100,'Normalization','probability');
axis([0 50 0 1])
hold on
    line(x,binopdf(0:1:50,N,Pi),'Color','red','MarkerSize',15);
hold off
title('Default repartition for Student 3 copula');
ylabel('Probability');
xlabel('Number of default');
subplot(1,3,3)
histogram(L_t10,100,'Normalization','probability');
axis([0 50 0 1])
hold on
    line(x,binopdf(0:1:50,N,Pi),'Color','red','MarkerSize',15);
hold off
title('Default repartition for Student 10 copula');
ylabel('Probability');
xlabel('Number of default');
set(gcf,'position',[10,10,1550,1200])
saveas(gcf,'Default_Repartition.png');
L_norm = sort(L_norm);
L_t3 = sort(L_t3);
L_t10 = sort(L_t10);
Q_norm = [L_norm(S*0.95) L_norm(S*0.99) L_norm(S*0.999)];
Q_t3 = [L_t3(S*0.95) L_t3(S*0.99) L_t3(S*0.999)];
Q_t10 = [L_t10(S*0.95) L_t10(S*0.99) L_t10(S*0.999)];


subplot(1,3,1)
histogram(Corr_norm,30,'Normalization','probability');
title('Default correlation for Gaussian copula');
ylabel('Number of occurence');
xlabel('Correlation');
subplot(1,3,2)
histogram(Corr_t3,30,'Normalization','probability');
title('Default correlation for Student 3 copula');
ylabel('Number of occurence');
xlabel('Correlation');
subplot(1,3,3)
histogram(Corr_t10,30,'Normalization','probability');
title('Default correlation for Student 10 copula');
ylabel('Number of occurence');
xlabel('Correlation');
set(gcf,'position',[10,10,1550,1200])
saveas(gcf,'Correlation.png');


