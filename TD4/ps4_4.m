rng(1); 
N = 10000; 
A = [ 1 0 0 0;
    1 1 0 0;
    -1 2 3 0;
    1 -1 1 1]; 
x = trnd(5,N,4); 
X = (A*x')';
alpha = 0.95;

L = sum(X,2);
VaR = quantile(L,alpha)
%VaR = 10.41

%b) PCA
Q = cov(X);
[V,D] = eig(Q);
first_PCA = V(:,4);
second_PCA = V(:,3);
mu = mean(X);
Y = (V'*((X-mu)'))';
X_approx = zeros(N,4);
for i=1:N
    X_approx(i,:) = mu + Y(i,4)*first_PCA' + Y(i,3)*second_PCA';
end

%we compute the approx Var using the first 2 EV associate to the biggest
%Eigen values
L_approx = sum(X_approx,2);
VaR_approx = quantile(L_approx,alpha)
