function [sig_t1] = sigmat1(args,X0, Xt,sig0)
l = length(Xt);
alpha0 = args(1);
alpha1 = args(2);
beta1 = args(3);

sig_t = zeros(l+1,1);

sig_t(1) = sqrt(alpha0 + alpha1*X0^2 + beta1*sig0^2);

for i = 1:l
    sig_t(i+1) = sqrt(alpha0 + alpha1*Xt(i)^2 + beta1*sig_t(i)^2);
end

sig_t1 = sig_t(end);
end