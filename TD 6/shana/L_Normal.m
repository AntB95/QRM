function [ logl ] = L_Normal( X0, Xt, sig0, args )
l = length(Xt);

alpha0 = args(1);
alpha1 = args(2);
beta1 = args(3);

sig_t = zeros(l+1,1);
sig_t(1) = sqrt(alpha0 + alpha1*X0^2 + beta1*sig0^2);

logl = 0;
for i = 1:l
    logl = logl - log((1/sig_t(i)) * normpdf(Xt(i) / sig_t(i)));
    sig_t(i+1) = sqrt(alpha0 + alpha1*Xt(i)^2 + beta1*sig_t(i)^2);
end

end

