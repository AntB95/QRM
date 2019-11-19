function [ logl ] = L_Student(X0, Xt, sigma0, args)
l = length(Xt);
alpha0 = args(1);
alpha1 = args(2);
beta1 = args(3);

v = args(4);
scale = sqrt(v/(v - 2));

sig_t = zeros(l+1,1);
sig_t(1) = sqrt(alpha0 + alpha1*X0^2 + beta1*sigma0^2);

logl = 0;
for i = 1:l
    logl = logl - log((1/sig_t(i)) * scale * tpdf(scale * Xt(i) / sig_t(i), v));
    sig_t(i+1) = sqrt(alpha0 + alpha1*Xt(i)^2 + beta1*sig_t(i)^2);
end

end