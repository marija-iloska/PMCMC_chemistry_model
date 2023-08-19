function [Y] = pertpdf(x, a, b, c)

r = c - a;

alpha = 1 + 4*(b - a)./r;
beta0 =  1 + 4*(c - b)./r;

top = (x - a).^(alpha - 1).* (c - x).^(beta0 - 1);
bottom = beta(alpha, beta0).* (c - a).^(alpha + beta0 + 1);

Y = top./bottom;

end