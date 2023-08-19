function [Y] = pertpdf(x, a, mu, c)

r = c - a;

b = (6*mu - a - c)/4;

alpha = 1 + 4*(b - a)./r;
beta0 =  1 + 4*(c - b)./r;

top = (x - a).^(alpha - 1).* (c - x).^(beta0 - 1);
bottom = beta(alpha, beta0).* (c - a).^(alpha + beta0 + 1);

Y = top./bottom;

end