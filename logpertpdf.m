function [ln_Y, alpha, beta0] = logpertpdf(x, a, mu ,c)

r = c - a;

b = (6*mu - a - c)/4;

alpha = 1 + 4*(b - a)./r;
beta0 =  1 + 4*(c - b)./r;

top = (alpha - 1).*log(x - a) + (beta0 - 1).*log(c - x);
bottom =  (alpha + beta0 + 1).*log(c - a);

ln_Y = top - bottom;



end