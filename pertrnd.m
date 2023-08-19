function [Y, alpha, beta] = pertrnd(a, mu ,c)

r = c - a;

b = (6*mu - a - c)/4;

alpha = 1 + 4*(b - a)./r;
beta =  1 + 4*(c - b)./r;


X = betarnd(alpha, beta);

Y = r.*X + a;


end