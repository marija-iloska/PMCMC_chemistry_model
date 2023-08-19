function [Y, alpha, beta] = pertrnd(a,b,c)

r = c - a;

alpha = 1 + 4*(b - a)./r;
beta =  1 + 4*(c - b)./r;

X = betarnd(alpha, beta);

Y = r.*X + a;


end