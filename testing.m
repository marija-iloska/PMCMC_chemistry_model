clear all
close all
clc

x = -10:0.01:10;
y1 = pertpdf(x, -10, 0, 10);
y2 = pertpdf(x, -10, -0, 10);
y3 = pertpdf(x, -10, 0, 10);


plot(x, y1)
hold on
plot(x, y2)
hold on
plot(x, y3)

