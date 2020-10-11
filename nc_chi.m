
%generates non-central chi-squared distribution
clear all;
clc;


x = (0:0.1:10)';
ncx2 = ncx2cdf(x,1,2);

del=0.5;
val=ncx2rnd(1,del)

% plot(p,I_ncx2,'o')