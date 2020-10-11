clear all;
clc;


fileID = fopen('param_vec_MKay.txt','r');
param=fscanf(fileID,'%f')';
fclose(fileID);
 
param2=sort(param);
length(param2)

param2(1)

param2(length(param2)-100)

vals=ncx2rnd(1,param2);

fileID = fopen('ncx2_vec_MKay.txt','w');
fprintf(fileID,'%f ',vals(1:74600));
fclose(fileID);

fileID = fopen('param_vec_37AB.txt','w');
fprintf(fileID,'%f ',param(1:74600));
fclose(fileID);