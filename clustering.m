
%from L=2 or 3 terminally lifted AB SC matrix, construct a larger L matrix

clear all;
clc;

% id = fopen('mval.txt','r'); 
% m=fscanf(id,'%d')'; 
% fclose(id); 
% 
% id2 = fopen('nval.txt','r'); 
% n=fscanf(id2,'%d')'; 
% fclose(id2);

% m=84; n=196; %for AB
m=48; n=96; %for non AB

rng('shuffle')

z=6
u=m/z
vec=randperm(m)-1 %for random case
% vec=(1:m)-1 %for contiguous case

for i=1:u
    strt=(i-1)*z+1;
    stp=strt+z-1;
    clus_ind_mat(i,:)=vec(strt:stp);
end

clus_ind_mat


%fileID = fopen('cls_ind_mat_contg_AB.txt','w');
% fileID = fopen('cls_ind_mat_ran_AB.txt','w');

% fileID = fopen('cls_ind_mat_contg.txt','w');
fileID = fopen('cls_ind_mat_ran.txt','w');

fprintf(fileID,'%d ',clus_ind_mat');
fclose(fileID);
