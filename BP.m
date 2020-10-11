clc; 
clear all;

% prompt='enter file no. ';
% fl=input(prompt);

iter=50; % maximum number of iterations
tic
tot_trans=1e1; %takes 24hrs

id = fopen('mval_kk.txt','r'); 
m=fscanf(id,'%d')'; 
fclose(id); 

id2 = fopen('nval_kk.txt','r'); 
n=fscanf(id2,'%d')'; 
fclose(id2);

% H=zeros(m,n);
% id = fopen('H_1_4_11_3_5_unc3.txt','r');
id = fopen('H_3_11_3_5_sc.txt','r');
% id = fopen('mat_sc2.txt','r'); 
tmp=fscanf(id,'%d'); 
fclose(id); 

m
n

for i=1:m 
    H(i,:)=tmp((i-1)*n+1:i*n); 
end

% H

 for i=1:m
     rw(i)=sum(H(i,:));
 end
 for i=1:n
     cw(i)=sum(H(:,i));
 end    
    cw
%     rw

% cnt=0;
% rng('shuffle');
% for p=0.05; %linspace(0.01,0.1,10) %BSC crossover
%     err=0;
%     fprintf(' %d\n', cnt);
%     for trans = 1:tot_trans
%         rx=rand(1,n)<=p;
%     %     rx=zeros(1,n);
% 
%         out=logical(decoder_BSCMEX(double(rx), double(H), p, iter));
% 
%         if sum(out)>0
%            err=err+1;
%         end
% 
%     end 
%     cnt=cnt+1;
%     Pe(cnt)=err/tot_trans;
%     
% end
% 
% Pe

% filename = sprintf('%s%d.txt','Pe',fl);
% fileID = fopen(filename,'w');
% fprintf(fileID,'%f ',Pe');
% fclose(fileID);

% fileID = fopen('mat_sc1.txt','w');
% fprintf(fileID,'%d ',H');
% fclose(fileID);

toc




