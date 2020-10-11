
%from L=2 or 3 terminally lifted AB SC matrix, construct a larger L matrx

clear all;
clc;

% id = fopen('mval.txt','r'); 
% m=fscanf(id,'%d')'; 
% fclose(id); 
% 
% id2 = fopen('nval.txt','r'); 
% n=fscanf(id2,'%d')'; 
% fclose(id2);

m=49; n=98;

id = fopen('mat_n98_AB.txt','r'); 
tmp=fscanf(id,'%d'); 
fclose(id);

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
    rw


% fileID = fopen('win_3_sc.txt','w');
% fprintf(fileID,'%d ',H3');
% fclose(fileID);
% 
% fileID = fopen('m_win.txt','w');
% fprintf(fileID,'%d ',m3);
% fclose(fileID);
% 
% fileID = fopen('W.txt','w');
% fprintf(fileID,'%d ',L);
% fclose(fileID);
