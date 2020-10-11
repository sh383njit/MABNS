
%majority rule based decoding taken from "Reed-muller error correcting codes", Ben Cooke

clear all;
clc;

tic


m=7+1;
n=2^(m-1)
x=zeros(1,n);
r=3; %change 0<r<=m-1 to adjust the rate of the RM code
trans=3e3;

j=0;
for i=1:r
    j=j+nchoosek(m-1,i);
end
k=j+1

R=k/n

mat_cvec=[]; %stores all characteristic vectors of each row of G matrix (except the 1st row)
G=zeros(k,n); %encoder matrix
msg=zeros(1,k);
[G,mat_cvec]=enc(r,m,n,x,G,mat_cvec);

snr=6; %3:1:5; %
% snr=[0 0.6 1.2 1.4 1.8]; %eb/no values in db

for i=1:length(snr)
    sigma(i)=sqrt(1/(2*R*10^(0.1*snr(i))));
% sigma(i)=sqrt(1/(2*R*snr(i)));
end

x=zeros(1,n);
for j=1:length(sigma)
    p=qfunc(1/sigma(j)); %AWGN approximated as BSC
    err=0;
    for i=1:trans
%         msg=rand(1,k)<=0.5;
%         x=msg*G;
%         x=mod(x,2);
        
        y=x+(rand(1,n)<=p);
        y=mod(y,2); %erronesous 
        
        x_hat=dec(r,m,n,G,mat_cvec,k,y);
        dist=sum(xor(x,x_hat));

        if(dist>0) 
            err=err+dist; %dist; %total no. of bits in error
        end
    end

    BER(j)=err/(n*trans);
end

snr
BER
 
toc


function[G,mat_cvec]=enc(r,m,n,x,G,mat_cvec)
    G(1,:)=ones(1,n);
    row=1:m-1; %row indices of G matrix exlcuding the 1st one
    
    num=nchoosek(m-1,1);
    for j=1:num
        div=2^j;
        val=n/div;
        for i=1:div
            strt=val*(i-1)+1;
            stp=strt+n/div-1;
            if mod(i,2)
                x(j,strt:stp)=1;
            else
                x(j,strt:stp)=0;
            end
        end
    end
    % x

    ridx=2;
    ridx2=1;
    for k=1:r
        indx=nchoosek(1:m-1,k);
        
        for i=1:size(indx,1)
            tmp(ridx2,1:size(indx,2))=indx(i,:);
            A(ridx2,1:length(setdiff(row,tmp(ridx2,:))))=setdiff(row,tmp(ridx2,:)); %a row of A has characteristic x vector indices of the correponding row in G 
            ridx2=ridx2+1;
        end
        
        for j=1:size(indx,1)
            for i=1:k
                if i==1
                    xvec=x(indx(j,i),:);
                else
                    xvec=xvec.*x(indx(j,i),:);
                end
            end
            G(ridx,:)=xvec;
            ridx=ridx+1;
        end
    end
    
    %finding characteristic vectors 
%     A 
    
    for i=1:size(A,1)
        num_x=0;
        for j=1:size(A,2) %j is a row of A
            if A(i,j)
                num_x=num_x+1; %num_x is the no. of x vectors not involved in a row of G
            end
        end
        
        if num_x
            combs=nchoosek(1:2*num_x,num_x); 
            for j=1:size(combs,1)
                for j2=1:num_x
                    if combs(j,j2)>num_x
                        col=combs(j,j2)-num_x;
                        temp=~x(A(i,col),:); %taking the complement of a x vector
                    else
                        col=combs(j,j2);
                        temp=x(A(i,col),:); %dont' take complement
                    end

                    if j2==1
                        xvec=temp;
                    else
                        xvec=xvec.*temp;
                    end
                end
                mat_cvec(j,:,i)=xvec; %storing all characteristic vectors for row i of G
            end
        else
            mat_cvec(j,:,i)=zeros(1,2^(m-1));
        end
    end
end



function x_hat=dec(r,m,n,G,mat_cvec,k,y)

    rows=k-1; %no. of rows of G excluding the 1st one
    
    h=size(mat_cvec,1); %no. of rows of a matrix in mat_vec
    

    temp=[];
    vec2=zeros(1,n);
    for i=rows:-1:1
        cnt=1;
        for j=1:h
            if sum(mat_cvec(j,:,i))
                val=mat_cvec(j,:,i)*y'; % y is the erroneous word
                temp(cnt)=mod(val,2);
                cnt=cnt+1;
            end
        end
        
        %t
        sum1=0;
        sum0=0;
        for j=1:length(temp)
            if temp(j) %if temp has more 1 than zeros
                sum1=sum1+1;
            else
                sum0=sum0+1;
            end
        end
        
%         temp
        if sum1>sum0
            coeff=1;
        else
            coeff=0;
        end
        
        vec=coeff*G(i+1,:);
        vec2=vec2+vec;
        
    end
   

    My=mod(vec2,2);
    
    My2=My+y;
    My2=mod(My2,2);
    
    sum1=0;
    sum0=0;
    for j=1:n
        if My2(j) %if temp has more 1 than zeros
        	sum1=sum1+1;
        else
            sum0=sum0+1;
        end
    end
    
    if sum1>sum0
        x_hat=My+1;
        x_hat=mod(x_hat,2);
    else
    	x_hat=mod(My,2);
    end
    

end
