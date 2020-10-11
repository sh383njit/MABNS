
//calls CS programs

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <time.h>
#include <omp.h>
#include <fstream>
#include <sstream>
#include <random>
#include <chrono>
//#include <Python.h>

using namespace std;

float *H,*H2,P,*LR,*pLR,*E_c_v,*E_c_v_mu,*E_c_v_old,*E_v_c,*E_v_c_mu,*H_D,*y_D,*LR_D,*res_c_v,*res_c_v_srtd,err,err2,err3,err4,err5,err6,biterr,biterr2,biterr3,biterr4,biterr5,biterr6,
varn,ncv,ncv2,ncv3,ncv4,ncv5,ncv6,ncv7,*ncv_vec,*ncv_vec2,*ncv_vec3,*ncv_vec4,*ncv_vec5,*ncv_vec6,*Q,*Q2,*Q3,*Q_cls,*Q_temp,*Q_temp2,*syn2,qstrt,gap,*Q_cnt,**G,*Gtemp,*Gmax,*mu,*mu_L,
*rep,*diff,*tmp2,*tmp2M,*x,*x_hat,*y,*res_cnt,*res_avg; //n is source length, m is code length
int m,n,n2,*vns,*cns,*vns_D,*cns_D,*excl_cw,BC,S,*indx,*indx2,*indx3,*vn_indx,*cn_indx,*indx_cls_CN,*indx_cls,*syn,*syn_cls,**cls_ind_mat_contg,**cls_ind_mat_ran,**cls_ind_mat_opt,*mu_cnt,*mu_cnt_L;
int col_wt,row_wt,num,dcw,dcw2,dcw3,dcw4,dcw5,dcw6,hshft,vshft,mem,W,call,meth,num_cls,cls_idx,cls_sz,*GI_arm,M,j_old; //row_wt is the max, row weight in case of SC codes
int *sv_ind; //coupling length (multiple of mem+1 and >=W) 

#define CW 1 //no. of codeword streams per transmission
#define I 2 //no. of iters
//#define ev 500 //min. no. of error events

string func(int n) {
	stringstream result;
	result << n;
	return result.str();
}

#include "max.cpp"
#include "min.cpp"
//#include "ipa4.cpp" //flooding	
//#include "ipa5.cpp" //residual based CN scheduling
//#include "ipa6.cpp" //LBP CN scheduling
//#include "MPCS.cpp" //MP based CS
#include "MPCS2.cpp" //sequential MP based CS


int main() {	
	srand(time(0));	
	clock_t tStart = clock();	
	int i,i2,i3,j,kk,k,AB,sum2,flg; 
	int fn,num_dat,cnt,gama,p,J,test=0,num,l1,l2,l3; 
	float trans=3,errprob=0,*acc,sum;

	float *blk_err,*blk_err2,*blk_err3,*blk_err4,*blk_err5,*blk_err6,*bit_err,*bit_err2,*bit_err3,*bit_err4,*bit_err5,*bit_err6,R;  
	//cout<<'\n'<<" device no: "; cin>>fn; ///////////////////////////////// 
	//cout<<'\n'<<" kk: "; cin>>kk; 

	BC=1;
	//construct a random generator engine:
	std::random_device rd;
    std::mt19937 e2(rd());

	ifstream inf,inf1,inf2,inf3,inf4,inf5,inf6,inf7,inf8,inf9,inf10,inf11,inf12,inf13,inf14,inf15; 

	//inf1.open("mval.txt"); inf1>>m; 
	//inf2.open("nval.txt"); inf2>>n;
	inf3.open("cls_sz.txt"); inf3>>cls_sz;
	//inf6.open("qstrt.txt"); inf6>>qstrt;
	//inf7.open("gap.txt"); inf7>>gap;
	inf8.open("M.txt"); inf8>>M;
	
	indx3= new int[M];
	rep= new float[M];
	diff= new float[M];

	qstrt=0; gap=8; 

	varn=0.01;
	std::normal_distribution<float> dist(0,sqrt(varn));

	m=98; n=196; AB=0; inf.open("mat_amp.txt");

	//m=98; n=196; AB=0; inf.open("mat_n196_1.txt"); inf12.open("Q1_1.txt"); inf14.open("Q2_1.txt"); inf15.open("Q3_1.txt"); rep[0]=-8.081; rep[1]=-0.423; rep[2]=1.27; rep[3]=9.82;
	//m=98; n=196; AB=0; inf.open("mat_n196_2.txt"); inf12.open("Q1_2.txt"); inf14.open("Q2_2.txt"); inf15.open("Q3_2.txt"); rep[0]=-8.081; rep[1]=-0.423; rep[2]=1.27; rep[3]=9.82;
	//m=98; n=196; AB=0; inf.open("mat_n196_3.txt"); inf12.open("Q1.txt"); inf14.open("Q2.txt"); inf15.open("Q3.txt"); rep[0]=-8.081; rep[1]=-0.423; rep[2]=1.27; rep[3]=9.82;

	//m=84; n=196; AB=1; inf.open("mat_n196_AB_1.txt"); inf12.open("Q4_1.txt"); inf14.open("Q5_1.txt"); inf15.open("Q6_1.txt"); rep[0]=-9.805; rep[1]=-0.189; rep[2]=1.77; rep[3]=13.36;
	//m=84; n=196; AB=1; inf.open("mat_n196_AB_2.txt"); inf12.open("Q4_2.txt"); inf14.open("Q5_2.txt"); inf15.open("Q6_2.txt"); rep[0]=-9.805; rep[1]=-0.189; rep[2]=1.77; rep[3]=13.36;
	//m=84; n=196; AB=1; inf.open("mat_n196_AB_3.txt"); inf12.open("Q4.txt"); inf14.open("Q5.txt"); inf15.open("Q6.txt"); rep[0]=-9.805; rep[1]=-0.189; rep[2]=1.77; rep[3]=13.36;

	//m=98; n=196; AB=0; inf.open("mat_n196_3.txt"); inf12.open("Q_0_1.8n196_opt_3.txt"); inf14.open("Q_0_1.8n196_opt_3.txt"); inf15.open("Q_0_1.8n196_opt_3.txt"); qstrt=4.8; gap=5.6; 
	//m=84; n=196; AB=1; inf.open("mat_n196_AB_1.txt"); inf12.open("Q_0_1.8n196_AB_opt_1.txt"); inf14.open("Q_0_1.8n196_AB_opt_1.txt"); inf15.open("Q_0_1.8n196_AB_opt_1.txt"); qstrt=0; gap=8;  

	if(test) {
		m=4; 
		n=5;
		col_wt=3; 
		row_wt=3;
		mem=1;
		gama=col_wt;
	}
	
	S=pow(M,cls_sz); //no. of syndromes for each cluster
	num_cls=m/cls_sz; //number of clusters

	R=1-float(m)/float(n); //code rate

	H= new float[m*n]; //H is the window. Since windows are identical, only the window Tanner graph is needed. The Lvalues "move" across the window
	H2= new float[m*n];
	ofstream outf,out_file2; string filename;

	Q=new float[S*num_cls*cls_sz]; 
	Q2=new float[S*num_cls*cls_sz];
	Q3=new float[S*num_cls*cls_sz];  
	Q_cnt=new float[S*num_cls*cls_sz];
	Q_temp=new float[cls_sz];
	Q_temp2=new float[m];
	Q_cls=new float[num_cls];
	
	//G=new float*[M]; for(i=0;i<M;i++) G[i]=new float[m];
	GI_arm=new int[I];
	mu=new float[I];
	mu_cnt=new int[I];
	mu_L=new float[I];
	mu_cnt_L=new int[I];

	sv_ind=new int[n];

	//inputs
	for(j=0;j<m;j++) for(i=0;i<n;i++) inf >> H[j*n+i];
	//for(j=0;j<m;j++) for(i=0;i<n;i++) H[j*n+i]=dist(e2); //for AMP

	for(i=0;i<S*num_cls;i++) for(j=0;j<cls_sz;j++) inf12 >> Q[i*cls_sz+j];  //load Q-table
	for(i=0;i<S*num_cls;i++) for(j=0;j<cls_sz;j++) inf14 >> Q2[i*cls_sz+j];
	for(i=0;i<S*num_cls;i++) for(j=0;j<cls_sz;j++) inf15 >> Q3[i*cls_sz+j];  
	
	//for(i=0;i<M;i++) for(j=0;j<m;j++) inf9 >> G[i][j];
	for(i=0;i<I;i++) inf9 >> GI_arm[i];  

	
	cls_ind_mat_contg=new int*[num_cls]; for(i=0;i<num_cls;i++) cls_ind_mat_contg[i]=new int[cls_sz];
	cls_ind_mat_ran=new int*[num_cls]; for(i=0;i<num_cls;i++) cls_ind_mat_ran[i]=new int[cls_sz];
	cls_ind_mat_opt=new int*[num_cls]; for(i=0;i<num_cls;i++) cls_ind_mat_opt[i]=new int[cls_sz];
	if(AB) {
		inf13.open("cls_ind_mat_contg_AB.txt"); for(i=0;i<num_cls;i++) for(j=0;j<cls_sz;j++) inf13 >> cls_ind_mat_contg[i][j]; //contiguous clustering
		inf10.open("cls_ind_mat_ran_AB.txt"); for(i=0;i<num_cls;i++) for(j=0;j<cls_sz;j++) inf10 >> cls_ind_mat_ran[i][j]; //random clustering
		inf4.open("cls_ind_mat_opt_AB_1.txt"); for(i=0;i<num_cls;i++) for(j=0;j<cls_sz;j++) inf4 >> cls_ind_mat_opt[i][j]; //optimized clustering
	}
	else {
		inf13.open("cls_ind_mat_contg.txt"); for(i=0;i<num_cls;i++) for(j=0;j<cls_sz;j++) inf13 >> cls_ind_mat_contg[i][j];
		inf10.open("cls_ind_mat_ran.txt"); for(i=0;i<num_cls;i++) for(j=0;j<cls_sz;j++) inf10 >> cls_ind_mat_ran[i][j]; 
		inf4.open("cls_ind_mat_opt_1.txt"); for(i=0;i<num_cls;i++) for(j=0;j<cls_sz;j++) inf4 >> cls_ind_mat_opt[i][j];
	}

	//cout<<'\n'<<"cls_ind_mat: "<<'\n'; for(i=0;i<num_cls;i++) {for(j=0;j<cls_sz;j++) cout<<cls_ind_mat_opt[i][j]<<" "; cout<<'\n';}

	/*for(i=0;i<S*num_cls;i++) for(j=0;j<cls_sz;j++) inf9 >> Q_cnt[i*cls_sz+j];
	for(i=0;i<S*num_cls;i++) 
		for(j=0;j<cls_sz;j++)
			if(Q_cnt[i*cls_sz+j]!=0)
				Q[i*cls_sz+j]/=Q_cnt[i*cls_sz+j];*/

	//cout<<'\n'<<"H: "<<'\n'; for(j=0;j<m;j++){for(i=0;i<n;i++) cout<<H[j*n+i]<<" "; cout<<'\n';}
	//cout<<'\n'<<"Q1: "<<'\n'; for(i=0;i<S*num_cls;i++) for(j=0;j<cls_sz;j++) cout<<Q[i*cls_sz+j]<<" ";
	//cout<<'\n'<<"G: "<<'\n'; for(i=0;i<M;i++) for(j=0;j<m;j++) cout<<G[i][j]<<" ";
	//cout<<'\n'<<"GI_arm: "<<'\n'; for(i=0;i<I;i++) cout<<GI_arm[i]<<" ";

	if(!test) {
		col_wt=0;
		for(j=0;j<m;j++) if(H[j*n+0]) col_wt++;
		row_wt=0;
		for(i=0;i<n;i++) if(H[0*n+i]) row_wt++;
		//row_wt=7;
	}

	cout<<'\n'<<"m: "<<m<<" n: "<<n<<" R: "<<R<<" col_wt: "<<col_wt<<" row_wt: "<<row_wt<<" W: "<<W;  
	cout<<'\n'<<"cls_sz: "<<cls_sz<<" num_cls: "<<num_cls<<" S: "<<S<<" qstrt: "<<qstrt<<" gap: "<<gap<<" M: "<<M;

	cout<<'\n'<<"hshft: "<<hshft<<" vshft: "<<vshft<<" n2: "<<n2<<endl; 
	
	//2D arrays
	x= new float[n]; 
	x_hat= new float[n]; 
	y=new float[m]; 
	tmp2=new float[col_wt]; 
	tmp2M=new float[col_wt]; 
	LR= new float[CW*n];  //Lvalues in the window
	pLR= new float[CW*n]; 
	syn=new int[m];
	syn_cls=new int[cls_sz]; 
	syn2=new float[m];
	excl_cw=new int[CW];

	Gtemp=new float[M];
	Gmax=new float[m];
	
	vn_indx=new int[m];
	cn_indx=new int[m];
	indx=new int[cls_sz];
	indx_cls=new int[num_cls];
	indx_cls_CN=new int[num_cls];

	indx2=new int[m];

	vns=new int[m*row_wt]; //all VNs of CN j in a row
	
	//3D arrays
	E_v_c= new float[m*row_wt*CW];
	E_v_c_mu= new float[m*row_wt*CW];
	E_c_v= new float[n*col_wt*CW];
    E_c_v_mu= new float[n*col_wt*CW];
	E_c_v_old= new float[n*col_wt*CW];
	res_c_v= new float[m*row_wt*CW];
	res_c_v_srtd= new float[m*row_wt*CW];
	res_cnt=new float[m];
	res_avg=new float[m];
	
	cns=new int[n*col_wt]; //all CNs of a VN in a row

	//float spar[]={0.05}; //,0.025,0.05,0.075,0.1,0.125,0.15};
	float spar[]={0.01}; //,0.02,0.03,0.04,0.05,0.06,0.075};
	num_dat=sizeof(spar)/sizeof(spar[0]);
	acc=new float[num_dat];

	blk_err= new float[num_dat]; 
	blk_err2= new float[num_dat];
	blk_err3= new float[num_dat]; 
	blk_err4= new float[num_dat]; 
	blk_err5= new float[num_dat];
	blk_err6= new float[num_dat];

	bit_err= new float[num_dat]; 
	bit_err2= new float[num_dat]; 	
	bit_err3= new float[num_dat]; 
	bit_err4= new float[num_dat]; 
	bit_err5= new float[num_dat];
	bit_err6= new float[num_dat];

	ncv_vec=new float[num_dat];
	ncv_vec2=new float[num_dat];
	ncv_vec3=new float[num_dat];
	ncv_vec4=new float[num_dat];
	ncv_vec5=new float[num_dat];
	ncv_vec6=new float[num_dat];

	//initializing Tanner graph information
	for(i=0;i<m;i++) for(j=0;j<row_wt;j++) vns[i*row_wt+j]=-1;
	for(i=0;i<n;i++) for(j=0;j<col_wt;j++) cns[i*col_wt+j]=-1;

	for(j=0;j<m;j++) {
		cnt=0; 
		for(i=0;i<n;i++) if(H[j*n+i]) {vns[j*row_wt+cnt]=i; cnt++;} 
	}
	//cout<<'\n'<<"vns "<<'\n'; for(i=0;i<m;i++) {for(j=0;j<row_wt;j++) cout<<vns[i*row_wt+j]<<" "; cout<<'\n';}
	for(i=0;i<n;i++) {
		cnt=0; 
		for(j=0;j<m;j++) if(H[j*n+i]) {cns[i*col_wt+cnt]=j; cnt++;} 
	}
	//cout<<'\n'<<"cns "<<'\n'; for(i=0;i<n;i++) {for(j=0;j<col_wt;j++) cout<<cns[i*col_wt+j]<<" "; cout<<'\n';}

	for(k=0;k<num_dat;k++) {
		errprob=0;
		cout<<"k: "<<k<<endl;
		for(i3=0;i3<trans;i3++) {
			//generating sparse x
			sum2=0; flg=1; 
			for(i=0;i<n;i++) x[i]=x_hat[i]=0; for(i=0;i<m;i++) y[i]=0; //refresh

			while(sum2<spar[k]*n) {
				i2=rand()%n; 
				for(j=0;j<sum2;j++) if(sv_ind[j]==i2) {flg=0; break;} else flg=1;
				if(flg)	{
					//x[i2]=dist(e2);
					x[i2]=1; 
					sv_ind[sum2]=i2; 
					sum2++;
				} //rand()%2
			}
			//x[0]=x[2]=x[6]=x[34]=x[56]=x[23]=x[11]=x[8]=x[56]=x[90]=1;
			//cout<<'\n'<<"x: "; for(i=0;i<n;i++) cout<<x[i]<<" "; cout<<endl;

			for(j=0;j<m;j++) {sum=0; for(i=0;i<n;i++) sum+=x[i]*H[j*n+i]; y[j]=sum;}
			//cout<<'\n'<<"y: "; for(i=0;i<m;i++) cout<<y[i]<<" ";
		
			//sipa();
			MPCS();
			//cout<<'\n'<<"x_h: "; for(i=0;i<n;i++) cout<<x_hat[i]<<" "; cout<<endl;
			
			if(err>0) errprob++;
		}
		acc[k]=1-errprob/trans;
	}
	cout<<'\n'<<"acc: "; for(i=0;i<num_dat;i++) cout<<acc[i]<<" ";
	cout<<'\n'<<"spar: "; for(i=0;i<num_dat;i++) cout<<spar[i]<<" ";
	cout<<"res_avg: "; for(i=0;i<m;i++) cout<<res_cnt[i]/res_avg[i]<<" "; cout<<'\n';

		//trans=(call*n2/hshft)*CW; //there are n2/hshft cws in a stream. Each cw has length hshft bits in a stream of length n2 bits
		/*blk_err[i2]=err/dcw;
		bit_err[i2]=biterr/(hshft*dcw);
		ncv_vec[i2]=ncv/(1*dcw); //no. of CN to VN msg. updates

		blk_err2[i2]=err2/dcw2;		
		bit_err2[i2]=biterr2/(hshft*dcw2);
		ncv_vec2[i2]=ncv2/(1*dcw2); 

		blk_err3[i2]=err3/dcw3;		
		bit_err3[i2]=biterr3/(hshft*dcw3);
		ncv_vec3[i2]=ncv3/(1*dcw3); 

		blk_err4[i2]=err4/dcw4;		
		bit_err4[i2]=biterr4/(hshft*dcw4);
		ncv_vec4[i2]=ncv4/(1*dcw4);

		blk_err5[i2]=err5/dcw5;		
		bit_err5[i2]=biterr5/(hshft*dcw5);
		ncv_vec5[i2]=ncv5/(1*dcw5); 

		blk_err6[i2]=err6/dcw6;		
		bit_err6[i2]=biterr6/(hshft*dcw6);
		ncv_vec6[i2]=ncv6/(1*dcw6); 		
		 
	}
	}

	cout<<'\n'<<"EbNo_BC: "; for(i=0;i<num_dat;i++) cout<<EbNo_BC[i]<<" "; cout<<'\n'<<endl;
	//cout<<'\n'<<"err: "<<err<<" dcw: "<<dcw<<'\n';
	//cout<<"blk_err: "; for(i=0;i<num_dat;i++) cout<<blk_err[i]<<" "; cout<<'\n';
	cout<<"bit_err: "; for(i=0;i<num_dat;i++) cout<<bit_err[i]<function [ S, M, T, H, Z, MAI, stopnorm ] = AMP_real_tau_fun(y, x, A, N, delta, varx, I, msetol, stoptol, sigma2, x0)
        

	cout<<"ncv_vec: "; for(i=0;i<num_dat;i++) cout<<ncv_vec[i]<<" "; cout<<'\n'<<endl;

	//cout<<'\n'<<"err2: "<<err2<<" dcw2: "<<dcw2<<'\n';
	//cout<<"blk_err2: "; for(i=0;i<num_dat;i++) cout<<blk_err2[i]<<" "; cout<<'\n';
	cout<<"bit_err2: "; for(i=0;i<num_dat;i++) cout<<bit_err2[i]<<" "; cout<<'\n';
	cout<<"ncv_vec2: "; for(i=0;i<num_dat;i++) cout<<ncv_vec2[i]<<" "; cout<<'\n'<<endl;

	//cout<<'\n'<<"err3: "<<err3<<" dcw3: "<<dcw3<<'\n';
	//cout<<"blk_err3: "; for(i=0;i<num_dat;i++) cout<<blk_err3[i]<<" "; cout<<'\n';
	cout<<"bit_err3: "; for(i=0;i<num_dat;i++) cout<<bit_err3[i]<<" "; cout<<'\n';
	cout<<"ncv_vec3: "; for(i=0;i<num_dat;i++) cout<<ncv_vec3[i]<<" "; cout<<'\n'<<endl;

	//cout<<'\n'<<"err4: "<<err4<<" dcw4: "<<dcw4<<'\n';
	//cout<<"blk_err3: "; for(i=0;i<num_dat;i++) cout<<blk_err3[i]<<" "; cout<<'\n';
	cout<<"bit_err4: "; for(i=0;i<num_dat;i++) cout<<bit_err4[i]<<" "; cout<<'\n';
	cout<<"ncv_vec4: "; for(i=0;i<num_dat;i++) cout<<ncv_vec4[i]<<" "; cout<<'\n'<<endl;

	//cout<<'\n'<<"err5: "<<err5<<" dcw5: "<<dcw5<<'\n';
	//cout<<"blk_err5: "; for(i=0;i<num_dat;i++) cout<<blk_err5[i]<<" "; cout<<'\n';
	cout<<"bit_err5: "; for(i=0;i<num_dat;i++) cout<<bit_err5[i]<<" "; cout<<'\n';
	cout<<"ncv_vec5: "; for(i=0;i<num_dat;i++) cout<<ncv_vec5[i]<<" "; cout<<'\n'<<endl;

	cout<<"bit_err6: "; for(i=0;i<num_dat;i++) cout<<bit_err6[i]<<" "; cout<<'\n';
	cout<<"ncv_vec6: "; for(i=0;i<num_dat;i++) cout<<ncv_vec6[i]<<" "; cout<<'\n'<<endl;

	for(i=0;i<I;i++) mu[i]/=mu_cnt[i]; 
	//cout<<"mean m_{c_v}: "; for(i=0;i<I;i++) cout<<mu[i]<<" "; cout<<'\n'<<endl;

	for(i=0;i<I;i++) mu_L[i]/=mu_cnt_L[i]; 
	//cout<<"mean mu_{Li_hat}: "; for(i=0;i<I;i++) cout<<mu_L[i]+col_wt*mu[i]<<" "; cout<<'\n'<<endl;

	cout<<"mean K*mu_{Li_hat}: "; for(i=0;i<I;i++) cout<<row_wt*(mu_L[i]+col_wt*mu[i])<<" "; cout<<'\n'<<endl;
	cout<<"std. dev. 2*K*mu_{Li_hat}: "; for(i=0;i<I;i++) cout<<sqrt(2*row_wt*(mu_L[i]+col_wt*mu[i]))<<" "; cout<<'\n'<<endl;*/

	//ofstream outf3; filename="Pe"+func(fn)+".txt"; outf3.open(filename.c_str()/*,fstream::app*/); for(i2=0;i2<num_dat;i2++) outf3<<blk_err[i2]<<" ";  outf3<<std::endl; outf3.close();
	//ofstream outf3; outf3.open("mat_amp.txt"); for(i=0;i<m;i++) for(j=0;j<n;j++) outf3<<H[i*n+j]<<" ";  outf3<<std::endl; outf3.close();

	cout<<'\n';
    	printf("executed in: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
	cout<<'\n';


}





