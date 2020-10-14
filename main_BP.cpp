
//for BC decoding L=W=3
//$ sudo apt-get install python3-dev (if dev not installed)
//$ g++ -o main2 main_BP.cpp -I/usr/include/python3.7 -L/usr/lib/python3.7/config-3.7m-x86_64-linux-gnu -lpython3.7
//g++ -o main2 -std=c++11 -I /home/salman/extern/include/ -L /home/salman/extern/bin/glnxa64/ -pthread main_BP.cpp -lMatlabDataArray -lMatlabEngine

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
//#include <fdeep/fdeep.hpp>
//#include "MatlabEngine.hpp"
//#include "MatlabDataArray.hpp"

using namespace std;

// Pass vector containing MATLAB data array scalar
/*using namespace matlab::engine;
// Start MATLAB engine synchronously
std::unique_ptr<MATLABEngine> matlabPtr = startMATLAB();
// Create MATLAB data array factory
matlab::data::ArrayFactory factory;*/


float *H,*H2,*x,*x_hat,*y,P,*LR,*pLR,*E_c_v,*E_c_v_old,*E_v_c,*E_c_v_mu,*E_v_c_mu,*H_D,*y_D,*LR_D,*res_c_v,*res_c_v_srtd,err,err2,err3,err4,err5,err6,err7,err8,biterr,biterr2,biterr3,biterr4,biterr5,biterr6,biterr7,biterr8,
varn,ncv,ncv2,ncv3,ncv4,ncv5,ncv6,ncv7,ncv8,*ncv_vec,*ncv_vec2,*ncv_vec3,*ncv_vec4,*ncv_vec5,*ncv_vec6,*ncv_vec7,*ncv_vec8,*Q,*Q2,*Q3,*Q4,*Q5,*Q6,*Q_cls,*Q_temp,*Q_temp2,*syn2,qstrt,gap,*Q_cnt,**G,*Gtemp,*Gmax,*mu,*mu_L,
*rep,*diff,*res_mean,*res_cnt,*res_mean2,*res_cnt2,*mean1_vec,*mean2_vec,*param_vec,*val_vec,*ncx2_vec; //n is source length, m is code length
long m,n,n2,*vns,*cns,*vns_D,*cns_D,*excl_cw,BC,S,*indx,*indx2,*indx3,*vn_indx,*cn_indx,*indx_cls_CN,*indx_cls,*cls_s_crnt,*cls_s_prev,*pick_cls,*syn,*syn_cls,**cls_ind_mat_contg,**cls_ind_mat_ran,**cls_ind_mat_opt,*mu_cnt,*mu_cnt_L;
long col_wt,row_wt,num,dcw,dcw2,dcw3,dcw4,dcw5,dcw6,dcw7,dcw8,hshft,vshft,mem,W,call,meth,num_cls,cls_idx,cls_sz,*GI_arm,M,cnt2=0,j_old,**Oc,dv,dc,DeepRL,boxplus=0; 
long L=99; //coupling length (multiple of mem+1 and >=W) 

#define CW 1 //no. of codeword streams per transmission
#define I 25 //no. of NS iters
#define Ifl 25 //no. of flooding iters
#define ev 500 //min. no. of error events
#define avg 10

string func(long n) {
	stringstream result;
	result << n;
	return result.str();
}

std::vector<float> inp;

#include "BP.cpp"	


int main() {	
	srand(time(0));	
	clock_t tStart = clock();	
	long i,i2,j,kk,AB,mac,bch; 
	long fn,num_dat,cnt,cw,gama,p,J,test=0,num,l1,l2,l3,matnum; 

	//cout<<'\n'<<"DeepRL ? "; cin>>DeepRL; //1 yes, 0 no
	DeepRL=0;
	float *blk_err,*blk_err2,*blk_err3,*blk_err4,*blk_err5,*blk_err6,*blk_err7,*blk_err8,*bit_err,*bit_err2,*bit_err3,*bit_err4,*bit_err5,*bit_err6,*bit_err7,*bit_err8,R;  
	cout<<'\n'<<" file no: "; cin>>fn; ///////////////////////////////// 	
		
	//indx3= new long[M];
	//diff= new float[M];
	 
	BC=1;
	//construct a random generator engine:
	std::random_device rd;
    std::mt19937 e2(rd());

	ifstream inf,inf1,inf2,inf3,inf4,inf5,inf6,inf7,inf8,inf9,inf10,inf11,inf12,inf13,inf14,inf15,inf16,inf17,inf18; 

	//m=500; n=1000; inf.open("mat_nonAB.txt"); dv=3; dc=6; 
	//m=420; n=980; AB=1; inf.open("mat_AB.txt"); dv=3; dc=7; 

	//****************************** the b (resp. c) version was trained with 5*10^7 (resp. 1.25*10^6) training examples

	//********************************************************** non-AB *****************************************************************
	
	if(fn==1 || fn==2 || fn==3) {
		m=98; n=196; AB=0; mac=0; bch=0; dv=3; dc=6; col_wt=3; row_wt=6; 
		rep= new float[M];
		if(!boxplus) {rep[0]=-8.081; rep[1]=-0.423; rep[2]=1.27; rep[3]=9.82;} 
		else {rep[0]=-2.901; rep[1]=-0.376; rep[2]=0.317; rep[3]= 1.12336;} //for box-plus syndrome
		cls_sz=7; M=4;
	}

	if(fn==1) {
		inf.open("matrices/mat_n196_1.txt"); matnum=1; 
		inf12.open("Qtables/Q1_1c.txt"); inf14.open("Qtables/Q2_1c.txt"); inf15.open("Qtables/Q3_1c.txt"); //model-free matrix 1
		inf16.open("QMtables/Q1_1_M.txt"); inf17.open("QMtables/Q2_1_M.txt"); inf18.open("QMtables/Q3_1_M.txt"); //model-based matrix 1
	}
	else if(fn==2) {
		inf.open("matrices/mat_n196_2.txt"); matnum=2; 
		inf12.open("Qtables/Q1_2c.txt"); inf14.open("Qtables/Q2_2c.txt"); inf15.open("Qtables/Q3_2c.txt"); //model-free matrix 2
		inf16.open("QMtables/Q1_2_M.txt"); inf17.open("QMtables/Q2_2_M.txt"); inf18.open("QMtables/Q3_2_M.txt"); //model-based matrix 2
	}
	else if(fn==3) {
		inf.open("matrices/mat_n196_3.txt"); matnum=3; 
		inf12.open("Qtables/Q1_3c.txt"); inf14.open("Qtables/Q2_3c.txt"); inf15.open("Qtables/Q3_3c.txt"); //model-free matrix 3
		inf16.open("QMtables/Q1_3_M.txt"); inf17.open("QMtables/Q2_3_M.txt"); inf18.open("QMtables/Q3_3_M.txt"); //model-based matrix 3
	}

	//********************************************************** AB *****************************************************************
	if(fn==4 || fn==5 || fn==6) {
		m=84; n=196; AB=1; mac=0; bch=0; dv=3; dc=7; col_wt=3; row_wt=7; cls_sz=7; M=4;
		rep= new float[M];
		rep[0]=-9.805; rep[1]=-0.189; rep[2]=1.77; rep[3]=13.36;
	}

	if(fn==4) {
		inf.open("matrices/mat_n196_AB_1.txt"); matnum=1; 
		inf12.open("Qtables/Q4_1c.txt"); inf14.open("Qtables/Q5_1c.txt"); inf15.open("Qtables/Q6_1c.txt"); //model-free matrix 1
		inf16.open("QMtables/Q4_1_M.txt"); inf17.open("QMtables/Q5_1_M.txt"); inf18.open("QMtables/Q6_1_M.txt"); //model-based matrix 1
	}
	else if(fn==5) {
		inf.open("matrices/mat_n196_AB_2.txt"); matnum=2; 
		inf12.open("Qtables/Q4_2c.txt"); inf14.open("Qtables/Q5_2c.txt"); inf15.open("Qtables/Q6_2c.txt"); //model-free matrix 2
		inf16.open("QMtables/Q4_2_M.txt"); inf17.open("QMtables/Q5_2_M.txt"); inf18.open("QMtables/Q6_2_M.txt"); //model-based matrix 2
	}	
	else if(fn==6) {
		inf.open("matrices/mat_n196_AB_3.txt"); matnum=3; 
		inf12.open("Qtables/Q4_3c.txt"); inf14.open("Qtables/Q5_3c.txt"); inf15.open("Qtables/Q6_3c.txt"); //model-free matrix 3 
		inf16.open("QMtables/Q4_3_M.txt"); inf17.open("QMtables/Q5_3_M.txt"); inf18.open("QMtables/Q6_3_M.txt"); //model-based matrix 3
	}

	//Mackay
	if(fn==7 || fn==8 || fn==9) {
		m=48; n=96; AB=0; mac=1; bch=0; col_wt=3; row_wt=6; cls_sz=6; M=4;
		rep= new float[M];
		if(M==8) {rep[0]=-3.236; rep[1]=9.65; rep[2]=22.2; rep[3]=34.4; rep[4]=46.8; rep[5]=61; rep[6]=76.3; rep[7]=96.8;}
		else {rep[0]=-3.236; rep[1]=13.7; rep[2]=32.4; rep[3]=57.3;}
	}
	
	if(fn==7) {
		inf.open("matrices/mackay_96_48.txt"); matnum=1; 
		inf14.open("Qtables/Q9_1.txt"); inf15.open("Qtables/Q10_1.txt"); //model-free opt-clus
		inf17.open("QMtables/Q9_1_M.txt"); inf18.open("QMtables/Q10_1_M.txt"); //model-based
	}
	else if(fn==8) {
		inf.open("matrices/mackay_96_48_2.txt"); matnum=2; 
		inf14.open("Qtables/Q9_2.txt"); inf15.open("Qtables/Q10_2.txt"); //model-free opt-clus
		inf17.open("QMtables/Q9_2_M.txt"); inf18.open("QMtables/Q10_2_M.txt"); //model-based
	}
	else if(fn==9) {
		inf.open("matrices/mackay_96_48_3.txt"); matnum=3; 
		inf14.open("Qtables/Q9_3.txt"); inf15.open("Qtables/Q10_3.txt"); //model-free opt-clus
		inf17.open("QMtables/Q9_3_M.txt"); inf18.open("QMtables/Q10_3_M.txt"); //model-based
	}
	
	//BCH
	if(fn==10) {
		m=12; n=63; AB=0; mac=0; bch=1; col_wt=8; row_wt=28; cls_sz=6; M=8;
		rep= new float[M];
		if(M==4) {rep[0]=74.61; rep[1]=104.6; rep[2]=119.58; rep[3]=131.56;}
		else {rep[0]=20.8; rep[1]=44.2; rep[2]=64; rep[3]=78.8; rep[4]=92.7; rep[5]=109; rep[6]=124; rep[7]=142;}
		inf.open("matrices/mat_BCH_63_51.txt"); matnum=1; 
		inf14.open("Qtables/Q7_1.txt"); inf15.open("Qtables/Q8_1.txt"); //model-free opt-clus
	}

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

	if(BC || test) {L=3;W=3;}

	Q=new float[S*num_cls*cls_sz]; 
	Q2=new float[S*num_cls*cls_sz];
	Q3=new float[S*num_cls*cls_sz]; 
	Q4=new float[S*num_cls*cls_sz]; 
	Q5=new float[S*num_cls*cls_sz];
	Q6=new float[S*num_cls*cls_sz]; 
 
	Q_cnt=new float[S*num_cls*cls_sz];
	Q_temp=new float[cls_sz];
	Q_temp2=new float[m];
	Q_cls=new float[num_cls];
	
	//G=new float*[M]; for(i=0;i<M;i++) G[i]=new float[m];
	GI_arm=new long[I];
	mu=new float[I];
	mu_cnt=new long[I];
	mu_L=new float[I];
	mu_cnt_L=new long[I];
	
	param_vec= new float[74600];
	ncx2_vec= new float[74600];
	val_vec= new float[100];
	
	if(!AB && !mac && !bch) {
		inf1.open("param_vec_36nonAB.txt"); for(i=0;i<74600;i++) inf1 >> param_vec[i]; //param is generated using main_RL and ncx2_vec is generated using non_cntrl_chi.m 
		inf3.open("ncx2_vec_36nonAB.txt"); for(i=0;i<74600;i++) inf3 >> ncx2_vec[i];
	}
	else if(mac) {
		inf1.open("param_vec_MKay.txt"); for(i=0;i<74600;i++) inf1 >> param_vec[i]; //param is generated using main_RL and ncx2_vec is generated using non_cntrl_chi.m 
		inf3.open("ncx2_vec_MKay.txt"); for(i=0;i<74600;i++) inf3 >> ncx2_vec[i];
	}
	else if(AB) {
		inf1.open("param_vec_37AB.txt"); for(i=0;i<74600;i++) inf1 >> param_vec[i];
		inf3.open("ncx2_vec_37AB.txt"); for(i=0;i<74600;i++) inf3 >> ncx2_vec[i];
	}
	else if(bch) {
		inf1.open("param_vec_BCH.txt"); for(i=0;i<74600;i++) inf1 >> param_vec[i];
		inf3.open("ncx2_vec_BCH.txt"); for(i=0;i<74600;i++) inf3 >> ncx2_vec[i];
	}
	
	for(j=0;j<m;j++) for(i=0;i<n;i++) inf >> H[j*n+i];
	
	//load Q-table
	for(i=0;i<S*num_cls;i++) for(j=0;j<cls_sz;j++) inf12 >> Q[i*cls_sz+j]; //contg. model-free
	for(i=0;i<S*num_cls;i++) for(j=0;j<cls_sz;j++) inf14 >> Q2[i*cls_sz+j]; //rand. model-free
	for(i=0;i<S*num_cls;i++) for(j=0;j<cls_sz;j++) inf15 >> Q3[i*cls_sz+j]; //opt. model-free
	for(i=0;i<S*num_cls;i++) for(j=0;j<cls_sz;j++) inf16 >> Q4[i*cls_sz+j]; //contg. model-based
	for(i=0;i<S*num_cls;i++) for(j=0;j<cls_sz;j++) inf17 >> Q5[i*cls_sz+j]; //rand. model-based
	for(i=0;i<S*num_cls;i++) for(j=0;j<cls_sz;j++) inf18 >> Q6[i*cls_sz+j]; //opt. model-based
	
	//for(i=0;i<M;i++) for(j=0;j<m;j++) inf9 >> G[i][j];
	for(i=0;i<I;i++) inf9 >> GI_arm[i];  

	cls_ind_mat_contg=new long*[num_cls]; for(i=0;i<num_cls;i++) cls_ind_mat_contg[i]=new long[cls_sz];
	cls_ind_mat_ran=new long*[num_cls]; for(i=0;i<num_cls;i++) cls_ind_mat_ran[i]=new long[cls_sz];
	cls_ind_mat_opt=new long*[num_cls]; for(i=0;i<num_cls;i++) cls_ind_mat_opt[i]=new long[cls_sz];

	if(AB) {
		inf13.open("clusters/cls_ind_mat_contg_AB.txt"); 
		for(i=0;i<num_cls;i++) for(j=0;j<cls_sz;j++) inf13 >> cls_ind_mat_contg[i][j]; //contiguous clustering
		
		inf10.open("clusters/cls_ind_mat_ran_AB.txt"); 
		for(i=0;i<num_cls;i++) for(j=0;j<cls_sz;j++) inf10 >> cls_ind_mat_ran[i][j]; //random clustering
		
		if(matnum==1) inf4.open("clusters/cls_ind_mat_opt_AB_1.txt"); 
		else if(matnum==2) inf4.open("clusters/cls_ind_mat_opt_AB_2.txt");
		else if(matnum==3) inf4.open("clusters/cls_ind_mat_opt_AB_3.txt");
		for(i=0;i<num_cls;i++) for(j=0;j<cls_sz;j++) inf4 >> cls_ind_mat_opt[i][j]; //optimized clustering
	}

	else if(!AB && !mac && !bch){
		inf13.open("clusters/cls_ind_mat_contg.txt"); 
		for(i=0;i<num_cls;i++) for(j=0;j<cls_sz;j++) inf13 >> cls_ind_mat_contg[i][j]; //for (3,6) LDPC
		
		inf10.open("clusters/cls_ind_mat_ran.txt"); //for (3,6) LDPC
		for(i=0;i<num_cls;i++) for(j=0;j<cls_sz;j++) inf10 >> cls_ind_mat_ran[i][j]; 

		if(matnum==1) inf4.open("clusters/cls_ind_mat_opt_1.txt"); 
		else if(matnum==2) inf4.open("clusters/cls_ind_mat_opt_2.txt");
		else if(matnum==3) inf4.open("clusters/cls_ind_mat_opt_3.txt");
		for(i=0;i<num_cls;i++) for(j=0;j<cls_sz;j++) inf4 >> cls_ind_mat_opt[i][j];
	}
	
	else if(mac) {		 
		inf10.open("clusters/cls_ind_mat_ran_mac_96_48.txt");
		for(i=0;i<num_cls;i++) for(j=0;j<cls_sz;j++) inf10 >> cls_ind_mat_ran[i][j]; 
		
		if(matnum==1) inf4.open("clusters/cls_ind_mat_opt_mac_96_48.txt");
		else if(matnum==2) inf4.open("clusters/cls_ind_mat_opt_mac_96_48_2.txt");
		else if(matnum==3) inf4.open("clusters/cls_ind_mat_opt_mac_96_48_3.txt");
		for(i=0;i<num_cls;i++) for(j=0;j<cls_sz;j++) inf4 >> cls_ind_mat_opt[i][j];
	}
	
	else if(bch) {		 
		inf10.open("clusters/cls_ind_mat_ran_BCH_63_51.txt");
		for(i=0;i<num_cls;i++) for(j=0;j<cls_sz;j++) inf10 >> cls_ind_mat_ran[i][j]; 

		inf4.open("clusters/cls_ind_mat_opt_BCH_63_51.txt");
		for(i=0;i<num_cls;i++) for(j=0;j<cls_sz;j++) inf4 >> cls_ind_mat_opt[i][j];
	}
	
	cout<<'\n'<<"cls_ind_mat_ran: "<<'\n'; for(i=0;i<num_cls;i++) {for(j=0;j<cls_sz;j++) cout<<cls_ind_mat_ran[i][j]<<" "; cout<<'\n';}
	//cout<<'\n'<<"cls_ind_mat_opt: "<<'\n'; for(i=0;i<num_cls;i++) {for(j=0;j<cls_sz;j++) cout<<cls_ind_mat_opt[i][j]<<" "; cout<<'\n';}
	
	/*for(i=0;i<S*num_cls;i++) for(j=0;j<cls_sz;j++) inf9 >> Q_cnt[i*cls_sz+j];
	for(i=0;i<S*num_cls;i++) 
		for(j=0;j<cls_sz;j++)
			if(Q_cnt[i*cls_sz+j]!=0)
				Q[i*cls_sz+j]/=Q_cnt[i*cls_sz+j];*/

	//cout<<'\n'<<"H: "<<'\n'; for(j=0;j<m;j++){for(i=0;i<n;i++) cout<<H[j*n+i]<<" "; cout<<'\n';}
	//cout<<'\n'<<"Q1: "<<'\n'; for(i=0;i<S*num_cls;i++) for(j=0;j<cls_sz;j++) cout<<Q[i*cls_sz+j]<<" ";
	//cout<<'\n'<<"G: "<<'\n'; for(i=0;i<M;i++) for(j=0;j<m;j++) cout<<G[i][j]<<" ";
	//cout<<'\n'<<"GI_arm: "<<'\n'; for(i=0;i<I;i++) cout<<GI_arm[i]<<" ";

	cout<<'\n'<<"m: "<<m<<" n: "<<n<<" R: "<<R<<" col_wt: "<<col_wt<<" row_wt: "<<row_wt<<" W: "<<W;  
	cout<<'\n'<<"cls_sz: "<<cls_sz<<" num_cls: "<<num_cls<<" S: "<<S<<" qstrt: "<<qstrt<<" gap: "<<gap<<" M: "<<M;

	//here n=WJp*p
	if(!BC) {
		n2=L*J*p*p; //length of codestream
		hshft=J*p*p*(mem+1); //amount of right shift of the window (length of a CW in the stream)
		vshft=gama*p*J*(mem+1); //amount of downward shift
	}
	else {
		n2=n;
		hshft=n;
		vshft=m;
	}

	cout<<'\n'<<"hshft: "<<hshft<<" vshft: "<<vshft<<" n2: "<<n2<<endl; 
	
	//2D arrays
	//x= new float[CW*n2]; //each row has a different CW stream
	x_hat= new float[CW*n2];
	y= new float[CW*n2];  
	LR= new float[CW*n];  //Lvalues in the window
	pLR= new float[CW*n]; 
	syn=new long[m];
	syn_cls=new long[cls_sz]; 
	syn2=new float[m];
	excl_cw=new long[CW];

	Gtemp=new float[M];
	Gmax=new float[m];
	

	vn_indx=new long[m];
	cn_indx=new long[m];
	indx=new long[cls_sz];
	indx_cls=new long[num_cls];
	indx_cls_CN=new long[num_cls];
	cls_s_crnt=new long[num_cls];
	cls_s_prev=new long[num_cls];
	pick_cls=new long[num_cls];

	indx2=new long[m];

	vns=new long[m*row_wt]; //all VNs of CN j in a row
	
	//3D arrays
	E_v_c= new float[m*row_wt*CW];
	E_c_v= new float[n*col_wt*CW];
	E_v_c_mu= new float[m*row_wt*CW];
	E_c_v_mu= new float[n*col_wt*CW];
	E_c_v_old= new float[n*col_wt*CW];
	res_c_v= new float[m*row_wt*CW];
	res_c_v_srtd= new float[m*row_wt*CW];
	res_mean= new float[I];
	res_cnt= new float[I];
	res_mean2= new float[I];
	res_cnt2= new float[I];
	
	cns=new long[n*col_wt]; //all CNs of a VN in a row

	//Eb/No values in dB (make sure SNR values match with main_RL
	//float EbNo_BC[]={0,0.6,1.2,1.4,1.8,2,2.2};
	//float EbNo_BC[]={0.5,0.75,1,1.25,1.5,1.75}; //for BCH
	float EbNo_BC[]={0.5,1,1.5,2,2.25,2.5}; //for (3,7) AB, (3,6), MKay, 
	
	num_dat=sizeof(EbNo_BC)/sizeof(EbNo_BC[0]);

	//Oc=new long*[num_dat]; for(i=0;i<num_dat;i++) Oc[i]=new long[I];
	//inf6.open("Oc.txt"); 
	//inf6.open("Oc_samp_5000.txt"); //for sample based residual
	//inf6.open("Oc_samp.txt"); 
	//for(i=0;i<num_dat;i++) for(j=0;j<I;j++) inf6 >> Oc[i][j];

	blk_err= new float[num_dat]; 
	blk_err2= new float[num_dat];
	blk_err3= new float[num_dat]; 
	blk_err4= new float[num_dat]; 
	blk_err5= new float[num_dat];
	blk_err6= new float[num_dat];
	blk_err7= new float[num_dat];
	blk_err8= new float[num_dat];

	bit_err= new float[num_dat]; 
	bit_err2= new float[num_dat]; 	
	bit_err3= new float[num_dat]; 
	bit_err4= new float[num_dat]; 
	bit_err5= new float[num_dat];
	bit_err6= new float[num_dat];
	bit_err7= new float[num_dat];
	bit_err8= new float[num_dat];

	ncv_vec=new float[num_dat];
	ncv_vec2=new float[num_dat];
	ncv_vec3=new float[num_dat];
	ncv_vec4=new float[num_dat];
	ncv_vec5=new float[num_dat];
	ncv_vec6=new float[num_dat];
	ncv_vec7=new float[num_dat];
	ncv_vec8=new float[num_dat];
	
	mean1_vec=new float[avg];
	mean2_vec=new float[avg];

	//***********************************************    test   ****************************************
	/*if(test) {
		meth=4;
		y[0*n2+0]=0.3; y[0*n2+1]=-0.37; y[0*n2+2]=0.5; y[0*n2+3]=-0.35; y[0*n2+4]=0.2; //1st noisy CW
	
		inf5.open("mat_test.txt"); for(j=0;j<m;j++) for(i=0;i<n;i++) inf5 >> H[j*n+i];
	
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
	
		varn=0.1;

		BP();	
	}*/
	//**************************************************************************************************


	//else {	
	if(BC) {
		//initializing
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
	}

	for(i=0;i<I;i++) mu[i]=mu_cnt[i]=mu_L[i]=mu_cnt_L[i]=0;
	for(i2=0;i2<num_dat;i2++) {
		
		varn=1/(2*R*pow(10,0.1*EbNo_BC[i2]));
		cout<<'\n'<<"i2, varn: "<<i2<<", "<<varn<<endl;
	
		//else varn=1/(2*R*pow(10,0.1*EbNo_kk[i2])); //if EbNo not in dB, then varn=1/(2*R*EbNo_kk[i2]);
		
		std::normal_distribution<float> dist(0,sqrt(varn));
		ncv=ncv2=ncv3=ncv4=ncv5=ncv6=ncv7=ncv8=0;
		dcw=dcw2=dcw3=dcw4=dcw5=dcw6=dcw7=dcw8=0;
		err=err2=err3=err4=err5=err6=err7=err8=biterr=biterr2=biterr3=biterr4=biterr5=biterr6=biterr7=biterr8=l1=l2=l3=0;
		while(/*(err<ev && i2<num_dat) || err2<ev ||*/ err3<ev || err4<ev || err5<ev /*|| err6<ev || err7<ev || err8<ev*/) { //take avg. of min. 3 blk. error events
			for(cw=0;cw<CW;cw++) 
				for(i=0;i<n2;i++) y[cw*n2+i]=1+dist(e2); //adding gaussian noise to all-zero CW
			
			//inf11.open("y.txt"); for(cw=0;cw<CW;cw++) for(i=0;i<n2;i++) inf11 >> y[cw*n2+i];
			//cw=0; cout<<'\n'<<"y: "; for(i=0;i<n2;i++) cout<<y[cw*n2+i]<<" ";

			//cw=0; ofstream outf; outf.open("y.txt"); for(i=0;i<n2;i++) outf<<y[cw*n2+i]<<" "; outf<<std::endl; outf.close();

			for(meth=1;meth<=8;meth++) {
				if(/*(meth==1 && err<ev && i2<num_dat) || (meth==2 && err2<ev) ||*/ (meth==3 && err3<ev) || (meth==4 && err4<ev)
					|| (meth==5 && err5<ev) /*|| (meth==6 && err6<ev) || (meth==7 && err7<ev) || (meth==8 && err8<ev)*/) {	

					num=BP(i2); 

					//else if(meth==2) l2+=num;
					//else if(meth==3) l3+=num;
					//cout<<'\n'<<"err: "<<err<<" err2: "<<err2<<" err3: "<<err3<<" err4: "<<err4<<" meth: "<<meth<<endl;
				}
			}
		}

		//trans=(call*n2/hshft)*CW; //there are n2/hshft cws in a stream. Each cw has length hshft bits in a stream of length n2 bits
		blk_err[i2]=err/dcw;
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

		blk_err7[i2]=err7/dcw7;		
		bit_err7[i2]=biterr7/(hshft*dcw7);
		ncv_vec7[i2]=ncv7/(1*dcw7); 
		//cout<<'\n'<<"ncv7: "<<ncv7<<" dcw7: "<<dcw7<<endl;
		
		blk_err8[i2]=err8/dcw8;		
		bit_err8[i2]=biterr8/(hshft*dcw8);
		ncv_vec8[i2]=ncv8/(1*dcw8); 	
		 
	}
	//}

	cout<<'\n'<<"EbNo_BC: "; for(i=0;i<num_dat;i++) cout<<EbNo_BC[i]<<" "; cout<<'\n'<<endl;
	//cout<<'\n'<<"err: "<<err<<" dcw: "<<dcw<<'\n';
	
	//cout<<"blk_err: "; for(i=0;i<num_dat;i++) cout<<blk_err[i]<<" "; cout<<'\n';
	cout<<"bit_err: "; for(i=0;i<num_dat;i++) cout<<bit_err[i]<<" "; cout<<'\n';
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
	//cout<<"blk_err3: "; for(i=0;i<num_dat;i+rep[0]=-1; rep[1]=10; rep[2]=20; rep[3]=30;+) cout<<blk_err3[i]<<" "; cout<<'\n';
	cout<<"bit_err4: "; for(i=0;i<num_dat;i++) cout<<bit_err4[i]<<" "; cout<<'\n';
	cout<<"ncv_vec4: "; for(i=0;i<num_dat;i++) cout<<ncv_vec4[i]<<" "; cout<<'\n'<<endl;

	//cout<<'\n'<<"err5: "<<err5<<" dcw5: "<<dcw5<<'\n';
	//cout<<"blk_err5: "; for(i=0;i<num_dat;i++) cout<<blk_err5[i]<<" "; cout<<'\n';
	cout<<"bit_err5: "; for(i=0;i<num_dat;i++) cout<<bit_err5[i]<<" "; cout<<'\n';
	cout<<"ncv_vec5: "; for(i=0;i<num_dat;i++) cout<<ncv_vec5[i]<<" "; cout<<'\n'<<endl;

	//cout<<"blk_err6: "; for(i=0;i<num_dat;i++) cout<<blk_err6[i]<<" "; cout<<'\n';
	cout<<"bit_err6: "; for(i=0;i<num_dat;i++) cout<<bit_err6[i]<<" "; cout<<'\n';
	cout<<"ncv_vec6: "; for(i=0;i<num_dat;i++) cout<<ncv_vec6[i]<<" "; cout<<'\n'<<endl;

	cout<<"bit_err7: "; for(i=0;i<num_dat;i++) cout<<bit_err7[i]<<" "; cout<<'\n';
	cout<<"ncv_vec7: "; for(i=0;i<num_dat;i++) cout<<ncv_vec7[i]<<" "; cout<<'\n'<<endl;
	
	cout<<"bit_err8: "; for(i=0;i<num_dat;i++) cout<<bit_err8[i]<<" "; cout<<'\n';
	cout<<"ncv_vec8: "; for(i=0;i<num_dat;i++) cout<<ncv_vec8[i]<<" "; cout<<'\n'<<endl;

	//for(i=0;i<I;i++) {res_mean[i]/=res_cnt[i]; res_mean2[i]/=res_cnt2[i];} 
	//cout<<"res_mean: "; for(i=0;i<I;i++) cout<<res_mean[i]<<" "; cout<<'\n'<<endl;
	//cout<<"res_mean2: "; for(i=0;i<I;i++) cout<<res_mean2[i]<<" "; cout<<'\n'<<endl;


	//for(i=0;i<I;i++) mu[i]/=mu_cnt[i]; 
	//cout<<"mean m_{c_v}: "; for(i=0;i<I;i++) cout<<mu[i]<<" "; cout<<'\n'<<endl;

	//for(i=0;i<I;i++) mu_L[i]/=mu_cnt_L[i]; 
	//cout<<"mean mu_{Li_hat}: "; for(i=0;i<I;i++) cout<<mu_L[i]+col_wt*mu[i]<<" "; cout<<'\n'<<endl;

	//cout<<"mean K*mu_{Li_hat}: "; for(i=0;i<I;i++) cout<<row_wt*(mu_L[i]+col_wt*mu[i])<<" "; cout<<'\n'<<endl;
	//cout<<"std. dev. 2*K*mu_{Li_hat}: "; for(i=0;i<I;i++) cout<<sqrt(2*row_wt*(mu_L[i]+col_wt*mu[i]))<<" "; cout<<'\n'<<endl;

	//ofstream outf3; filename="Pe"+func(fn)+".txt"; outf3.open(filename.c_str()/*,fstream::app*/); for(i2=0;i2<num_dat;i2++) outf3<<blk_err[i2]<<rep[0]=-1; rep[1]=10; rep[2]=20; rep[3]=30;" ";  outf3<<std::endl; outf3.close();
	cout<<'\n';
    	printf("executed in: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
	cout<<'\n';


}





