
//does RL
//for Matlab stuff, first do: $ nano /home/salman/.bashrc
//then paste: $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/salman/extern/bin/glnxa64/:matlabroot/sys/os/glnxa64
//to compile: $ g++ -o main -std=c++11 -I /home/salman/extern/include/ -L /home/salman/extern/bin/glnxa64/ -pthread main_RL.cpp -lMatlabDataArray -lMatlabEngine
//the matlabroot is /home/salman

//for python stuff (first check python version):
//g++ -o main main_RL.cpp -I/usr/include/python3.8 -L/usr/lib/python3.8/config-3.8-x86_64-linux-gnu -lpython3.8 

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
//#include "MatlabEngine.hpp"
//#include "MatlabDataArray.hpp"

using namespace std;


// Pass vector containing MATLAB data array scalar
/*using namespace matlab::engine;
// Start MATLAB engine synchronously
std::unique_ptr<MATLABEngine> matlabPtr = startMATLAB();
// Create MATLAB data array factory
matlab::data::ArrayFactory factory;*/

float *H,*H2,*x,*x_hat,*y,P,*LR,*pLR,*E_c_v,*E_v_c,*E_c_v_mu,*E_v_c_mu,*H_D,*y_D,*LR_D,*res_c_v,*res_c_v_srtd,err,err2,err3,biterr,biterr2,biterr3,varn,ncv,ncv2,ncv3,*ncv_vec,*ncv_vec2,*ncv_vec3,*Q,*Q_temp,**Q2,**Q2_inv,
**rw_vec,**rw_vec_srtd,*d_vec,*b_vec,**tran_prob,**G,*GI2,*Gmax,*divv,qstrt,gap,*Q_cnt,*zeta,*Gtemp,*rep,*diff,**R_sv,*param_vec,*val_vec,*ncx2_vec;
long m,n,n2,*vns,*cns,*vns_D,*cns_D,*excl_cw,BC,S,*indx,*indx2,*states,*syn,*syn_cls,*syn_old,***syn_sv,***syn_sv2,*batch_sz,**a_sv,**rw_cnt,*C_set,*GI_arm,*lmax_vec,num_lmax,**cls_ind_mat;
long col_wt,row_wt,num,dcw,dcw2,dcw3,hshft,vshft,mem,W,call,flg,num_cls,L=99,cnt_st=0,**Oc,dc,dv,*cn_indx,bch_sz_max,DeepRL,M,cls_sz,boxplus=0; //row_wt is the max, row weight in case of SC codes

#define CW 1 //1 for model based, 25 for model free schemes
#define alpha 0.1
#define gamma 0.9
#define eps 0.6
#define rnd_max 25 //(25 for model-free schemes)

string func(long n) {
	stringstream result;
	result << n;
	return result.str();
}


//#include "RL.cpp"
//#include "RL2.cpp" //clustered Q-learning
//#include "RL3.cpp" //for GIs	
#include "RL4.cpp" //for optimimzed or random clustering	
//#include "MP_RL3.cpp" //RL for MAB-NS-TS-1
//#include "recursion_Git.cpp"

int main() {	
	srand(time(0));	
	clock_t tStart = clock();	
	int modl; long i,i2,k,j,i3,cw_cnt=0,matnum;
	long fn,num_dat,cnt,cw,gama,p,J,test=0,num,l1,l2,l3,load,tmp2,samps=1000000,AB; //with 5 SNR points, CW=25, samps=10000 creates 1.25*10^6 training examples
	bch_sz_max=samps; 
	
	float *blk_err,*blk_err2,*blk_err3,*bit_err,*bit_err2,*bit_err3,R,tmp;  
	//cout<<'\n'<<" device no: "; cin>>fn; ///////////////////////////////// 
	
	cout<<'\n'<<"fn: "; cin>>fn;
	cout<<'\n'<<"matnum?: "; cin>>matnum;
	//cout<<'\n'<<"load Q table: "; cin>>load; //1 loads 0 doesn't
	load=0;
	cout<<'\n'<<"model based?: "; cin>>modl;  //1 yes, 0 no
	//cout<<'\n'<<"DeepRL?: "; cin>>DeepRL;  //1 yes, 0 no
	DeepRL=0;
	  
	
	BC=1;

	//construct a random generator engine:
	std::random_device rd;
    std::mt19937 e2(rd());

	ifstream inf,inf1,inf2,inf3,inf4,inf5,inf6,inf7,inf8,inf9,inf10,inf11;  


	if(fn==1 || fn==2 || fn==3) { //fn=1,2,3 for contg, rand, optimized
		if(matnum==1) inf.open("matrices/mat_n196_1.txt"); 
		else if(matnum==2) inf.open("matrices/mat_n196_2.txt");
		else if(matnum==3) inf.open("matrices/mat_n196_3.txt"); 
		AB=0;
		m=98; n=196;
		cls_sz=7; M=4;
		col_wt=3; //max. col. wt.
		row_wt=6; //max. row. wt.
		dv=3; dc=6;
	} //(3,6) regular

	if(fn==4 || fn==5 || fn==6) { //fn=4,5,6 for contg, rand, optimized
		if(matnum==1) inf.open("matrices/mat_n196_AB_1.txt"); 
		if(matnum==2) inf.open("matrices/mat_n196_AB_2.txt"); 
		if(matnum==3) inf.open("matrices/mat_n196_AB_3.txt"); 
		AB=1;
		m=84; n=196;
		cls_sz=7; M=4;
		col_wt=3; //max. col. wt.
		row_wt=7; //max. row. wt.
		dv=3; dc=7;
	} //(3,7) AB

	if(fn==7|| fn==8) {
		if(matnum==1) inf.open("matrices/mat_BCH_63_51.txt"); 
		m=12; n=63;
		cls_sz=6; M=4;
		col_wt=8; //max. col. wt.
		row_wt=28; //max. row. wt.
		dv=col_wt; dc=row_wt;
	} //(63,51) BCH

	if(fn==9|| fn==10) {

		if(matnum==1) inf.open("matrices/mackay_96_48.txt");
		else if(matnum==2) inf.open("matrices/mackay_96_48_2.txt");
		else if(matnum==3) inf.open("matrices/mackay_96_48_3.txt"); 
		//inf.open("mat.txt"); matnum=1;
		m=48; n=96;
		cls_sz=6; M=4;
		col_wt=3; //max. col. wt.
		row_wt=6; //max. row. wt.
		dv=col_wt; dc=row_wt;
	} //(96,48) Mackay LDPC



	if(test) {
		m=4; n=5; col_wt=3; row_wt=3;
		mem=1; gama=col_wt;
	}

	if(BC || test) {L=3;W=3;}

	S=pow(M,cls_sz); //no. of syndromes for each cluster
	num_cls=m/cls_sz; //number of clusters

	R=1-float(m)/float(n); //code rate

	H= new float[m*n]; //H is the window. Since windows are identical, only the window Tanner graph is needed. The Lvalues "move" across the window
	H2= new float[m*n];
	ofstream outf,outf2,outf3,outf4,outf5,outf6,outf7,outf8,outf9,outf10; string filename;

	//inf.open("H_3_7_24_5_sc_rnd1.txt"); for(j=0;j<m;j++) for(i=0;i<n;i++) inf >> H[j*n+i]; //upload the parity-check matrix
	//inf.open("mat_p5_sc.txt"); for(j=0;j<m;j++) for(i=0;i<n;i++) inf >> H[j*n+i];


	for(j=0;j<m;j++) for(i=0;i<n;i++) inf >> H[j*n+i];
	//cout<<'\n'<<"H: "<<'\n'; for(j=0;j<m;j++){for(i=0;i<n;i++) {cout<<H[j*n+i]<<" "; if(i>0 && (i+1)%(row_wt)==0) cout<<"  ";} cout<<'\n'; if(j>0 && (j+1)%(row_wt)==0) cout<<'\n';}

	//col_wt=0;
	//for(j=0;j<m;j++) if(H[j*n+0]) col_wt++;

	cout<<'\n'<<"m: "<<m<<" n: "<<n<<" R: "<<R<<" col_wt: "<<col_wt<<" row_wt: "<<row_wt<<" W: "<<W;  
	cout<<'\n'<<"cls_sz: "<<cls_sz<<" num_cls: "<<num_cls<<" S: "<<S<<" M: "<<M;

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
	syn= new long[m];
	syn_old= new long[m];
	syn_cls= new long[cls_sz];
	
	if(DeepRL) {
		syn_sv=new long**[num_cls]; 
		for(i=0;i<num_cls;i++) {
			syn_sv[i]=new long*[bch_sz_max];
			for(j=0;j<bch_sz_max;j++)
				syn_sv[i][j]=new long[cls_sz];
		}
	
		syn_sv2=new long**[num_cls]; 
		for(i=0;i<num_cls;i++) {
			syn_sv2[i]=new long*[bch_sz_max];
			for(j=0;j<bch_sz_max;j++)
				syn_sv2[i][j]=new long[cls_sz];
		}
	
		a_sv=new long*[num_cls]; for(i=0;i<num_cls;i++) a_sv[i]=new long[bch_sz_max];
		R_sv=new float*[num_cls]; for(i=0;i<num_cls;i++) R_sv[i]=new float[bch_sz_max];
	}
	
	batch_sz=new long[num_cls];
	excl_cw= new long[CW];
	indx= new long[n];
	indx2= new long[M];
	rep= new float[M];
	diff= new float[M];
	param_vec= new float[74600];
	ncx2_vec= new float[74600];
	val_vec= new float[100];
	

	if(!AB && (fn==1 || fn==2 || fn==3)) {
		inf1.open("param_vec_36nonAB.txt"); for(i=0;i<74600;i++) inf1 >> param_vec[i]; //param is generated using main_RL and ncx2_vec is generated using non_cntrl_chi.m 
		inf3.open("ncx2_vec_36nonAB.txt"); for(i=0;i<74600;i++) inf3 >> ncx2_vec[i];
	}
	else if(!AB && (fn==9 || fn==10)) {
		inf1.open("param_vec_MKay.txt"); for(i=0;i<74600;i++) inf1 >> param_vec[i]; //param is generated using main_RL and ncx2_vec is generated using non_cntrl_chi.m 
		inf3.open("ncx2_vec_MKay.txt"); for(i=0;i<74600;i++) inf3 >> ncx2_vec[i];
	}
	else if(AB) {
		inf1.open("param_vec_37AB.txt"); for(i=0;i<74600;i++) inf1 >> param_vec[i];
		inf3.open("ncx2_vec_37AB.txt"); for(i=0;i<74600;i++) inf3 >> ncx2_vec[i];
	}

	vns=new long[m*row_wt]; //all VNs of CN j in a row
	
	//3D arrays
	E_v_c= new float[m*row_wt*CW];
	E_c_v= new float[n*col_wt*CW];
	E_v_c_mu= new float[m*row_wt*CW];
	E_c_v_mu= new float[n*col_wt*CW];
	res_c_v= new float[m*row_wt*CW];
	res_c_v_srtd= new float[m*row_wt*CW];
	
	cns=new long[n*col_wt]; //all CNs of a VN in a row

	if(!DeepRL) {
		Q=new float[S*num_cls*cls_sz]; 
		Q_cnt=new float[S*num_cls*cls_sz];
		Q_temp=new float[cls_sz];
	}
	
	Q2=new float*[M]; for(i=0;i<M;i++) Q2[i]=new float[M]; 
	tran_prob=new float*[M]; for(i=0;i<M;i++) tran_prob[i]=new float[M]; 
	Q2_inv=new float*[M]; for(i=0;i<M;i++) Q2_inv[i]=new float[2*M]; 
	
	rw_vec=new float*[m]; for(i=0;i<m;i++) rw_vec[i]=new float[M];
	rw_vec_srtd=new float*[m]; for(i=0;i<m;i++) rw_vec_srtd[i]=new float[M];
	G=new float*[M]; for(i=0;i<M;i++) G[i]=new float[m];
	cls_ind_mat=new long*[num_cls]; for(i=0;i<num_cls;i++) cls_ind_mat[i]=new long[cls_sz];

	Gmax=new float[m];
	GI2=new float[m*M];
	GI_arm=new long[m*M];
	rw_cnt=new long*[m]; for(i=0;i<m;i++) rw_cnt[i]=new long[M];
	d_vec=new float[M];
	b_vec=new float[M];
	divv=new float[M];
	states= new long[M];
	C_set= new long[M];
	Gtemp=new float[M];
	cn_indx=new long[m];


	if(load && !DeepRL) {
		filename="Q_M_"+func(fn)+".txt"; inf10.open(filename.c_str()); for(i=0;i<S*num_cls;i++) for(j=0;j<cls_sz;j++) inf10 >> Q[i*cls_sz+j];
		//filename="Q_8a.txt"; inf10.open(filename.c_str()); for(i=0;i<S*num_cls;i++) for(j=0;j<cls_sz;j++) inf10 >> Q[i*cls_sz+j];
	} //load Q table from prev. training

	if(fn==1) inf9.open("clusters/cls_ind_mat_contg.txt"); 
	else if(fn==2) inf9.open("clusters/cls_ind_mat_ran.txt"); 
	else if(fn==3 && matnum==1) inf9.open("clusters/cls_ind_mat_opt_1.txt"); //for 3,6 LDPC
	else if(fn==3 && matnum==2) inf9.open("clusters/cls_ind_mat_opt_2.txt"); //for 3,6 LDPC
	else if(fn==3 && matnum==3) inf9.open("clusters/cls_ind_mat_opt_3.txt"); //for 3,6 LDPC

	if(fn==1 || fn==2 || fn==3) {
		if(!boxplus) {rep[0]=-8.081; rep[1]=-0.423; rep[2]=1.27; rep[3]=9.82;}
		else {rep[0]=-2.901; rep[1]=-0.376; rep[2]=0.317; rep[3]= 1.12336;} //for box-plus syndrome
	} 

	if(fn==4) inf9.open("clusters/cls_ind_mat_contg_AB.txt"); 
	else if(fn==5) inf9.open("clusters/cls_ind_mat_ran_AB.txt"); 
	else if(fn==6 && matnum==1) inf9.open("clusters/cls_ind_mat_opt_AB_1.txt"); //for 3,7 AB 
	else if(fn==6 && matnum==2) inf9.open("clusters/cls_ind_mat_opt_AB_2.txt"); //for 3,7 AB
	else if(fn==6 && matnum==3) inf9.open("clusters/cls_ind_mat_opt_AB_3.txt"); //for 3,7 AB
	
	if(fn==4 || fn==5 || fn==6) {rep[0]=-9.805; rep[1]=-0.189; rep[2]=1.77; rep[3]=13.36;}  

	//BCH
	if(fn==7) inf9.open("clusters/cls_ind_mat_ran_BCH_63_51.txt");
	else if(fn==8) inf9.open("clusters/cls_ind_mat_opt_BCH_63_51.txt");
	if(fn==7 || fn==8) {
		if(M==4) {rep[0]=74.61; rep[1]=104.6; rep[2]=119.58; rep[3]=131.56;}
		else {rep[0]=20.8; rep[1]=44.2; rep[2]=64; rep[3]=78.8; rep[4]=92.7; rep[5]=109; rep[6]=124; rep[7]=142;}
	} 

	//Mackay
	if(fn==9) inf9.open("clusters/cls_ind_mat_ran_mac_96_48.txt");
	else if(fn==10 && matnum==1) inf9.open("clusters/cls_ind_mat_opt_mac_96_48.txt");
	else if(fn==10 && matnum==2) inf9.open("clusters/cls_ind_mat_opt_mac_96_48_2.txt");
	else if(fn==10 && matnum==3) inf9.open("clusters/cls_ind_mat_opt_mac_96_48_3.txt");
	if(fn==9 || fn==10) {
		if(M==4) {rep[0]=-3.236; rep[1]=13.7; rep[2]=32.4; rep[3]=57.3;} 
		//else{inf2.open("quantization/rep_mac.txt"); for(i=0;i<M;i++) inf2>>rep[i];}
		else{rep[0]=-3.236; rep[1]=9.65; rep[2]=22.2; rep[3]=34.4; rep[4]=46.8; rep[5]=61; rep[6]=76.3; rep[7]=96.8;}
	}  
	
	cout<<'\n'<<"rep: "; for(j=0;j<M;j++) cout<<rep[j]<<" "; cout<<'\n';

	for(i=0;i<num_cls;i++) for(j=0;j<cls_sz;j++) inf9 >> cls_ind_mat[i][j];
	cout<<'\n'<<"cls_ind_mat: "<<'\n'; for(i=0;i<num_cls;i++) {for(j=0;j<cls_sz;j++) cout<<cls_ind_mat[i][j]<<" "; cout<<'\n';}

	float EbNo_BC[]={0,0.6,1.2,1.4,1.8}; //,2,2.2};
	//float EbNo_BC[]={4}; float EbNo_BC[]={1.8};

	num_dat=sizeof(EbNo_BC)/sizeof(EbNo_BC[0]);
	
	Oc=new long*[num_dat]; for(i=0;i<num_dat;i++) Oc[i]=new long[rnd_max];

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

	for(i=0;i<m;i++) 
		for(j=0;j<M;j++) 
			rw_vec[i][j]=rw_cnt[i][j]=0; //refresh



lmax_vec=new long[M];
lmax_vec[0]=5; lmax_vec[1]=10; lmax_vec[2]=15; lmax_vec[3]=20;
num_lmax=sizeof(lmax_vec)/sizeof(lmax_vec[0]);
zeta=new float[num_lmax];


cw_cnt=0; 


//for(i3=0;i3<batches;i3++) {
	while(cw_cnt<samps) {
		for(i2=0;i2<num_dat;i2++) {
			//cout<<'\n'<<"i2: "<<i2<<endl;
			varn=1/(2*R*pow(10,0.1*EbNo_BC[i2])); //if EbNo not in dB, then varn=1/(2*R*EbNo[i2]);

			//cout<<'\n'<<"std_dev: "<<sqrt(varn)<<endl;
			std::normal_distribution<float> dist(0,sqrt(varn));
			for(cw=0;cw<CW;cw++) for(i=0;i<n2;i++) y[cw*n2+i]=1+dist(e2); //adding gaussian noise to all-zero CW

			//inf11.open("y.txt"); for(cw=0;cw<CW;cw++) for(i=0;i<n2;i++) inf11 >> y[cw*n2+i];
			//cw=0; cout<<'\n'<<"y: "; for(i=0;i<n2;i++) cout<<y[cw*n2+i]<<" ";
		
			RL(modl); 
		}
		cw_cnt++;
		cout<<'\n'<<"cw_cnt: "<<cw_cnt<<endl;

		//if(cw_cnt==10000000 || cw_cnt==20000000) {
			//filename="Q_"+func(fn)+".txt"; 
			//outf.open(filename.c_str()); for(i=0;i<S*num_cls;i++) for(j=0;j<cls_sz;j++) outf<<Q[i*cls_sz+j]<<" ";  outf<<std::endl; outf.close();
		//}
	}
	
	//create training dataset
	if(DeepRL) { 
		outf4.open("DeepRL/batch_sz.txt"); 
		for(i=0;i<num_cls;i++) {
			filename="DeepRL/syn_cls_"+func(i)+".txt"; 
			outf3.open(filename.c_str()); 
			filename="DeepRL/syn_cls_new_"+func(i)+".txt"; 
			outf5.open(filename.c_str()); 
			filename="DeepRL/R_"+func(i)+".txt"; 
			outf6.open(filename.c_str()); 
			filename="DeepRL/a_"+func(i)+".txt"; 
			outf7.open(filename.c_str()); 
			for(k=0;k<batch_sz[i];k++) {
				for(j=0;j<cls_sz;j++) {
					outf3<<syn_sv[i][k][j]<<" ";  
					outf5<<syn_sv2[i][k][j]<<" ";  
				}
				outf6<<R_sv[i][k]<<" "; 
				outf7<<a_sv[i][k]<<" "; 
			}
			outf3<<std::endl; outf3.close();
			outf5<<std::endl; outf5.close();
			outf6<<std::endl; outf6.close();
			outf7<<std::endl; outf7.close();
			
			outf4<<batch_sz[i]<<" ";  
		}
		outf4<<std::endl; outf4.close();
	}
//}


	//deciding which arms to play via GIs
	/*for(j=0;j<m;j++) 
		for(i=0;i<M;i++) {
			GI2[j*M+i]=G[i][j];
			GI_arm[j*M+i]=j;
		}

	//cout<<'\n'<<"GI_arm: "<<'\n'; for(i=0;i<m*M;i++) cout<<GI_arm[i]<<" ";
	for(i=0;i<m*M;i++) 
		for(i2=i+1;i2<m*M;i2++) 
			if(GI2[i]<GI2[i2]) {   
				tmp=GI2[i];
				GI2[i]=GI2[i2];
				GI2[i2]=tmp;
	
				tmp2=GI_arm[i];
				GI_arm[i]=GI_arm[i2];
				GI_arm[i2]=tmp2;
			}*/

	//cout<<'\n'<<"GI_arm: "<<'\n'; for(i=0;i<m*M;i++) cout<<GI_arm[i]<<" ";


	if(!DeepRL) {
		if(!modl) filename="Q"+func(fn)+"_"+func(matnum)+".txt";
		else filename="Q"+func(fn)+"_"+func(matnum)+"_M.txt"; 
		outf.open(filename.c_str()); for(i=0;i<S*num_cls;i++) for(j=0;j<cls_sz;j++) outf<<Q[i*cls_sz+j]<<" ";  outf<<std::endl; outf.close();
	
		//outf2.open("Oc.txt"); 
		//outf2.open("Oc_samp.txt"); //for residual sampling case
		//for(i=0;i<num_dat;i++) for(j=0;j<rnd_max;j++) outf2<<Oc[i][j]<<" "; outf2<<std::endl; outf2.close();
	}
	
	/*outf3.open("GI_arm.txt"); for(i=0;i<M*m;i++) outf3<<GI_arm[i]<<" "; outf3<<std::endl; outf3.close();
	outf4.open("qstrt.txt"); outf4<<qstrt;  outf4<<std::endl; outf4.close();
	outf5.open("gap.txt"); outf5<<gap;  outf5<<std::endl; outf5.close();
	outf6.open("M.txt"); outf6<<M;  outf6<<std::endl; outf6.close();*/


	cout<<'\n';
    	printf("executed in: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
	cout<<'\n'<<"fn "<<fn;


}




/*//computing Gittins index of each arm/CN (Varaiya, Walrand, Buyui3oc)
	//estimate of state trans_prob
	for(i=0;i<M;i++) 
		for(j=0;j<M;j++) 
			tran_prob[i][j]/=cnt_st;

	//estimate of reward function for each arm (CN) state
	for(i=0;i<m;i++) 
		for(j=0;j<M;j++) 
			rw_vec[i][j]/=rw_cnt[i][j];
	
	cout<<"\n\tran_prob\n\n";
   	for(i=0;i<M;i++) {
     		for(j=0;j<M;j++)
        		cout<<"\t"<<tran_prob[i][j];
      			cout<<"\n";
    		}

	for(j=0;j<m;j++) {
		cnt=0;
		//update reward vector for each CN
		for(i=0;i<M;i++) //M is no. of states of a bandit process
			rw_vec_srtd[j][i]=rw_vec[j][i]; 
		//sort reward vec
		for(i=0;i<M;i++) states[i]=i;
		for(i=0;i<M;i++) {
			for(i2=i+1;i2<M;i2++) {
				if(rw_vec_srtd[j][i]<rw_vec_srtd[j][i2]) {   
					tmp=rw_vec_srtd[j][i];
					rw_vec_srtd[j][i]=rw_vec_srtd[j][i2];
					rw_vec_srtd[j][i2]=tmp;

					tmp2=states[i];
					states[i]=states[i2];
					states[i2]=tmp2;
				}			
			}
		}

		//rw_vec_srtd[j][0] highest reward
		//states[0] state with highest reward

		G[j][states[0]]=rw_vec_srtd[j][0]; //state of a bandit process with highest reward will have G.I.=reward
		C_set[cnt]=states[0]; //continuation set
		cnt++; //no. of states whose G.I. has been found
		
		recursion_Git(j,cnt); //for other states G.I. is found via recursion

	}
		

	//cout<<"\n\trw_vec\n\n";
   		//for(i=0;i<m;i++) {
     			//for(j=0;j<M;j++)
        			//cout<<"\t"<<rw_vec[i][j];
      			//cout<<"\n";
    	//}

	//cout<<"\n\trw_cnt\n\n";
   		//for(i=0;i<m;i++) {
     			//for(j=0;j<M;j++)
        			//cout<<"\t"<<rw_cnt[i][j];
      			//cout<<"\n";
    	//}

				


	//cout<<"\n\tGI2\n\n"; for(i=0;i<m*M;i++) cout<<"\t"<<GI2[i];
	//cout<<"\n\tGI2_arm\n\n"; */


