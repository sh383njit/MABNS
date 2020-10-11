//matrix generator: L must be equal to mem+1 for generating terminally lifted replicas
//for BC generation, L=1 and J=(mem+1)*J_sc. No need to call matlab
//Generates terminally lifted AB matrix from H(3,p)

#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <time.h>
#include <omp.h>
#include <fstream>
#include <sstream>
using namespace std;

int **H_proto,**H,**Hsc,**Hsc_sv,**Hsc_proto,**H0,**H1,**I,**I2,**I3,**B,**sigma,**sigma2,**sigma3,**tau2,**Hsc_2,**Hsc_3,**lsft_mat,*E; 
int gama=3,p=7,q,m2,n2,J=20,m,n,L=1,opt,randm;  //H has total m rows and n cols

//params for find_cycl4
long par_num=0,chld_num=0,tcyc=pow(10,5),revs=pow(10,5),abs_64=0; //max. no. of reverses in a rnd can be 12
int *snk_r,*snk_c,**fbd_r,**fbd_c,*par_r,*par_c,**chld_r,**chld_c,*col_wt,**cyc_sv_r,**cyc_sv_c,**cyc_sv_rg,**cyc_sv_cg,*cyc_sel_r,*cyc_sel_c,**cyc_param,**cyc_cg,**cyc_r_nel,**cyc_c_nel; 
int strt_r,strt_c,stp_r,stp_c,row_strt,col_strt,row_stp,col_stp,*edge_r,*edge_c; //H has total m rows and n cols
long while_tot=0,rev_tot=0,tot_r,tot_c,sum=0,ccnt=0,ccnt_old=0,num_edgs=pow(10,5),ge=1,mem,s,cyc_len=6; 

#include "H_matb.cpp"
#include "Hsc2_AB.cpp"
#include "Hsc_matb.cpp"
#include "Hsc_mat2.cpp"
#include "main_cyc4.cpp"

string func(int n) 
{
    stringstream result;
    result << n;
    return result.str();
}

int main() {	
	
	srand(time(0));	
	double t2,t1= time(NULL);	
	long i,i5,i2,i3,i4,j,j2,k,l,count=0,A,*outp; 
	int f=1,totf=1,param,ccnt2=0,flg=0,rr,cnt,ncyc=0,Allison=1,unc=0,blks,ers_r=-1,rf,glob=0;   
	
	ifstream inf4,inf5,inf6,inf7,inf10,inf11,inf12; 

	/*cout<<'\n'<<" f: "; cin>>f; 
	if(f==1){
		cout<<'\n'<<" p: "; 
		cin>>p; cout<<'\n'<<" gama: "; 
		cin>>gama; cout<<'\n'<<" sub-code: "; 
		cin>>s; cout<<'\n'<<" J: "; 
		cin>>J; 
		L=1;
	} 
	else {inf6.open("p.txt"); inf6>>p; inf10.open("gama.txt"); inf10>>gama; inf12.open("s.txt"); inf12>>s;}*/

	q=p; m2=gama*q; n2=p*q;	

	//if(s>0) rf=1; else {cout<<'\n'<<" read file: "; cin>>rf;} //1 will read matrix from mat.txt, 0 won't
	//1 for generating the BC, 2 for generating the SC code	
	if(f==2) {
		cout<<'\n'<<" uncoupled? "; cin>>unc; //1 means uncoupled, 0 is not 
		cout<<'\n'<<" mem: "; cin>>mem; 
		cout<<'\n'<<" L: "; cin>>L; 
		cout<<'\n'<<" J: "; cin>>J;
		if(J>1) {
			cout<<'\n'<<" optimize t-lift?: "; cin>>opt; //1 optimizes for second lift, 0 does not
			if(!opt) {cout<<'\n'<<" randm t-lift?: "; cin>>randm;} //one does random terminal lift, 0 lifts using previously saved circulant shifts
		}
		else opt=0;
	} 

	//inf11.open("ge.txt"); inf11>>ge;

	m=m2; n=n2;
	
	if(!ge && f==2) {if(Allison) {L++; tot_r=L*m;} else tot_r=(L+1)*m;}
	else if(f==2) {tot_r=(L+mem)*m; tot_c=L*n;} //for SC protograph
	else if(f==1) {tot_r=m; tot_c=n;}

	E=new int[gama]; 
	H_proto=new int*[m2]; for(i5=0;i5<m2;i5++) H_proto[i5]=new int[n2]; //m rows n cols
	H=new int*[m]; for(i5=0;i5<m;i5++) H[i5]=new int[n]; //m rows n cols
	H0=new int*[m]; for(i5=0;i5<m;i5++) H0[i5]=new int[n]; H1=new int*[m]; for(i5=0;i5<m;i5++) H1[i5]=new int[n]; 
	Hsc_proto=new int*[tot_r]; for(i5=0;i5<tot_r;i5++) Hsc_proto[i5]=new int[tot_c]; 
	Hsc=new int*[tot_r*J]; for(i5=0;i5<tot_r*J;i5++) Hsc[i5]=new int[tot_c*J]; 
	Hsc_sv=new int*[tot_r*J]; for(i5=0;i5<tot_r*J;i5++) Hsc_sv[i5]=new int[tot_c*J]; 
	Hsc_2=new int*[tot_r]; for(i5=0;i5<tot_r;i5++) Hsc_2[i5]=new int[tot_c]; 
	Hsc_3=new int*[tot_r]; for(i5=0;i5<tot_r;i5++) Hsc_3[i5]=new int[tot_c];
	I=new int*[q]; for(i5=0;i5<q;i5++) I[i5]=new int[q];
	I2=new int*[L]; for(i5=0;i5<L;i5++) I2[i5]=new int[L];
	tau2=new int*[L]; for(i5=0;i5<L;i5++) tau2[i5]=new int[L];
	sigma2=new int*[L]; for(i4=0;i4<L;i4++) sigma2[i4]=new int[L];
	B=new int*[gama]; for(i5=0;i5<gama;i5++) B[i5]=new int[p];
	sigma=new int*[q]; for(i5=0;i5<q;i5++) sigma[i5]=new int[q];
	sigma3=new int*[J]; for(i4=0;i4<J;i4++) sigma3[i4]=new int[J];
	lsft_mat=new int*[tot_r]; for(i5=0;i5<tot_r;i5++) lsft_mat[i5]=new int[tot_c];

	edge_r=new int[num_edgs]; edge_c=new int[num_edgs]; outp=new long[12]; 
	par_r=new int[cyc_len]; par_c=new int[cyc_len]; 
	snk_r= new int[cyc_len]; snk_c= new int[cyc_len]; 
	chld_r=new int*[tcyc]; for(i4=0;i4<tcyc;i4++) chld_r[i4]=new int[cyc_len]; 
	chld_c=new int*[tcyc]; for(i4=0;i4<tcyc;i4++) chld_c[i4]=new int[cyc_len]; 
	cyc_sv_r=new int*[tcyc]; for(i4=0;i4<tcyc;i4++) cyc_sv_r[i4]=new int[cyc_len]; 
	cyc_sv_c=new int*[tcyc]; for(i4=0;i4<tcyc;i4++) cyc_sv_c[i4]=new int[cyc_len];
	cyc_sv_rg=new int*[tcyc]; for(i4=0;i4<tcyc;i4++) cyc_sv_rg[i4]=new int[cyc_len]; 
	cyc_sv_cg=new int*[tcyc]; for(i4=0;i4<tcyc;i4++) cyc_sv_cg[i4]=new int[cyc_len];
	cyc_cg=new int*[tcyc]; for(i4=0;i4<tcyc;i4++) cyc_cg[i4]=new int[cyc_len];
	col_wt=new int[tot_c];
	fbd_r=new int*[revs]; for(i5=0;i5<revs;i5++) fbd_r[i5]=new int[cyc_len]; 
	fbd_c=new int*[revs]; for(i5=0;i5<revs;i5++) fbd_c[i5]=new int[cyc_len]; 
	for(i2=0;i2<revs;i2++) for(i3=0;i3<cyc_len;i3++) fbd_r[i2][i3]=fbd_c[i2][i3]=-1; //initializing

	ifstream inf1,inf2,inf3;
	ofstream outf,outf2,outf3,outf4,outf5,outf6,outf7,outf8,outf9,outf10,outf11,outf12,outf13,outf14,outf15,outf16; string filename; 

	//if(f==2) {cout<<'\n'<<"lsft_mat: "<<'\n'; for(i3=0;i3<tot_r;i3++){for(j=0;j<tot_c;j++) {cout<<lsft_mat[i3][j]<<""; if(j>0 && (j+1)%p==0) cout<<" ";} cout<<'\n'; if(i3>0 && (i3+1)%p==0) cout<<'\n';}}

	H_mat(); //produces the m by n BC protograph matrix
	//else H2(); //produces the m by n AB-BC from protograph
	for(i=0;i<m;i++) for(j=0;j<n;j++) Hsc_proto[i][j]=H_proto[i][j];	

		//ifstream in_file; in_file.open("B1p5.txt"); for(i5=0;i5<gama;i5++) for(j=0;j<p;j++) in_file >> B[i5][j];
		//else for(i5=0;i5<gama;i5++) for(j=0;j<p;j++) B[i5][j]=0; //uncoupled case
		//cout<<'\n'<<"B: "<<'\n'; for(i=0;i<gama;i++) {for(j=0;j<p;j++) cout<<B[i][j]<<" "; cout<<'\n';}

	//Hsc_mat2(); //1st lift
	Hsc2_AB(); //do terminal lift

		//cout<<"here"<<endl;proto

		/*if(!ge) tot_c-=n; //if not generalized edge spreading
		if(J==1) {for(i=0;i<tot_r;i++) for(j=0;j<tot_c;j++) Hsc[i][j]=Hsc_proto[i][j];}
		else if (J>1 && !opt) Hsc2_AB(); //do terminal lift based on saved shifts
		else if(J>1 && (opt || randm)) { //tries to minimize 6-cycles in Hsc via the 2nd lift
			ccnt=ccnt_old=1000000; 
			while(ccnt>0) {
				Hsc2_AB(); //do terminal lift for optimization
				//for(i=0;i<tot_r*J;i++) for(j=0;j<tot_c*J;j++) Hsc_sv[i][j]=Hsc[i][j];

				if(randm) break; //random terminal lift

				ccnt=0;
				main_cyc4();
				//for(i=0;i<tot_r*J;i++) for(j=0;j<tot_c*J;j++) Hsc[i][j]=Hsc_sv[i][j];
				//cout<<'\n'<<"cycles: "<<ccnt;

				if(ccnt<ccnt_old) { //save the 
					ccnt_old=ccnt;
					//for(i=0;i<tot_r;i++) for(j=0;j<tot_c;j++) outf11<<lsft_mat[i][j]<<" "; outf11<<std::endl; outf11.close();
					//outf14<<J<<" "; outf14<<std::endl; outf14.close();
					outf15.open("ccnt.txt"); outf15<<ccnt<<" "; outf15<<std::endl; outf15.close();
				}
			}
		}*/

	outf4.open("mat.txt"); for(i=0;i<tot_r*J;i++) for(j=0;j<tot_c*J;j++) outf4<<Hsc[i][j]<<" "; outf4<<std::endl; outf4.close(); 

	//cout<<'\n'<<"H_proto: "<<'\n'; for(i3=0;i3<m2;i3++){for(j=0;j<n2;j++) {cout<<H_proto[i3][j]<<""; if(j>0 && (j+1)%p==0) cout<<" ";} cout<<'\n'; if(i3>0 && (i3+1)%p==0) cout<<'\n';}
	//cout<<'\n'<<"H: "<<'\n'; for(i3=0;i3<m;i3++){for(j=0;j<n;j++) {cout<<H[i3][j]<<" "; if(j>0 && (j+1)%(p)==0) cout<<"  ";} cout<<'\n'; if(i3>0 && (i3+1)%(p)==0) cout<<'\n';}	
	//cout<<'\n'<<"H0: "<<'\n'; for(i3=0;i3<m;i3++){for(j=0;j<n;j++) {cout<<H0[i3][j]<<""; if(j>0 && (j+1)%(p*J)==0) cout<<" ";} cout<<'\n'; if(i3>0 && (i3+1)%(p*J)==0) cout<<'\n';}
	//cout<<'\n'<<"H1: "<<'\n'; for(i3=0;i3<m;i3++){for(j=0;j<n;j++) {cout<<H1[i3][j]<<""; if(j>0 && (j+1)%(p*J)==0) cout<<" ";} cout<<'\n'; if(i3>0 && (i3+1)%(p*J)==0) cout<<'\n';}
	//cout<<'\n'<<"Hsc: "<<'\n'; for(i3=0;i3<tot_r*J;i3++){for(j=0;j<tot_c*J;j++) {cout<<Hsc[i3][j]<<""; if(j>0 && (j+1)%(p*J)==0) cout<<" ";} cout<<'\n'; if(i3>0 && (i3+1)%(p*J)==0) cout<<'\n';}	

	for(j=0;j<tot_c*J;j++) {
		cnt=0; 
		for(i=0;i<tot_r*J;i++) if(Hsc[i][j]) cnt++;
		if(cnt!=gama) {flg=1; cout<<'\n'<<"col: "<<j; break;}
		else flg=0;
	}
	if(!flg) cout<<'\n'<<"valid"; else cout<<'\n'<<"invalid";
	
	if(!unc) {cout<<'\n'<<"tot_r "<<tot_r*J<<" tot_c "<<tot_c*J;}
	else {cout<<'\n'<<"tot_r "<<gama*p*J<<" tot_c "<<tot_c*J;} 
 

		
	//f=1 for BC
	//f=2 for SC
	//outf6.open("nval.txt"); outf6<<tot_c*J<<" ";  outf6<<std::endl; outf6.close();	
	//outf7.open("fval.txt"); outf7<<f<<" "; outf7<<std::endl; outf7.close();	
	//outf8.open("Lval.txt"); outf8<<L<<" "; outf8<<std::endl; outf8.close();	
	//outf9.open("J.txt"); outf9<<J<<" "; outf9<<std::endl; outf9.close();
	//outf16.open("mem.txt"); outf16<<mem<<" "; outf16<<std::endl; outf16.close();


	//outf5.open("mval.txt"); 
	if(f==1) outf5<<tot_r*J<<" "; 
	else if(f==2 && !unc) outf5<<tot_r*J<<" ";
	else if(f==2 && unc) outf5<<gama*p*J<<" "; //for BC
	outf5<<std::endl; outf5.close();

	t2= time(NULL); cout<<'\n'<<"Executed in "<<t2-t1<<"s"<<'\n';
	cout<<'\n';
}
