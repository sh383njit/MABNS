//enumerating cycles in cluster induced sub-graph

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

long par_num=0,chld_num=0,tcyc=pow(10,5),revs=pow(10,6),abs_64=0; //max. no. of reverses in a rnd can be 12
int *snk_r,*snk_c,*E,**E_vec,**fbd_r,**fbd_c,**H,**Hsc,**Hsc2,**H0,**H1,**Hsc_2,**Hsc_3,**I,**I2,**sigma,**sigma2,**sigma3,**tau2,**B,*par_r,*par_c,**chld_r,**chld_c,*col_wt,
**cyc_sv_r,**cyc_sv_c,**cyc_sv_rg,**cyc_sv_cg,*cyc_sel_r,*cyc_sel_c,**cyc_param,**cyc_cg,**cyc_r_nel,**cyc_c_nel,**regn; 
int p,q,gama,J,m,n,L,cyc_len,strt_r,strt_c,stp_r,stp_c,row_strt,col_strt,row_stp,col_stp; //H has total m rows and n cols
long while_tot=0,rev_tot=0,tot_r,tot_c,sum=0,ccnt=0; 

/************************************** the eqn for abs33 (p^2)(p-1) holds only for odd values of p. It is upper-bound for other p values ***********************************************************/

#include "find_cycl4.cpp"

string func(int n) 
{
    stringstream result;
    result << n;
    return result.str();
}

int main() {	
	
	srand(time(0));	
	double t2,t1= time(NULL);	
	long i,i5,i2,i3,i4,j,j2,k,l,count=0,A,b2,*outp; 
	int f,fl=1,totf=1,*edge_r,*edge_c,num_edgs=pow(10,6),param,ccnt2=0,flg=0,rr,cnt,ncyc=0,c,BC,*sset_r,*sset_c,mu1,mu2,mu3,*cyc_rows,num_cls,cls_sz,**cls_ind_mat;   
	
	ifstream inf1,inf2,inf3,inf4,inf5,inf6,inf7,inf10; 
	//inf1.open("m_win.txt"); inf1>>m; 
	//inf2.open("nval.txt"); inf2>>n; 

	ifstream inf; 
	//inf.open("matrices/mat_n196_AB_2.txt"); m=84; n=196; inf2.open("clusters/cls_ind_mat_opt_AB_2.txt"); cyc_len=6;
	inf.open("matrices/mat_n196_1.txt"); m=98; n=196; inf2.open("clusters/cls_ind_mat_opt_1.txt"); cyc_len=6;
	//inf.open("matrices/mat_n196_1.txt"); m=98; n=196; inf2.open("clusters/cls_ind_mat_contg.txt"); cyc_len=4;
	
	//inf.open("matrices/mat_BCH_63_51.txt"); m=12; n=63; inf2.open("clusters/cls_ind_mat_ran_BCH_63_51.txt"); cyc_len=6;
	//inf.open("matrices/mat_BCH_63_51.txt"); m=12; n=63; inf2.open("clusters/cls_ind_mat_opt_BCH_63_51.txt"); cyc_len=6;

	//inf.open("matrices/mackay_96_48.txt"); m=48; n=96; inf2.open("clusters/cls_ind_mat_ran_mac_96_48.txt"); cyc_len=6;

	cls_sz=6;

	num_cls=m/cls_sz;
	cout<<'\n'<<"m: "<<m<<" n: "<<n<<" cls_sz: "<<cls_sz<<" num_cls: "<<num_cls;

	cls_ind_mat=new int*[num_cls]; for(i=0;i<num_cls;i++) cls_ind_mat[i]=new int[cls_sz];

	tot_r=m; tot_c=n; 

	edge_r=new int[num_edgs]; edge_c=new int[num_edgs]; outp=new long[12]; 
	par_r=new int[cyc_len]; par_c=new int[cyc_len]; cyc_rows=new int[m]; 
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

	regn=new int*[gama]; for(i=0;i<gama;i++) regn[i]=new int[p];

	H=new int*[m]; for(i5=0;i5<m;i5++) H[i5]=new int[n]; //m rows n cols
	Hsc=new int*[tot_r]; for(i5=0;i5<tot_r;i5++) Hsc[i5]=new int[tot_c];
	Hsc2=new int*[tot_r]; for(i5=0;i5<tot_r;i5++) Hsc2[i5]=new int[tot_c];

	
	ofstream outf,outf2,outf3,outf4,outf5,outf6; string filename; 

	for(i=0;i<m;i++) for(j=0;j<n;j++) inf >> Hsc2[i][j];

	for(i=0;i<num_cls;i++) for(j=0;j<cls_sz;j++) inf2 >> cls_ind_mat[i][j]; //for opt. clusters
	cout<<'\n'<<"cls_ind_mat: "<<'\n'; for(i=0;i<num_cls;i++) {for(j=0;j<cls_sz;j++) cout<<cls_ind_mat[i][j]<<" "; cout<<'\n';}

	
	//cout<<'\n'<<"Hsc2: "<<'\n'; for(i3=0;i3<m;i3++){for(j=0;j<n;j++) cout<<Hsc2[i3][j]<<""; cout<<'\n';} 

	//cout<<'\n'<<"row_stp "<<row_stp<<" col_stp "<<col_stp<<" cyc_len "<<cyc_len<<std::endl;

for(i4=0;i4<num_cls;i4++) {

	row_strt=col_strt=0; row_stp=tot_r; col_stp=tot_c; 
	
	//choosing the cluster induced Tanner graph
	for(i=0;i<cls_sz;i++)
		for(k=0;k<n;k++)
			Hsc[cls_ind_mat[i4][i]][k]=Hsc2[cls_ind_mat[i4][i]][k];

	sum=count=ccnt=0; for(i2=0;i2<tot_c;i2++) col_wt[i2]=0; //refresh
	for(i=row_strt;i<row_stp;i++) for(i2=col_strt;i2<col_stp;i2++) if(Hsc[i][i2]) {edge_r[count]=i; edge_c[count]=i2; count++;} 
	//**************imp for standard dec. testing, col_wt=3 always**********************/
	for(i2=col_strt;i2<col_stp;i2++) for(i=row_strt;i<row_stp;i++) if(Hsc[i][i2]) col_wt[i2]=3/*++*/; 
	//cout<<'\n'<<"count: "<<count<<std::endl;
	if(totf>1) {A=(fl-1)*(count/(totf-1)); if(fl<totf) b2=A+count/(totf-1); else b2=count;}
	else {A=0; b2=count;}
	for(i=0;i<A;i++) Hsc[edge_r[i]][edge_c[i]]=0; //decimating unnecessary edges
	//cout<<'\n'<<"A: "<<A<<" b2-1: "<<b2-1<<std::endl;

	for(i=A;i<b2;i++) { 
		//while_tot=0
		par_num=chld_num=0; for(i3=0;i3<cyc_len;i3++) snk_r[i3]=snk_c[i3]=par_r[i3]=par_c[i3]=-1; //refresh
		strt_r=edge_r[i]; strt_c=edge_c[i];
		//cout<<'\n'<<"start: "<<strt_r<<","<<strt_c;	
		find_cycl();
		Hsc[strt_r][strt_c]=0; //decimating unnecessary edges

		//cout<<'\n'<<"while_tot: "<<while_tot<<" rev_tot: "<<rev_tot;
		sum+=par_num+chld_num;
		//cout<<'\n'<<i<<std::endl;
		//if(count==30) break;
	}
	cout<<'\n'<<"no. of "<<cyc_len<<" cycles: "<<ccnt;
	row_strt=row_stp=col_strt=col_stp=0;

	//cout<<'\n'; for(j=0;j<ccnt;j++) {for(i3=0;i3<cyc_len;i3++) cout<<cyc_sv_r[j][i3]<<","<<cyc_sv_c[j][i3]<<" "; cout<<'\n';} //all cycles
	//cout<<'\n'; for(j=0;j<ccnt;j++) {for(i3=0;i3<cyc_len;i3++) if(!(i3%2)) cout<<cyc_sv_r[j][i3]<<" "; cout<<'\n';} //view rows only

	//cout<<'\n'<<"cyc_rows: "; for(i=0;i<m;i++) cout<<cyc_rows[i]<<" "; 

}
	//then apply read2 to get trg and tcg etc

	t2= time(NULL); cout<<'\n'<<"Executed in "<<t2-t1<<"s"<<'\n';
	cout<<'\n';
}
