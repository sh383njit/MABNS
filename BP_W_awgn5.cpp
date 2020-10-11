
//GPU implementation of sliding window BP alg.
//for BC decoding L=W=3

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
using namespace std;

float *H,*H2,*x,*x_hat,*y,P,*LR,*pLR,*E_c_v,*E_v_c,*H_D,*y_D,*LR_D,*pLR_D,*E_c_v_D,*E_v_c_D,err,err2,biterr,biterr2,varn; //n is source length, m is code length
long m,n,n2,*vns,*cns,*vns_D,*cns_D,BC;
long col_wt,row_wt,num,dcw,dcw2,hshft,vshft,mem,W,call,flg; //row_wt is the max, row weight in case of SC codes
long L=99; //coupling length (multiple of mem+1 and >=W) 

#define CW 256 //no. of codeword streams per transmission
#define I 20 //no. of iters
#define ev 10 //min. no. of error events

string func(long n) {
	stringstream result;
	result << n;
	return result.str();
}






void init1() {

	long i,cw,k,strt,stp;
	
	//printf("cw: %d\n", cw);

	for(cw=0;cw<CW;cw++) {
		if(num==0) //if no shift
			for(i=0;i<n;i++) {
				//LR[cw*n+i]=(0.5-y[cw*n+i])/varn; //AWGN case, LLr= ln [exp{-(y-0)^2/2sig^2)}/exp{-(y-1)^2/2sig^2)}]
				LR[cw*n+i]=2*y[cw*n+i]/varn; //AWGN case, LLr= ln(Pr[x=0|y]/Pr[x=1|y])= ln [exp{-(y-1)^2/2sig^2)}/exp{-(y-(-1))^2/2sig^2)}], x=0,1 means y=1,-1
				pLR[cw*n+i]=0;
			}
		else {
			//shifting the msgs (from bottom to top)
			for(i=0;i<n-hshft;i++) {
				LR[cw*n+i]=LR[cw*n+i+hshft]; //hshft is the no. of new VNs entering window
				pLR[cw*n+i]=pLR[cw*n+i+hshft];
				for(k=0;k<col_wt;k++)
					//E_c_v[i][k][cw]=E_c_v[i+hshft][k][cw]; 
					E_c_v[CW*(col_wt*i+k)+cw]=E_c_v[CW*(col_wt*(i+hshft)+k)+cw];
			}

			strt=num*hshft;
			stp=strt+hshft;

			for(i=strt;i<stp;i++) 
				LR[cw*n+(i-strt)+(n-hshft)]=2*y[cw*n+i]/varn; //new LLR values entering the window

		} 
	}

	//cout<<'\n'<<"LR: "; for(i=0;i<n;i++) for(cw=0;cw<CW;cw++) cout<<LR[cw*n+i]<<" "; 

	//printf("LR: \n"); for(i=0;i<n;i++) printf(" %f",LR[cw*n+i]);
	//printf("y: \n"); for(i=0;i<n;i++) printf(" %f",y[cw*n+i]);
}





void init2() {

	long cw,j,k;

	//printf("cw= %d\n", cw);

	for(cw=0;cw<CW;cw++) {
		if(num==0) //if no shifts
			//for(j=0;j<m;j++) for(k=0;k<row_wt;k++) E_v_c[j][k][cw]=LR[cw][vns[j][k]];
			for(j=0;j<m;j++) 
				for(k=0;k<row_wt;k++) 
					E_v_c[CW*(row_wt*j+k)+cw]=LR[cw*n+vns[j*row_wt+k]]; //CN accumulating initial msgs

		else {
			//shifting the msgs (from bottom to top)
			for(j=0;j<m-vshft;j++) {
				for(k=0;k<row_wt;k++)
					//E_v_c[j][k][cw]=E_v_c[j+vshft][k][cw]; //vshft is the no. of new CNs entering window
					E_v_c[CW*(row_wt*j+k)+cw]=E_v_c[CW*(row_wt*(j+vshft)+k)+cw];
			}

			for(j=0;j<vshft;j++) 
				for(k=0;k<row_wt;k++)
					//E_v_c[m-vshft+j][k][cw]=LR[cw*n+vns[(m-vshft+j)*row_wt+k]]; 
					E_v_c[CW*(row_wt*(m-vshft+j)+k)+cw]=LR[cw*n+vns[(m-vshft+j)*row_wt+k]]; //CNs of window initialized with new LLR values 

		}
	}

	//printf("E_v_c: \n"); for(j=0;j<n;j++) printf(" %.1f",E_v_c[j]);
	//cw=0;
	//cout<<'\n'<<"E_v_c: "<<'\n'; for(j=0;j<m;j++){for(k=0;k<row_wt;k++) cout<<E_v_c[CW*(row_wt*j+k)+cw]<<" "; cout<<'\n';}
	//printf("num= %d\n", num);
}






void horz() { 
	long i,cw,j,k,vidx;
	float tmp;

	//printf("cw,j=%d,%d\n", cw,j);
	//printf("\n");

	for(cw=0;cw<CW;cw++) {
		for(j=0;j<m;j++) { //m is no. of CNs
			for(i=0;i<row_wt;i++) { //row_wt is the no. of neighboring VNs of CN j

				vidx=vns[j*row_wt+i]; //index of ith neighboring VN of CN j

				if(vidx>-1) {
					tmp=1; 
					for(k=0;k<row_wt;k++) 
						if(vns[j*row_wt+k]>-1 && k!=i) tmp*=tanh(0.5*E_v_c[CW*(row_wt*j+k)+cw]); //jth CN accumulating msgs from all neighboring VNs except i
					
					//cout<<'\n'<<"j: "<<j<<" vidx: "<<vidx<<" tmp: "<<tmp<<endl;
		
					for(k=0;k<col_wt;k++) 
						if(cns[vidx*col_wt+k]==j) {E_c_v[CW*(vidx*col_wt+k)+cw]=2*atanhf(tmp); break;} //msg sent by jth CN to ith VN				
				}
			} 
		}
	}

	//cw=0; cout<<'\n'<<"E_c_v: "<<'\n'; for(j=0;j<n;j++){for(k=0;k<col_wt;k++) cout<<E_c_v[CW*(col_wt*j+k)+cw]<<" "; cout<<'\n';}
	
}






void vert() {
	long i,cw,j,k,cidx;
	float tmp;

	//printf("cw,i=%d,%d\n", cw,i);

	for(cw=0;cw<CW;cw++) {
		for(i=0;i<n;i++) { //n is no. of VNs
			for(j=0;j<col_wt;j++) { //col_wt is the no. of neighboring CNs of VN i

				cidx=cns[i*col_wt+j];  //index of jth neighboring CN of VN i

				if(cidx>-1) {
					tmp=0; 
					for(k=0;k<col_wt;k++) 
						if(cns[i*col_wt+k]>-1 && k!=j) tmp+=E_c_v[CW*(col_wt*i+k)+cw]; //ith VN accumulating msg from all neighboring CNs except j 
	
					for(k=0;k<row_wt;k++) 
						if(vns[cidx*row_wt+k]==i) {E_v_c[CW*(row_wt*cidx+k)+cw]=tmp+LR[cw*n+i]; break;} //msg sent by ith VN to jth CN
						//printf("tmp2=%f\n",tmp2);					
				}
			} 

			//updating the aposteriori LLR

			tmp=0; 
			for(k=0;k<col_wt;k++) 
				if(cns[col_wt*i+k]>-1) tmp+=E_c_v[CW*(col_wt*i+k)+cw];
			pLR[cw*n+i]=LR[cw*n+i]+tmp; 

		}
	}

	//cw=0; cout<<'\n'<<"E_v_c: "<<'\n'; for(j=0;j<m;j++){for(k=0;k<row_wt;k++) cout<<E_v_c[CW*(row_wt*j+k)+cw]<<" "; cout<<'\n';}
	//cout<<'\n'<<"pLR: "; for(i=0;i<n;i++) for(cw=0;cw<CW;cw++) cout<<pLR[cw*n+i]<<" ";
}





void BP() {
	long i,cw,j,l,strt=0,stp,STP,flg3;

	//E_v_c[j][i][cw]: msg accumulated by CN j from VN i for codeword stream cw, E_v_c[j][i][cw]=E_v_c[CW*(no. of cols *j + i)+cw]
	//E_c_v[i][j][cw]: msg accumulated by VN i from CN j for codeword stream cw

	num=0;
	while(num<=(L-W)/(mem+1)) {

		//cout<<'\n'<<"num: "<<num<<endl;
		//initializing the window
		init1();
		init2(); 

		for(l=0;l<I;l++) {
			//horizontal step
			horz();
			
			//vertical step
			vert();			
		
		} 

		//cout<<'\n'<<"pLR: "; for(cw=0;cw<CW;cw++){ for(i=0;i<n;i++) cout<<pLR[cw*n+i]<<" "; cout<<'\n'<<'\n';}
		
		//hard decision
		if(num<(L-W)/(mem+1) || BC) STP=1;
		else STP=W/(mem+1); //decision in the last window position

		for(j=0;j<STP;j++) {
			if(num<(L-W)/(mem+1))
				strt=num*hshft;
			else if(num==0 && j==0) 
				strt=0;
			else strt+=hshft;

			stp=strt+hshft;
			//cout<<'\n'<<"strt: "<<strt<<" stp: "<<stp<<endl;
			
			for(cw=0;cw<CW;cw++) {
				flg3=0;
				for(i=strt;i<stp;i++) {
					if(pLR[cw*n+i-strt]<0) {
						if(!flg3) {
							if(!flg) err++; //for blk err
							else err2++;
							flg3=1;
						}
						if(!flg) biterr++; //for bit err
						else biterr2++;
						//break;
					}
 					//if(pLR[cw*n+i-strt]>=0) 
						//x_hat[cw*n+i]=0; 
					//else 
						//x_hat[cw*n+i]=1;
				}
				if(!flg) dcw++; //no. of decoded cws
				else dcw2++;

				//for(i=strt;i<stp;i++) 
					//if(x[cw*n+i]!=x_hat[cw*n+i]) {
						//err++; 
						//break;
				//}	
			}
		
			//cout<<'\n'<<"pLR: "; for(cw=0;cw<CW;cw++){ for(i=strt;i<stp;i++) cout<<pLR[cw*n+i-strt]<<" "; cout<<'\n'<<'\n';}

		}

		num++; //no. of window shifts

		if(!flg && err>=ev) break; //to avoid waiting until the last window position
		else if(flg && err2>=ev) break;
		
	}
		
}	


int main() {	
	srand(time(0));	
	clock_t tStart = clock();	
	long i,i2,j,kk; 
	long fn,num_dat,cnt,cw,gama,p,J,test=0; 

	float *blk_err,*blk_err2,*bit_err,*bit_err2,R;  
	//cout<<'\n'<<" device no: "; cin>>fn; ///////////////////////////////// 
	//cout<<'\n'<<" kk: "; cin>>kk; 

	kk=1;	

	if(!kk) { //if not Kelley Kliewer code
		cout<<'\n'<<" col_wt: "; cin>>col_wt; 
		cout<<'\n'<<" BC?: "; cin>>BC; //1 for BC, 0 for SC
	}
	else BC=1;


	//construct a random generator engine:
	std::random_device rd;
    	std::mt19937 e2(rd());

	ifstream inf,inf1,inf2,inf3,inf4,inf5,inf6,inf7,inf8,inf9,inf10; 
	if(!BC && !kk) {
		inf1.open("mval.txt"); inf1>>m; 
		inf2.open("nval.txt"); inf2>>n;
	} 
	else if(BC && !kk) {
		inf1.open("mval_unc.txt"); inf1>>m; 
		inf2.open("nval_unc.txt"); inf2>>n;
	} 
	else if(kk) {
		inf1.open("mval_kk.txt"); inf1>>m; 
		inf2.open("nval_kk.txt"); inf2>>n;
	} 
	
	R=1-float(m)/float(n); //code rate

	inf3.open("J.txt"); inf3>>J; 
	inf4.open("p.txt"); inf4>>p; 
	inf5.open("mem.txt"); inf5>>mem;
	inf6.open("gama.txt"); inf6>>gama;
	inf7.open("W.txt"); inf7>>W; //no. of col. blks of window (multiple of mem+1). width of col. blk is Jp^2
	if(BC) {L=3;W=3;}
	row_wt=p;

	if(test) {
		inf1.open("mval_test.txt"); inf1>>m; 
		inf2.open("nval_test.txt"); inf2>>n;
		//maximum values
		inf9.open("col_wt_test.txt"); inf9 >> col_wt; 
		inf10.open("row_wt_test.txt"); inf10 >> row_wt;
	}

	H= new float[m*n]; //H is the window. Since windows are identical, only the window Tanner graph is needed. The Lvalues "move" across the window
	H2= new float[m*n];
	ofstream outf,out_file2; string filename;

	//inf.open("H_3_7_24_5_sc_rnd1.txt"); for(j=0;j<m;j++) for(i=0;i<n;i++) inf >> H[j*n+i]; //upload the parity-check matrix
	//inf.open("mat_p5_sc.txt"); for(j=0;j<m;j++) for(i=0;i<n;i++) inf >> H[j*n+i];

	if(kk) { //Kelley Kliewer random lifted BC
		inf.open("mat_kk.txt"); for(j=0;j<m;j++) for(i=0;i<n;i++) inf >> H[j*n+i];
		inf9.open("col_wt_kk.txt"); inf9 >> col_wt; //maximum values
		inf10.open("row_wt_kk.txt"); inf10 >> row_wt;
	}
	else if(col_wt==3 && !BC) {
		inf.open("win_3_sc.txt"); for(j=0;j<m;j++) for(i=0;i<n;i++) inf >> H[j*n+i];
		inf8.open("win_3_sc_rnd1.txt"); for(j=0;j<m;j++) for(i=0;i<n;i++) inf8 >> H2[j*n+i];
	}
	else if(col_wt==3 && BC) {
		inf.open("win_3_unc.txt"); for(j=0;j<m;j++) for(i=0;i<n;i++) inf >> H[j*n+i];
	}
	else if(col_wt==4 && !BC) {
		inf.open("win_4_sc.txt"); for(j=0;j<m;j++) for(i=0;i<n;i++) inf >> H[j*n+i];
		inf8.open("win_4_sc_rnd1.txt"); for(j=0;j<m;j++) for(i=0;i<n;i++) inf8 >> H2[j*n+i];
	}
	else if(col_wt==4 && BC) {
		inf.open("win_4_unc.txt"); for(j=0;j<m;j++) for(i=0;i<n;i++) inf >> H[j*n+i];
	}
	//cout<<'\n'<<"H: "<<'\n'; for(j=0;j<m;j++){for(i=0;i<n;i++) {cout<<H[j*n+i]<<" "; if(i>0 && (i+1)%(row_wt)==0) cout<<"  ";} cout<<'\n'; if(j>0 && (j+1)%(row_wt)==0) cout<<'\n';}

	//col_wt=0;
	//for(j=0;j<m;j++) if(H[j*n+0]) col_wt++;

	cout<<'\n'<<"m: "<<m<<" n: "<<n<<" R: "<<R<<" col_wt: "<<col_wt<<" row_wt: "<<row_wt<<" gama: "<<gama<<" W: "<<W<<endl;  

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
	//x_hat= new float[CW*n2];
	y= new float[CW*n2];  
	LR= new float[CW*n];  //Lvalues in the window
	pLR= new float[CW*n];  

	vns=new long[m*row_wt]; //all VNs of CN j in a row
	
	//3D arrays
	E_v_c= new float[m*row_wt*CW];
	E_c_v= new float[n*col_wt*CW];

	cns=new long[n*col_wt]; //all CNs of a VN in a row


	//*********************col_wt 3****************************
	//Eb/No values in dB
	float EbNo[]={1.42,1.47,1.52,1.54,1.56,1.58,1.6,1.62,1.64,1.66,1.67}; //p=7, SC 
	float EbNo_BC[]={1.55,1.65,1.72,1.76,1.78,1.8,1.82,1.84,1.85}; //p=7, BC 
	
	//*********************col_wt 4****************************
	float EbNo2[]={1,1.1,1.15,1.2,1.25,1.3,1.32,1.33,1.34,1.35,1.36}; //p=7, SC
	float EbNo2_BC[]={1.44,1.47,1.5,1.523,1.55,1.58,1.61,1.63,1.64}; //p=7, BC


	//float EbNo_kk[]={0.5,1,1.1,1.15,1.2,1.22,1.24,1.26,1.28,1.3,1.32}; //for Kelley Kliewer random sub-code1
	//float EbNo_kk[]={0,0.4,0.8,1.15,1.45,1.8,2}; 
	float EbNo_kk[]={0,0.6,1.2,1.4,1.8,2,2.4,2.8}; //,3.2};

	if(!kk) {
		if(col_wt==3 && !BC) num_dat=sizeof(EbNo)/sizeof(EbNo[0]);
		else if(col_wt==3 && BC) num_dat=sizeof(EbNo_BC)/sizeof(EbNo_BC[0]);
		else if(col_wt==4 && !BC) num_dat=sizeof(EbNo2)/sizeof(EbNo2[0]);
		else if(col_wt==4 && BC) num_dat=sizeof(EbNo2_BC)/sizeof(EbNo2_BC[0]);
	}
	else num_dat=sizeof(EbNo_kk)/sizeof(EbNo_kk[0]);

	blk_err= new float[num_dat]; blk_err2= new float[num_dat];
	bit_err= new float[num_dat]; bit_err2= new float[num_dat]; 





	


	//***********************************************    test   ****************************************
	if(test) {
		y[0*n2+0]=0.322; y[0*n2+1]=-0.037; y[0*n2+2]=0; y[0*n2+3]=0; y[0*n2+4]=0; //1st noisy CW
	
		inf.open("mat_test.txt"); for(j=0;j<m;j++) for(i=0;i<n;i++) inf >> H[j*n+i];
	
		for(i=0;i<m;i++) for(j=0;j<row_wt;j++) vns[i*row_wt+j]=-1;
		for(i=0;i<n;i++) for(j=0;j<col_wt;j++) cns[i*col_wt+j]=-1;
	
		for(j=0;j<m;j++) {
			cnt=0; 
			for(i=0;i<n;i++) if(H[j*n+i]) {vns[j*row_wt+cnt]=i; cnt++;} 
		}
		cout<<'\n'<<"vns "<<'\n'; for(i=0;i<m;i++) {for(j=0;j<row_wt;j++) cout<<vns[i*row_wt+j]<<" "; cout<<'\n';}
		
		for(i=0;i<n;i++) {
			cnt=0; 
			for(j=0;j<m;j++) if(H[j*n+i]) {cns[i*col_wt+cnt]=j; cnt++;} 
		}
		cout<<'\n'<<"cns "<<'\n'; for(i=0;i<n;i++) {for(j=0;j<col_wt;j++) cout<<cns[i*col_wt+j]<<" "; cout<<'\n';}
	
		varn=0.1;
		
		init1();
		init2(); 
		for(i=0;i<I;i++) {
	 		horz(); //horizontal step
			vert();	//vertical step		
		}
	}
	//**************************************************************************************************


		
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


	for(i2=0;i2<num_dat;i2++) {
		if(!kk) {
			if(col_wt==3 && !BC) varn=1/(2*R*pow(10,0.1*EbNo[i2]));
			else if(col_wt==3 && BC) varn=1/(2*R*pow(10,0.1*EbNo_BC[i2]));
			else if(col_wt==4 && !BC) varn=1/(2*R*pow(10,0.1*EbNo2[i2]));
			else if(col_wt==4 && BC) varn=1/(2*R*pow(10,0.1*EbNo2_BC[i2]));
		}
		else varn=1/(2*R*pow(10,0.1*EbNo_kk[i2])); //if EbNo not in dB, then varn=1/(2*R*EbNo_kk[i2]);

		//std::normal_distribution<float> dist(0,varn);
		std::normal_distribution<float> dist(0,sqrt(varn));

		err=err2=biterr=biterr2=dcw=dcw2=call=0;
		while((!BC && (err<ev || err2<ev)) || (BC && err<ev)) { //take avg. of min. 3 blk. error events
			for(cw=0;cw<CW;cw++) for(i=0;i<n2;i++) y[cw*n2+i]=1+dist(e2); //adding gaussian noise to all-zero CW
			//cw=0; cout<<'\n'<<"y: "; for(i=0;i<n2;i++) cout<<y[cw*n2+i]<<" ";
	
			//ofstream outf; outf.open("y.txt"); for(i=0;i<n2;i++) outf<<y[cw*n2+i]<<" "; outf<<std::endl; outf.close();

			flg=0;
			while((!BC && flg<2) || (BC && !flg)) {
				
				if(!flg && !BC) {
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
				else if(flg && !BC) {
					//initializing
					for(i=0;i<m;i++) for(j=0;j<row_wt;j++) vns[i*row_wt+j]=-1;
					for(i=0;i<n;i++) for(j=0;j<col_wt;j++) cns[i*col_wt+j]=-1;

					for(j=0;j<m;j++) {
						cnt=0; 
						for(i=0;i<n;i++) if(H2[j*n+i]) {vns[j*row_wt+cnt]=i; cnt++;} 
					}
	
					for(i=0;i<n;i++) {
						cnt=0; 
						for(j=0;j<m;j++) if(H2[j*n+i]) {cns[i*col_wt+cnt]=j; cnt++;} 
					}
				}

				BP(); //call wind. dec.

				//if(!flg) {cout<<"err: "<<err<<" dcw: "<<dcw<<endl; }
				//else {cout<<"err2: "<<err2<<" dcw2: "<<dcw2<<endl; cout<<'\n';}

				//if(!flg) call++;	
			
				flg++;
				
			}

		}

		cout<<'\n'<<"i2: "<<i2<<endl;

		//trans=(call*n2/hshft)*CW; //there are n2/hshft cws in a stream. Each cw has length hshft bits in a stream of length n2 bits
		blk_err[i2]=err/dcw;
		bit_err[i2]=biterr/(hshft*dcw);

		blk_err2[i2]=err2/dcw2;		
		bit_err2[i2]=biterr2/(hshft*dcw2);

		cout<<"EbNo: "; 
		for(i=0;i<=i2;i++) 
			if(!kk) {
				if(col_wt==3 && !BC) cout<<EbNo[i]<<" ";
				else if(col_wt==3 && BC) cout<<EbNo_BC[i]<<" "; 
				else if(col_wt==4 && !BC) cout<<EbNo2[i]<<" "; 
				else if(col_wt==4 && BC) cout<<EbNo2_BC[i]<<" "; 
			}
			else cout<<EbNo_kk[i]<<" ";
		cout<<'\n'<<endl;

		cout<<"blk_err: "; for(i=0;i<=i2;i++) cout<<blk_err[i]<<" "; cout<<'\n'<<endl;
		cout<<"bit_err: "; for(i=0;i<=i2;i++) cout<<bit_err[i]<<" "; cout<<'\n'<<endl;
		if(!BC) {
			cout<<"blk_err2: "; for(i=0;i<=i2;i++) cout<<blk_err2[i]<<" "; cout<<'\n'<<endl;
			cout<<"bit_err2: "; for(i=0;i<=i2;i++) cout<<bit_err2[i]<<" "; cout<<'\n'<<endl;
		}
		 
	}

	//cout<<'\n'<<"blk_err: "; for(i=0;i<num_dat;i++) cout<<blk_err[i]<<" ";
	//ofstream outf3; filename="Pe"+func(fn)+".txt"; outf3.open(filename.c_str()/*,fstream::app*/); for(i2=0;i2<num_dat;i2++) outf3<<blk_err[i2]<<" ";  outf3<<std::endl; outf3.close();

	cout<<'\n';
    	printf("executed in: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
	cout<<'\n';


}





