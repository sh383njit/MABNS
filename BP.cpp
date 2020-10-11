
//called by main_BP

#include "init.cpp"
#include "MP0.cpp"
#include "MP1.cpp"
#include "MP1b.cpp"
#include "MP2.cpp"
#include "MP3.cpp" //for contigious, optimized or random clustering
//#include "MP4.cpp" //for Deep-RL
//#include "MP5.cpp" //for model-based (learned CN ordering)
#include "MP5c.cpp" //for model-based using Matlab


//int BP(PyObject *pFunc) {
int BP(int cnt2) {
	long i,k,cw,j,l,strt=0,stp,iter,STP,flg2,flg3,flg4,cnt,niter;
	float cnt1;

	//E_v_c[j][i][cw]: msg accumulated by CN j from VN i for codeword stream cw, E_v_c[j][i][cw]=E_v_c[CW*(no. of cols *j + i)+cw]
	//E_c_v[i][j][cw]: msg accumulated by VN i from CN j for codeword stream cw
	
	//refresh
	for(i=0;i<CW;i++) excl_cw[i]=-1;

	num=0;
	while(num<=(L-W)/(mem+1)) {

		//cout<<'\n'<<"num: "<<num<<endl;
		//initializing the window
		init1();
		init2(); 
		//initialization
		for(cw=0;cw<CW;cw++) for(j=0;j<m;j++) for(k=0;k<row_wt;k++) {
			//E_v_c_mu[CW*(row_wt*j+k)+cw]=mu_LR; //mean of LLR 
			E_v_c_mu[CW*(row_wt*j+k)+cw]=0;
			//std::normal_distribution<float> dist(mu_LR,sqrt(2*mu_LR));
			//E_v_c[CW*(row_wt*j+k)+cw]=dist(e2);
		}

		cnt=0; flg4=1;

		for(cw=0;cw<CW;cw++) for(i=0;i<m;i++) for(j=0;j<row_wt;j++) res_c_v[CW*(row_wt*i+j)+cw]=0; //refresh residual

		if(meth==1) iter=Ifl; //for flooding
		else iter=I; 
		
		//cout<<'\n'<<"iter: "<<iter<<endl;
		cnt1=0;
		for(l=0;l<iter;l++) {
			if(meth==1) MP0(cnt); //flooding
			else if(meth==2) MP1(cnt,l); //GI scheduling
			//else if(meth==2) MP1b(cnt,l); //random scheduling
			
			else if(meth==3) MP2(cnt,l); //NS

			else if(meth==4) MP3(cnt,l); //random model-free
			else if(meth==5) MP3(cnt,l); //optimized model-free
			else if(meth==6) MP3(cnt,l); //random model-based for JSAIT, contiguous model-free for NIPS
			else if(meth==7) MP3(cnt,l); //optimized model-based
			else if(meth==8) MP5c(cnt,cnt1,l); //NSTS

			//else if(meth==7) MP4(cnt,l,pFunc); //Deep RL based
			
			//checking syndromes	
			for(cw=0;cw<CW;cw++) {
				for(j=0;j<m;j++) syn[j]=0; //refresh
				for(k=0;k<cnt;k++) 
					if(cw!=excl_cw[k]) flg4=1;
					else {flg4=0; break;}

				if(flg4) { //if cw was not recovered previously
					flg2=0;
					for(j=0;j<m;j++) {
 						for(i=0;i<n;i++) {
							if(pLR[cw*n+i]>=0) x_hat[cw*n+i]=0; 
							else x_hat[cw*n+i]=1;
						}
						for(i=0;i<n;i++) 
							syn[j]+=H[j*n+i]*x_hat[cw*n+i]; //finding cw syndrome after each iter.
						syn[j]=syn[j]%2;
						if(syn[j]) {flg2=1; break;} //CN is not satisfied hence error
					}	
					if(!flg2) { //CW has been recovered
						excl_cw[cnt]=cw;
						cnt++;
					}
				}	
			}
	
			//cout<<'\n'<<"excl_cw: "; for(i=0;i<cnt;i++) cout<<excl_cw[i]<<" "; cout<<'\n';
			//cout<<'\n'<<"l: "<<l<<" cnt: "<<cnt;			
			if(cnt==CW) break; //all cws were recovered
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
							if(meth==1) err++; 
							else if(meth==2) err2++;
							else if(meth==3) err3++;
							else if(meth==4) err4++;
							else if(meth==5) err5++;
							else if(meth==6) err6++;
							else if(meth==7) err7++;
							else if(meth==8) err8++;
							flg3=1;
						}
						if(meth==1) biterr++; //for bit err
						else if(meth==2) biterr2++;
						else if(meth==3) biterr3++;
						else if(meth==4) biterr4++;
						else if(meth==5) biterr5++;
						else if(meth==6) biterr6++;
						else if(meth==7) biterr7++;
						else if(meth==8) biterr8++;
						//break;
					}
 					//if(pLR[cw*n+i-strt]>=0) 
						//x_hat[cw*n+i]=0; 
					//else 
						//x_hat[cw*n+i]=1;
				}
				if(meth==1) dcw++; //no. of decoded cws
				else if(meth==2) dcw2++;
				else if(meth==3) dcw3++;
				else if(meth==4) dcw4++;
				else if(meth==5) dcw5++;
				else if(meth==6) dcw6++;
				else if(meth==7) dcw7++;
				else if(meth==8) dcw8++;
	
				//for(i=strt;i<stp;i++) 
					//if(x[cw*n+i]!=x_hat[cw*n+i]) {
						//err++; 
						//break;
				//}	
			}
		
			//cout<<'\n'<<"pLR: "; for(cw=0;cw<CW;cw++){ for(i=strt;i<stp;i++) cout<<pLR[cw*n+i-strt]<<" "; cout<<'\n'<<'\n';}

		}

		num++; //no. of window shifts
		
	}
 

	return l;
}
