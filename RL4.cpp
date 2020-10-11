
//called my main_RL
//for any type of clustering

#include "init.cpp"
#include "MP_RL.cpp" //direct residual calculation of the pulled CN
//#include "MP_RL2.cpp" //calculates residual according to the NS scheme
#include "MP_RL4.cpp" //RL for model-based scheme

int RL(int modl) {
	long a,s=0,s_new,i,i2,k,cw,j,j2,l,strt=0,stp,STP,flg,flg2,flg3,flg4,cnt,niter,rnd,a_new,tmp2,row,cls_idx,cls_CN;
	float tmp,Qmax,R,val,res_max,mu_LR;

	//E_v_c[j][i][cw]: msg accumulated by CN j from VN i for codeword stream cw, E_v_c[j][i][cw]=E_v_c[CW*(no. of cols *j + i)+cw]
	//E_c_v[i][j][cw]: msg accumulated by VN i from CN j for codeword stream cw
	
	//for(i=0;i<CW;i++) excl_cw[i]=-1; //refresh

	//ofstream outf; outf.open("syn.txt",fstream::app);  

	num=0;
	//initialize BP decoder
	init1();
	init2(); 

	//for model based case
	if(modl) {
		mu_LR=0;
		for(cw=0;cw<CW;cw++) for(i=0;i<n;i++) 
			mu_LR+=LR[cw*n+i];
		mu_LR/=n;
		//cout<<'\n'<<"mu_LR: "<<mu_LR;

		//initialization
		for(cw=0;cw<CW;cw++) for(j=0;j<m;j++) for(k=0;k<row_wt;k++) {
			E_v_c_mu[CW*(row_wt*j+k)+cw]=mu_LR; //mean of LLR 
			//std::normal_distribution<float> dist(mu_LR,sqrt(2*mu_LR));
			//E_v_c[CW*(row_wt*j+k)+cw]=dist(e2);
		} 
		for(cw=0;cw<CW;cw++) for(j=0;j<n;j++) for(k=0;k<col_wt;k++) {
			E_c_v_mu[CW*(col_wt*j+k)+cw]=E_c_v[CW*(col_wt*j+k)+cw]=0;
		}
	}

	//cout<<'\n'<<"S: "<<S<<endl;
		
	for(cw=0;cw<CW;cw++) {
		//cout<<'\n'<<"cw: "<<cw;
		//compute observed syndrome
		for(j2=0;j2<m;j2++) {
 			//for(i=0;i<n;i++) 
				//if(LR[cw*n+i]>=0) x_hat[cw*n+i]=0; 
				//else x_hat[cw*n+i]=1;

			if(!boxplus) {
				val=0;
				for(i=0;i<n;i++) 
					val+=H[j2*n+i]*LR[cw*n+i]; //syndrome
			}
			else {
				val=1;
				for(i=0;i<n;i++) 
					if(H[j2*n+i]) 
						val*=tanh(0.5*LR[cw*n+i]); 
				val=2*atanhf(val); //for box-plus 
			}
				
			if(M==2) {
				if(val>=0) syn[j2]=0; 
				else syn[j2]=1;
			}	

			//cout<<'\n'<<"qstrt+2*gap: "<<qstrt+2*gap;

			else if(M==4) {
				if(val<=(rep[0]+rep[1])/2) syn[j2]=0; 
				else if(val>(rep[0]+rep[1])/2 && val<=(rep[1]+rep[2])/2) syn[j2]=1;
				else if(val>(rep[1]+rep[2])/2 && val<=(rep[2]+rep[3])/2) syn[j2]=2;
				else if(val>=(rep[2]+rep[3])/2) syn[j2]=3;
			}
				
			else if(M>4) {
				if(val<=(rep[0]+rep[1])/2) syn[j2]=0;
				else if(val>=(rep[M-2]+rep[M-1])/2) syn[j2]=M-1;
				else 
					for(i=1;i<M-1;i++) 
						if(val>rep[i] && val<=rep[i+1]) syn[j2]=i;							
			}
	
		}

		//cout<<'\n'<<"syn_obs: "; for(i=0;i<m;i++) cout<<syn[i]<<" ";
		//cout<<'\n'<<"s_obs: "<<s;

		rnd=R=0;
		//Q learning
		while(rnd<rnd_max) { //rnd_max rounds
			//choose action a
			if(!rnd) a=rand()%m;  
			else a=a_new;

			//cout<<'\n'<<"s,a: "<<s<<","<<a<<endl;

			//finding the current cluster index cls_idx of CN 'a'
			flg=0;
			for(j=0;j<num_cls;j++) {
				for(i=0;i<cls_sz;i++) {
					if(a==cls_ind_mat[j][i]) { //each row of cls_ind_mat represents a cluster with CN indices in each col.
						cls_idx=j;
						cls_CN=i;
						//cls_sz=cls_sz_vec[j];
						flg=1;
						break;
					}
					if(flg) break;
				}
			}

			for(i=0;i<cls_sz;i++) 
				syn_cls[i]=syn[cls_ind_mat[cls_idx][i]]; 
						
			if(DeepRL==1 && batch_sz[cls_idx]<bch_sz_max) {
				for(i=0;i<cls_sz;i++) {
					//cout<<'\n'<<"cls_idx, batch_sz[cls_idx], i: "<<cls_idx<<", "<<batch_sz[cls_idx]<<", "<<i;
					syn_sv[cls_idx][batch_sz[cls_idx]][i]=syn_cls[i]; //saving current cluster syndromes
				}
				//batch_sz[cls_idx]++; //done later
			}
							
			s=0; 
			for(j=0;j<cls_sz;j++) 
				s+=syn_cls[j]*pow(M,cls_sz-j-1);

			if(!modl) MP_RL(a,cw); //execute action a (schedule CN a) for model-free case
			else MP_RL4(a,cw,mu_LR); //for model-based

			//Reward computation
			R=res_c_v_srtd[CW*(row_wt*a+0)+cw]; //pick highest residual of the scheduled CN
			if(DeepRL==1 && batch_sz[cls_idx]<bch_sz_max) {
				a_sv[cls_idx][batch_sz[cls_idx]]=cls_CN; //saving the action taken (index of CN 'a' in the cluster)
				R_sv[cls_idx][batch_sz[cls_idx]]=R; //saving the corresponding reward earned
			}
			//cout<<'\n'<<"R: "<<R;

			//compute new main syndrome
			for(j2=0;j2<m;j2++) {
 				//for(i=0;i<n;i++) 
					//if(pLR[cw*n+i]>=0) x_hat[cw*n+i]=0; 
					//else x_hat[cw*n+i]=1;
				syn_old[j2]=syn[j2];
				
				if(!boxplus) {
					val=0;
					for(i=0;i<n;i++) 
						val+=H[j2*n+i]*pLR[cw*n+i];
				}
				else {
					val=1;
					for(i=0;i<n;i++) 
						if(H[j2*n+i]) 
							val*=tanh(0.5*LR[cw*n+i]); 
					val=2*atanhf(val); //box-plus 
				}
				
				if(M==2) {
					if(val>=0) syn[j2]=0; 
					else syn[j2]=1;
				}	

				else if(M==4) {
					if(val<=(rep[0]+rep[1])/2) syn[j2]=0; 
					else if(val>(rep[0]+rep[1])/2 && val<=(rep[1]+rep[2])/2) syn[j2]=1;
					else if(val>(rep[1]+rep[2])/2 && val<=(rep[2]+rep[3])/2) syn[j2]=2;
					else if(val>=(rep[2]+rep[3])/2) syn[j2]=3;
			 	}
				
				else if(M>4) {
					if(val<=(rep[0]+rep[1])/2) syn[j2]=0;
					else if(val>=(rep[M-2]+rep[M-1])/2) syn[j2]=M-1;
					else 
						for(i=1;i<M-1;i++) 
							if(val>rep[i] && val<=rep[i+1]) syn[j2]=i;							
				}

				//outf<<val<<" "; //for finding quantization levels
			}

			//cout<<'\n'<<"syn_old: "; for(i=0;i<m;i++) cout<<syn_old[i]<<" ";
			//cout<<'\n'<<"syn: "; for(i=0;i<m;i++) cout<<syn[i]<<" ";

			//picking cluster synd. from main synd.
			//for(i=cls_idx*cls_sz;i<(cls_idx+1)*cls_sz;i++) 
				//syn_cls[i-cls_idx*cls_sz]=syn[i]; //cls_sz is no. of CNs in a cluster
				
			for(i=0;i<cls_sz;i++) 
				syn_cls[i]=syn[cls_ind_mat[cls_idx][i]]; 
						
			if(DeepRL==1 && batch_sz[cls_idx]<bch_sz_max) {
				for(i=0;i<cls_sz;i++) {
					//cout<<'\n'<<"cls_idx, batch_sz[cls_idx], i: "<<cls_idx<<", "<<batch_sz[cls_idx]<<", "<<i;
					syn_sv2[cls_idx][batch_sz[cls_idx]][i]=syn_cls[i]; //saving current cluster syndromes
				}
				batch_sz[cls_idx]++;
			}
					
			s_new=0; 
			for(j=0;j<cls_sz;j++) 
				s_new+=syn_cls[j]*pow(M,cls_sz-j-1);

			//cout<<'\n'<<"syn_new : "; for(i=0;i<cls_sz;i++) cout<<syn[i]<<" "; cout<<endl;
			//cout<<'\n'<<"xhat_new: "; for(i=0;i<n;i++) cout<<x_hat[i]<<" ";

			//cout<<'\n';		
			//cout<<'\n'<<"s_new: "<<s;

			//finding action that gives max. Q value in state s_new
			if(!DeepRL) {
				for(i=0;i<cls_sz;i++) Q_temp[i]=0; //refresh			
				row=cls_idx*S+s_new; //row of the Q table, S is no. of syndromes per cluster
				for(i=0;i<cls_sz;i++) 
					if(i!=a%cls_sz) 
						Q_temp[i]=Q[row*cls_sz+i];

				//cout<<'\n'<<"Q: "; for(i=0;i<cls_sz;i++) cout<<Q_temp[i]<<" "; cout<<'\n';
	
				//sort descending
				for(i=0;i<cls_sz;i++) indx[i]=i; //refresh
				for(i=0;i<cls_sz;i++) {
					for(i2=i+1;i2<cls_sz;i2++) {
						if(Q_temp[i]<Q_temp[i2]) {   
							tmp=Q_temp[i];
							Q_temp[i]=Q_temp[i2];
							Q_temp[i2]=tmp;
			
							tmp2=indx[i];
						indx[i]=indx[i2];
							indx[i2]=tmp2;
						}			
					}
				}
				Qmax=Q_temp[0];
	
				//cout<<'\n'<<"Q: "; for(i=0;i<m;i++) cout<<Q_temp[i]<<" ";
				//cout<<'\n'<<"Qmax: "<<Qmax<<" ";
	
				row=cls_idx*S+s; 
				Q[row*cls_sz+a%cls_sz]=(1-alpha)*Q[row*cls_sz+a%cls_sz] + alpha*(R+gamma*Qmax);
	
				//cout<<'\n'<<"Q: "; cout<<Q[row*cls_sz+a%cls_sz]<<" ";
				//cout<<'\n'<<"Q_cnt: "; cout<<Q_cnt[row*cls_sz+a%cls_sz]<<" ";
			}
			
			//s=s_new;
			num=rand()%100;
			if(num<=eps*100) a_new=rand()%m;
			else a_new=indx[0]; //index of the best new action for a fixed s_new

			//if(!Qmax) a_new=rand()%m;
			//else a_new=indx[0]; 

			//cout<<'\n'<<"a_new: "<<a_new;

			rnd++;
			//cout<<'\n'<<"rnd: "<<rnd<<endl;

		}

		//cout<<'\n'<<"y: "; for(i=0;i<n2;i++) cout<<y[cw*n2+i]<<" ";
			
	}	

	//cout<<'\n'<<"Q: "; for(i=0;i<m;i++) cout<<Q[s*m+i]<<" "; cout<<'\n';
	//outf<<std::endl; outf.close();
}
