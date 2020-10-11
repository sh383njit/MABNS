
//called my main_RL
//uses CN clustering

#include "init.cpp"
#include "MP_RL.cpp" //direct residual calculation of the pulled CN
//#include "MP_RL2.cpp" //calculates residual according to the NS scheme

//reinforcement learning
int RL() {
	long a,s,s_new,i,i2,k,cw,j,j2,l,strt=0,stp,STP,flg2,flg3,flg4,cnt,niter,rnd,a_new,tmp2,row,cls_idx;
	float tmp,Qmax,R,val,res_max;

	ofstream outf,outf2,outf3,outf4; 

	//E_v_c[j][i][cw]: msg accumulated by CN j from VN i for codeword stream cw, E_v_c[j][i][cw]=E_v_c[CW*(no. of cols *j + i)+cw]
	//E_c_v[i][j][cw]: msg accumulated by VN i from CN j for codeword stream cw
	
	//for(i=0;i<CW;i++) excl_cw[i]=-1; //refresh

	num=0;
	//initialize BP decoder
	init1();
	init2(); 

	//cout<<'\n'<<"S: "<<S<<endl;
		
	for(cw=0;cw<CW;cw++) {
		//cout<<'\n'<<"cw: "<<cw;
		//compute observed syndrome
		for(j2=0;j2<m;j2++) {
 			//for(i=0;i<n;i++) 
				//if(LR[cw*n+i]>=0) x_hat[cw*n+i]=0; 
				//else x_hat[cw*n+i]=1;

			val=0;
			for(i=0;i<n;i++) 
				val+=H[j2*n+i]*LR[cw*n+i]; //syndrome
				//syn[j2]+=H[j2*n+i]*x_hat[cw*n+i];

			//syn[j2]=syn[j2]%2;

			if(M==2) {
				if(val>=0) syn[j2]=0; else syn[j2]=1;
			}	

			if(M==4) {
				if(val<=qstrt) syn[j2]=0; 
				else if(val>qstrt && val<qstrt+gap) syn[j2]=1;
				else if(val>=qstrt+gap && val<qstrt+2*gap) syn[j2]=2;
				else if(val>=qstrt+2*gap) syn[j2]=3;
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

			//cout<<'\n'<<endl;
			//cout<<'\n'<<"s,a: "<<s<<","<<a;

			//finding the current cluster index cls_idx of CN 'a'
			for(j=0;j<num_cls;j++) 
				if(a>=j*cls_sz && a<(j+1)*cls_sz) {
					cls_idx=j;
					break;
				}
			for(i=cls_idx*cls_sz;i<(cls_idx+1)*cls_sz;i++) 
					syn_cls[i-cls_idx*cls_sz]=syn[i]; 	
			s=0; 
			for(j=0;j<cls_sz;j++) 
				s+=syn_cls[j]*pow(M,cls_sz-j-1);


			MP_RL(a,cw); //execute action a (schedule CN a)

			outf2.open("NN_a.txt",fstream::app); outf2<<a<<" ";  outf2<<std::endl; outf2.close(); //save

			//Reward computation
			R=res_c_v_srtd[CW*(row_wt*a+0)+cw];

			outf3.open("NN_R.txt",fstream::app); outf3<<R<<" "; outf3<<std::endl; outf3.close(); //save
			//cout<<'\n'<<"R: "<<R;

			outf.open("NN_S.txt",fstream::app); for(i=0;i<m;i++) outf<<syn[i]<<" "; outf<<std::endl; outf.close(); //save

			//compute new main syndrome
			for(j2=0;j2<m;j2++) {
 				//for(i=0;i<n;i++) 
					//if(pLR[cw*n+i]>=0) x_hat[cw*n+i]=0; 
					//else x_hat[cw*n+i]=1;

				syn_old[j2]=syn[j2];

				val=0;
				for(i=0;i<n;i++) 
					val+=H[j2*n+i]*pLR[cw*n+i];
					//syn[j2]+=H[j2*n+i]*x_hat[cw*n+i];

				if(M==2) {
					if(val>=0) syn[j2]=0; else syn[j2]=1;
				}	

				if(M==4) {
					if(val<=qstrt) syn[j2]=0; 
					else if(val>qstrt && val<qstrt+gap) syn[j2]=1;
					else if(val>=qstrt+gap && val<qstrt+2*gap) syn[j2]=2;
					else if(val>=qstrt+2*gap) syn[j2]=3;
			 	}

				//syn[j2]=syn[j2]%2;
			}

			//cout<<'\n'<<"syn_old: "; for(i=0;i<m;i++) cout<<syn_old[i]<<" ";
			//cout<<'\n'<<"syn: "; for(i=0;i<m;i++) cout<<syn[i]<<" ";

			outf4.open("NN_S2.txt",fstream::app); for(i=0;i<m;i++) outf4<<syn[i]<<" ";  outf4<<std::endl; outf4.close(); //save

			//picking cluster synd. from main synd.
			for(i=cls_idx*cls_sz;i<(cls_idx+1)*cls_sz;i++) 
				syn_cls[i-cls_idx*cls_sz]=syn[i]; //cls_sz is no. of CNs in a cluster
			s_new=0; 
			for(j=0;j<cls_sz;j++) 
				s_new+=syn_cls[j]*pow(M,cls_sz-j-1);
			
			//cout<<'\n'<<"syn_new : "; for(i=0;i<cls_sz;i++) cout<<syn[i]<<" "; cout<<endl;
			//cout<<'\n'<<"xhat_new: "; for(i=0;i<n;i++) cout<<x_hat[i]<<" ";

			//cout<<'\n';		
			//cout<<'\n'<<"s_new: "<<s;

			//finding action that gives max. Q value in state s_new
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

			//Q[row*cls_sz+a%cls_sz]+=R;	
			//Q_cnt[row*cls_sz+a%cls_sz]++;

			//cout<<'\n'<<"Q: "; cout<<Q[row*cls_sz+a%cls_sz]<<" ";
			//cout<<'\n'<<"Q_cnt: "; cout<<Q_cnt[row*cls_sz+a%cls_sz]<<" ";


			//s=s_new;
			num=rand()%100;
			if(num<=eps*100) a_new=rand()%m;
			else a_new=indx[0]; //index of the best new action for a fixed s_new

			//if(!Qmax) a_new=rand()%m;
			//else a_new=indx[0]; 

			//cout<<'\n'<<"a_new: "<<a_new;

			rnd++;

		}

		//cout<<'\n'<<"y: "; for(i=0;i<n2;i++) cout<<y[cw*n2+i]<<" ";
			
	}	

	//cout<<'\n'<<"Q: "; for(i=0;i<m;i++) cout<<Q[s*m+i]<<" "; cout<<'\n';

}
