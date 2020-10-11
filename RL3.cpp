
//called my main_RL
//for GIs

#include "init.cpp"
#include "MP_RL.cpp" //direct residual calculation of the pulled CN
//#include "MP_RL2.cpp" //calculates residual according to the NS scheme

//reinforcement learning
int RL() {
	long a,s,s_new,i,i2,k,cw,j,j2,l,strt=0,stp,STP,flg2,flg3,flg4,cnt,niter,a_new,tmp2,row,cls_idx;
	float tmp,Qmax,R,val,res_max,V1,V2;

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
		//compute syndrome s_0
		for(j2=0;j2<m;j2++) {
 			//for(i=0;i<n;i++) 
				//if(LR[cw*n+i]>=0) x_hat[cw*n+i]=0; 
				//else x_hat[cw*n+i]=1;
			val=0;
			for(i=0;i<n;i++) 
				val+=H[j2*n+i]*LR[cw*n+i]; //syndrome
				//syn[j2]+=H[j2*n+i]*x_hat[cw*n+i];

			//syn[j2]=syn[j2]%2;

			//quantizing
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
	
		//G estimation starts
		//finding a_new
		num=rand()%100;
		if(num<=eps*100) 
			a_new=rand()%m;
		else {
			for(a=0;a<m;a++) {
				for(j=0;j<M;j++)
					Gtemp[j]=G[j][i];

				//sort descending
				for(i=0;i<M;i++) 
					for(i2=i+1;i2<M;i2++) 
						if(Gtemp[i]<Gtemp[i2]) {   
							tmp=Gtemp[i];
							Gtemp[i]=Gtemp[i2];
							Gtemp[i2]=tmp;
						}			
				Gmax[a]=Gtemp[0]; //highest G for action a
			}

			//finding the best a
			for(i=0;i<m;i++) indx[i]=i; //refresh
			//sort descending
			for(i=0;i<m;i++) 
				for(i2=i+1;i2<m;i2++) 
					if(Gmax[i]<Gmax[i2]) {   
						tmp=Gmax[i];
						Gmax[i]=Gmax[i2];
						Gmax[i2]=tmp;
			
						tmp2=indx[i];
						indx[i]=indx[i2];
						indx[i2]=tmp2;
					}			
			a_new=indx[0];
		}
		s=syn[a_new]; //observed state
		//cout<<'\n'<<"s,a_new: "<<s<<","<<a_new<<endl;

		for(i=0;i<num_lmax;i++) {
			V1=V2=l=0;
			//cout<<'\n'<<"lmax_vec: "<<lmax_vec[i]<<endl;

			while(l<lmax_vec[i]) { //lmax rounds
				//choose action a
				a=a_new;
				//cout<<'\n'<<endl;
				//cout<<'\n'<<"s,a: "<<s<<","<<a<<endl;
	
				MP_RL(a,cw); //execute action a (schedule CN a)
	
				//Reward computation
				R=res_c_v_srtd[CW*(row_wt*a+0)+cw];
				V1=V1+pow(gamma,l)*R;
				V2=V2+pow(gamma,l);
				l++;
			}

			V1/=lmax_vec[i];
			V2/=lmax_vec[i];
			zeta[i]=V1/V2;
	
		}

		//sort descending
		for(i=0;i<num_lmax;i++) {
			for(i2=i+1;i2<num_lmax;i2++) {
				if(zeta[i]<zeta[i2]) {   
					tmp=zeta[i];
					zeta[i]=zeta[i2];
					zeta[i2]=zeta[i];
		
					tmp2=indx[i];
					indx[i]=indx[i2];
					indx[i2]=tmp2;
				}			
			}
		}

		//G[s][a]=(1-alpha)*G[s][a] + alpha*zeta[indx[0]];
		G[s][a]+= alpha*zeta[indx[0]];
			
	}	

	for(j=0;j<M;j++) 
     		for(i=0;i<m;i++)
			G[j][i]/=CW;

	cout<<"G"<<'\n';
	for(j=0;j<M;j++) {
     		for(i=0;i<m;i++)
        		cout<<G[j][i]<<" ";
      		cout<<"\n";
    	}

}


			//cout<<'\n'<<"syn_old: "; for(i=0;i<m;i++) cout<<syn_old[i]<<" ";
			//cout<<'\n'<<"syn: "; for(i=0;i<m;i++) cout<<syn[i]<<" ";

			//picking cluster synd. from main synd.
			/*for(i=cls_idx*cls_sz;i<(cls_idx+1)*cls_sz;i++) 
				syn_cls[i-cls_idx*cls_sz]=syn[i];  //cls_sz is no. of CNs in a cluster
			s_new=0; 
			for(j=0;j<cls_sz;j++) 
				s_new+=syn_cls[j]*pow(M,cls_sz-j-1);*/

			
			//cout<<'\n'<<"syn_new : "; for(i=0;i<cls_sz;i++) cout<<syn[i]<<" "; cout<<endl;
			//cout<<'\n'<<"xhat_new: "; for(i=0;i<n;i++) cout<<x_hat[i]<<" ";

			//cout<<'\n';		
			//cout<<'\n'<<"s_new: "<<s;

			//cout<<'\n'<<"Q: "; for(i=0;i<m;i++) cout<<Q_temp[i]<<" ";
			//cout<<'\n'<<"Qmax: "<<Qmax<<" ";

		//cout<<'\n'<<"y: "; for(i=0;i<n2;i++) cout<<y[cw*n2+i]<<" ";


