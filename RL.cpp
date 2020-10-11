
//called my main_RL

#include "init.cpp"
#include "MP_RL.cpp"

//reinforcement learning
int RL() {
	long a,s,s_new,i,i2,k,cw,j,l,strt=0,stp,STP,flg2,flg3,flg4,cnt,niter,rnd,a_new,tmp2;
	float tmp,Qmax,R,val,res_max;

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
		for(j=0;j<m;j++) {
 			for(i=0;i<n;i++) 
				if(LR[cw*n+i]>=0) x_hat[cw*n+i]=0; 
				else x_hat[cw*n+i]=1;
			for(i=0;i<n;i++) 
				syn[j]+=H[j*n+i]*LR[cw*n+i]; //finding cw syndrome after each iter.
			//syn[j]=syn[j]%2;
			if(syn[j]>=0) syn[j]=0;
			else syn[j]=1;
		}

		//mapping syndrome vector to state value
		s=0;
		for(j=0;j<m;j++) 
			s+=syn[j]*pow(2,m-j-1); 

		//cout<<'\n'<<"syn_obs: "; for(i=0;i<m;i++) cout<<syn[i]<<" ";
		//cout<<'\n'<<"s_obs: "<<s;

		rnd=R=0;
		//Q learning
		while(rnd<rnd_max) { //each CW is allowed 100 schedules
			//choose action a
			if(!rnd) a=rand()%m;  
			else a=a_new;

			//cout<<'\n'<<endl;
			//cout<<'\n'<<"s,a: "<<s<<","<<a;

			MP_RL(a,cw); //execute action a (schedule CN a)

			//Reward computation
			R=res_c_v_srtd[CW*(row_wt*a+0)+cw];

			//cout<<'\n'<<"R: "<<R;

			//compute new syndrome
			for(j=0;j<m;j++) {
 				/*for(i=0;i<n;i++) 
					if(pLR[cw*n+i]>=0) x_hat[cw*n+i]=0; 
					else x_hat[cw*n+i]=1;*/
				for(i=0;i<n;i++) 
					syn[j]+=H[j*n+i]*pLR[cw*n+i]; //finding cw syndrome after each iter.
				//syn[j]=syn[j]%2;
				if(syn[j]>=0) syn[j]=0;
				else syn[j]=1;
			}

			//cout<<'\n'<<"syn_new: "; for(i=0;i<m;i++) cout<<syn[i]<<" "; cout<<endl;
			//cout<<'\n'<<"xhat_new: "; for(i=0;i<n;i++) cout<<x_hat[i]<<" ";	

			//observe the new state
			s_new=0; 
			for(j=0;j<m;j++) 
				s_new+=syn[j]*pow(2,m-j-1);

			 //cout<<'\n';		
			//cout<<'\n'<<"s_new: "<<s;
	
			//finding action that gives max. Q value in state s_new
			for(i=0;i<m;i++) Q_temp[i]=0; //refresh
			for(i=0;i<m;i++) 
				if(i!=a) 
					Q_temp[i]=Q[s_new*m+i];

			//cout<<'\n'<<"Q: "; for(i=0;i<m;i++) cout<<Q_temp[i]<<" "; cout<<'\n';

			//sort descending
			for(i=0;i<m;i++) indx[i]=i; //refresh
			for(i=0;i<m;i++) {
				for(i2=i+1;i2<m;i2++) {
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

			Q[s*m+a]=(1-alpha)*Q[s*m+a] + alpha*(R+gamma*Qmax);
			//Q[s*m+a]=R+gamma*Qmax;	
			//cout<<'\n'<<"Q(s,a): "; cout<<Q[s*m+a]<<" ";
			
			s=s_new;	
			num=rand()%100;
			if(num>=eps*100) a_new=rand()%m;
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
