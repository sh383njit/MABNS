

void recursion_Git(int cn, int cnt) {

	long i,i2,j,k,flg,tmp2;
	float t,tmp;
	
	/*Q2[0][2]=0.8; 
	Q2[1][2]=0.1;
	rw_vec[0]=16; rw_vec[1]=19; rw_vec[2]=30; rw_vec[3]=4;*/

	while(cnt<num_st) {
	
		//updating the Q2 matrix
		for(i=0;i<num_st;i++) 
			for(j=0;j<cnt;j++) 
				Q2[i][C_set[j]]=tran_prob[i][C_set[j]];
	
		for(i=0;i<num_st;i++)
     			for(j=0;j<num_st;j++) 
        	  		if(i==j) Q2_inv[i][j]=1-gamma*Q2[i][j];
        	 		else Q2_inv[i][j]=-1*gamma*Q2[i][j];
			
		for(i=0;i<num_st;i++)
     			for(j=num_st;j<2*num_st;j++) 
        	  		if(i==j-num_st) Q2_inv[i][j]=1;
        	 		else Q2_inv[i][j]=0;
	
	
		/*cout<<"\n\nQ2\n\n";
   		for(i=0;i<num_st;i++) {
     			for(j=0;j<2*num_st;j++)
        	 		cout<<"\t"<<Q2_inv[i][j];
      			cout<<"\n";
    		}*/
	
		//computing the inverse
   		for(i=0;i<num_st;i++) {
      			t=Q2_inv[i][i];
      			for(j=i;j<2*num_st;j++)
        	  		Q2_inv[i][j]=Q2_inv[i][j]/t;
      			for(j=0;j<num_st;j++) 
        	 		if(i!=j) {
        	    			t=Q2_inv[j][i];
        	    			for(k=0;k<2*num_st;k++)
        	        			Q2_inv[j][k]=Q2_inv[j][k]-t*Q2_inv[i][k];
        	  		}
   		}
	
  		/*cout<<"\n\nQ2_inv\n\n";
   		for(i=0;i<num_st;i++) {
     			for(j=num_st;j<2*num_st;j++)
        			cout<<"\t"<<Q2_inv[i][j];
      			cout<<"\n";
    		}*/
	
		for(i=0;i<num_st;i++) d_vec[i]=b_vec[i]=0; //refresh
		//computing the d and b vectors
		for(i=0;i<num_st;i++) 
			for(j=num_st;j<2*num_st;j++) {
				d_vec[i]+=Q2_inv[i][j]*rw_vec[cn][j-num_st];
				b_vec[i]+=Q2_inv[i][j];
			}
	
		//cout<<"\n\nd_vec\n\n"; for(i=0;i<num_st;i++) cout<<"\t"<<d_vec[i];
		//cout<<"\n\nb_vec\n\n"; for(i=0;i<num_st;i++) cout<<"\t"<<b_vec[i];
	
		for(i=0;i<num_st;i++)
			divv[i]=d_vec[i]/b_vec[i];
		
		//sorting for divv descending
		for(i=0;i<num_st;i++) states[i]=i;
		for(i=0;i<num_st;i++) 
			for(i2=i+1;i2<num_st;i2++) 
				if(divv[i]<divv[i2]) {   
					tmp=divv[i];
					divv[i]=divv[i2];
					divv[i2]=tmp;
	
					tmp2=states[i];
					states[i]=states[i2];
					states[i2]=tmp2;
				}			
	
		flg=0;
		for(i=0;i<num_st;i++) {
			for(j=0;j<cnt;j++) 
				if(C_set[j]==states[i]) {flg=0; break;}
				else flg=1;
			if(flg) {
				GI[cn][states[i]]=divv[i];
				C_set[cnt]=states[i];
				break;
			}
		}
		cnt++;
	}	
	
	
}	
	

