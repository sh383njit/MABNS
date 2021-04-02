
//called by RL4.cpp

void MP_RL(long j, long cw) { 
	long i,i2,j3,j2,k,vidx,vidx2,cidx,flg4,stp,tmp2;
	float tmp,val;

	//for(cw=0;cw<CW;cw++) for(i=0;i<m;i++) res_c_v[CW*(row_wt*i+0)+cw]=0; //refresh

	//j is index of scheduled CN
	//cout<<'\n'<<"sh. CN: "<<j;
	
	for(i=0;i<row_wt;i++) { 
	
		vidx=vns[j*row_wt+i]; 
		//cout<<'\n'<<"vidx: "<<vidx<<endl;

		if(vidx>-1) {
			tmp=1; 
			for(k=0;k<row_wt;k++) {
				if(vns[j*row_wt+k]>-1 && k!=i) tmp*=tanh(0.5*E_v_c[CW*(row_wt*j+k)+cw]); 
				//cout<<'\n'<<"j: "<<j<<" vidx: "<<vidx<<" tmp: "<<tmp<<endl;
				//cout<<'\n'<<"E_v_c1: "<<E_v_c[CW*(row_wt*j+k)+cw];
			}
			for(k=0;k<col_wt;k++) 
				if(cns[vidx*col_wt+k]==j) {
					//cout<<'\n'<<"tmp: "<<tmp;
					if(tmp<-0.999) tmp=-0.999;
					else if(tmp>0.999) tmp=0.999;
					val=2*atanhf(tmp); //new msg. 
					//cout<<'\n'<<"val: "<<val;
					res_c_v[CW*(row_wt*j+i)+cw]=abs(val-E_c_v[CW*(vidx*col_wt+k)+cw]); 
					//cout<<'\n'<<"res_c_v: "<<res_c_v[CW*(row_wt*j+i)+cw];
					E_c_v[CW*(vidx*col_wt+k)+cw]=val; //msg sent by jth CN to VN vidx
					break;											
				}	
									
	
			///////////////////////////////////////////
			for(j2=0;j2<col_wt;j2++) { 
		
				cidx=cns[vidx*col_wt+j2]; 
			
				if(cidx>-1 && cidx!=j) {
	
					//cout<<'\n'<<"cidx: "<<cidx<<endl;
					tmp=0; 
					for(k=0;k<col_wt;k++) 
						if(cns[vidx*col_wt+k]>-1 && cns[vidx*col_wt+k]!=cidx) tmp+=E_c_v[CW*(col_wt*vidx+k)+cw]; 			
					for(k=0;k<row_wt;k++) {
						if(vns[cidx*row_wt+k]==vidx) {E_v_c[CW*(row_wt*cidx+k)+cw]=tmp+LR[cw*n+vidx]; break;} 
							//printf("tmp2=%f\n",tmp2);
							//cout<<'\n'<<"E_v_c2: "<<E_v_c[CW*(row_wt*cidx+k)+cw];	
					}
				}
			}
				
			//updating the aposteriori LLR of VN vidx
			tmp=0; 
			for(k=0;k<col_wt;k++) 
				if(cns[col_wt*vidx+k]>-1) 
					tmp+=E_c_v[CW*(col_wt*vidx+k)+cw]; 
			pLR[cw*n+vidx]=LR[cw*n+vidx]+tmp; 
		}
	}

	//finding highest residual of CN j
	//cout<<'\n'<<"res:"<<'\n'; for(i=0;i<row_wt;i++) cout<<res_c_v[CW*(row_wt*j+i)+cw]<<" "; 
					
	//update residual
	for(i=0;i<row_wt;i++) res_c_v_srtd[CW*(row_wt*j+i)+cw]=res_c_v[CW*(row_wt*j+i)+cw]; 

	//sort residual descending
	for(i=0;i<row_wt;i++) {
		for(i2=i+1;i2<row_wt;i2++) {
			if(res_c_v_srtd[CW*(row_wt*j+i)+cw]<res_c_v_srtd[CW*(row_wt*j+i2)+cw]) {   
				tmp=res_c_v_srtd[CW*(row_wt*j+i)+cw];
				res_c_v_srtd[CW*(row_wt*j+i)+cw]=res_c_v_srtd[CW*(row_wt*j+i2)+cw];
				res_c_v_srtd[CW*(row_wt*j+i2)+cw]=tmp;
			}			
		}
	}
		
	//cout<<'\n'<<"res_srtd: "; for(i=0;i<row_wt;i++) cout<<res_c_v_srtd[CW*(row_wt*j+i)+cw]<<" "; 

	//cout<<'\n'<<"indx: "; for(i=0;i<m;i++) cout<<indx[i]<<" "; 

	//cout<<'\n'<<"E_c_v: "<<'\n'; for(j2=0;j2<n;j2++){for(k=0;k<col_wt;k++) cout<<E_c_v[CW*(col_wt*j2+k)+cw]<<" "; cout<<'\n';}
	//cout<<'\n'<<"E_v_c: "<<'\n'; for(j2=0;j2<m;j2++){for(k=0;k<row_wt;k++) cout<<E_v_c[CW*(row_wt*j2+k)+cw]<<" "; cout<<'\n';}

	//cout<<'\n'<<"pLR: "; for(i=0;i<n;i++) for(cw=0;cw<CW;cw++) cout<<pLR[cw*n+i]<<" ";
	
}
