
//called by BP
//MAB-NS for optimized, contiguous or random clusters

void MP3(long cnt, long l) { 
	long i,i2,i3,cw,j,j3,j2,k,vidx,vidx2,cidx,flg4,stp,tmp2,s,row,num,cnt2;
	float tmp,val;
	
	//for(cw=0;cw<CW;cw++) for(i=0;i<m;i++) res_c_v[CW*(row_wt*i+0)+cw]=0; //refresh

	for(cw=0;cw<CW;cw++) {
		//cout<<'\n'<<"cw: "<<cw;
		flg4=1;
		for(k=0;k<cnt;k++) 
			if(cw==excl_cw[k]) {flg4=0; break;} //if flg4=0, CN is excluded from scheduling

		//cout<<'\n'<<'\n'<<"l: "<<l;
		//cout<<'\n'<<"val:";
		if(flg4) {
			//compute main syndrome
			for(i=0;i<m;i++) syn[i]=syn2[i]=0; //refresh
			for(j2=0;j2<m;j2++) {
				if(!boxplus) val=0;
				else val=1; //for box-plus
 				for(i=0;i<n;i++) 
					if(!l) { 
						if(!boxplus) val+=H[j2*n+i]*LR[cw*n+i];
						else 
							if(H[j2*n+i]) val*=tanh(0.5*LR[cw*n+i]);
					}
					else { 
						if(!boxplus) val+=H[j2*n+i]*pLR[cw*n+i];
						else  
							if(H[j2*n+i]) val*=tanh(0.5*pLR[cw*n+i]);
					}
				if(boxplus) val=2*atanhf(val); //for box-plus
			 
				//cout<<val<<" "<<endl;

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

				/*for(i=0;i<M;i++) 
					diff[i]=abs(val-rep[i]);
				for(i=0;i<M;i++) indx3[i]=i; //refresh
					for(i=0;i<M;i++) {
					for(i2=i+1;i2<M;i2++) {
						if(diff[i]>diff[i2]) {   
							tmp=diff[i];
							diff[i]=diff[i2];
							diff[i2]=tmp;
		
							tmp2=indx3[i];
							indx3[i]=indx3[i2];
							indx3[i2]=tmp2;
						}			
					}
				}
				syn[j2]=indx3[0];*/ //pick the closest quantization representation point index

			}

			//cout<<'\n'<<"syn: "; for(i=0;i<m;i++) cout<<syn[i]<<" "<<endl;

			//ofstream outf; outf.open("syn.txt",fstream::app); for(i=0;i<m;i++) outf<<syn2[i]<<" "; outf<<std::endl; outf.close();
		
			//mapping cluster syndrome to state value
			for(cls_idx=0;cls_idx<num_cls;cls_idx++) {
				//picking cluster synd. from main synd.
				//for(i=cls_idx*cls_sz;i<(cls_idx+1)*cls_sz;i++) //cls_sz is no. of CNs in a cluster
					//syn_cls[i-cls_idx*cls_sz]=syn[i]; 

				for(i=0;i<cls_sz;i++) 
					if(meth==4) syn_cls[i]=syn[cls_ind_mat_ran[cls_idx][i]];
					else if(meth==5) syn_cls[i]=syn[cls_ind_mat_opt[cls_idx][i]]; 
					//else if(meth==6) syn_cls[i]=syn[cls_ind_mat_contg[cls_idx][i]]; //for NIPS
					else if(meth==6) syn_cls[i]=syn[cls_ind_mat_ran[cls_idx][i]]; //for JSAIT
					else if(meth==7) syn_cls[i]=syn[cls_ind_mat_opt[cls_idx][i]];
					

				//cls_sz=cls_sz_vec[cls_idx];
				//cout<<'\n'<<"syn_cls: "; for(i=0;i<cls_sz;i++) cout<<syn_cls[i]<<" ";
				//cout<<'\n'<<"cls: "; for(i=0;i<cls_sz;i++) cout<<cls_ind_mat_opt[cls_idx][i]<<" ";

				//cluster synd. to state value
				s=0;
				for(j2=0;j2<cls_sz;j2++) 
					s+=syn_cls[j2]*pow(M,cls_sz-j2-1);
				//cout<<'\n'<<"s: "<<s<<endl;

				cls_s_crnt[cls_idx]=s;

				S=pow(M,cls_sz);

				//corresponding row of cluster in Q table
				row=cls_idx*S+s; 
				for(i=0;i<cls_sz;i++) 
					/*if(meth==4) Q_temp[i]=Q[row*cls_sz+i];
					else if(meth==5) Q_temp[i]=Q2[row*cls_sz+i];
					else if(meth==6) Q_temp[i]=Q3[row*cls_sz+i];*/
					
					//if(!DeepRL) {
						if(meth==4) Q_temp[i]=Q2[row*cls_sz+i]; //random model-free
						else if(meth==5) Q_temp[i]=Q3[row*cls_sz+i]; //optimized model-free
						//else if(meth==6) Q_temp[i]=Q[row*cls_sz+i]; //contiguous model-free	//for NIPS					
						else if(meth==6) Q_temp[i]=Q5[row*cls_sz+i]; //random model-based //for JSAIT
						else if(meth==7) Q_temp[i]=Q6[row*cls_sz+i]; //optimized model-based
					//}
					
					/*else 
						if(meth==5) { //optimized model-free
							inp.assign(syn_cls, syn_cls+cls_sz);
							const auto model_0 = fdeep::load_model("DeepRL/fdeep_model.json");
							const fdeep::tensor inp_tnsr(fdeep::tensor_shape(static_cast<std::size_t>(cls_sz)),inp); //converting vector to tensor
							const auto result = model_0.predict({inp_tnsr});
							const auto y=result[0];
    						const std::vector<float> vec = y.to_vector();
							Q_temp[i]=vec[i]; 
						}*/
				//cout<<'\n'<<"s: "<<s<<" row: "; cout<<row;
				//cout<<'\n'<<"Q_temp: "; for(i=0;i<cls_sz;i++) cout<<Q_temp[i]<<" "<<endl;
			
				//sorting descending
				for(i=0;i<cnt2;i++) indx[i]=i; //refresh
				for(i=0;i<cnt2;i++) {  
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

				//cout<<'\n'<<"cls_idx: "<<cls_idx;
				Q_cls[cls_idx]=Q_temp[0]; //best Q value of a cluster
				indx_cls_CN[cls_idx]=indx[0]; //best CN index of a cluster

				//cout<<'\n'<<"Q_temp[0]: "; cout<<Q_temp[0]<<" ";
			}

			//cout<<'\n'<<"Q_cls: "; for(i=0;i<num_cls;i++) cout<<Q_cls[i]<<" ";
			//cout<<'\n'<<"indx_cls_CN: "; for(i=0;i<num_cls;i++) cout<<indx_cls_CN[i]<<" "<<endl;
			

			//finding which cluster to pick the CN from
			/*if(!l) 
				for(i=0;i<num_cls;i++) 
					cls_s_prev[i]=-1;

			//cout<<'\n'<<"cls_s_prev: "; for(i=0;i<num_cls;i++) cout<<cls_s_prev[i]<<" ";
			//cout<<'\n'<<"cls_s_crnt: "; for(i=0;i<num_cls;i++) cout<<cls_s_crnt[i]<<" ";
			
			for(i=0;i<cnt2;i++) pick_cls[i]=-1; //refresh
			cnt2=0;
			for(i=0;i<num_cls;i++) 
				if(cls_s_crnt[i]!=cls_s_prev[i]) { //check current vs prev. cluster states
					pick_cls[cnt2]=i;
					cnt2++;
				}
			//cout<<'\n'<<"pick_cls: "; for(i=0;i<cnt2;i++) cout<<pick_cls[i]<<" ";

			for(i=0;i<num_cls;i++)
				cls_s_prev[i]=cls_s_crnt[i];

			//finding best CN of all valid clusters
			for(i=0;i<num_cls;i++) indx_cls[i]=i; //refresh
			for(i=0;i<cnt2;i++) {
				for(i2=i+1;i2<cnt2;i2++) {
					if(Q_cls[pick_cls[i]]<Q_cls[pick_cls[i2]]) {   
						tmp=Q_cls[pick_cls[i]];
						Q_cls[pick_cls[i]]=Q_cls[pick_cls[i2]];
						Q_cls[pick_cls[i2]]=tmp;
		
						tmp2=indx_cls_CN[pick_cls[i]]; //best CN index of a cluster
						indx_cls_CN[pick_cls[i]]=indx_cls_CN[pick_cls[i2]];
						indx_cls_CN[pick_cls[i2]]=tmp2;

						tmp2=indx_cls[pick_cls[i]]; //cluster index
						indx_cls[pick_cls[i]]=indx_cls[pick_cls[i2]];
						indx_cls[pick_cls[i2]]=tmp2;
					}			
				}
			}*/
			
			
			//finding best CN of all clusters
			for(i=0;i<num_cls;i++) indx_cls[i]=i; //refresh
			for(i=0;i<num_cls;i++) {
				for(i2=i+1;i2<num_cls;i2++) {
					if(Q_cls[i]<Q_cls[i2]) {   
						tmp=Q_cls[i];
						Q_cls[i]=Q_cls[i2];
						Q_cls[i2]=tmp;
		
						tmp2=indx_cls_CN[i]; //best CN index of a cluster
						indx_cls_CN[i]=indx_cls_CN[i2];
						indx_cls_CN[i2]=tmp2;

						tmp2=indx_cls[i]; //cluster index
						indx_cls[i]=indx_cls[i2];
						indx_cls[i2]=tmp2;
					}			
				}
			}
		

			if(meth==4) j=cls_ind_mat_ran[indx_cls[0]][indx_cls_CN[0]];
			else if(meth==5) j=cls_ind_mat_opt[indx_cls[0]][indx_cls_CN[0]];
			//else if(meth==6) j=cls_ind_mat_contg[indx_cls[0]][indx_cls_CN[0]]; //for NIPS 
			else if(meth==6) j=cls_ind_mat_ran[indx_cls[0]][indx_cls_CN[0]]; //for JSAIT
			else if(meth==7) j=cls_ind_mat_opt[indx_cls[0]][indx_cls_CN[0]];
			
			//cout<<'\n'<<"sh. CN RL: "<<j<<endl;
			

			for(i=0;i<row_wt;i++) { //row_wt is the no. of neighboring VNs of CN j
				vidx=vns[j*row_wt+i]; //index of ith neighboring VN of CN j
	
				//cout<<'\n'<<"vidx: "<<vidx<<endl;

				if(vidx>-1) {
					tmp=1; 
					for(k=0;k<row_wt;k++) 
						if(vns[j*row_wt+k]>-1 && k!=i) tmp*=tanh(0.5*E_v_c[CW*(row_wt*j+k)+cw]); //jth CN accumulating msgs from all neighboring VNs except vidx
						//cout<<'\n'<<"j: "<<j<<" vidx: "<<vidx<<" tmp: "<<tmp<<endl;
			
					for(k=0;k<col_wt;k++) 
						if(cns[vidx*col_wt+k]==j) {
							E_c_v[CW*(vidx*col_wt+k)+cw]=2*atanhf(tmp); //msg sent by jth CN to VN vidx
							break;
						}	
															
					///////////////////////////////////////////
					for(j2=0;j2<col_wt;j2++) { //col_wt is the no. of neighboring CNs of VN vidx
						cidx=cns[vidx*col_wt+j2]; //index of j2 th neighboring CN of VN vidx
			
						if(cidx>-1 && cidx!=j) {
							//if(meth==7) {cout<<'\n'<<"cidx: "<<cidx<<endl;}
							tmp=0; 
							for(k=0;k<col_wt;k++) 
								if(cns[vidx*col_wt+k]>-1 && cns[vidx*col_wt+k]!=cidx) {
									tmp+=E_c_v[CW*(col_wt*vidx+k)+cw]; 
									if(meth==4) ncv4++;
									else if(meth==5) ncv5++;
									else if(meth==6) ncv6++;
									else if(meth==7) ncv7++;
									else if(meth==8) ncv8++;
								} //VN vidx accumulating msg from all neighboring CNs except cidx 
							//if(meth==7) {cout<<'\n'<<"ncv7: "<<ncv7<<" dcw7: "<<dcw7<<endl;}	
										
							for(k=0;k<row_wt;k++) 
								if(vns[cidx*row_wt+k]==vidx) {E_v_c[CW*(row_wt*cidx+k)+cw]=tmp+LR[cw*n+vidx]; break;} //msg sent by VN vidx to CN cidx
									//printf("tmp2=%f\n",tmp2);	
						}
					} 
					//updating the aposteriori LLR of vidx
					tmp=0; 
					for(k=0;k<col_wt;k++) 
						if(cns[col_wt*vidx+k]>-1) {
						//if(!E_c_v[CW*(col_wt*vidx+k)+cw]) cout<<"zero ";
						tmp+=E_c_v[CW*(col_wt*vidx+k)+cw]; 
					}
					pLR[cw*n+vidx]=LR[cw*n+vidx]+tmp;
					//cout<<'\n'<<"pLR: "; for(i3=0;i3<n;i3++) cout<<pLR[cw*n+i3]<<" ";
				}
			}

		}
		//cout<<'\n'<<"E_c_v: "<<'\n'; for(j2=0;j2<n;j2++){for(k=0;k<col_wt;k++) cout<<E_c_v[CW*(col_wt*j2+k)+cw]<<" "; cout<<'\n';}
		//cout<<'\n'<<"E_v_c: "<<'\n'; for(j2=0;j2<m;j2++){for(k=0;k<row_wt;k++) cout<<E_v_c[CW*(row_wt*j2+k)+cw]<<" "; cout<<'\n';}
		
	}
	//cout<<'\n'<<"pLR: "; for(i=0;i<n;i++) for(cw=0;cw<CW;cw++) cout<<pLR[cw*n+i]<<" ";

	//find mean of m_{c_v} and LLR
	/*for(cw=0;cw<CW;cw++)
		for(j2=0;j2<n;j2++)
			for(k=0;k<col_wt;k++) {
				mu[l]+=E_c_v[CW*(col_wt*j2+k)+cw];
				mu_cnt[l]++;
			}

	for(cw=0;cw<CW;cw++)
		for(i=0;i<n;i++) {
			mu_L[l]+=LR[cw*n+i];
			mu_cnt_L[l]++;
		}*/
	
}