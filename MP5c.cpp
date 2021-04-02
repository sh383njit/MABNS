//called by BP.cpp
//NS with residual sampling

void MP5c(int cnt, int cnt1, int l) { 
	long c,i,i2,i3,cw,j,j3,j2,k,k2,vidx,vidx2,cidx,cidx2,flg4,stp,tmp2,cnt2,flg;
	float tmp,val,mean,mean1,mean2,var,param;

	//printf("cw,j=%d,%d\n", cw,j);
	//printf("\n");

	//for(cw=0;cw<CW;cw++) for(i=0;i<m;i++) res_c_v[CW*(row_wt*i+0)+cw]=0; //refresh

	for(cw=0;cw<CW;cw++) {
		flg4=1;
		for(k=0;k<cnt;k++) 
			if(cw==excl_cw[k]) {flg4=0; break;} //if flg4=0, CN is excluded from scheduling

		if(flg4) { 
			//finding index of CN with highest residual
			//cout<<'\n'<<"res:"<<'\n'; for(i=0;i<m;i++) {for(j=0;j<row_wt;j++) cout<<res_c_v[CW*(row_wt*i+j)+cw]<<" "; cout<<'\n';}

			for(i=0;i<m;i++) for(j=0;j<row_wt;j++) res_c_v_srtd[CW*(row_wt*i+j)+cw]=res_c_v[CW*(row_wt*i+j)+cw]; 
			//cout<<'\n'<<"res_c_v_srtd:"<<'\n'; for(i=0;i<m;i++) {for(j=0;j<row_wt;j++) cout<<res_c_v_srtd[CW*(row_wt*i+j)+cw]<<" "; cout<<'\n';}

			for(j=0;j<m;j++) {
				for(i=0;i<row_wt;i++) {
					for(i2=i+1;i2<row_wt;i2++) {
						if(res_c_v_srtd[CW*(row_wt*j+i)+cw]<res_c_v_srtd[CW*(row_wt*j+i2)+cw]) {   
							tmp=res_c_v_srtd[CW*(row_wt*j+i)+cw];
							res_c_v_srtd[CW*(row_wt*j+i)+cw]=res_c_v_srtd[CW*(row_wt*j+i2)+cw];
							res_c_v_srtd[CW*(row_wt*j+i2)+cw]=tmp;
						}			
					}
				}
			}

			//cout<<'\n'<<"vn_indx: "; for(i=0;i<m;i++) cout<<vn_indx[i]<<" "; 

			for(i=0;i<m;i++) cn_indx[i]=i; //refresh
			for(i=0;i<m;i++) {
				for(i2=i+1;i2<m;i2++) {
					if(res_c_v_srtd[CW*(row_wt*i+0)+cw]<res_c_v_srtd[CW*(row_wt*i2+0)+cw]) {   
						tmp=res_c_v_srtd[CW*(row_wt*i+0)+cw];
						res_c_v_srtd[CW*(row_wt*i+0)+cw]=res_c_v_srtd[CW*(row_wt*i2+0)+cw];
						res_c_v_srtd[CW*(row_wt*i2+0)+cw]=tmp;
					
						tmp2=cn_indx[i];
						cn_indx[i]=cn_indx[i2];
						cn_indx[i2]=tmp2;
					}			
				}
			}

			//cout<<'\n'<<'\n'<<"res_srtd: "; for(i=0;i<m;i++) cout<<res_c_v_srtd[CW*(row_wt*i+0)+cw]<<" "; 
			//cout<<'\n'<<"cn_indx: "; for(i=0;i<m;i++) cout<<cn_indx[i]<<" "; 

			//if(!l) stp=m; else stp=1;
			//stp=1;

			if(!l) j=rand()%m;
			else j=cn_indx[0]; //index of scheduled CN
			
			//cout<<'\n'<<endl;
			//cout<<'\n'<<"sh. CN res: "<<j;
	
			for(i=0;i<row_wt;i++) { //row_wt is the no. of neighboring VNs of CN j
	
				vidx=vns[j*row_wt+i]; //index of ith neighboring VN of CN j
	
				//cout<<'\n'<<"vidx: "<<vidx<<endl;

				if(vidx>-1) {
					tmp=1; 
					for(k=0;k<row_wt;k++) 
						if(vns[j*row_wt+k]>-1 && k!=i) tmp*=tanh(0.5*E_v_c[CW*(row_wt*j+k)+cw]); 
										//cout<<'\n'<<"j: "<<j<<" vidx: "<<vidx<<" tmp: "<<tmp<<endl;
			
					for(k=0;k<col_wt;k++) 
						if(cns[vidx*col_wt+k]==j) {
							E_c_v[CW*(col_wt*vidx+k)+cw]=2*atanhf(tmp); //msg sent by jth CN to VN vidx
							break;
						}	

					//cout<<'\n'<<"E_c_v: "<<'\n'; for(j2=0;j2<n;j2++){for(k=0;k<col_wt;k++) cout<<E_c_v[CW*(col_wt*j2+k)+cw]<<" "; cout<<'\n';}
							
					for(i2=0;i2<row_wt;i2++) res_c_v[CW*(row_wt*j+i2)+cw]=0; 

					///////////////////////////////////////////
					for(j2=0;j2<col_wt;j2++) { 
		
						cidx=cns[vidx*col_wt+j2];
			
						if(cidx>-1 && cidx!=j) {
	
							//cout<<'\n'<<"cidx: "<<cidx;
	
							tmp=0; 
							for(k=0;k<col_wt;k++) 
								if(cns[vidx*col_wt+k]>-1 && cns[vidx*col_wt+k]!=cidx) {
									tmp+=E_c_v[CW*(col_wt*vidx+k)+cw]; 
									if(meth==1) ncv++;
									else if(meth==2) ncv2++;
									else if(meth==3) ncv3++;
									else if(meth==4) ncv4++;
									else if(meth==5) ncv5++;
									else if(meth==6) ncv6++;
									else if(meth==7) ncv7++;
									else if(meth==8) ncv8++;
								} 
			
							for(k=0;k<row_wt;k++) 
								if(vns[cidx*row_wt+k]==vidx) {E_v_c[CW*(row_wt*cidx+k)+cw]=tmp+LR[cw*n+vidx]; break;} 
									//printf("tmp2=%f\n",tmp2);	

							//cout<<'\n'<<"E_v_c: "<<'\n'; for(j3=0;j3<m;j3++){for(k=0;k<row_wt;k++) cout<<E_v_c[CW*(row_wt*j3+k)+cw]<<" "; cout<<'\n';}



							//residual update
							///////////////////////////////////////////////
							for(i2=0;i2<row_wt;i2++) { 
	
								vidx2=vns[cidx*row_wt+i2]; 
	
								//cout<<'\n'<<"vidx2: "<<vidx2;
	
								if(vidx2>-1) {
									tmp=1; 
									for(k=0;k<row_wt;k++) 
										if(vns[cidx*row_wt+k]>-1 && vns[cidx*row_wt+k]!=vidx2) 
											tmp*=tanh(0.5*E_v_c[CW*(row_wt*cidx+k)+cw]); 
						
									//cout<<'\n'<<"tmp: "<<tmp;
			
									for(k=0;k<col_wt;k++) 
										if(cns[vidx2*col_wt+k]==cidx) {
											val=2*atanhf(tmp); //new msg. CN cidx to VN vidx2
											//cout<<'\n'<<"val: "<<val;
											//res_c_v[CW*(row_wt*cidx+i2)+cw]=abs(val-E_c_v[CW*(vidx2*col_wt+k)+cw]); //residual update of msg. CN cidx to VN vidx2
											//res_c_v[CW*(row_wt*cidx+i2)+cw]=pow(val-E_c_v[CW*(vidx2*col_wt+k)+cw],2); //residual update (square method)

											if(cnt1==avg-1) {
												for(c=0;c<cnt1;c++) {
													mean1_vec[c]=mean1_vec[c+1];
													mean2_vec[c]=mean2_vec[c+1];
												}
											}	

											mean1_vec[cnt1]=val;
											mean2_vec[cnt1]=E_c_v[CW*(vidx2*col_wt+k)+cw];
											cnt1++;

											mean1=mean2=0;
											for(c=0;c<cnt1;c++) {
												mean1+=mean1_vec[c];
												mean2+=mean2_vec[c]; 
											}
											mean1/=cnt1;
											mean2/=cnt1;
											
											if(cnt1>avg-1) cnt1--;											

											//mean=val-E_c_v[CW*(vidx2*col_wt+k)+cw]; 
											mean=mean1-mean2; 

											//cout<<'\n'<<"mean, cnt1: "<<mean<<" "<<cnt1<<endl;

											var=2*abs(mean1)+2*abs(mean2); //variance of residual= sum of the variance of messages
											param=pow(mean/sqrt(var),2); //non-centrality parameter (mean/s.d.)^2
											/*std::vector<matlab::data::Array> args({factory.createScalar<float>(1),factory.createScalar<float>(param)});
											matlab::data::TypedArray<float> result = matlabPtr->feval(u"ncx2rnd",args); //sampling from non-central chi^2 distr.
											res_c_v[CW*(row_wt*cidx+i2)+cw]=result[0];*/ //sampled residual of msg. CN cidx to VN vidx2
											
											cnt2=0; flg=0;
											for(k2=0;k2<74600;k2++) {
												if(param>0.9*param_vec[k2] && param<1.1*param_vec[k2]) {
													val_vec[cnt2]=ncx2_vec[k2];
													cnt2++;
													flg=1;
												}
												else 
													flg=0;
												if((cnt2 && !flg) || cnt2>99)
													break;
											}	
									
											//cout<<'\n'<<"cnt2:"; cout<<cnt2<<endl;
										
											res_c_v[CW*(row_wt*cidx+i2)+cw]=val_vec[rand()%100]; //pseudo sampling

											//res_mean2[l]+=res_c_v[CW*(row_wt*cidx+i2)+cw]; res_cnt2[l]++;	//used for plotting residual histogram
										
											break;
										}
									//cout<<'\n'<<"res_c_v:"; cout<<res_c_v[CW*(row_wt*cidx+i2)+cw];
									//cout<<'\n'<<"res_c_v:"<<'\n'; for(i3=0;i3<m;i3++) {for(j3=0;j3<row_wt;j3++) cout<<res_c_v[CW*(row_wt*i3+j3)+cw]<<" "; cout<<'\n';}
								}
							}
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


			//updating the aposteriori LLR of all VNs of CN j
			/*for(i=0;i<row_wt;i++) { //row_wt is the no. of neighboring VNs of CN j
				vidx=vns[j*row_wt+i]; //index of ith neighboring VN of CN j
				//cout<<'\n'<<"vidx: "<<vidx<<endl;
				if(vidx>-1) {
					tmp=0; 
					for(k=0;k<col_wt;k++) 
						if(cns[col_wt*vidx+k]>-1) {
							//if(!E_c_v[CW*(col_wt*vidx+k)+cw]) cout<<"zero ";
							tmp+=E_c_v[CW*(col_wt*vidx+k)+cw]; 
						}

					pLR[cw*n+vidx]=LR[cw*n+vidx]+tmp;
				}
			}*/

		}
		
		 //cout<<'\n'<<"pLR: "; for(i=0;i<n;i++) cout<<pLR[cw*n+i]<<" ";	
	}

	//cw=0; cout<<'\n'<<"E_c_v: "<<'\n'; for(j2=0;j2<n;j2++){for(k=0;k<col_wt;k++) cout<<E_c_v[CW*(col_wt*j2+k)+cw]<<" "; cout<<'\n';}
	//cw=0; cout<<'\n'<<"E_v_c: "<<'\n'; for(j2=0;j2<m;j2++){for(k=0;k<row_wt;k++) cout<<E_v_c[CW*(row_wt*j2+k)+cw]<<" "; cout<<'\n';}
	//cout<<'\n'<<"res_c_v:"<<'\n'; for(i=0;i<m;i++) {for(j=0;j<row_wt;j++) cout<<res_c_v[CW*(row_wt*i+j)+cw]<<" "; cout<<'\n';}
	
}
