//called by BP.cpp
//MAB-NS-TS-1 new

void MP5b(int cnt) { 
	long i,i2,i3,cw,j,j3,j2,k,k2,l,vidx,vidx2,cidx,cidx2,flg4,stp,tmp2;
	float tmp,tmp4,tmp3,phi,mu_c_v,mu_v_c,mean,val,var,param,mu_Lv=0;

	std::random_device rd;
    std::mt19937 e2(rd());

	for(cw=0;cw<CW;cw++) {
		flg4=1;
		for(k=0;k<cnt;k++) 
			if(cw==excl_cw[k]) {flg4=0; break;} //if flg4=0, CN is excluded from scheduling

		if(flg4) {
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
	
			for(i=0;i<row_wt;i++) { 
	
				vidx=vns[j*row_wt+i];
	
				//cout<<'\n'<<"vidx: "<<vidx<<endl;
				if(vidx>-1) {
					//tmp=1; 
					for(k=0;k<row_wt;k++) 
						if(vns[j*row_wt+k]>-1 && k==i) {
							tmp=E_v_c_mu[CW*(row_wt*j+k)+cw]; 
							//cout<<'\n'<<"tmp: "<<tmp; 
							break;
						} 
						//cout<<'\n'<<"j: "<<j<<" vidx: "<<vidx<<" tmp: "<<tmp<<endl;
					
					for(k=0;k<col_wt;k++) 
						if(cns[vidx*col_wt+k]==j) {
							phi=exp(-0.4527*pow(tmp,0.86)+0.0218);
							tmp4= 1-pow(1-phi,dc-1);
							tmp3=(log(tmp4)-0.0218)/(-1*0.4527);
							mu_c_v=pow(tmp3,1/0.86); 
							std::normal_distribution<float> dist(mu_c_v,sqrt(2*mu_c_v));
							E_c_v[CW*(col_wt*vidx+k)+cw]=dist(e2); 
							E_c_v_mu[CW*(col_wt*vidx+k)+cw]=mu_c_v;
							//cout<<'\n'<<"mu_c_v: "<<mu_c_v;
							break;
						}	

					//cout<<'\n'<<"E_c_v_mu: "<<'\n'; for(j2=0;j2<n;j2++){for(k=0;k<col_wt;k++) cout<<E_c_v_mu[CW*(col_wt*j2+k)+cw]<<" "; cout<<'\n';}
					for(i2=0;i2<row_wt;i2++) res_c_v[CW*(row_wt*j+i2)+cw]=0; 

					///////////////////////////////////////////
					for(j2=0;j2<col_wt;j2++) { 
						cidx=cns[vidx*col_wt+j2]; 
			
						if(cidx>-1 && cidx!=j) {
							//cout<<'\n'<<"cidx: "<<cidx;
							//tmp=0; 
							for(k=0;k<col_wt;k++) 
								if(cns[vidx*col_wt+k]>-1 && cns[vidx*col_wt+k]==cidx) {tmp=E_c_v_mu[CW*(col_wt*vidx+k)+cw]; break;} 
			
							for(k=0;k<row_wt;k++) 
								if(vns[cidx*row_wt+k]==vidx) {
									mu_v_c=mu_Lv+(dv-1)*tmp; 
									//std::normal_distribution<float> dist(mu_v_c,sqrt(2*mu_v_c));
									//E_v_c[CW*(row_wt*cidx+k)+cw]=dist(e2);
									E_v_c_mu[CW*(row_wt*cidx+k)+cw]=mu_v_c; //msg sent by VN vidx to CN cidx
								 	break;
								}

							//cout<<'\n'<<"E_v_c: "<<'\n'; for(j3=0;j3<m;j3++){for(k=0;k<row_wt;k++) cout<<E_v_c[CW*(row_wt*j3+k)+cw]<<" "; cout<<'\n';}



							//residual update
							///////////////////////////////////////////////
							for(i2=0;i2<row_wt;i2++) { 
	
								vidx2=vns[cidx*row_wt+i2];
								//cout<<'\n'<<"vidx2: "<<vidx2;
	
								if(vidx2>-1) {
									//tmp=1; 
									for(k=0;k<row_wt;k++) 
										if(vns[cidx*row_wt+k]>-1 && vns[cidx*row_wt+k]==vidx2) {
											tmp=E_v_c_mu[CW*(row_wt*cidx+k)+cw];
											break;
										}
						
									//cout<<'\n'<<"tmp: "<<tmp;
			
									for(k=0;k<col_wt;k++) 
										if(cns[vidx2*col_wt+k]==cidx) {
											//val=2*atanhf(tmp); //new msg. CN cidx to VN vidx2
											//cout<<'\n'<<"val: "<<val;									
											phi=exp(-0.4527*pow(tmp,0.86)+0.0218);
											tmp4= 1-pow(1-phi,dc-1);
											tmp3=(log(tmp4)-0.0218)/(-1*0.4527);
											mu_c_v=pow(tmp3,1/0.86); 

											//cout<<'\n'<<"mu_c_v: "<<mu_c_v;
											std::normal_distribution<float> dist(mu_c_v,sqrt(2*mu_c_v));
											val=dist(e2); //msg sent by CN cidx to VN vidx2
											//res_c_v[CW*(row_wt*cidx+i2)+cw]=abs(val-E_c_v[CW*(vidx2*col_wt+k)+cw]); 	

											//mean=mu_c_v-E_c_v_mu[CW*(vidx2*col_wt+k)+cw]; //residual mean= the difference between means of new and current message 

											mean=val-E_c_v[CW*(vidx2*col_wt+k)+cw]; 

											var=2*mu_c_v+2*E_c_v_mu[CW*(vidx2*col_wt+k)+cw]; //variance of residual= sum of the variance of messages
											param=pow(mean/sqrt(var),2); //non-centrality parameter (mean/s.d.)^2
											std::vector<matlab::data::Array> args({factory.createScalar<float>(1),factory.createScalar<float>(param)});
											matlab::data::TypedArray<float> result = matlabPtr->feval(u"ncx2rnd",args); //sampling from non-central chi^2 distr.
											res_c_v[CW*(row_wt*cidx+i2)+cw]=result[0]; //sampled residual of msg. CN cidx to VN vidx2	

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
		}
			
	}
 
	//cw=0; cout<<'\n'<<"E_c_v: "<<'\n'; for(j2=0;j2<n;j2++){for(k=0;k<col_wt;k++) cout<<E_c_v[CW*(col_wt*j2+k)+cw]<<" "; cout<<'\n';}
	//cw=0; cout<<'\n'<<"E_v_c: "<<'\n'; for(j2=0;j2<m;j2++){for(k=0;k<row_wt;k++) cout<<E_v_c[CW*(row_wt*j2+k)+cw]<<" "; cout<<'\n';}
	//cout<<'\n'<<"res_c_v:"<<'\n'; for(i=0;i<m;i++) {for(j=0;j<row_wt;j++) cout<<res_c_v[CW*(row_wt*i+j)+cw]<<" "; cout<<'\n';}
	
}
