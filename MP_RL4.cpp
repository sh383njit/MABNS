//called by RL4.cpp
//RL for MAB-NS-TS-2

void MP_RL4(long j, long cw, float mu_LR) { 
	long i,i2,i3,j3,j2,k,k2,l,vidx,vidx2,cidx,cidx2,flg4,stp,tmp2,cnt,flg;
	float tmp,tmp4,tmp3,phi,mu_c_v,mu_v_c,val,var,param;

	std::random_device rd;
    std::mt19937 e2(rd());
    
    //ofstream outf; 
	//outf.open("param.txt",fstream::app);

	//for(cw=0;cw<CW;cw++) for(i=0;i<m;i++) res_c_v[CW*(row_wt*i+0)+cw]=0; //refresh
	for(i=0;i<m;i++) for(j3=0;j3<row_wt;j3++) res_c_v_srtd[CW*(row_wt*i+j3)+cw]=res_c_v[CW*(row_wt*i+j3)+cw]; 
	//cout<<'\n'<<"res_c_v_srtd:"<<'\n'; for(i=0;i<m;i++) {for(j3=0;j3<row_wt;j3++) cout<<res_c_v_srtd[CW*(row_wt*i+j3)+cw]<<" "; cout<<'\n';}

	for(j3=0;j3<m;j3++) {
		for(i=0;i<row_wt;i++) {
			for(i2=i+1;i2<row_wt;i2++) {
				if(res_c_v_srtd[CW*(row_wt*j3+i)+cw]<res_c_v_srtd[CW*(row_wt*j3+i2)+cw]) {   
					tmp=res_c_v_srtd[CW*(row_wt*j3+i)+cw];
					res_c_v_srtd[CW*(row_wt*j3+i)+cw]=res_c_v_srtd[CW*(row_wt*j3+i2)+cw];
					res_c_v_srtd[CW*(row_wt*j3+i2)+cw]=tmp;
				}			
			}
		}
	}

	//j is index of scheduled CN
	//cout<<'\n'<<"sh. CN: "<<j;

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
					E_c_v[CW*(col_wt*vidx+k)+cw]=dist(e2); //sampling msg from jth CN to VN vidx
					E_c_v_mu[CW*(col_wt*vidx+k)+cw]=mu_c_v;
					//cout<<'\n'<<"mu_c_v: "<<mu_c_v;
					break;
				}	

			//cout<<'\n'<<"E_c_v_mu: "<<'\n'; for(j2=0;j2<n;j2++){for(k=0;k<col_wt;k++) cout<<E_c_v_mu[CW*(col_wt*j2+k)+cw]<<" "; cout<<'\n';}
							
			//for(i2=0;i2<row_wt;i2++) res_c_v[CW*(row_wt*j+i2)+cw]=0; //erasing residuals of scheduled CN

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
							mu_v_c=mu_LR+(dv-1)*tmp; 
							//std::normal_distribution<float> dist(mu_v_c,sqrt(2*mu_v_c));
							//E_v_c[CW*(row_wt*cidx+k)+cw]=dist(e2);
							E_v_c_mu[CW*(row_wt*cidx+k)+cw]=mu_v_c; //mean of msg sent by VN vidx to CN cidx
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

									/*std::normal_distribution<float> dist(mu_c_v,sqrt(2*mu_c_v));
									val=dist(e2); //msg sent by CN cidx to VN vidx2
									res_c_v[CW*(row_wt*cidx+i2)+cw]=abs(val-E_c_v[CW*(vidx2*col_wt+k)+cw]);*/ //residual update of msg. CN cidx to VN vidx2	

									val=mu_c_v-E_c_v_mu[CW*(vidx2*col_wt+k)+cw]; //residual mean= the difference between means of new and current message 
									var=2*mu_c_v+2*E_c_v_mu[CW*(vidx2*col_wt+k)+cw]; //variance of residual= sum of the variance of messages
									param=pow(val/sqrt(var),2); //non-centrality parameter (mean/s.d.)^2
									
									//outf<<param<<" ";  
									
									/*std::vector<matlab::data::Array> args({factory.createScalar<float>(1),factory.createScalar<float>(param)});
									matlab::data::TypedArray<float> sample = matlabPtr->feval(u"ncx2rnd", args); //sampling from non-central chi^2 distr.
									res_c_v[CW*(row_wt*cidx+i2)+cw]=sample[0]; //sampled residual of msg. CN cidx to VN vidx2	
									*/
									
									cnt=0; flg=0;
									for(k2=0;k2<74600;k2++) {
										if(param>0.9*param_vec[k2] && param<1.1*param_vec[k2]) {
											val_vec[cnt]=ncx2_vec[k2];
											cnt++;
											flg=1;
										}
										else 
											flg=0;
										if((cnt && !flg) || cnt>99)
											break;
									}	
									
									//cout<<'\n'<<"cnt:"; cout<<cnt<<endl;
										
									res_c_v[CW*(row_wt*cidx+i2)+cw]=val_vec[rand()%100]; //pseudo sampling
									//cout<<'\n'<<"res_c_v:"; cout<<res_c_v[CW*(row_wt*cidx+i2)+cw]<<endl;
									
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
	
	//outf<<std::endl; outf.close();
	
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

	//cw=0; cout<<'\n'<<"E_c_v: "<<'\n'; for(j2=0;j2<n;j2++){for(k=0;k<col_wt;k++) cout<<E_c_v[CW*(col_wt*j2+k)+cw]<<" "; cout<<'\n';}
	//cw=0; cout<<'\n'<<"E_v_c: "<<'\n'; for(j2=0;j2<m;j2++){for(k=0;k<row_wt;k++) cout<<E_v_c[CW*(row_wt*j2+k)+cw]<<" "; cout<<'\n';}
	//cout<<'\n'<<"res_c_v_srtd:"<<'\n'; for(i=0;i<m;i++) {for(j=0;j<row_wt;j++) cout<<res_c_v_srtd[CW*(row_wt*i+j)+cw]<<" "; cout<<'\n';}
	
}
