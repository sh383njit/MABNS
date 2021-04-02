
//Alg. 3 of Casado paper Q function based scheduling
//for Deep-RL

void MP4(long cnt, long l, PyObject *pFunc) { 
	long i,i2,i3,cw,j,j3,j2,k,vidx,vidx2,cidx,flg4,stp,tmp2,s,row,num;
	float tmp,val;

	ofstream outf; ifstream inf; 

	for(cw=0;cw<CW;cw++) {
		//cout<<'\n'<<"cw: "<<cw;
		flg4=1;
		for(k=0;k<cnt;k++) 
			if(cw==excl_cw[k]) {flg4=0; break;} //if flg4=0, CN is excluded from scheduling

		if(flg4) {
			//compute main syndrome
			for(i=0;i<m;i++) syn[i]=syn2[i]=0; //refresh
			for(j2=0;j2<m;j2++) {
				val=0;
 				for(i=0;i<n;i++) 
					if(!l) { 
						//if(y[cw*n+i]>=0) x_hat[cw*n+i]=0; else x_hat[cw*n+i]=1;
						val+=H[j2*n+i]*LR[cw*n+i];
					}
					else { 
						//if(pLR[cw*n+i]>=0) x_hat[cw*n+i]=0; else x_hat[cw*n+i]=1;
						val+=H[j2*n+i]*pLR[cw*n+i];
					}

				//for(i=0;i<n;i++) 
					//syn[j2]+=H[j2*n+i]*x_hat[cw*n+i]; //finding cw syndrome 
				//syn[j2]=syn[j2]%2;

				//syn[j2]=val;

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

			//cout<<'\n'<<"pLR: "; for(i=0;i<m;i++) cout<<pLR[i]<<" ";
			//cout<<'\n'<<"syn: "; for(i=0;i<m;i++) cout<<syn[i]<<" ";

			outf.open("NN_S_pred.txt"/*,fstream::app*/); for(i=0;i<m;i++) outf<<syn[i]<<" "; outf<<std::endl; outf.close(); 
			PyObject_CallObject(pFunc, NULL); //call the NN for inference
			inf.open("out.txt"); for(j=0;j<m;j++) inf >> Q_temp2[j];
			//cout<<'\n'<<"Q_temp2: "; for(i=0;i<m;i++) cout<<Q_temp2[i]<<" ";

			//sorting descending
			for(i=0;i<m;i++) indx2[i]=i; //refresh
			for(i=0;i<m;i++) {
				for(i2=i+1;i2<m;i2++) {
					if(Q_temp2[i]<Q_temp2[i2]) {   
						tmp=Q_temp2[i];
						Q_temp2[i]=Q_temp2[i2];
						Q_temp2[i2]=tmp;
		
						tmp2=indx2[i];
						indx2[i]=indx2[i2];
						indx2[i2]=tmp2;
					}			
				}
			}

			j=indx2[0]; //index of scheduled CN
			
			/*if(j==j_old) 
				cnt2++;

			//cout<<'\n'<<"l: "<<l<<" cnt2: "<<cnt2;

			if(cnt2>2) {
				cnt2=0;
				j=rand()%m;
			}*/

			cout<<'\n'<<"sh. CN RL: "<<j<<endl;

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
							E_c_v[CW*(vidx*col_wt+k)+cw]=2*atanhf(tmp); 
							break;
						}	
	
														
					///////////////////////////////////////////
					for(j2=0;j2<col_wt;j2++) { 
		
						cidx=cns[vidx*col_wt+j2];
			
						if(cidx>-1 && cidx!=j) {
	
							//cout<<'\n'<<"cidx: "<<cidx<<endl;
	
							tmp=0; 
							for(k=0;k<col_wt;k++) 
								if(cns[vidx*col_wt+k]>-1 && cns[vidx*col_wt+k]!=cidx) {tmp+=E_c_v[CW*(col_wt*vidx+k)+cw]; ncv7++;} 			
							for(k=0;k<row_wt;k++) 
								if(vns[cidx*row_wt+k]==vidx) {E_v_c[CW*(row_wt*cidx+k)+cw]=tmp+LR[cw*n+vidx]; break;}
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
		//cout<<'\n'<<"E_c_v: "<<'\n'; for(j2=0;j2<n;j2++){for(k=0;k<col_wt;k++) cout<<E_c_v[CW*(col_wt*j2+k)+cw]<<" "; cout<<'\n';}
		//cout<<'\n'<<"E_v_c: "<<'\n'; for(j2=0;j2<m;j2++){for(k=0;k<row_wt;k++) cout<<E_v_c[CW*(row_wt*j2+k)+cw]<<" "; cout<<'\n';}
		
	}
	//cout<<'\n'<<"pLR: "; for(i=0;i<n;i++) for(cw=0;cw<CW;cw++) cout<<pLR[cw*n+i]<<" ";
	
}
