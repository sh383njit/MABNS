
//Alg. 3 of Casado paper with GI scheduling
void MP1(int cnt, int l) { 
	long i,cw,j,j2,k,vidx,cidx,flg4;
	float tmp;

	//printf("cw,j=%d,%d\n", cw,j);
	//printf("\n");

	for(cw=0;cw<CW;cw++) {
		flg4=1;
		for(k=0;k<cnt;k++) 
			if(cw==excl_cw[k]) {flg4=0; break;} //if flg4=0, CN is excluded from scheduling

		if(flg4) { 
			if(meth==2) j=rand()%m; //index of randomly scheduled CN
			else if(meth==3) j=GI_arm[l];
			//cout<<'\n'<<"sh. CN: "<<j<<endl;

			for(i=0;i<row_wt;i++) { 

				vidx=vns[j*row_wt+i]; 

				//cout<<'\n'<<"vidx: "<<vidx<<endl;

				if(vidx>-1) {
					tmp=1; 
					for(k=0;k<row_wt;k++) 
						if(vns[j*row_wt+k]>-1 && k!=i) tmp*=tanh(0.5*E_v_c[CW*(row_wt*j+k)+cw]); 
					
					//cout<<'\n'<<"j: "<<j<<" vidx: "<<vidx<<" tmp: "<<tmp<<endl;
		
					for(k=0;k<col_wt;k++) 
						if(cns[vidx*col_wt+k]==j) {E_c_v[CW*(col_wt*vidx+k)+cw]=2*atanhf(tmp); break;} 	
												

					///////////////////////////////////////////
					for(j2=0;j2<col_wt;j2++) { 
	
						cidx=cns[vidx*col_wt+j2]; 
		
						if(cidx>-1 && cidx!=j) {

							//cout<<'\n'<<"cidx: "<<cidx<<endl;

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
								}			
							for(k=0;k<row_wt;k++) 
								if(vns[cidx*row_wt+k]==vidx) {E_v_c[CW*(row_wt*cidx+k)+cw]=tmp+LR[cw*n+vidx]; break;} 
								//printf("tmp2=%f\n",tmp2);					
						}
					} 

					//updating the aposteriori LLR of VN vidx	
					tmp=0; 
					for(k=0;k<col_wt;k++) 
						if(cns[col_wt*vidx+k]>-1) {
							//if(!E_c_v[CW*(col_wt*vidx+k)+cw]) cout<<"zero ";
							tmp+=E_c_v[CW*(col_wt*vidx+k)+cw]; 
						}

					pLR[cw*n+vidx]=LR[cw*n+vidx]+tmp;
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



