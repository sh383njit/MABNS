
//Alg. 3 of Casado paper with residual CN scheduling
void MP_RL(long j, long cw) { 
	long i,i2,i3,j3,j2,k,k2,vidx,vidx2,cidx,cidx2,flg4,stp,tmp2;
	float tmp,val;

	//finding index of CN with highest residual
	//cout<<'\n'<<"res:"<<'\n'; for(i=0;i<m;i++) {for(j=0;j<row_wt;j++) cout<<res_c_v[CW*(row_wt*i+j)+cw]<<" "; cout<<'\n';}

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

	//cout<<'\n'<<"vn_indx: "; for(i=0;i<m;i++) cout<<vn_indx[i]<<" "; 

	/*for(i=0;i<m;i++) cn_indx[i]=i; //refresh
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
	}*/

	//cout<<'\n'<<'\n'<<"res_srtd: "; for(i=0;i<m;i++) cout<<res_c_v_srtd[CW*(row_wt*i+0)+cw]<<" "; 
	//cout<<'\n'<<"cn_indx: "; for(i=0;i<m;i++) cout<<cn_indx[i]<<" "; 

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
					if(tmp<-0.999) tmp=-0.999;
					else if(tmp>0.999) tmp=0.999;
					E_c_v[CW*(col_wt*vidx+k)+cw]=2*atanhf(tmp); //msg sent by jth CN to VN vidx
					break;
				}	

			//cout<<'\n'<<"E_c_v: "<<'\n'; for(j2=0;j2<n;j2++){for(k=0;k<col_wt;k++) cout<<E_c_v[CW*(col_wt*j2+k)+cw]<<" "; cout<<'\n';}
							
			for(i2=0;i2<row_wt;i2++) res_c_v[CW*(row_wt*j+i2)+cw]=0; //erasing residuals of updated CN

			///////////////////////////////////////////
			for(j2=0;j2<col_wt;j2++) { //col_wt is the no. of neighboring CNs of VN vidx
		
				cidx=cns[vidx*col_wt+j2]; //index of j2 th neighboring CN of VN vidx
			
				if(cidx>-1 && cidx!=j) {
	
					//cout<<'\n'<<"cidx: "<<cidx;
	
					tmp=0; 
					for(k=0;k<col_wt;k++) 
						if(cns[vidx*col_wt+k]>-1 && cns[vidx*col_wt+k]!=cidx) tmp+=E_c_v[CW*(col_wt*vidx+k)+cw]; //VN vidx accumulating msg from all neighboring CNs except cidx 			
					for(k=0;k<row_wt;k++) 
						if(vns[cidx*row_wt+k]==vidx) {E_v_c[CW*(row_wt*cidx+k)+cw]=tmp+LR[cw*n+vidx]; break;} //msg sent by VN vidx to CN cidx
							//printf("tmp2=%f\n",tmp2);	

					//cout<<'\n'<<"E_v_c: "<<'\n'; for(j3=0;j3<m;j3++){for(k=0;k<row_wt;k++) cout<<E_v_c[CW*(row_wt*j3+k)+cw]<<" "; cout<<'\n';}



					//residual update
					///////////////////////////////////////////////
					for(i2=0;i2<row_wt;i2++) { //row_wt is the no. of neighboring VNs of CN cidx
	
						vidx2=vns[cidx*row_wt+i2]; //index of i2 th neighboring VN of CN cidx
	
						//cout<<'\n'<<"vidx2: "<<vidx2;
	
						if(vidx2>-1) {
							tmp=1; 
							for(k=0;k<row_wt;k++) 
								if(vns[cidx*row_wt+k]>-1 && vns[cidx*row_wt+k]!=vidx2) 
									tmp*=tanh(0.5*E_v_c[CW*(row_wt*cidx+k)+cw]); //CN cidx accumulating msgs from all neighboring VNs except vidx2
						
							//cout<<'\n'<<"tmp: "<<tmp;
			
							for(k=0;k<col_wt;k++) 
								if(cns[vidx2*col_wt+k]==cidx) {
									val=2*atanhf(tmp); //new msg. CN cidx to VN vidx2
									//cout<<'\n'<<"val: "<<val;
									res_c_v[CW*(row_wt*cidx+i2)+cw]=abs(val-E_c_v[CW*(vidx2*col_wt+k)+cw]); //residual update of msg. CN cidx to VN vidx2		
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


		//}
		
		 //cout<<'\n'<<"pLR: "; for(i=0;i<n;i++) cout<<pLR[cw*n+i]<<" ";	
	//}

	//cw=0; cout<<'\n'<<"E_c_v: "<<'\n'; for(j2=0;j2<n;j2++){for(k=0;k<col_wt;k++) cout<<E_c_v[CW*(col_wt*j2+k)+cw]<<" "; cout<<'\n';}
	//cw=0; cout<<'\n'<<"E_v_c: "<<'\n'; for(j2=0;j2<m;j2++){for(k=0;k<row_wt;k++) cout<<E_v_c[CW*(row_wt*j2+k)+cw]<<" "; cout<<'\n';}
	//cout<<'\n'<<"res_c_v:"<<'\n'; for(i=0;i<m;i++) {for(j=0;j<row_wt;j++) cout<<res_c_v[CW*(row_wt*i+j)+cw]<<" "; cout<<'\n';}
	
}
