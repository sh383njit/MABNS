
//flooding schedule
void MP0(int cnt) { 
	long i,cw,j,k,vidx,cidx,flg4;
	float tmp;

	//printf("cw,j=%d,%d\n", cw,j);
	//printf("\n");

	for(cw=0;cw<CW;cw++) {
		flg4=1;
		for(k=0;k<cnt;k++) 
			if(cw==excl_cw[k]) {flg4=0; break;} //if flg4=0, cw has been recovered already
		
		if(flg4) {

			//horizontal step
			for(j=0;j<m;j++) { //m is no. of CNs
				for(i=0;i<row_wt;i++) { //row_wt is the no. of neighboring VNs of CN j
	
					vidx=vns[j*row_wt+i]; //index of ith neighboring VN of CN j
	
					if(vidx>-1) {
						tmp=1; 
						for(k=0;k<row_wt;k++) 
							if(vns[j*row_wt+k]>-1 && k!=i) tmp*=tanh(0.5*E_v_c[CW*(row_wt*j+k)+cw]); //jth CN accumulating msgs from all neighboring VNs except i
						
						//cout<<'\n'<<"j: "<<j<<" vidx: "<<vidx<<" tmp: "<<tmp<<endl;
			
						for(k=0;k<col_wt;k++) 
							if(cns[vidx*col_wt+k]==j) {E_c_v[CW*(vidx*col_wt+k)+cw]=2*atanhf(tmp); break;} //msg sent by jth CN to ith VN				
					}
				} 
			}

			//vertical step
			for(i=0;i<n;i++) { //n is no. of VNs
				for(j=0;j<col_wt;j++) { //col_wt is the no. of neighboring CNs of VN i
	
					cidx=cns[i*col_wt+j];  //index of jth neighboring CN of VN i
	
					if(cidx>-1) {
						tmp=0; 
						for(k=0;k<col_wt;k++) 
							if(cns[i*col_wt+k]>-1 && k!=j) {
								tmp+=E_c_v[CW*(col_wt*i+k)+cw]; 
								if(meth==1) ncv++;
								else if(meth==2) ncv2++;
								else if(meth==3) ncv3++;
								else if(meth==4) ncv4++;
								else if(meth==5) ncv5++;
								else if(meth==6) ncv6++;
							}//ith VN accumulating msg from all neighboring CNs except j 
		
						for(k=0;k<row_wt;k++) 
							if(vns[cidx*row_wt+k]==i) {E_v_c[CW*(row_wt*cidx+k)+cw]=tmp+LR[cw*n+i]; break;} //msg sent by ith VN to jth CN
							//printf("tmp2=%f\n",tmp2);					
					}
				} 

				//updating the aposteriori LLR
				tmp=0; 
				for(k=0;k<col_wt;k++) 
					if(cns[col_wt*i+k]>-1) tmp+=E_c_v[CW*(col_wt*i+k)+cw]; 
				pLR[cw*n+i]=LR[cw*n+i]+tmp; 

			}

		}
	}

	//cw=0; cout<<'\n'<<"E_c_v: "<<'\n'; for(j=0;j<n;j++){for(k=0;k<col_wt;k++) cout<<E_c_v[CW*(col_wt*j+k)+cw]<<" "; cout<<'\n';}
	
}


/*void horz() { 
	long i,cw,j,k,vidx;
	float tmp;

	//printf("cw,j=%d,%d\n", cw,j);
	//printf("\n");


	for(cw=0;cw<CW;cw++) {
		for(j=0;j<m;j++) { //m is no. of CNs
			for(i=0;i<row_wt;i++) { //row_wt is the no. of neighboring VNs of CN j

				vidx=vns[j*row_wt+i]; //index of ith neighboring VN of CN j

				if(vidx>-1) {
					tmp=1; 
					for(k=0;k<row_wt;k++) 
						if(vns[j*row_wt+k]>-1 && k!=i) tmp*=tanh(0.5*E_v_c[CW*(row_wt*j+k)+cw]); //jth CN accumulating msgs from all neighboring VNs except i
					
					//cout<<'\n'<<"j: "<<j<<" vidx: "<<vidx<<" tmp: "<<tmp<<endl;
		
					for(k=0;k<col_wt;k++) 
						if(cns[vidx*col_wt+k]==j) {E_c_v[CW*(vidx*col_wt+k)+cw]=2*atanhf(tmp); break;} //msg sent by jth CN to ith VN				
				}
			} 
		}
	}

	//cw=0; cout<<'\n'<<"E_c_v: "<<'\n'; for(j=0;j<n;j++){for(k=0;k<col_wt;k++) cout<<E_c_v[CW*(col_wt*j+k)+cw]<<" "; cout<<'\n';}
	
}






void vert() {
	long i,cw,j,k,cidx;
	float tmp;

	//printf("cw,i=%d,%d\n", cw,i);

	for(cw=0;cw<CW;cw++) {
		for(i=0;i<n;i++) { //n is no. of VNs
			for(j=0;j<col_wt;j++) { //col_wt is the no. of neighboring CNs of VN i

				cidx=cns[i*col_wt+j];  //index of jth neighboring CN of VN i

				if(cidx>-1) {
					tmp=0; 
					for(k=0;k<col_wt;k++) 
						if(cns[i*col_wt+k]>-1 && k!=j) tmp+=E_c_v[CW*(col_wt*i+k)+cw]; //ith VN accumulating msg from all neighboring CNs except j 
	
					for(k=0;k<row_wt;k++) 
						if(vns[cidx*row_wt+k]==i) {E_v_c[CW*(row_wt*cidx+k)+cw]=tmp+LR[cw*n+i]; break;} //msg sent by ith VN to jth CN
						//printf("tmp2=%f\n",tmp2);					
				}
			} 

			//updating the aposteriori LLR

			tmp=0; 
			for(k=0;k<col_wt;k++) 
				if(cns[col_wt*i+k]>-1) tmp+=E_c_v[CW*(col_wt*i+k)+cw];
			pLR[cw*n+i]=LR[cw*n+i]+tmp; 

		}
	}

	//cw=0; cout<<'\n'<<"E_v_c: "<<'\n'; for(j=0;j<m;j++){for(k=0;k<row_wt;k++) cout<<E_v_c[CW*(row_wt*j+k)+cw]<<" "; cout<<'\n';}
	//cout<<'\n'<<"pLR: "; for(i=0;i<n;i++) for(cw=0;cw<CW;cw++) cout<<pLR[cw*n+i]<<" ";
}*/



