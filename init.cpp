

void init1() {
	long i,cw,k,strt,stp;
	
	//printf("cw: %d\n", cw);

	for(cw=0;cw<CW;cw++) {
		if(num==0) //if no shift
			for(i=0;i<n;i++) {
				LR[cw*n+i]=2*y[cw*n+i]/varn; //AWGN case, LLr= ln(Pr[x=0|y]/Pr[x=1|y])= ln [exp{-(y-1)^2/2sig^2)}/exp{-(y-(-1))^2/2sig^2)}], x=0,1 means y=1,-1
				pLR[cw*n+i]=0;
			}
		else {
			//shifting the msgs (from bottom to top)
			for(i=0;i<n-hshft;i++) {
				LR[cw*n+i]=LR[cw*n+i+hshft]; //hshft is the no. of new VNs entering window
				pLR[cw*n+i]=pLR[cw*n+i+hshft];
				for(k=0;k<col_wt;k++)
					//E_c_v[i][k][cw]=E_c_v[i+hshft][k][cw]; 
					E_c_v[CW*(col_wt*i+k)+cw]=E_c_v[CW*(col_wt*(i+hshft)+k)+cw];
			}

			strt=num*hshft;
			stp=strt+hshft;

			for(i=strt;i<stp;i++) 
				LR[cw*n+(i-strt)+(n-hshft)]=2*y[cw*n+i]/varn; //new LLR values entering the window

		} 
	}

	//cout<<'\n'<<"LR: "; for(i=0;i<n;i++) for(cw=0;cw<CW;cw++) cout<<LR[cw*n+i]<<" "; 

	//printf("LR: \n"); for(i=0;i<n;i++) printf(" %f",LR[cw*n+i]);
	//printf("y: \n"); for(i=0;i<n;i++) printf(" %f",y[cw*n+i]);
}





void init2() {
	long cw,j,k;

	//printf("cw= %d\n", cw);

	for(cw=0;cw<CW;cw++) {
		if(num==0) //if no shifts
			//for(j=0;j<m;j++) for(k=0;k<row_wt;k++) E_v_c[j][k][cw]=LR[cw][vns[j][k]];
			for(j=0;j<m;j++) 
				for(k=0;k<row_wt;k++) 
					E_v_c[CW*(row_wt*j+k)+cw]=LR[cw*n+vns[j*row_wt+k]]; //CN accumulating initial msgs

		else {
			//shifting the msgs (from bottom to top)
			for(j=0;j<m-vshft;j++) {
				for(k=0;k<row_wt;k++)
					//E_v_c[j][k][cw]=E_v_c[j+vshft][k][cw]; //vshft is the no. of new CNs entering window
					E_v_c[CW*(row_wt*j+k)+cw]=E_v_c[CW*(row_wt*(j+vshft)+k)+cw];
			}

			for(j=0;j<vshft;j++) 
				for(k=0;k<row_wt;k++)
					//E_v_c[m-vshft+j][k][cw]=LR[cw*n+vns[(m-vshft+j)*row_wt+k]]; 
					E_v_c[CW*(row_wt*(m-vshft+j)+k)+cw]=LR[cw*n+vns[(m-vshft+j)*row_wt+k]]; //CNs of window initialized with new LLR values 

		}

		for(j=0;j<n;j++) for(k=0;k<col_wt;k++) E_c_v[CW*(col_wt*j+k)+cw]=0; //erasing all prev. CN msgs.

	}


	//printf("E_v_c: \n"); for(j=0;j<n;j++) printf(" %.1f",E_v_c[j]);
	//cw=0;
	//cout<<'\n'<<"E_v_c: "<<'\n'; for(j=0;j<m;j++){for(k=0;k<row_wt;k++) cout<<E_v_c[CW*(row_wt*j+k)+cw]<<" "; cout<<'\n';}
	//printf("num= %d\n", num);
}

