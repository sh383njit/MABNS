
//terminal lift

int Hsc2_AB() { 
	//lifts the H2 matrix
	int i,j,i2,i3,j2,lsft,sum,sum2=0,flg=1;
	
	//generating I matrix
	//for(i=0;i<J;i++) for(i2=0;i2<J;i2++) if(i==i2) I3[i][i2]=1;

	for(i=0;i<tot_r;i++) {//check the permutation assignment
		for(j=0;j<tot_c;j++) {
			//generate sigma
			lsft=0;
			if(Hsc_proto[i][j]==1) {
				while(!lsft) {
					for(i3=0;i3<J;i3++) for(i2=0;i2<J;i2++) sigma3[i3][i2]=0; //refresh
					/*if(opt || randm) {
						if(s==1) lsft=rand()%J; //for first sub-matrix 1
						if(s==2) { //for sub. mat. 2
							if(i<3*p) lsft=lsft_mat[i][j]; //read from the sub. mat. 1 case
							else lsft=rand()%J;
						}
					}
					else lsft=lsft_mat[i][j];*/
				
					lsft=rand()%J;

					for(i3=1;i3<=lsft;i3++) sigma3[lsft-i3][J-i3]=1;
					for(i3=lsft;i3<J;i3++) for(i2=0;i2<J-lsft;i2++) if(i3==i2+lsft) sigma3[i3][i2]=1;  
				}
				//cout<<'\n'<<"sigma3 "<<'\n'; for(i=0;i<J;i++){for(j=0;j<J;j++) cout<<sigma3[i][j]<<" "; cout<<'\n';}
				//replace 1 by sigma3
				for(i2=i*J;i2<(i+1)*J;i2++) 
					for(j2=j*J;j2<(j+1)*J;j2++) 
						Hsc[i2][j2]=sigma3[i2-i*J][j2-j*J]; 	

				if(opt) lsft_mat[i][j]=lsft; //storing the shift factors for future use
			}
		}
	}

	

}
