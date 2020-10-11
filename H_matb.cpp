

void H_mat() {
	int i,r,j,i2,shft;
	
	//generating I matrix
	for(i=0;i<q;i++) for(i2=0;i2<q;i2++) if(i==i2) I[i][i2]=1;

	//filling H_proto with I column wise. j is col group index
	for(j=0;j<p;j++) {for(i=0;i<q;i++) for(i2=j*q;i2<(j+1)*q;i2++) H_proto[i][i2]=I[i][i2-j*q];}

	//filling H_proto with I row wise. r is row group index
	for(r=1;r<gama;r++) {for(i=r*q;i<(r+1)*q;i++) for(i2=0;i2<q;i2++) H_proto[i][i2]=I[i-r*q][i2];}

	for(r=1;r<gama;r++)
		for(j=1;j<p;j++) {
			shft=(r*j)%q;
			for(i=0;i<q;i++) for(i2=0;i2<q;i2++) sigma[i][i2]=0; //refresh sigma matrix
			for(i=0;i<q;i++) 
				for(i2=0;i2<q;i2++) 
					if(i==i2) {
						if(i2-shft>=0) sigma[i][i2-shft]=1;
						else sigma[i][q+i2-shft]=1;
					}
			for(i=r*q;i<(r+1)*q;i++) for(i2=j*q;i2<(j+1)*q;i2++) H_proto[i][i2]=sigma[i-r*q][i2-j*q];
		}

}


