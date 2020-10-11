

void Hsc_mat() {
	int i,j,k;
	
	//generating H0
	for(i=0;i<gama;i++) 
		for(j=i*p*J;j<(i+1)*p*J;j++) 
			for(k=0;k<E[i]*p*J;k++)
				H0[j][k]=H[j][k]; 

	/*for(i=0;i<cut_el;i++) {
		for(j=(m/cut_el)*i;j<(m/cut_el)*(i+1);j++) {
			for(k=0;k<E[i]*q;k++) {
				H0[j][k]=H[j][k]; 
	
			}
		}
	}*/

	//generating H1
	for(i=0;i<m;i++) for(j=0;j<n;j++) H1[i][j]=H[i][j]-H0[i][j]; 
	
	//putting H0 & H1 in Hsc
	for(i=0;i<L;i++) 
		for(j=i*m;j<(i+1)*m;j++) 
			for(k=i*n;k<(i+1)*n;k++) {
				Hsc[j][k]=H0[j-i*m][k-i*n]; 
				Hsc[j+m][k]=H1[j-i*m][k-i*n]; 
			}
			

}


