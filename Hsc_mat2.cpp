

void Hsc_mat2() {
	int i,r,j,i2,j2,i3,j3,k,shft=1,temp;

	if(!ge) {
		//generating I matrix
		for(i=0;i<L;i++) for(i2=0;i2<L;i2++) if(i==i2) I2[i][i2]=1;
		//generate sigma (which is Allison tau)
		sigma2[0][L-1]=1;
		for(i=1;i<L;i++)  for(i2=0;i2<L-1;i2++)  if(i==i2+1) sigma2[i][i2]=1; //sigma1
		//generate sigma (which is Allison tau)
		tau2[0][L-2]=1; tau2[1][L-1]=1;
		for(i=2;i<L;i++)  for(i2=0;i2<L-2;i2++)  if(i==i2+2) tau2[i][i2]=1; //sigma1
		//cout<<'\n'<<"tau2 "<<'\n'; for(i=0;i<L;i++) {for(j=0;j<L;j++) cout<<tau2[i][j]<<" ";  cout<<'\n';}
	}

	for(i=0;i<p;i++) {//check the permutation assignment
		for(j=0;j<gama;j++) {
 			//enter a block of H
			for(i2=i*p;i2<(i+1)*p;i2++) {
				for(j2=j*p;j2<(j+1)*p;j2++) {
					//Hsc_proto[j2][i2]=H[j2][i2]; //needed for generating original code from reordered one
					if(H[j2][i2] && !ge) { //for Allison method 
						//enter Hsc //L is the Allison L
						for(i3=i2*L;i3<(i2+1)*L;i3++) for(j3=j2*L;j3<(j2+1)*L;j3++) 
							if(B[j][i]==1) Hsc_2[j3][i3]=sigma2[j3-j2*L][i3-i2*L]; //tau and I are L by L
							else if(B[j][i]==2) Hsc_2[j3][i3]=tau2[j3-j2*L][i3-i2*L];
							else Hsc_2[j3][i3]=I2[j3-j2*L][i3-i2*L];
					}
					else if(ge && mem==1) { //generalized edge spreading
						for(k=0;k<L;k++) {
							if(!B[j][i]) Hsc_proto[j2+k*m][i2+k*n]=H[j2][i2];							
							else if(B[j][i]==1) Hsc_proto[j2+(k+1)*m][i2+k*n]=H[j2][i2];
						} 
					}	
					else if(ge && mem==2) { //generalized edge spreading
						for(k=0;k<L;k++) {
							//cout<<j2+k*m<<","<<i2+k*n<<endl;
							//else if(i2+k*n>=tot_c) cout<<"inv. col"<<endl;
							if(!B[j][i]) Hsc_proto[j2+k*m][i2+k*n]=H[j2][i2];	
							else if(B[j][i]==1) Hsc_proto[j2+(k+1)*m][i2+k*n]=H[j2][i2];
							else if(B[j][i]==2) Hsc_proto[j2+(k+2)*m][i2+k*n]=H[j2][i2];
						}
					}		
				}
			}
		}
	}

	//cout<<'\n'<<"Hsc2 "<<'\n'; for(i=0;i<tot_r;i++){for(j=0;j<tot_c;j++) {cout<<Hsc_2[i][j]<<" "; if(!((j+1)%L)) cout<<" ";} cout<<'\n'; if(!((i+1)%L)) cout<<'\n';}

	//for testing
	//for(i=0;i<tot_r;i++) for(j=0;j<tot_c;j++) Hsc_proto[i][j]=Hsc_2[i][j];

	//for getting back orginal from reordered
	/*ifstream in_file; in_file.open("H0.txt"); for(i=0;i<gama*p;i++) for(j=0;j<p*p;j++) in_file >> Hsc_2[i][j];  
	ifstream in_file2; in_file2.open("H0.txt"); for(i=gama*p;i<2*gama*p;i++) for(j=p*p;j<2*p*p;j++) in_file2 >> Hsc_2[i][j];
	ifstream in_file3; in_file3.open("H1.txt"); for(i=0;i<gama*p;i++) for(j=p*p;j<2*p*p;j++) in_file3 >> Hsc_2[i][j];  
	ifstream in_file4; in_file4.open("H1.txt"); for(i=gama*p;i<2*gama*p;i++) for(j=0;j<p*p;j++) in_file4 >> Hsc_2[i][j];*/
	
	if(!ge) {
		//col reordering
		for(i=0;i<tot_r;i++)
			for(j=0;j<L;j++)
				for(k=j*n;k<(j+1)*n;k++)
					Hsc_3[i][k]=Hsc_2[i][(k%n)*L+j];
					//Hsc_3[i][(k%n)*L+j]=Hsc_2[i][k]; //to get back original from reordered

		//cout<<'\n'<<"Hsc3 "<<'\n'; for(i=0;i<tot_r;i++){for(j=0;j<tot_c;j++) {cout<<Hsc_3[i][j]<<" "; if(!((j+1)%L)) cout<<" ";} cout<<'\n'; if(!((i+1)%L)) cout<<'\n';}	

		//row reordering
		for(i=0;i<tot_c;i++)
			for(j=0;j<L;j++)
				for(k=j*m;k<(j+1)*m;k++)
					Hsc[k][i]=Hsc_3[(k%m)*L+j][i];
					//Hsc_proto[(k%m)*L+j][i]=Hsc_3[k][i]; //to get back original from reordered
		}
		//cout<<'\n'<<"Hsc "<<'\n'; for(i=0;i<L*m;i++){for(j=0;j<L*n;j++) {cout<<Hsc_proto[i][j]<<""; if(!((j+1)%p)) cout<<" "; if(!((j+1)%n)) cout<<" ";} cout<<'\n'; if(!((i+1)%p)) cout<<'\n'; else if(!((i+1)%m)) cout<<'\n';}
		
	}

