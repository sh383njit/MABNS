

int minn(float *in, int len) {

	int i,j; float temp;
	for(i=0;i<len;i++) 
		for(j=i+1;j<len;j++) 
			if(in[i]>in[j]) {temp=in[i]; in[i]=in[j]; in[j]=temp;}

	return in[0]; 
}
