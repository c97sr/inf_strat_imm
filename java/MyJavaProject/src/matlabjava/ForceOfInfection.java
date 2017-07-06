package matlabjava;

public class ForceOfInfection {
    //fields
	double[][] matFOI; //force of infection matrix
	Parameters pars;
	
	public void ForceOfInfection(){
		
	}
	
	public void addParameters(Parameters p){
		pars = p;
		matFOI = new double[pars.maxX][pars.maxa];
	}
	
	
	public void calmatFOI(double[] y){
		for (int X=0;X<pars.maxX;X++) {
			for (int a=0;a<pars.maxa;a++) {	
				double effective_contact = 0; //matM x (f x I)
				for (int b=0;b<pars.maxa;b++) {	
					double infectiousness = 0; //f x I
					for (int i=0;i<pars.maxi;i++) {
						for (int j=0;j<pars.maxj;j++) {
							for (int k=0;k<pars.maxk;k++) {
								infectiousness = infectiousness + pars.arrf[b][i][j][k]*y[pars.arrIlu[X][b][i][j][k]];
							}			
						}			
					}
					effective_contact = effective_contact + infectiousness*pars.matM[a][b];//b infects a
				}
				matFOI[X][a] = effective_contact*pars.beta; //beta x matM x (f x I)
			}
		}
	}
	
	public double[][] getmatFOI() {
		return matFOI;
	}

		
}
