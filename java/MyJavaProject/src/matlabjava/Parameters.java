package matlabjava;

public class Parameters {
//todo: set beta
	//Serological parameters
	int maxX = 1; //number of strains
	public int maxa = 4; //number of age groups
	int maxi = 10; //Ab level to strain A
	int maxj = 1; //Ab level to strain B
	int maxk = 1; //Ab level to strain C
	int[][][][] arrSlu;  //S lookup [age][titres1][titres2][titres3]
	int[][][][][] arrIlu;  //I lookup [strain][age][titres1][titres2][titres3]
	int[][][][] arrRlu;  //R lookup [age][titres1][titres2][titres3]
	int[][][][][] arrCIlu; //CI lookup [strain][age][titres1][titres2][titres3]
	int[][][][][] arrCINaivelu;
	int[][][][][] arrCIImmlu; 
	
	//Basic epidemiological parameters
	double beta = 0.5;
	double Tg = 3.3;
	double[] Tga = new double[]{3.3, 3.3, 3.3, 3.3}; //age dependent 
	double wan = 1/(3*7); //21 days for high CTL level 
	double alpha = 0/365;
	double N = 1E+7;
	double[] s0_imm = new double[]{0, 0, 0, 0};
	double seed = 10;
	double trickle = 0;
	boolean age_flag = false;
	boolean semi_mechanistic = false;

	//Antibody boosting level
	double AbB;
	double AbB1;
	double AbB2;
	double AbB3;
	double AbB4;

	double[][] matM;
	double[][][][] arrf;//
	double[][][][][][][][] arrg;//(maxX,maxa,maxl,maxm,maxn,maxl,maxm,maxn)
	double[][][][][] arrh;
	double immune_alpha = 0;
	double immune_beta = 2.102;
	
	
	//pa.novars = max(max(max(pa.arrCIlu)));
	int novars; 
	double[][] ages;


public Parameters(){
	//arrSlu = new int[4][9][1][1];
	make_SIRlookup();
	make_M_china(); //age mixing array
	make_f_simple(); //infectivity array
	make_g_simple(); //boosting array 
	make_h_simple(); //susceptibility array
	make_age_range();
	setBoosting(4, 4, 4, 4);
}

public Parameters(int a, int i, int j, int k){
	//arrSlu = new int[4][9][1][1];
	maxa = a;
	maxi = i;
	maxj = j;
	maxk = k;
	
	make_SIRlookup();
	make_M_china(); //age mixing array
	make_f_simple(); //infectivity array
	make_g_simple(); //boosting array 
	make_h_simple(); //susceptibility array
	make_age_range();
	setBoosting(4, 4, 4, 4);
}


public void make_age_range() {
	ages = new double[][]{{0,18}, {19,39}, {40,65}, {65,100}};
}

public void setSlu(int s[][][][]){
  arrSlu = s;
}

public int[][][][] getSlu(){
	return arrSlu;
}

void make_SIlookup(){
	System.out.print("make SI lookup");
	int count = 0;//for java, array index starts at zero
	arrSlu = new int[maxa][maxi][maxj][maxk];
	for (int a = 0; a < maxa; a++){
		for(int i = 0; i < maxi; i++){
			for (int j = 0; j < maxj; j++){
				for (int k = 0; k < maxk; k++){
					arrSlu[a][i][j][k] = count;
					count++;
				}	
			}
		}
	}
	
	arrIlu = new int[maxX][maxa][maxi][maxj][maxk];
	for (int x = 0; x < maxX; x++){
		for (int a = 0; a < maxa; a++){
			for(int i = 0; i < maxi; i++){
				for (int j = 0; j < maxj; j++){
					for (int k = 0; k < maxk; k++){
						arrIlu[x][a][i][j][k] = count;
						count++;
					}
				}	
			}
		}
	}
	
	arrCIlu = new int[maxX][maxa][maxi][maxj][maxk];
	for (int x = 0; x < maxX; x++){
		for (int a = 0; a < maxa; a++){
			for(int i = 0; i < maxi; i++){
				for (int j = 0; j < maxj; j++){
					for (int k = 0; k < maxk; k++){
						arrCIlu[x][a][i][j][k] = count;
						count++;
					}
				}	
			}
		}
	}
	
	novars = count;
}

void make_SIRlookup(){
	System.out.print("make SIR lookup");
	int count = 0;//for java, array index starts at zero
	arrSlu = new int[maxa][maxi][maxj][maxk];
	for (int a = 0; a < maxa; a++){
		for(int i = 0; i < maxi; i++){
			for (int j = 0; j < maxj; j++){
				for (int k = 0; k < maxk; k++){
					arrSlu[a][i][j][k] = count;
					count++;
				}	
			}
		}
	}
	
	arrIlu = new int[maxX][maxa][maxi][maxj][maxk];
	for (int x = 0; x < maxX; x++){
		for (int a = 0; a < maxa; a++){
			for(int i = 0; i < maxi; i++){
				for (int j = 0; j < maxj; j++){
					for (int k = 0; k < maxk; k++){
						arrIlu[x][a][i][j][k] = count;
						count++;
					}
				}	
			}
		}
	}
	
	arrRlu = new int[maxa][maxi][maxj][maxk];
	for (int a = 0; a < maxa; a++){
		for(int i = 0; i < maxi; i++){
			for (int j = 0; j < maxj; j++){
				for (int k = 0; k < maxk; k++){
					arrRlu[a][i][j][k] = count;
					count++;
				}	
			}
		}
	}
	
	
	arrCIlu = new int[maxX][maxa][maxi][maxj][maxk];
	for (int x = 0; x < maxX; x++){
		for (int a = 0; a < maxa; a++){
			for(int i = 0; i < maxi; i++){
				for (int j = 0; j < maxj; j++){
					for (int k = 0; k < maxk; k++){
						arrCIlu[x][a][i][j][k] = count;
						count++;
					}
				}	
			}
		}
	}
	
	novars = count;
}

void make_M_china(){
	matM = new double[][]{
		    {10.92,0.6,0.8,0.4},
			{0.5,1.6,1,0.91},
			{0.6,0.79,1.13,0.94},
			{0.1,0.56,0.7,2.33}
			};
}


void make_f_simple(){
	arrf = new double[maxa][maxi][maxj][maxk];
	for (int a = 0; a < maxa; a++){
		for(int i = 0; i < maxi; i++){
			for (int j = 0; j < maxj; j++){
				for (int k = 0; k < maxk; k++){
					arrf[a][i][j][k] = 1;		
				}	
			}
		}
	}
}

void make_g_simple(){
	arrg = new double[maxX][maxa][maxi][maxj][maxk][maxi][maxj][maxk];
	//RandomEngine engine = new DRand();
	//Poisson pois = new Poisson(10, engine);
}

void make_h_simple() {
	arrh =new double[maxX][maxa][maxi][maxj][maxk];
}

public void setGArray(double g[][][][][][][][]){
	arrg = g;
}

public double[][][][][][][][] getGArray(){
	return arrg;
}

public void setHArray(double h[][][][][]){
	arrh = h;
}

public double[][][][][] getHArray(){
	return arrh;
}

public void setBeta(double b){
	beta = b;
}

public double getBeta(){
	return beta;
}

public void setmaxi(int i){
	maxi = i;
}

public int getmaxi(){
	return maxi;
}

public void setmaxa(int a){
	maxa = a;
}

public int getmaxa(){
	return maxa;
}

public void setAgeFlag(int i) {
    if (i == 1) {
       age_flag = true;
    } else {
    	age_flag = false;
    }
}

public boolean getAgeFlag() {
	return age_flag;
}

public void setBoosting(double ab1, double ab2, double ab3, double ab4){
	double AbB1 = ab1;
	double AbB2 = ab2;
	double AbB3 = ab3;
	double AbB4 = ab4;
}

public double[] getBoosting(){
	double[] AbBList = new double[]{AbB1,AbB2,AbB3,AbB4};
	return AbBList;
}

}

