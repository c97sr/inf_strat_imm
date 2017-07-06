package matlabjava;

import java.util.Random;

public class Serology {

	double[] x; //state 
	double[] xdot; //state derivative
	Parameters pars;
    ImmuneBoosting imboost;
    ForceOfInfection foi;
	double[][][][] InBoost;
	double[][][][][] OutBoost;
	double[][][][] ImmDecayS;
	double[][][][][] ImmDecayI;
	double[][] matFOI;
	double gamma_basic = 0;

public Serology() {
}

public void setStates(double[] st){
	x = st;
}

public void setParameters(Parameters p){
	pars = p;
}

public void initParameters() {
	pars = new Parameters();
}

// update parameters
public void updateParametersG(double [][][][][][][][] g){
	pars.arrg = g;
}
public void updateParametersH(double [][][][][] h){
	pars.arrh = h;
}
public void updateParametersBeta(double b){
	pars.beta = b;
}
public void updateParametersM(double [][] m){
	pars.matM = m;
}

public void updateParameters(String pname, double pval){
	if (pname.equalsIgnoreCase("wan")) {
		pars.wan = pval;
		//System.out.println("update period 1/wan " + 1/pars.wan);
	}
	if (pname.equalsIgnoreCase("seed")) {
		pars.seed = pval;
	}
	
	if (pname.equalsIgnoreCase("maxi")) {
		pars.maxi = (int)pval;
		pars.make_SIRlookup();
		System.out.println("update total titres levels" + pars.maxi);
	}
}

public void updateParameters(String pname, double [] pval){

	if (pname.equalsIgnoreCase("s0_imm")) {
		pars.s0_imm = pval;
		//System.out.println("update s0_imm[1] " + pars.s0_imm[1]);
	}
	
}

public double[] getStatesDeriv() {
	//pars.make_SIlookup();
	xdot = new double[pars.novars];
	foi = new ForceOfInfection();
	foi.addParameters(pars);
	imboost = new ImmuneBoosting();
	imboost.addParameters(pars);
	calMOF();
	calImmBoost();
	double[] xd = updateStateDeriv();
	return xd;
}

//Using this subfunction in matlab
public double[] getStatesSIRDeriv() {
	//pars.make_SIRlookup();
	xdot = new double[pars.novars];
	foi = new ForceOfInfection();
	foi.addParameters(pars);
	imboost = new ImmuneBoosting();
	imboost.addParameters(pars);
	calMOF();
	calImmBoost();
	double[] xd = updateSIRDeriv();
	//  for(int i = 0; i<=xd.length-1; i++ ) {
	//        System.out.println(xd[i]);
	//    }
	return xd;
}

//public void updateGArray(double [][][][][][][][] arrg) {
//	pars.setGArray(arrg);
//}
//public void updateGArray(double [][][][][][][][] arrg) {
//	pars.setGArray(arrg);
//}

private void calMOF(){
	foi.calmatFOI(x);
	matFOI = foi.getmatFOI();
}

private void calImmBoost(){
	//imboost.addParameters(pars);
	imboost.calImmBoosting(x);
	InBoost = imboost.getInBoost();
	OutBoost = imboost.getOutBoost();
	ImmDecayS = imboost.getImmDecayS();
	ImmDecayI = imboost.getImmDecayI();
}

//update state derivative
private double[] updateStateDeriv(){
	double NewInfect;
	double DecayInfect;
	double RecoverInfect;
	double[] xd = new double[x.length];
	
	for (int X=0;X<pars.maxX;X++) {
		for (int a=0;a<pars.maxa;a++) {	
			for (int i=0;i<pars.maxi;i++) {
				for (int j=0;j<pars.maxj;j++) {
					for (int k=0;k<pars.maxk;k++) {
						NewInfect = 0;
						DecayInfect = 0;
						RecoverInfect = 0;
						NewInfect = x[pars.arrSlu[a][i][j][k]]*pars.arrh[X][a][i][j][k]*(matFOI[X][a]+gamma_basic);
						//System.out.println(""+pars.arrh[0][0][0][0][0]);
						//System.out.println(NewInfect);
						RecoverInfect = (1/pars.Tg)*OutBoost[X][a][i][j][k];
					    xd[ pars.arrIlu[X][a][i][j][k] ]= NewInfect + DecayInfect - RecoverInfect;
					    xd[ pars.arrCIlu[X][a][i][j][k] ] = NewInfect;
	                    xd[ pars.arrSlu[a][i][j][k] ] = xd[ pars.arrSlu[a][i][j][k] ] - NewInfect + (1/pars.Tg)*InBoost[a][i][j][k]*(1/pars.maxX) + pars.alpha*ImmDecayS[a][i][j][k];
					}
				}
			}
		}
	}
	return xd;
}

//update state derivative 
//The function is replaced by getStatesSIRDeriv() 
private double[] updateSIRDeriv(){
	double NewInfect;
	double DecayInfect;
	double RecoverInfect;
	double[] xd = new double[x.length];
	
	//for (int X=0;X<pars.maxX;X++) {
	int X = 0;
	    for (int a=0;a<pars.maxa;a++) {	
			for (int i=0;i<pars.maxi;i++) {
				for (int j=0;j<pars.maxj;j++) {
					for (int k=0;k<pars.maxk;k++) {
						NewInfect = 0;
						DecayInfect = 0;
						RecoverInfect = 0;
						
						double S = x[pars.arrSlu[a][i][j][k]];
						if (i==0 & j==0 & k==0) { // immune but undetectable individuals
							S = S - pars.s0_imm[a];
						}
						NewInfect = S*pars.arrh[X][a][i][j][k]*(matFOI[X][a]+gamma_basic);
						//NewInfect = x[pars.arrSlu[a][i][j][k]]*pars.arrh[X][a][i][j][k]*(matFOI[X][a]+gamma_basic);
						//System.out.println(""+pars.arrh[0][0][0][0][0]);
						//System.out.println(NewInfect);
						RecoverInfect = (1/pars.Tga[a])*OutBoost[X][a][i][j][k];
					    xd[ pars.arrIlu[X][a][i][j][k] ]= NewInfect + DecayInfect - RecoverInfect;
					    double Waning = x[ pars.arrRlu[a][i][j][k] ]*pars.wan; 
					    xd[ pars.arrRlu[a][i][j][k] ] = (1/pars.Tga[a])*InBoost[a][i][j][k] - Waning;
					    xd[ pars.arrSlu[a][i][j][k] ] = - NewInfect + Waning;
					    //xd[ pars.arrSlu[a][i][j][k] ] = xd[ pars.arrSlu[a][i][j][k] ] - NewInfect + x[ pars.arrRlu[a][i][j][k] ]*pars.wan;
	                    //System.out.println(pars.arrCIlu[X][a][i][j][k] );
					    xd[ pars.arrCIlu[X][a][i][j][k] ] = NewInfect;
	                    
					}
				}
			}
		}
	//}
	return xd;
}

public double getBeta(){
	return pars.beta;
}

public double[][] getFOI(){
	return matFOI;
}

public double[][][][][][][][] getGArray(){
	return pars.getGArray();
}

public double[][][][] getInBoost(){
	return InBoost;
}

public double[][][][][] getOutBoost(){
	return OutBoost;
}

public static void main(String[] args) {
	// TODO Auto-generated method stub
	//Testing
    Serology meser = new Serology();
    Parameters mepar = new Parameters();
    meser.setParameters(mepar);
    meser.updateParametersBeta(0.4);
    double[] xini = {0.146953112028553,0,0.003061868,0.006123735,0,0.006123735,0,0.003061868,0,0.26261229,0,0.002824116,0.019768812,0.008472348,0.008472348,0,0.002824116,0,0.368792732,0.002426529,0.006066322,0.007279587,0.006066322,0.001213264,0.002426529,0.002426529,0,0.122984648,0,0.003967675,0.003967675,0,0.001983838,0,0,0,1.65341E-05,0,0,0,0,0,0,0,0,3.05005E-05,0,0,0,0,0,0,0,0,3.96737E-05,0,0,0,0,0,0,0,0,1.32917E-05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    meser.setStates(xini);
    ForceOfInfection foi = new ForceOfInfection();
    foi.addParameters(mepar);
	foi.calmatFOI(xini);
	double[][] matfoi = foi.getmatFOI();
	double[] xder = meser.getStatesSIRDeriv();
	  for(int i = 0; i<=xder.length-1; i++ ) {
	        System.out.println(xder[i]);
	    }
	System.out.println("complete");

} 
}