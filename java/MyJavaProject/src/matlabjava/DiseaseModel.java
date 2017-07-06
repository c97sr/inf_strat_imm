package matlabjava;

import java.util.Random;

public class DiseaseModel {

	double[] x; 		//state S I R CI 
	double[] xdot; 		//state derivative
	Parameters pars; 	//parameters
	double FOI;
	double[][][][] InBoost;
	double[][][][][] OutBoost;
	double[][][][] ImmDecayS;
	double[][][][][] ImmDecayI;
	double gamma_basic = 0;

public DiseaseModel() {
}



// set parameters
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
		System.out.println("update period 1/wan " + 1/pars.wan);
	}
	
	if (pname.equalsIgnoreCase("maxa")) {
		pars.maxa = (int)pval;
		pars.make_SIRlookup();
		System.out.println("update total age groups" + pars.maxa);
	}
	
	if (pname.equalsIgnoreCase("maxi")) {
		pars.maxi = (int)pval;
		pars.make_SIRlookup();
		System.out.println("update total titres levels" + pars.maxi);
	}
	
	
	if (pname.equalsIgnoreCase("seed")) {
		pars.seed = pval;
		//System.out.println("update seed number" + pars.maxi);
	}
}

public void updateParameters(String pname, int pval){
	if (pname.equalsIgnoreCase("maxa")) {
		pars.maxa = pval;
		pars.make_SIRlookup();
		System.out.println("update total age groups" + pars.maxa);
	}
	
	if (pname.equalsIgnoreCase("maxi")) {
		pars.maxa = pval;
		pars.make_SIRlookup();
		System.out.println("update total titres levels" + pars.maxi);
	}
}


public void updateParameters(String pname, double [] pval){
	if (pname.equalsIgnoreCase("s0_imm")) {
		pars.s0_imm = pval;
		System.out.println("update s0_imm[1] " + pars.s0_imm[1]);
	}
}

// set states
public void setStates(double[] st){
	x = st;
}

private void calFOI(){
	FOI = pars.beta*x[2];
}

//update state derivative
private double[] getStateDeriv(){
	//double S = x[1];
	//double I = x[2];
	//double R = x[3];
	double[] xd = new double[x.length];
	//xd[1] = pars.wan*x[3] - FOI*x[1];
	//xd[2] = pars.beta*x[1]*x[2] - (1/pars.Tg)*x[2];
	//xd[3] = (1/pars.Tg)*x[2] - pars.wan*x[3];
	//return xd;
	return xd;
}


public double getBeta(){
	return pars.beta;
}

public double getFOI(){
	return FOI;
}

public static void main(String[] args) {
	// TODO Auto-generated method stub
	//Testing
    Serology meser = new Serology();
    Parameters mepar = new Parameters();
    meser.setParameters(mepar);
    meser.updateParametersBeta(0.4);
    double[] xini = {0.146953112028553,0,0.003061868,0.006123735,0,0.006123735,0,0.003061868,0,0.26261229,0,0.002824116,0.019768812,0.008472348,0.008472348,0,0.002824116,0,0.368792732,0.002426529,0.006066322,0.007279587,0.006066322,0.001213264,0.002426529,0.002426529,0,0.122984648,0,0.003967675,0.003967675,0,0.001983838,0,0,0,1.65341E-05,0,0,0,0,0,0,0,0,3.05005E-05,0,0,0,0,0,0,0,0,3.96737E-05,0,0,0,0,0,0,0,0,1.32917E-05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    meser.setStates(xini);
    ForceOfInfection foi = new ForceOfInfection();
    foi.addParameters(mepar);
	foi.calmatFOI(xini);
	double[][] matfoi = foi.getmatFOI();
	double[] xder = meser.getStatesSIRDeriv();
	System.out.println("complete");

} 
}