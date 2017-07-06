package matlabjava;

import java.util.Random;

public class basicSIR extends DiseaseModel {



public basicSIR() {
}

public void setStates(double[] st){
	x = st;
}


public void initParameters() {
	pars = new Parameters();
}

// update parameters
public void updateParametersBeta(double b){
	pars.beta = b;
}


private void calFOI(){
	FOI = pars.beta* x[pars.arrIlu[0][0][1][0][0]];
}

//update state derivative
public double[] getStatesSIRDeriv(){
	double[] xd = new double[x.length];
	double S = x[pars.arrSlu[0][0][0][0]];
	double I = x[pars.arrIlu[0][0][1][0][0]];
	double R = x[pars.arrRlu[0][1][0][0]];
	calFOI();
	xd[pars.arrSlu[0][0][0][0]] = pars.wan*R - FOI*S;
	xd[pars.arrIlu[0][0][1][0][0]] = FOI*S - (1/pars.Tg)*I;
	xd[pars.arrRlu[0][1][0][0]] = (1/pars.Tg)*I - pars.wan*R;
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
    basicSIR sir = new basicSIR();
    Parameters mepar = new Parameters(1,2,1,1); //new Parameters(a,i,j,k)
    sir.setParameters(mepar);
    sir.updateParametersBeta(0.4);
    sir.updateParameters("wan", 1/2*365);
    //sir.updateParameters("maxa", 1);
    //sir.updateParameters("maxi", 2);
    double[] xini = {0.99,0,0,0.01,0,0};
    sir.setStates(xini);
	double[] xder = sir.getStatesSIRDeriv();
	System.out.println("complete");
} 
}