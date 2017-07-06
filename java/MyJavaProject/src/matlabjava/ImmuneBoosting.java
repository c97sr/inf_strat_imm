package matlabjava;
import java.util.*;

public class ImmuneBoosting {

	//Fields
	double[][][][] InBoost;
	double[][][][][] OutBoost;
	double[][][][] ImmDecayS;
	double[][][][][] ImmDecayI;
	
	Parameters pars;
	int maxl; //Ab level to strain A
	int maxm; //Ab level to strain B
	int maxn; //Ab level to strain C

    
	
	public ImmuneBoosting() {
	}
	
	public void addParameters(Parameters p){
		pars = p;
		maxl = p.maxi;
		maxm = p.maxj;
		maxn = p.maxk;
		InBoost = new double[pars.maxa][pars.maxi][pars.maxj][pars.maxk];
		OutBoost = new double[pars.maxX][pars.maxa][pars.maxi][pars.maxj][pars.maxk];
		ImmDecayS = new double[pars.maxa][pars.maxi][pars.maxj][pars.maxk];
		ImmDecayI = new double[pars.maxX][pars.maxa][pars.maxi][pars.maxj][pars.maxk];
	}
	
	public void calImmBoosting(double[] y){
		InBoost = new double[pars.maxa][pars.maxi][pars.maxj][pars.maxk];
		OutBoost = new double[pars.maxX][pars.maxa][pars.maxi][pars.maxj][pars.maxk];
		for (int a=0;a<pars.maxa;a++) {
			for (int i=0;i<pars.maxi;i++) {
				for (int j=0;j<pars.maxj;j++) {
					for (int k=0;k<pars.maxk;k++) {
						//Immune boosting
						for (int X=0;X<pars.maxX;X++) {
							for (int l=0;l<maxl;l++) {
								for (int m=0;m<maxm;m++) {
									for (int n=0;n<maxn;n++) {
										double boosting_toS= pars.arrg[X][a][l][m][n][i][j][k];
										InBoost[a][i][j][k] = InBoost[a][i][j][k] + y[pars.arrIlu[X][a][l][m][n]]*boosting_toS;
										double boosting_fromI= pars.arrg[X][a][i][j][k][l][m][n];
										OutBoost[X][a][i][j][k] = OutBoost[X][a][i][j][k] + y[pars.arrIlu[X][a][i][j][k]]*boosting_fromI;
									}
								}
							}
						}
						//Immune decay for I
						//Skip
						//Immune decay of S
						//Skip
					}
				}	
			}
		}
	}
	
	public double[][][][] getInBoost(){
		return InBoost;
	}
	
	public double[][][][][] getOutBoost(){
		return OutBoost;
	}
	
	public double[][][][] getImmDecayS(){
		return ImmDecayS;
	}
	
	public double[][][][][] getImmDecayI(){
		return ImmDecayI;
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		//Testing
        ImmuneBoosting imboost = new ImmuneBoosting();
        Parameters mypar = new Parameters();
        imboost.addParameters(mypar);
        Random t = new Random();
        double[] x = new double[128];
        for (int c = 1; c < 128; c++) {
            x[c] = t.nextInt(100);
        }
        imboost.calImmBoosting(x);
	}

}
