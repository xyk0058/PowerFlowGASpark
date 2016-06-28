package com.dhcc.powerflow;

import java.io.Serializable;

import com.dhcc.Global.Variable;
import com.dhcc.Model.Gene;
import com.dhcc.Model.Info;
import com.dhcc.Model.Load;
import com.dhcc.util.IOUtil;

public class NewtonPowerFlow implements Serializable {
	private double jacob[][];
	private Variable variable = null;
	
	public NewtonPowerFlow(Variable variable) {
		this.variable = variable;
	}
	
	public Variable getVariable() {
		return variable;
	}

	public void CalcJacobian() {
		Info info = variable.getPf_info();
		double B[][] = variable.getB();
		double G[][] = variable.getG();
		double Um[] = variable.getOriU();
		double Ua[] = variable.getOriTheta();
		double Pi[] = variable.getP();
		double Qi[] = variable.getQ();
		int n =  info.getN();
		int nu = 2*n+1, n2 = 2*n;
		
		//System.out.println(jacob.length+" "+jacob[0].length);
		for (int i=0; i<info.getN(); ++i) {
			double vi = Um[i], di = Ua[i], dp = 0.0, dq = 0.0;
			for (int j=0; j<info.getN(); ++j) {
				if (i==j) continue;
				double g = G[i][j], b = B[i][j], dj = Ua[j];
				//double dij = (di-dj) * Math.PI / 180.0;
				double dij = di - dj;
				double hij = -Um[i]*Um[j]*(g*Math.sin(dij)-b*Math.cos(dij));
				jacob[i][j]=hij;
				jacob[i+n][j+n]=hij;
				hij = -Um[i]*Um[j]*(g*Math.cos(dij)+b*Math.sin(dij));
				jacob[i][j+n]=hij;
				jacob[i+n][j]=-hij;
				double p = Um[j]*(g*Math.cos(dij)+b*Math.sin(dij));
				double q = Um[j]*(g*Math.sin(dij)-b*Math.cos(dij));
				dp += p;
				dq += q;
			}
			double g = G[i][i], b = B[i][i];
			jacob[i][i]=vi*dq;
			jacob[i][i+n]=-vi*dp-2*vi*vi*g;
			jacob[i+n][i]=-vi*dp;
			jacob[i+n][i+n]=-vi*dq+2*vi*vi*b;
			
			jacob[i][nu-1] = -vi*(dp+vi*g);
			jacob[i+n][nu-1] = -vi*(dq-vi*b);
			Pi[i] = vi*(dp+vi*g);
			Qi[i] = vi*(dq-vi*b);
		}

		Load load[] = variable.getLoad();
		for (int i=0; i<info.getNl(); ++i) {
			int kk=load[i].getI();
			double lp=load[i].getP();
			double lq=load[i].getQ();
			jacob[kk][nu-1] = -lp+jacob[kk][nu-1];
			jacob[kk+n][nu-1] = -lq+jacob[kk+n][nu-1];
		}
		
		Gene gene[] = variable.getGenerator();
		for (int i=0; i<info.getNg(); ++i) {
			int kk=gene[i].getI();
			double gp=gene[i].getP();
			//double gq=gene[i].getQ();
			double gq=0;
			jacob[kk][nu-1] = gp+jacob[kk][nu-1];
			jacob[kk+n][nu-1] = gq+jacob[kk+n][nu-1];
		}
		
		for (int k=0; k<info.getNg(); ++k) {
			int ii=gene[k].getI();
			if (gene[k].getJ() == variable.REF) {
				for (int j=0; j<n2; ++j) {
					if(ii == j) continue;
					jacob[ii][j]=0;
					jacob[n+ii][j]=0;
					jacob[j][ii]=0;
					jacob[j][n+ii]=0;
				}
				jacob[ii][ii]=1.0;
				jacob[n+ii][n+ii]=1.0;
				jacob[ii][nu-1]=0;
				jacob[n+ii][nu-1]=0;
			} else if (gene[k].getJ() == variable.PV) {
				for (int j=0; j<n2; ++j) {
					jacob[n+ii][j]=0;
					jacob[j][n+ii]=0;
				}
				jacob[n+ii][n+ii] = 1;
				jacob[n+ii][nu-1] = 0;
			}
		}
		
		variable.setJacob(jacob);
	}

	public void SolvEqn() {
		Info info = variable.getPf_info();
		int n2 = 2*info.getN();
		int nu = n2+1;

		for (int i=0; i<n2; ++i) {
			double d = 1.0/jacob[i][i];
//			double d = jacob[i][i];
			for (int j=i+1; j<nu; ++j) {
				double e = jacob[i][j];
				if (e==0.0) continue;
				jacob[i][j] = e*d;
//				jacob[i][j] = e / d;
			}
			if (i==n2-1)continue;
			for (int j=i+1; j<n2; ++j) {
				double e = jacob[j][i];
				if (e==0.0) continue;
				for (int k=i+1; k<nu; ++k) 
					jacob[j][k] = jacob[j][k] - jacob[i][k] * e;
			}
		}
		for (int k=1; k<n2; ++k) {
			int i = n2-k-1;
			//System.out.println(i+" "+ nu);
			for (int j=i+1; j<n2; ++j){
				//System.out.println(j +" "+ (nu-1));
				jacob[i][nu-1] = jacob[i][nu-1] - jacob[i][j]*jacob[j][nu-1];
			}
		}
	}
	
	int k=0;
	
	public boolean checkPV() {
		Gene[] gene = variable.getGenerator();
		Load[] load = variable.getLoad();
		double Um[] = variable.getOriU();
		for (int i=0; i<gene.length; ++i) {
			if (Um[gene[i].getI()]>gene[i].getMax() || Um[gene[i].getI()]<gene[i].getMin())
				return false;
		}
		for (int i=0; i<load.length; ++i) {
			if (Um[load[i].getI()]>load[i].getMax() || Um[load[i].getI()]<load[i].getMin())
				return false;
		}
		return true;
	}
	
	public boolean Run() {
		Info info = variable.getPf_info();
		int n =  info.getN(), nu = 2*n+1, n2 = 2*n;
		jacob = new double[n2][nu];
		double Um[] = variable.getOriU();
		double Ua[] = variable.getOriTheta();

		double Pi[] = new double[info.getN()];
		double Qi[] = new double[info.getN()];
		variable.setP(Pi);
		variable.setQ(Qi);
		
		while (true) {
			if(k>1000) return false;
			CalcJacobian();
			
			double error=0.0;
			for (int i=0; i<2*n; ++i)
				if (Math.abs(jacob[i][2*n])>error)
					error = Math.abs(jacob[i][2*n]);
			//System.out.println("\n"+k +" "+ error+"\n");
			if (error<info.getEps()) {
				System.out.println("The End! " + k);
				break;
			}
			++k;
			SolvEqn();
			for (int i=0; i<n; ++i) {
				//System.out.println(i +" "+ (n2+1));
				double a = jacob[i][n2];
				Ua[i] = Ua[i] - a;
				a = jacob[n+i][n2];
				Um[i] = Um[i]-(Um[i]*a);
			}
		}

		return true;
	}

	
}
