package com.dhcc.powerflow;

import java.io.Serializable;

import com.dhcc.Global.Variable;
import com.dhcc.Model.Branch;
import com.dhcc.Model.Gene;
import com.dhcc.Model.Info;
import com.dhcc.Model.Load;
import com.dhcc.Model.Tran;
import com.dhcc.util.MatrixUtil;

public class ProcData implements Serializable {
	
	private Variable variable = null;
	
	public ProcData(Variable variable) {
		this.variable = variable;
	}
	public Variable getVariable() {
		return variable;
	}

	public void AdmtMatrix() {
		Info info = variable.getPf_info();
		Branch branch[] = variable.getBranch();
		Tran tran[] = variable.getTrans();
		Gene gene[] = variable.getGenerator();
		Load load[] = variable.getLoad();
		double G[][] = new double[info.getN()][info.getN()];
		double B[][] = new double[info.getN()][info.getN()];
		for (int k=0; k<gene.length; ++k) {
			G[gene[k].getI()][gene[k].getI()] = gene[k].getG();
			B[gene[k].getI()][gene[k].getI()] = gene[k].getB();
		}
		for (int k=0; k<load.length; ++k) {
			G[load[k].getI()][load[k].getI()] = load[k].getG();
			B[load[k].getI()][load[k].getI()] = load[k].getB();
		}
		double r,x,b,kt;
		int i,j;
		for (int k=0; k<info.getNb(); ++k) {
			i = branch[k].getFrom();
			j = branch[k].getTo();
			r = branch[k].getR();
			x = branch[k].getX();
			b = r*r + x*x;
			r = r/b;
			x = -x/b;
//			if (i==j) {
//				G[i][j] += r;
//				B[i][j] += x;
//				continue;
//			}
			b = branch[k].getY0();
			
			G[i][j] = G[i][j] - r;
			B[i][j] = B[i][j] - x;
			G[j][i] = G[j][i] - r;
			B[j][i] = B[j][i] - x;
			
			G[i][i] = G[i][i] + r;
			B[i][i] = B[i][i] + x + b/2;
			//B[i][i] = B[i][i] + x + b;
			G[j][j] = G[j][j] + r;
			B[j][j] = B[j][j] + x + b/2;
			//B[j][j] = B[j][j] + x + b;
		}
		
		for (int k=0; k<info.getNt(); ++k) {
			//TODO变压器支路计算有问题
			i = tran[k].getFrom();
			j = tran[k].getTo();
			r = tran[k].getR();
			x = tran[k].getX();
			b = r*r + x*x;
			r = r/b;
			x = -x/b;
			kt = tran[k].getK();
			G[j][j] += r;
			B[j][j] += x;
//			G[i][i] = G[i][i] + r;
//			B[i][i] = B[i][i] + x;
			G[i][j] = G[i][j] - r/kt;
			B[i][j] = B[i][j] - x/kt;
			G[j][i] = G[i][j];
			B[j][i] = B[i][j];
			r = r/kt/kt;x = x/kt/kt;
			//System.out.println("Lanlan " + B[j][j] + " " + r + " " + x);
			G[i][i] += r;
			B[i][i] += x;
//			G[j][j] += r;
//			B[j][j] += x;
		}
		variable.setB(B);
		variable.setG(G);
	}

	public void CalcFactor() {
		Info info = variable.getPf_info();
		double Bc[][] = variable.getB();
		
		double Bp[][] = new double[info.getN()-1][info.getN()-1];
		double Bpp[][] = new double[info.getNpq()][info.getNpq()];
		
		for (int i=0; i<info.getN()-1; ++i) 
			for (int j=0; j<info.getN()-1; ++j)
				Bp[i][j] = Bc[i][j];
		for (int i=0; i<info.getNpq(); ++i) 
			for (int j=0; j<info.getNpq(); ++j)
				Bpp[i][j] = Bc[i][j];
		double invBp[][] = MatrixUtil.Inverse(Bp);
		double invBpp[][] = MatrixUtil.Inverse(Bpp);
		variable.setInvBp(invBp);
		variable.setInvBpp(invBpp);
		variable.setBp(Bp);
		variable.setBpp(Bpp);
	}
	
	public void InitOri() {
		Info info = variable.getPf_info();
		Gene gene[] = variable.getGenerator();
		Load load[] = variable.getLoad();
		//double bus[][] = _mpc.getBus();
		double oriU[] = new double[info.getN()];
		double oriTheta[] = new double[info.getN()];
		//System.out.println("GENE: " + info.getNt() + " " + load.length);
		for (int i=0; i<info.getN(); ++i) {
			oriU[i] = 1.0; 
			oriTheta[i] = 0.0; 
		}
		for (int i=0; i<info.getNg(); ++i) {
			if (gene[i].getJ() == variable.REF) oriTheta[(int)gene[i].getI()] = variable.getRefTheta();
		}
		for (int i=0; i<info.getNg(); ++i)
			if (gene[i].getJ() == variable.PV || gene[i].getJ() == variable.REF)
				oriU[(int) gene[i].getI()] = gene[i].getV();
		for (int i=0; i<load.length; ++i)
			if (load[i].getJ() == variable.PV || load[i].getJ() == variable.REF)
				oriU[(int) load[i].getI()] = load[i].getV();
		variable.setOriTheta(oriTheta);
		variable.setOriU(oriU);
//		System.out.println("Um " + oriU.length );
//		for (int i=0; i<oriU.length; ++i) 
//			System.out.print(oriU[i] + " ");
//		System.out.println();
	}
		
	public void calcPQ() {
		Info info = variable.getPf_info();
		Gene gene[] = variable.getGenerator();
		Load load[] = variable.getLoad();
		double Pi[] = new double[info.getN()];
		double Qi[] = new double[info.getN()];
		for (int i=0; i<gene.length; ++i) {
			Pi[gene[i].getI()] = gene[i].getP();
			if (gene[i].getJ() == variable.PQ) {
				Qi[gene[i].getI()] = gene[i].getQ();
			}
		}
		for (int i=0; i<load.length; ++i) {
			Pi[load[i].getI()] = load[i].getP();
			if (load[i].getJ() == variable.PQ) {
				Qi[load[i].getI()] = load[i].getQ();
			}
		}
		System.out.println("P");
		for (int i=0; i<info.getN(); ++i) System.out.print(Pi[i] + " ");
		System.out.println("\r\nQ");
		for (int i=0; i<info.getN(); ++i) System.out.print(Qi[i] + " ");
		System.out.println();
		variable.setP(Pi);
		variable.setQ(Qi);
	}
	
	public void CalBusPQFlow() {
		Info info = variable.getPf_info();
		double B[][] = variable.getB();
		double G[][] = variable.getG();
		double Um[] = variable.getOriU();
		double Ua[] = variable.getOriTheta();
		double Pi[] = variable.getP();
		double Qi[] = variable.getQ();
		Load load[] = variable.getLoad();
		
		System.out.println("sf:" + info.getNpq() + " " + info.getN() + " " + load.length);
		
		for (int i=0;i<load.length; ++i) {
			Pi[load[i].getI()] = load[i].getP();
			Qi[load[i].getI()] = load[i].getQ();
		}
		
		for (int i=info.getNpq(); i<info.getN(); ++i) {
			double sump = 0.0, di = Ua[i];
			for (int j=0; j<info.getN(); ++j) {
				double g = G[i][j], b = B[i][j], dj = Ua[j];
				double dij = di-dj;
				sump += Um[i]*Um[j]*(g*Math.cos(dij * Math.PI / 180)+b*Math.sin(dij * Math.PI / 180));
				//sump += Um[i]*Um[j]*(g*Math.cos(dij)+b*Math.sin(dij));
			}
			Pi[i] = sump; 
		}
		
		for (int i=info.getNpq(); i<info.getN(); ++i) {
			double sumq = 0.0, di = Ua[i];
			for (int j=0; j<info.getN(); ++j) {
				double g = G[i][j], b = B[i][j], dj = Ua[j];
				double dij = di-dj;
				sumq += Um[i]*Um[j]*(g*Math.sin(dij * Math.PI / 180)-b*Math.cos(dij * Math.PI / 180));
				//sumq += Um[i]*Um[j]*(g*Math.sin(dij)-b*Math.cos(dij));
			}
			Qi[i] = sumq;
		}
		
		for (int i=0; i<info.getN(); ++i) {
			System.out.println("P"+i+": " + Pi[i] + "\tQ:" + Qi[i]);
		}
	}
	
	public void CalBusPQGFlow() {
		double[] Pi = variable.getP();
		double[] Qi = variable.getQ();
		Gene gene[] = variable.getGenerator();
		Load load[] = variable.getLoad();
		Info info = variable.getPf_info();
		for (int i=0; i<info.getN(); ++i) {
			for (int j=0; j<gene.length; ++j) {
				if (i!=gene[j].getI())continue;
				double pl = gene[j].getPl(),ql = gene[j].getQl();
				System.out.println(i+" "+pl+" "+ql);
				System.out.println((Pi[i] - pl)+" "+(Qi[i] - ql)+"\n");
				//System.out.println(i +" "+pl+ " " + (Pi[i] - pl) + " " + ql + " " + (Qi[i] - ql) );
				break;
			}
			for (int j=0; j<load.length; ++j) {
				if (i!=load[j].getI())continue;
				double pl = load[j].getPl(),ql = load[j].getQl();
				System.out.println(i+" "+pl+" "+ql);
				System.out.println((Pi[i]*100 - pl)+" "+(Qi[i]*100 - ql)+"\n");
				//System.out.println(i +" "+pl+ " " + (Pi[i] - pl) + " " + ql + " " + (Qi[i] - ql) );
				break;
			}
		}
		
	}
	
	public void CalBusFlow() {
		Info info = variable.getPf_info();
		Gene gene[] = variable.getGenerator();
		Load load[] = variable.getLoad();
		double Pi[] = variable.getP();
		double Qi[] = variable.getQ();
		double Ua[] = variable.getOriTheta();
		double Um[] = variable.getOriU();
		for (int i=0; i<info.getN(); ++i) {
			double b1=0.0,b2=0.0,c1=0.0,c2=0.0;
			for (int j=0; j<gene.length; ++j) {
				int ii = gene[j].getI();
				int kk = gene[j].getJ();
				if (i == ii && kk == variable.REF) {
					b1 = Pi[ii];
					b2 = Qi[ii];
					for (int k=0; k<load.length; ++k) {
						ii = load[k].getI();
						if (i == ii) {
							c1 = load[k].getP();
							c2 = load[k].getQ();
							b1 = b1 + c1;
							b2 = b2 + c2;
						}
					}
					break;
				}
				if (i == ii && kk == variable.PV) {
					b1 = gene[j].getP();
					b2 = Qi[ii];
					for (int k=0; k<load.length; ++k) {
						ii = load[k].getI();
						if (i == ii) {
							c1 = load[k].getP();
							c2 = load[k].getQ();
							b2 = b2 + c2;
						}
					}
				}
			}
			for (int j=0; j<load.length; ++j) {
				int ii = load[j].getI();
				if (i == ii) {
					c1 = load[j].getP();
					c2 = load[j].getQ();
					break;
				}
			}
			System.out.println(i + " " + Um[i] + " " + Ua[i] + " " + b1 + " " + b2 + " " + c1 + " " + c2);
		}
	}
	
	private double dp=0,dq=0;
	public void BranchFlowByMat() {
		Branch[] br = variable.getBranch();
		Tran[] trans = variable.getTrans();
		Info info = variable.getPf_info();
		int numbr = info.getNb(), numtr = info.getNt();
		double ph=0, qh=0;
		for (int i=0; i<numbr; ++i) {
			int from = br[i].getFrom(),to = br[i].getTo();
			double r = br[i].getR(), x = br[i].getX();
			MatPower(from,to,r,x,1);
			ph+=dp;
			qh+=dq;
		}
		for (int i=0; i<numtr; ++i) {
			int from = trans[i].getFrom(),to = trans[i].getTo();
			double r = trans[i].getR(), x = trans[i].getX();
			MatPower(from,to,r,x,trans[i].getK());
			ph+=dp;
			qh+=dq;
		}
		System.out.println("sum loss： "+ph*100 +" " + qh*100);
	}
	
	public void MatPower(int from, int to, double r, double x,double k) {
		double b = r*r+x*x;
		r = r/b;
		x = x/b;
		double[] Ua = variable.getOriTheta();
		double[] Um = variable.getOriU();
		double vi = Um[from], vj = Um[to];
		double vii = vi*vi, vjj = vj*vj, vij = vi*vj;
		double di = Ua[from]*Math.PI/180, dj = Ua[to]*Math.PI/180;
		double cdi = Math.cos(di)*Math.cos(di);
		double sdi = Math.sin(di)*Math.sin(di);
		double cdj = Math.cos(dj)*Math.cos(dj);
		double sdj = Math.sin(dj)*Math.sin(dj);
		double cij = Math.cos(di)*Math.cos(dj);
		double sij = Math.sin(di)*Math.sin(dj);
		double sici = Math.sin(di)*Math.cos(di);
		double sjcj = Math.sin(dj)*Math.cos(dj);
		

//		double real = vii*(cdi-sdi)/k/k
//				- 2*vij*(cij-sij)/k
//				+ vjj*(cdj-sdj);
//		
//		double imag = 2*vii*sici/k/k
//				- 2*vij*(Math.sin(dj)*Math.cos(di)+Math.sin(di)*Math.cos(dj))
//				+ 2*vjj*sjcj;
//		
//		dp = r*real+x*imag;
//		dq = r*imag-x*real;
		double hi = vii+k*k*vjj-2*k*vij*(cij+sij);
		dp = hi*r/k/k;
		dq = hi*x/k/k;
		System.out.println(from +" "+ to +" "+ dp*100 +" "+ dq*100);
	}
	
	public void BranchFlow() {
		Branch[] br = variable.getBranch();
		Tran[] trans = variable.getTrans();
		Info info = variable.getPf_info();
		double[] Ua = variable.getOriTheta();
		double[] Um = variable.getOriU();
		int numbr = info.getNb(), numtr = info.getNt();
		double ph=0, qh=0;
		
		for (int i=0; i<numbr; ++i) {
			int from = br[i].getFrom(),to = br[i].getTo();
			double r = br[i].getR(), x = br[i].getX();
			double b = r*r+x*x;
			//没有from=to
			if(from==to) System.out.println("from == to");
			r = r/b;
			x = -x/b;
			b = br[i].getY0();
			double dij = Ua[from] - Ua[to];
			double vi = Um[from], vj = Um[to];
			double vij = vi*vj;
			vi = vi*vi;
			vj = vj*vj;
			double cd = vij*Math.cos(dij*Math.PI/180);
			double sd = vij*Math.sin(dij*Math.PI/180);
			double pij = vi*r - r*cd - x*sd;
			double pji = vj*r - r*cd + x*sd;
			double dpb = pij+pji;
			ph = ph+dpb;
			//double qij = -vi * (b+x) + x*cd - r*sd;
			double qij = -vi * (b/2+x) + x*cd - r*sd;
			double qji = -vj * (b/2+x) + x*cd + r*sd;
			double dqb = qij+qji;
			qh = qh+dqb;
			System.out.println(from +" "+ to +" "+ pij*100 +" "+ qij*100 +" "+ pji*100 +" "+ qji*100 +" "+ dpb*100 +" "+ dqb*100);
		}
		for (int i=0; i<numtr; ++i) {
			int from = trans[i].getFrom(),to = trans[i].getTo();
			double r = trans[i].getR(), x = trans[i].getX(), t = trans[i].getK();
			double b = t*(r*r+x*x);
			r/=b;
			x/=-b;
			b=t-1.0;
			double ri = r*b,xi = x*b, rj=-ri/t, xj=-xi/t;
			double dij = Ua[from] - Ua[to];
			double vi = Um[from], vj = Um[to];
			double vij = vi*vj;
			vi = vi*vi;
			vj = vj*vj;
			double cd = vij*Math.cos(dij*Math.PI/180);
			double sd = vij*Math.sin(dij*Math.PI/180);
			double pij = vi*(r+ri) - r*cd - x*sd;
			double pji = vj*(r+rj) - r*cd + x*sd;
			double dpb = pij+pji;
			ph = ph+dpb;
			double qij = -vi * (xj+x) + x*cd - r*sd;
			double qji = -vj * (xi+x) + x*cd + r*sd;
			double dqb = qij+qji;
			qh = qh+dqb;
			System.out.println(from +" "+ to +" "+ pij*100 +" "+ qij*100 +" "+ pji*100 +" "+ qji*100 +" "+ dpb*100 +" "+ dqb*100);
		}
		System.out.println("sum loss： "+ph*100 +" " + qh*100);
	}
	
}
