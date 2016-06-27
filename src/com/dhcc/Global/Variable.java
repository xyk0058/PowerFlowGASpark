package com.dhcc.Global;

import com.dhcc.Model.Branch;
import com.dhcc.Model.Gene;
import com.dhcc.Model.Info;
import com.dhcc.Model.Load;
import com.dhcc.Model.Tran;

public class Variable {
	
	public int PQ = 1;
	public int PV = 2;
	public int REF = 3;
	private Info pf_info;
	private Branch[] branch;		//Nb
	private Gene[] generator;		//Ng
	private Load[] load;			//Nl
	private Tran[] trans;
	private double[][] G;
	private double[][] B;
	private double[][] Bp;
	private double[][] Bpp;
	private double[] oriU;
	private double[] oriTheta;
	private double[] P;
	private double[] Q;
	private double[] Ptemp;
	private double[] Qtemp;
	private double[][] invBp;
	private double[][] invBpp;
	private Load[] realLoad;
	private double[][] jacob;
	private double refTheta = 0.0;
	public int getPQ() {
		return PQ;
	}
	public void setPQ(int pQ) {
		PQ = pQ;
	}
	public int getPV() {
		return PV;
	}
	public void setPV(int pV) {
		PV = pV;
	}
	public int getREF() {
		return REF;
	}
	public void setREF(int rEF) {
		REF = rEF;
	}
	public Info getPf_info() {
		return pf_info;
	}
	public void setPf_info(Info pf_info) {
		this.pf_info = pf_info;
	}
	public Branch[] getBranch() {
		return branch;
	}
	public void setBranch(Branch[] branch) {
		this.branch = branch;
	}
	public Gene[] getGenerator() {
		return generator;
	}
	public void setGenerator(Gene[] generator) {
		this.generator = generator;
	}
	public Load[] getLoad() {
		return load;
	}
	public void setLoad(Load[] load) {
		this.load = load;
	}
	public Tran[] getTrans() {
		return trans;
	}
	public void setTrans(Tran[] trans) {
		this.trans = trans;
	}
	public double[][] getG() {
		return G;
	}
	public void setG(double[][] g) {
		G = g;
	}
	public double[][] getB() {
		return B;
	}
	public void setB(double[][] b) {
		B = b;
	}
	public double[][] getBp() {
		return Bp;
	}
	public void setBp(double[][] bp) {
		Bp = bp;
	}
	public double[][] getBpp() {
		return Bpp;
	}
	public void setBpp(double[][] bpp) {
		Bpp = bpp;
	}
	public double[] getOriU() {
		return oriU;
	}
	public void setOriU(double[] oriU) {
		this.oriU = oriU;
	}
	public double[] getOriTheta() {
		return oriTheta;
	}
	public void setOriTheta(double[] oriTheta) {
		this.oriTheta = oriTheta;
	}
	public double[] getP() {
		return P;
	}
	public void setP(double[] p) {
		P = p;
	}
	public double[] getQ() {
		return Q;
	}
	public void setQ(double[] q) {
		Q = q;
	}
	public double[] getPtemp() {
		return Ptemp;
	}
	public void setPtemp(double[] ptemp) {
		Ptemp = ptemp;
	}
	public double[] getQtemp() {
		return Qtemp;
	}
	public void setQtemp(double[] qtemp) {
		Qtemp = qtemp;
	}
	public double[][] getInvBp() {
		return invBp;
	}
	public void setInvBp(double[][] invBp) {
		this.invBp = invBp;
	}
	public double[][] getInvBpp() {
		return invBpp;
	}
	public void setInvBpp(double[][] invBpp) {
		this.invBpp = invBpp;
	}
	public Load[] getRealLoad() {
		return realLoad;
	}
	public void setRealLoad(Load[] realLoad) {
		this.realLoad = realLoad;
	}
	public double[][] getJacob() {
		return jacob;
	}
	public void setJacob(double[][] jacob) {
		this.jacob = jacob;
	}
	public double getRefTheta() {
		return refTheta;
	}
	public void setRefTheta(double refTheta) {
		this.refTheta = refTheta;
	}
	
}
