package com.dhcc.Model;

import java.io.Serializable;

public class Load implements Serializable {
	private int i;
	private double p;
	private double q;
	private int j;//type
	private double v;
	private double g;
	private double b;
	private double pl;
	private double ql;
	private double max;
	private double min;
	
	
	public double getMax() {
		return max;
	}

	public void setMax(double max) {
		this.max = max;
	}

	public double getMin() {
		return min;
	}

	public void setMin(double min) {
		this.min = min;
	}
	public double getPl() {
		return pl;
	}

	public void setPl(double pl) {
		this.pl = pl;
	}

	public double getQl() {
		return ql;
	}

	public void setQl(double ql) {
		this.ql = ql;
	}

	public double getG() {
		return g;
	}

	public void setG(double g) {
		this.g = g;
	}
	public double getB() {
		return b;
	}

	public void setB(double b) {
		this.b = b;
	}
	public Load(){}
	public Load(int i, int j,double p,double q,double v) {
		this.i=i;
		this.j=j;
		this.p=p;
		this.q=q;
		this.v=v;
		this.g=0;
		this.b=0;
	}

	public int getJ() {
		return j;
	}
	public void setJ(int j) {
		this.j = j;
	}
	public double getV() {
		return v;
	}
	public void setV(double v) {
		this.v = v;
	}
	public int getI() {
		return i;
	}
	public void setI(int i) {
		this.i = i;
	}
	public double getP() {
		return p;
	}
	public void setP(double p) {
		this.p = p;
	}
	public double getQ() {
		return q;
	}
	public void setQ(double q) {
		this.q = q;
	}
	
}
