package com.dhcc.util;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;

import com.dhcc.Global.Variable;
import com.dhcc.Model.Branch;
import com.dhcc.Model.Gene;
import com.dhcc.Model.Info;
import com.dhcc.Model.Load;
import com.dhcc.Model.MPC;
import com.dhcc.Model.Tran;


public class IOUtil {
	
	private Variable variable = null;
	public IOUtil(Variable variable) {
		this.variable = variable;
	}
	
	public void UsingAdmtMatrix() {
		double G[][] = {{6.02503,-4.99913,0,0,-1.0259,0,0,0,0,0,0,0,0,0},
				{-4.99913,9.52132,-1.13502,-1.68603,-1.70114,0,0,0,0,0,0,0,0,0},
				{0,-1.13502,3.12099,-1.98598,0,0,0,0,0,0,0,0,0,0},
				{0,-1.68603,-1.98598,10.513,-6.84098,0,0,0,0,0,0,0,0,0},
				{-1.0259,-1.70114,0,-6.84098,9.56802,0,0,0,0,0,0,0,0,0},
				{0,0,0,0,0,6.57992,0,0,0,0,-1.95503,-1.52597,-3.09893,0},
				{0,0,0,0,0,0,0,0,0,0,0,0,0,0},
				{0,0,0,0,0,0,0,0,0,0,0,0,0,0},
				{0,0,0,0,0,0,0,0,5.32606,-3.90205,0,0,0,-1.42401},
				{0,0,0,0,0,0,0,0,-3.90205,5.78293,-1.88088,0,0,0},
				{0,0,0,0,0,-1.95503,0,0,0,-1.88088,3.83591,0,0,0},
				{0,0,0,0,0,-1.52597,0,0,0,0,0,4.01499,-2.48902,0},
				{0,0,0,0,0,-3.09893,0,0,0,0,0,-2.48902,6.72495,-1.13699},
				{0,0,0,0,0,0,0,0,-1.42401,0,0,0,-1.13699,2.561}};
		double B[][] = {{-19.4471,15.2631,0,0,4.23498,0,0,0,0,0,0,0,0,0},
				{15.2631,-30.2707,4.78186,5.11584,5.19393,0,0,0,0,0,0,0,0,0},
				{0,4.78186,-9.81148,5.06882,0,0,0,0,0,0,0,0,0,0},
				{0,5.11584,5.06882,-38.6352,21.5786,0,4.88951,0,1.8555,0,0,0,0,0},
				{4.23498,5.19393,0,21.5786,-35.5275,4.25745,0,0,0,0,0,0,0,0},
				{0,0,0,0,4.25745,-17.3407,0,0,0,0,4.09407,3.17596,6.10276,0},
				{0,0,0,4.88951,0,0,-19.549,5.67698,9.09008,0,0,0,0,0},
				{0,0,0,0,0,0,5.67698,-5.67698,0,0,0,0,0,0},
				{0,0,0,1.8555,0,0,9.09008,0,-24.0925,10.3654,0,0,0,3.02905},
				{0,0,0,0,0,0,0,0,10.3654,-14.7683,4.40294,0,0,0},
				{0,0,0,0,0,4.09407,0,0,0,4.40294,-8.49702,0,0,0},
				{0,0,0,0,0,3.17596,0,0,0,0,0,-5.42794,2.25197,0},
				{0,0,0,0,0,6.10276,0,0,0,0,0,2.25197,-10.6697,2.31496},
				{0,0,0,0,0,0,0,0,3.02905,0,0,0,2.31496,-5.34401}};

		variable.setB(B);
		variable.setG(G);
	}
	
	public void TestInfo() {
		Info info = new Info(4, 4,1, 2, 2, 0, 1, 3, 0.0001);
		info.setNpq(2);
		variable.setPf_info(info);
		Branch[] branch = new Branch[info.getNb()];
		Tran[] tran = new Tran[info.getNt()];
		Gene[] generator = new Gene[info.getNg()];
		Load[] load = new Load[info.getNl()];
		
		branch[0] = new Branch(0, 3, 0.173554, 0.330579, 0.017243);
		branch[1] = new Branch(1,1,0.0,-20.0,0.0);
		branch[2] = new Branch(2,0,.130165,0.247934,0.012932);
		branch[3] = new Branch(3, 2, 0.260331, 0.495868, 0.025864);

		tran[0] = new Tran(0,1,0.0,0.166667,1.128205);
		
		generator[0] = new Gene(3,variable.REF,0,0,1.05);
		generator[1] = new Gene(2,variable.PV,0.2,0,1.05);
		
		load[0] = new Load(1,variable.PQ,0.5,0.3,0.0);
		load[1] = new Load(3,variable.PQ,0.15,0.1,0.0);
		
		variable.setTrans(tran);
		variable.setBranch(branch);
		variable.setGenerator(generator);
		variable.setLoad(load);
	}
	
	public void readCDFData(String filename) {
		int n_bus = 0;
		int n_branch = 0;
		int nt = 0;
		int ng = 0;
		int nl = 0;
		int npv = 0;
		int npq = 0;
		int nb = 0;
		int nrl = 0;
		Branch[] branch;
		Tran[] tran;
		Gene[] generator;
		Load[] load;
		Load[] realLoad;
		
		try {
			InputStreamReader instrr = new InputStreamReader(new FileInputStream(filename));
			BufferedReader br = new BufferedReader(instrr);
			String row = null;
			String[] rowdata = null;
			row = br.readLine();n_bus = Integer.parseInt(row);
			generator = new Gene[n_bus];
			load = new Load[n_bus];
			realLoad = new Load[n_bus];
			
			int phIdx = 0;
			double php = 0, phq = 0, phv = 0, phtheta = 0, phpl = 0, phql = 0;
			
			int _npv=0;
			double[] max = new double[n_bus],min = new double[n_bus];
			for (int i = 0; i < n_bus; ++i) {
				row = br.readLine();
				rowdata = row.split(",");
				int type = Integer.parseInt(rowdata[6]);
				//System.out.println("type:" + type);
				if (type == 2) {		//generator
					int idx = Integer.parseInt(rowdata[0]);
					npv++;
					double p,q,v;
					v = Double.parseDouble(rowdata[7]);
					p = Double.parseDouble(rowdata[11]) - Double.parseDouble(rowdata[9]);
					q = Double.parseDouble(rowdata[12]) - Double.parseDouble(rowdata[10]);
					generator[ng++] = new Gene(idx,variable.PV,p/100.0,q/100.0,v);
					generator[ng-1].setPl(Double.parseDouble(rowdata[9]));
					generator[ng-1].setQl(Double.parseDouble(rowdata[10]));
					generator[ng-1].setG(Double.parseDouble(rowdata[17]));
					generator[ng-1].setB(Double.parseDouble(rowdata[18]));
					if (Double.parseDouble(rowdata[15]) == 0) generator[ng-1].setMax(2147483647);
					else generator[ng-1].setMax(Double.parseDouble(rowdata[15]));
					generator[ng-1].setMin(Double.parseDouble(rowdata[16]));
					_npv+=Double.parseDouble(rowdata[9]);
					max[ng-1] = 300;
					min[ng-1] = 0;
					//System.out.println("Gene:" + Double.parseDouble(rowdata[18]));
				} else if (type == 0 || type == 1) {
					int idx = Integer.parseInt(rowdata[0]);
					npq++;
					double p,q,v;
					v = Double.parseDouble(rowdata[7]);
					p = Double.parseDouble(rowdata[11]) - Double.parseDouble(rowdata[9]);
					q = Double.parseDouble(rowdata[12]) - Double.parseDouble(rowdata[10]);
					load[nl++] = new Load(idx,variable.PQ,p/100.0,q/100.0,v);
					if(type == 0) realLoad[nrl++] = new Load(idx,variable.PQ,p/100,q/100,v);
					load[nl-1].setPl(Double.parseDouble(rowdata[9]));
					load[nl-1].setQl(Double.parseDouble(rowdata[10]));
					load[nl-1].setG(Double.parseDouble(rowdata[17]));
					load[nl-1].setB(Double.parseDouble(rowdata[18]));
					if (Double.parseDouble(rowdata[15]) == 0) load[nl-1].setMax(2147483647);
					else load[nl-1].setMax(Double.parseDouble(rowdata[15]));
					load[nl-1].setMin(Double.parseDouble(rowdata[16]));
					//System.out.println("Load:" + Double.parseDouble(rowdata[18]));
				} else if (type == 3) {
					phIdx = Integer.parseInt(rowdata[0]);
					npv++;
					phv = Double.parseDouble(rowdata[7]);
					php = Double.parseDouble(rowdata[11]) - Double.parseDouble(rowdata[9]);
					phq = Double.parseDouble(rowdata[12]) - Double.parseDouble(rowdata[10]);
					phtheta = Double.parseDouble(rowdata[8]);
					phpl = Double.parseDouble(rowdata[9]);
					phql = Double.parseDouble(rowdata[10]);
				}
			}
			
			/*
			Opt.setNpq(_npv);
			double[] _max = new double[npv];
			System.arraycopy(max, 0, _max, 0, npv);
			Opt.setMax(_max);
			System.arraycopy(min, 0, _max, 0, npv);
			Opt.setMin(_max);
			*/
			
			generator[ng++] = new Gene(phIdx,variable.REF,php/100.0,phq/100.0,phv);
			generator[ng-1].setPl(phpl);
			generator[ng-1].setQl(phql);
			if (Double.parseDouble(rowdata[15]) == 0) generator[ng-1].setMax(2147483647);
			else generator[ng-1].setMax(Double.parseDouble(rowdata[15]));
			generator[ng-1].setMin(Double.parseDouble(rowdata[16]));
			
			
			int newIdx = 0;
			int[] newIndex = new int[n_bus + 1];
			for (int i = 0; i < nl; ++i) {
				newIndex[load[i].getI()] = newIdx;
				load[i].setI(newIdx);
				newIdx++;
			}
			for (int i = 0; i < ng - 1; ++i) {
				newIndex[generator[i].getI()] = newIdx;
				generator[i].setI(newIdx);
				newIdx++;
			}
			
			newIndex[generator[ng-1].getI()] = newIdx;
			generator[ng-1].setI(newIdx);

//			System.out.println("newidx");
//			for (int i = 1; i <= n_bus; ++i)
//				System.out.println(i + " " + newIndex[i]);
			
			row = br.readLine();n_branch = Integer.parseInt(row);
			tran = new Tran[n_branch];
			branch = new Branch[n_branch];

			for (int i = 0; i < n_branch; ++i) {
				row = br.readLine();
				rowdata = row.split(",");
				double k = Double.parseDouble(rowdata[14]);
				int from, to;
				double r, x, b;
				from = Integer.parseInt(rowdata[0]);
				from = newIndex[from];
				to = Integer.parseInt(rowdata[1]);
				to = newIndex[to];
				r = Double.parseDouble(rowdata[6]);
				x = Double.parseDouble(rowdata[7]);
				b = Double.parseDouble(rowdata[8]);
				if (k == 0) {
					branch[nb++] = new Branch(from,to,r,x,b);
				} else {
					tran[nt++] = new Tran(from,to,r,x,k);
				}
			}
			//System.out.println("nt" + nt);
			Info info = new Info(n_bus,nb,nt,ng,nl,1,npv,nrl,0.0001);
			variable.setPf_info(info);
			variable.setRefTheta(phtheta);
			variable.setTrans(Arrays.copyOf(tran, nt));
			variable.setBranch(Arrays.copyOf(branch, nb));
			variable.setGenerator(Arrays.copyOf(generator, ng));
			variable.setLoad(Arrays.copyOf(load, nl));
			variable.setRealLoad(Arrays.copyOf(realLoad, nrl));
			
			instrr.close();
		} catch (UnsupportedEncodingException e) {
			e.printStackTrace();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void readCDFDataWithOriIdx(String filename) {
		int n_bus = 0;
		int n_branch = 0;
		int nt = 0;
		int ng = 0;
		int nl = 0;
		int npv = 0;
		int npq = 0;
		int nb = 0;
		int nrl = 0;
		Branch[] branch;
		Tran[] tran;
		Gene[] generator;
		Load[] load;
		Load[] realLoad;
		
		try {
			InputStreamReader instrr = new InputStreamReader(new FileInputStream(filename));
			BufferedReader br = new BufferedReader(instrr);
			String row = null;
			String[] rowdata = null;
			row = br.readLine();n_bus = Integer.parseInt(row);
			generator = new Gene[n_bus];
			load = new Load[n_bus];
			realLoad = new Load[n_bus];
			
			int phIdx = 0;
			double php = 0, phq = 0, phv = 0, phg = 0, phb = 0, phpl = 0, phql = 0;
			
			int _npv=0;
			double[] max = new double[n_bus],min = new double[n_bus];
			for (int i = 0; i < n_bus; ++i) {
				row = br.readLine();
				rowdata = row.split(",");
				int type = Integer.parseInt(rowdata[6]);
				if (type == 2) {		//generator
					int idx = Integer.parseInt(rowdata[0]) - 1;
					npv++;
					double p,q,v;
					v = Double.parseDouble(rowdata[7]);
					p = Double.parseDouble(rowdata[11]) - Double.parseDouble(rowdata[9]);
					q = Double.parseDouble(rowdata[12]) - Double.parseDouble(rowdata[10]);
					p = -p;
					q = -q;
					generator[ng++] = new Gene(idx,variable.PV,p/100.0,q/100.0,v);
					generator[ng-1].setPl(Double.parseDouble(rowdata[9]));
					generator[ng-1].setQl(Double.parseDouble(rowdata[10]));
					generator[ng-1].setG(Double.parseDouble(rowdata[17]));
					generator[ng-1].setB(Double.parseDouble(rowdata[18]));
					if (Double.parseDouble(rowdata[15]) == 0) generator[ng-1].setMax(2147483647);
					else generator[ng-1].setMax(Double.parseDouble(rowdata[15]));
					generator[ng-1].setMin(Double.parseDouble(rowdata[16]));
					_npv+=Double.parseDouble(rowdata[9]);
					max[ng-1] = 300;
					min[ng-1] = 0;
//					System.out.println("Gene:" + idx);
				} else if (type == 0 || type == 1) {
					int idx = Integer.parseInt(rowdata[0]) - 1;
					npq++;
					double p,q,v;
					v = Double.parseDouble(rowdata[7]);
					p = Double.parseDouble(rowdata[11]) - Double.parseDouble(rowdata[9]);
					q = Double.parseDouble(rowdata[12]) - Double.parseDouble(rowdata[10]);
					p = -p;
					q = -q;
					load[nl++] = new Load(idx,variable.PQ,p/100.0,q/100.0,v);
					if(type == 0) realLoad[nrl++] = new Load(idx,variable.PQ,p/100,q/100,v);
					load[nl-1].setPl(Double.parseDouble(rowdata[11]));
					load[nl-1].setQl(Double.parseDouble(rowdata[12]));
					load[nl-1].setG(Double.parseDouble(rowdata[17]));
					load[nl-1].setB(Double.parseDouble(rowdata[18]));
					if (Double.parseDouble(rowdata[15]) == 0) load[nl-1].setMax(2147483647);
					else load[nl-1].setMax(Double.parseDouble(rowdata[15]));
					load[nl-1].setMin(Double.parseDouble(rowdata[16]));
				} else if (type == 3) {
					phIdx = Integer.parseInt(rowdata[0]) - 1;
					npv++;
					phv = Double.parseDouble(rowdata[7]);
					php = Double.parseDouble(rowdata[11]) - Double.parseDouble(rowdata[9]);
					phq = Double.parseDouble(rowdata[12]) - Double.parseDouble(rowdata[10]);
					php = -php;
					phq = -phq;
					phg = Double.parseDouble(rowdata[17]);
					phb = Double.parseDouble(rowdata[18]);
					phpl = Double.parseDouble(rowdata[11]);
					phql = Double.parseDouble(rowdata[12]);
;				}
			}
			
			/*
			Opt.setNpq(_npv);
			double[] _max = new double[npv];
			System.arraycopy(max, 0, _max, 0, npv);
			Opt.setMax(_max);
			System.arraycopy(min, 0, _max, 0, npv);
			Opt.setMin(_max);
			*/
			
			generator[ng] = new Gene(phIdx,variable.REF,php/100.0,phq/100.0,phv);
			generator[ng].setG(phg);
			generator[ng].setPl(phpl);
			generator[ng].setQl(phql);
			generator[ng++].setB(phb);
			if (Double.parseDouble(rowdata[15]) == 0) generator[ng-1].setMax(2147483647);
			else generator[ng-1].setMax(Double.parseDouble(rowdata[15]));
			generator[ng-1].setMin(Double.parseDouble(rowdata[16]));
			

			row = br.readLine();n_branch = Integer.parseInt(row);
			tran = new Tran[n_branch];
			branch = new Branch[n_branch];

			for (int i = 0; i < n_branch; ++i) {
				row = br.readLine();
				rowdata = row.split(",");
				double k = Double.parseDouble(rowdata[14]);
				int from, to;
				double r, x, b;
				from = Integer.parseInt(rowdata[0]) - 1;
				to = Integer.parseInt(rowdata[1]) - 1;
				r = Double.parseDouble(rowdata[6]);
				x = Double.parseDouble(rowdata[7]);
				b = Double.parseDouble(rowdata[8]);
				if (k == 0) {
					branch[nb++] = new Branch(from,to,r,x,b);
				} else {
					tran[nt++] = new Tran(from,to,r,x,k);
				}
			}
			
			Info info = new Info(n_bus,nb,nt,ng,nl,1,npv,nrl,0.0001);
			variable.setPf_info(info);
			variable.setTrans(Arrays.copyOf(tran, nt));
			variable.setBranch(Arrays.copyOf(branch, nb));
			variable.setGenerator(Arrays.copyOf(generator, ng));
			variable.setLoad(Arrays.copyOf(load, nl));
			variable.setRealLoad(Arrays.copyOf(realLoad, nrl));
			
			instrr.close();
		} catch (UnsupportedEncodingException e) {
			e.printStackTrace();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void ReadCase14(String filename) {
		ReadData(filename);
		InitData();
	}
	
	public void PrintInfo_b() {
		Info info = variable.getPf_info();
//		double bus[][] = _mpc.getBus();
//		double bra[][] = _mpc.getBranch();
		Branch[] branch = variable.getBranch();
		Tran[] tran = variable.getTrans();
		Gene[] gen = variable.getGenerator();
		Load[] load = variable.getLoad();	
		
//		System.out.println("Bus " + bus.length + " " + bus[0].length);
//		for (int i=0; i<info.getN(); ++i) {
//			for (int j=0; j<bus[i].length; ++j)
//				System.out.print(bus[i][j] + " ");
//			System.out.print("\n");
//		}
//		System.out.println("Bra " + bra.length + " " + bra[0].length);
//		for (int i=0; i<bra.length; ++i) {
//			for (int j=0; j<bra[i].length; ++j)
//				System.out.print(bra[i][j] + " ");
//			System.out.print("\n");
//		}
		
		System.out.println("Branch " + branch.length);
		for (int i=0; i<info.getNb(); ++i) {
			//System.out.println(branch[i].getFrom() + " " + branch[i].getTo() + " "
			//		+ branch[i].getR() + " "+ branch[i].getX() + " "+ branch[i].getY0() + " ");
			System.out.println(branch[i].getFrom() + " " + branch[i].getTo() + " "
							+ branch[i].getR() + " "+ branch[i].getX() + " "+  "1 ");
		}
		//System.out.println("Tran " + tran.length );
		for (int i=0; i<info.getNt(); ++i) {
			System.out.println(tran[i].getFrom() + " " + tran[i].getTo() + " "
					+ tran[i].getR() + " "+ tran[i].getX() + " "+ tran[i].getK() + " ");
		}
		System.out.println("Gen " + gen.length);
		for (int i=0; i<gen.length; ++i) {
			System.out.println(gen[i].getI() + " " + gen[i].getJ() 
						+ " " + gen[i].getP() + " " + gen[i].getQ() + " "
						+ gen[i].getV());
		}
		System.out.println("Load " + load.length);
		for (int i=0; i<load.length; ++i) {
			System.out.println(load[i].getI() + " " + load[i].getJ() + " " + load[i].getP() 
					+ " " + load[i].getQ()+ " "+ load[i].getV());
		}
	}

	public void PrintInfo(){
		Info info = variable.getPf_info();
		double B[][] = variable.getB();
		double Bp[][] = variable.getBp();
		double Bpp[][] = variable.getBpp();
		double G[][] = variable.getG();
		
		System.out.println("G " + G.length + " " + G[0].length);
		for (int i=0; i<info.getN(); ++i) {
			for (int j=0; j<G[i].length; ++j)
				System.out.print(G[i][j] + " ");
			System.out.print("\n");
		}
		System.out.println("B " + B.length + " " + B[0].length);
		for (int i=0; i<info.getN(); ++i) {
			for (int j=0; j<B[i].length; ++j)
				System.out.print(B[i][j] + " ");
			System.out.print("\n");
		}
//		System.out.println("BP " + Bp.length + " " + Bp[0].length);
//		for (int i=0; i<info.getN()-1; ++i) {
//			for (int j=0; j<Bp[i].length; ++j)
//				System.out.print(Bp[i][j] + " ");
//			System.out.print("\n");
//		}
//		System.out.println("BPP " + Bpp.length + " " + Bpp[0].length);
//		for (int i=0; i<info.getNpq(); ++i) {
//			for (int j=0; j<Bpp[i].length; ++j)
//				System.out.print(Bpp[i][j] + " ");
//			System.out.print("\n");
//		}
	}
	
	public void PrintInfo_iter(int k) {
		
		double Um[] = variable.getOriU();
		double Ua[] = variable.getOriTheta();
		double Ptemp[] = variable.getPtemp();
		double Qtemp[] = variable.getQtemp();
		
		System.out.println(k + "'s iterator.");
		
		System.out.println("ID UM UA");
		for (int i=0; i<Um.length; ++i) {
			System.out.println(i+" "+Um[i] + " "+Ua[i]);
		}
//		System.out.println("Um " + Um.length );
//		for (int i=0; i<Um.length; ++i) 
//			System.out.print(Um[i] + " ");
//		System.out.println();
//		System.out.println("Ua " + Ua.length );
//		for (int i=0; i<Ua.length; ++i) 
//			//System.out.print((Ua[i]*180/Math.PI) + " ");
//			System.out.print((Ua[i]) + " ");
//		System.out.println();
		
//		System.out.println("Ptemp " + Ptemp.length );
//		for (int i=0; i<Ptemp.length; ++i) 
//			System.out.print(Ptemp[i] + " ");
//		System.out.println();
//		System.out.println("Qtemp " + Qtemp.length );
//		for (int i=0; i<Qtemp.length; ++i) 
//			System.out.print(Qtemp[i] + " ");
//		System.out.println();
	}
	private MPC _mpc;
	
	public void ReadData(String filename) {
		
		try {
			InputStreamReader instrr = new InputStreamReader(new FileInputStream(filename));
			BufferedReader br = new BufferedReader(instrr);
			String row = null;
			String[] rowdata = null;
			
			row = br.readLine();int nbus = Integer.parseInt(row);
			row = br.readLine();int ngen = Integer.parseInt(row);
			row = br.readLine();int nbranch = Integer.parseInt(row);
			
			int N = nbus;
			int Nb = nbranch;
			int Ng = ngen;
			int Nl = nbus - ngen;
			int Npv = 0;
			int V0 = 0;
			
			_mpc = new MPC(nbus, ngen, nbranch);

			double[][] bus = _mpc.getBus();
			//System.out.println("lanlan");
			for (int i=0; i<nbus; ++i) {
				row = br.readLine();
				rowdata = row.split(",");
				for (int j=0; j<rowdata.length; ++j) {
					bus[i][j] = Double.parseDouble(rowdata[j]);
				}
				V0 += Math.abs(bus[i][7]);
				bus[i][0] = bus[i][0] - 1;
				bus[i][3] = bus[i][3] / 100;
				bus[i][2] = bus[i][2] / 100;
				bus[i][8] = bus[i][8] * Math.PI / 180;
				if(bus[i][1] == variable.PV) ++Npv;
			}
			V0 = V0 / nbus;
			
			double[][] gen = _mpc.getGen();
			for (int i=0; i<ngen; ++i) {
				row = br.readLine();
				rowdata = row.split(",");
				for (int j=0; j<rowdata.length; ++j) {
					gen[i][j] = Double.parseDouble(rowdata[j]);
				}
				gen[i][0] = gen[i][0] - 1;
				gen[i][1] = gen[i][1] / 100;
				gen[i][2] = gen[i][2] / 100;
			}
			int Nt=0;
			double[][] branch = _mpc.getBranch();
			for (int i=0; i<nbranch; ++i) {
				row = br.readLine();
				rowdata = row.split(",");
				//System.out.println(row);
				for (int j=0; j<rowdata.length; ++j) {
					branch[i][j] = Double.parseDouble(rowdata[j]);
				}
				branch[i][0] = branch[i][0] - 1;
				branch[i][1] = branch[i][1] - 1;
				if (branch[i][8]!=0) 
					++Nt;
			}
			
			//System.out.print(Nt);
			
			//index
			int index[] = new int[N];
			double busc[][] = new double[bus.length][bus[0].length];
			double branchc[][] = new double[branch.length][branch[0].length];
			int _ref=N-1,  _pv=N-Npv-1, _pq = 0;
			
			for (int i=0; i<Nb; ++i) 
				for (int j=0; j<branch[i].length; ++j)
					branchc[i][j]=branch[i][j];
			
			for (int i=0; i<N; ++i) {
				int id = (int)bus[i][0];
				if ((int)bus[i][1] == variable.REF) {
					index[(int)bus[i][0]] = _ref;
					busc[_ref] = bus[i];
					busc[_ref][0] = _ref;
					for (int j=0; j<Nb; ++j) {
						if ((int)branch[j][0] == id) 
							branchc[j][0] = _ref;
						if ((int)branch[j][1] == id)
							branchc[j][1] = _ref;
					}
				}else if ((int)bus[i][1] == variable.PQ) {
					index[(int)bus[i][0]] = _pq;
					busc[_pq] = bus[i];
					busc[_pq][0] = _pq;
					for (int j=0; j<Nb; ++j) {
						if ((int)branch[j][0] == id) 
							branchc[j][0] = _pq;
						if ((int)branch[j][1] == id)
							branchc[j][1] = _pq;
					}
					++_pq;
				}else if ((int)bus[i][1] == variable.PV) {
					index[(int)bus[i][0]] = _pv;
					busc[_pv] = bus[i];
					busc[_pv][0] = _pv;
					for (int j=0; j<Nb; ++j) {
						if ((int)branch[j][0] == id) 
							branchc[j][0] = _pv;
						if ((int)branch[j][1] == id)
							branchc[j][1] = _pv;
					}
					++_pv;
				}
			}
			for (int i=0; i<gen.length; ++i)
				gen[i][0] = index[(int)gen[i][0]];
			
			_mpc.setBus(busc);
			_mpc.setBranch(branchc);
			_mpc.setGen(gen);
			Nb = nbranch - Nt;
			Info pf_info = new Info(N, Nb,Nt, Ng, Nl, V0, Npv, Nl, 0.0001);
			variable.setPf_info(pf_info);
			
			
			br.close();
			instrr.close();
		} catch (UnsupportedEncodingException e) {
			e.printStackTrace();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return;
	}
	
	public void InitData() {
		Info info = variable.getPf_info();
		Branch[] branch = new Branch[info.getNb()];
		Tran[] tran = new Tran[info.getNt()];
		Gene[] generator = new Gene[info.getNg()];
		Load[] load = new Load[info.getNl()];
		
		ArrayList<Integer> genIdx = new ArrayList<Integer>();
		
		double[][] brch = _mpc.getBranch();
		double[][] gen = _mpc.getGen();
		double[][] bus = _mpc.getBus();
		
		//֧·
		int _nb = 0, _nt = 0;
		for (int i = 0; i < info.getNb() + info.getNt(); ++i) {
			if (brch[i][8] == 0) {
				branch[_nb] = new Branch();
				branch[_nb].setFrom((int) brch[i][0]);
				branch[_nb].setTo((int) brch[i][1]);
				branch[_nb].setR(brch[i][2]);
				branch[_nb].setX(brch[i][3]);
				branch[_nb++].setY0(brch[i][4]);
			} else {
				tran[_nt] = new Tran();
				tran[_nt].setFrom((int) brch[i][0]);
				tran[_nt].setTo((int) brch[i][1]);
				tran[_nt].setR(brch[i][2]);
				tran[_nt].setX(brch[i][3]);
				tran[_nt++].setK(brch[i][8]);
			}
		}
		
		for (int i = 0; i < info.getNg(); ++i) {
			generator[i] = new Gene();
			generator[i].setI((int) gen[i][0]);
			generator[i].setP(gen[i][1]);
			generator[i].setQ(gen[i][2]);
			generator[i].setV(gen[i][5]);
			//TODO
			if ((int)bus[(int) gen[i][0]][1] == variable.REF) {
				generator[i].setJ(variable.REF);
			}else if ((int)bus[(int) gen[i][0]][1] == variable.PV) {
				generator[i].setJ(variable.PV);
			}else if ((int)bus[(int) gen[i][0]][1] == variable.PQ) {
				generator[i].setJ(variable.PQ);
			}
			genIdx.add((int) gen[i][0]);
		}

		int j = 0;
		for (int i = 0; i < info.getN(); ++i) {
			if (genIdx.contains((int) bus[i][0])) continue;
			//if ((int)bus[i][1] == variable.REF) continue;
			//System.out.println(load.length + " " + j + " " + bus[i][0]);
			load[j] = new Load();
			load[j].setI((int) bus[i][0]);
			if ((int)bus[i][1] == variable.REF) {
				load[i].setJ(variable.REF);
			}else if ((int)bus[i][1] == variable.PV) {
				load[i].setJ(variable.PV);
			}else if ((int)bus[i][1] == variable.PQ) {
				load[i].setJ(variable.PQ);
			}
			load[j].setP(bus[i][2]);
			load[j].setQ(bus[i][3]);
			load[j].setV(bus[i][8]);
			j++;
		}
		variable.setTrans(tran);
		variable.setBranch(branch);
		variable.setGenerator(generator);
		variable.setLoad(load);
	}

	public Variable getVariable() {
		return variable;
	}

}
