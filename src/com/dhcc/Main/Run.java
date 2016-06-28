package com.dhcc.Main;

import java.util.ArrayList;
import java.util.List;

import com.dhcc.GA.GA;
import com.dhcc.GA.Result;
import com.dhcc.Global.Variable;
import com.dhcc.Model.Gene;
import com.dhcc.policy.OptFlow;
import com.dhcc.powerflow.ProcData;
import com.dhcc.util.IOUtil;

public class Run {
	
	private Variable variable = null;

	public void run() {
		
		List<Result> res = GA.run();
		ArrayList<Variable> list = new ArrayList<Variable>();
		
		variable = new Variable();
		IOUtil io = new IOUtil(variable);
		io.readCDFDataWithOriIdx("/Users/xyk0058/Git/PowerFlow/src/com/dhcc/casedata/ieee14cdf.txt");
		variable = io.getVariable();
		ProcData pd = new ProcData(variable);
		pd.AdmtMatrix();
		pd.InitOri();
		variable = pd.getVariable();
		list.add(variable);
		
		for (int i = 0; i < res.size(); ++i) {
			Result result = res.get(i);
			double[] x = result.getX();
			variable = new Variable();
			io = new IOUtil(variable);
			io.readCDFDataWithOriIdx("/Users/xyk0058/Git/PowerFlow/src/com/dhcc/casedata/ieee14cdf.txt");
			variable = io.getVariable();
			for (int j = 0; j < x.length; ++j) {
				Gene[] gene = variable.getGenerator();
				gene[j].setP(x[j]);
			}
			pd = new ProcData(variable);
			pd.AdmtMatrix();
			pd.InitOri();
			variable = pd.getVariable();
			list.add(variable);
		}
		
		OptFlow of = new OptFlow();
		of.run(list);
		
//		NewtonPowerFlow pf = new NewtonPowerFlow(variable);
//		pf.Run();
//		variable = pf.getVariable();
		//io.PrintInfo_iter(0);
		
//		double[] Pi = variable.getP();
//		double[] Qi = variable.getQ();
//		System.out.println("\nP Q");
//		for (int i=0;i<Pi.length;i++) {
//			System.out.println("P"+i+": " + Pi[i] + "\tQ"+i+":" + Qi[i]);
//		}
//		//pd.CalBusPQFlow();
//		System.out.println("\nPG QG PL QL");
//		pd.CalBusFlow();
//		System.out.println("\nloss");
//		pd.BranchFlow();
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Run run = new Run();
		run.run();
	}
	
	//spark-submit --class com.dhcc.Main.Run --master local[3] /Users/xyk0058/Desktop/powerflow.jar

}
