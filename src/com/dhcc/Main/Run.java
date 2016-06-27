package com.dhcc.Main;

import com.dhcc.Global.Variable;
import com.dhcc.powerflow.NewtonPowerFlow;
import com.dhcc.powerflow.ProcData;
import com.dhcc.util.IOUtil;

public class Run {
	
	private Variable variable = null;

	public void run() {
		variable = new Variable();
		IOUtil io = new IOUtil(variable);
		variable = io.getVariable();
		ProcData pd = new ProcData(variable);
		variable = pd.getVariable();
		io.readCDFDataWithOriIdx("/Users/xyk0058/Git/PowerFlow/src/com/dhcc/casedata/ieee14cdf.txt");
		pd.AdmtMatrix();
		pd.InitOri();
		NewtonPowerFlow pf = new NewtonPowerFlow(variable);
		pf.Run();
		variable = pf.getVariable();
		io.PrintInfo_iter(0);
		
		double[] Pi = variable.getP();
		double[] Qi = variable.getQ();
		System.out.println("\nP Q");
		for (int i=0;i<Pi.length;i++) {
			System.out.println("P"+i+": " + Pi[i] + "\tQ"+i+":" + Qi[i]);
		}
		//pd.CalBusPQFlow();
		System.out.println("\nPG QG PL QL");
		pd.CalBusFlow();
		System.out.println("\nloss");
		pd.BranchFlow();
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Run run = new Run();
		run.run();
	}

}
