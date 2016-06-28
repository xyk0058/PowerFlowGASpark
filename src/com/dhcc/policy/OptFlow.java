package com.dhcc.policy;

import java.io.Serializable;
import java.util.ArrayList;

import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.function.VoidFunction;

import com.dhcc.Global.GV;
import com.dhcc.Global.Variable;
import com.dhcc.powerflow.NewtonPowerFlow;
import com.dhcc.powerflow.ProcData;
import com.dhcc.util.IOUtil;

public class OptFlow implements Serializable{
	
	
	public void run(ArrayList<Variable> list) {
		
		JavaRDD<Variable> vars= GV.getSC().parallelize(list);
		
		JavaRDD<Variable> result = vars.map(new Function<Variable, Variable>() {
			@Override
			public Variable call(Variable variable) throws Exception {
				// TODO Auto-generated method stub
				NewtonPowerFlow pf = new NewtonPowerFlow(variable);
				pf.Run();
				variable = pf.getVariable();
				return variable;
			}
		});
		
		result.foreach(new VoidFunction<Variable>() {
			@Override
			public void call(Variable variable) throws Exception {
				// TODO Auto-generated method stub
				System.out.println("----------------------------------------------------------------");
				
				ProcData pd = new ProcData(variable);

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
				
				System.out.println("----------------------------------------------------------------");
			}
		});
		//result.saveAsTextFile("/Users/xyk0058/Desktop/RESULTALL");
	}
}
