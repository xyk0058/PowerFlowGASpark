package com.dhcc.Global;

import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.sql.SQLContext;

public class GV {
	
	private static SparkConf conf = null;
	private static JavaSparkContext sc = null;
	private static SQLContext sqlContext = null;

	public static SparkConf getConf() {
		if (conf == null) {
			conf = new SparkConf().setAppName("PowerFlow").setMaster("local");
		}
		return conf;
	}

	public static JavaSparkContext getSC() {
		if (sc == null) {
			sc = new JavaSparkContext(getConf());
		}
		return sc;
	}

	public static SQLContext getSqlContext() {
		if (sqlContext == null) {
			sqlContext = new org.apache.spark.sql.SQLContext(getSC());
		}
		return sqlContext;
	}
	
}
