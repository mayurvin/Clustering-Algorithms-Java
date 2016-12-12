package com.talem;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;
import java.util.Queue;
import java.util.Random;
import java.util.Scanner;
import java.util.Set;
import java.util.TreeMap;

public class DBScanImpl {
	CommonMethods objCommon;

	static HashMap<Integer, Integer> geneIdGroundTruthMap;
	static HashMap<Integer, ArrayList<Double>> geneInfoMap;
	
	static HashMap<Integer, ArrayList<Integer>> dbScanClusters;
	static ArrayList<HashSet<Integer>> allClusters = new ArrayList<HashSet<Integer>>();
	static ArrayList<Integer> corePoints = new ArrayList<Integer>();
	
	HashMap<Integer, ArrayList<Double>> centroidMap = new HashMap<Integer, ArrayList<Double>>();
	static HashMap<Integer, ArrayList<Integer>> clusters = new HashMap<Integer, ArrayList<Integer>>();
	
	static double[][] distMatrix;
	
	public DBScanImpl(HashMap<Integer, ArrayList<Double>> geneInfoMap,
			HashMap<Integer, Integer> geneIdGroundTruthMap){
		objCommon = new CommonMethods();
		DBScanImpl.dbScanClusters = new HashMap<Integer, ArrayList<Integer>>();
		DBScanImpl.geneInfoMap = geneInfoMap;
		DBScanImpl.geneIdGroundTruthMap = geneIdGroundTruthMap;
	}
	
	@SuppressWarnings("rawtypes")
	public void DBScanImplementation(int sampleSize) {
		
		Scanner scan = new Scanner(System.in);
		
		System.out.println("sampleSize:" + sampleSize);
		int randomPoint = 1;
		System.out.println("Enter eps:");
		double eps = scan.nextDouble();
		System.out.println("Enter minpts:");
		int minpts = scan.nextInt();
		scan.close();
		boolean isProcessed[] = new boolean[sampleSize + 1];
		for (int i = 1; i <= sampleSize; i++){
			isProcessed[i] = false;
		}
		
		distMatrix = new double[geneInfoMap.size() + 1][geneInfoMap.size() + 1];
		for (int i = 1; i <= geneInfoMap.size(); i++) {
			for (int j = 1; j <= geneInfoMap.size(); j++) {
				distMatrix[i][j] = objCommon.getEuclideanDistance(geneInfoMap.get(i),
						geneInfoMap.get(j));
			}
		}
		
		for (int i = 1; i <= geneInfoMap.size(); i++) {
			calculateEps(new ArrayList(), new Integer(geneInfoMap.size()));
		}
		Random rand = new Random();
		int point = rand.nextInt(sampleSize);
		while (!checkProcessed(sampleSize, isProcessed)) {
			if (point > 0) {
				HashSet<Integer> clusters = new HashSet<Integer>();
				if (!isProcessed[randomPoint]) {

					ArrayList<Integer> neighbors = new ArrayList<Integer>();

					neighbors = findNeighbors(randomPoint, eps, neighbors);
					if (neighbors.size() > minpts) {
						Queue<Integer> queue = new LinkedList<Integer>();
						for (int elem : neighbors) {
							clusters.add(elem);
							queue.add(elem);
						}
						while (!queue.isEmpty()) {
							int individualNeighbor = queue.poll();
							if (isProcessed[individualNeighbor] == false) {
								neighbors = new ArrayList<Integer>();
								neighbors = findNeighbors(individualNeighbor, eps, neighbors);
								if (neighbors.size() >= minpts) {
									for (int elem : neighbors) {
										clusters.add(elem);
										queue.add(elem);
									}
								}

								corePoints.add(individualNeighbor);
							}
							isProcessed[individualNeighbor] = true;

						}

						corePoints.add(randomPoint);
						isProcessed[randomPoint] = true;
					} else {
						isProcessed[randomPoint] = true;
					}

				}
				if (clusters.size() > 0)
					allClusters.add(clusters);
				isProcessed[randomPoint] = true;
				randomPoint++;
			}
		}
		
		Iterator<HashSet<Integer>> allClusterIterator = allClusters.iterator();

		while(!checkCentroidEquality(centroidMap, updateCentroid(centroidMap)))
		{
			centroidMap = updateCentroid(centroidMap);
			//updatedDistanceMatrix(centroidMap);
		}
		
		int i = 1;
		while (allClusterIterator.hasNext()) {
			HashSet<Integer> cluster = allClusterIterator.next();
			Iterator<Integer> iter = cluster.iterator();
			ArrayList<Integer> clusterList = new ArrayList<Integer>();
			while (iter.hasNext()) {
				clusterList.add(iter.next());
			}
			dbScanClusters.put(i, clusterList);
			i++;
			//System.out.println("cluster size is:" + cluster.size());
		}
		
		String dbScanFileName = "Gene-ClusterId-Map-DBScan.txt";
		HashMap<Integer, Integer> tempClusterMap = writeToFile(dbScanFileName);
		
		//System.out.println("clusters size:" + dbScanClusters.size());
		TreeMap<Integer, Integer> gtTreeMap = new TreeMap<Integer, Integer>(geneIdGroundTruthMap);
		Set<Integer> gtKeys = gtTreeMap.keySet();
		ArrayList<Integer> gtKeysList = new ArrayList<Integer>();
		Iterator<Integer> gtIter = gtKeys.iterator();
		while (gtIter.hasNext()){
			gtKeysList.add(gtIter.next());
		}

		System.out.println("Jaccard coefficient is: " + objCommon.calculateJaccardCoeff(gtKeysList, 
				tempClusterMap, geneIdGroundTruthMap));
		System.out.println("Silhoutte index is: " + objCommon.calculateSilhoutte(geneInfoMap, dbScanClusters));

		try {
			String PCAFilename = "DBScan_new_dataset_1.txt";
			objCommon.generatePCAFile(dbScanClusters, geneInfoMap, PCAFilename);
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
	
	
	public static boolean checkCentroidEquality(HashMap<Integer, ArrayList<Double>> centroidMap1,HashMap<Integer, ArrayList<Double>> centroidMap2)
	{
		if(centroidMap1.size()!=centroidMap2.size())
		{
			return false;
		}
		else
		{
			for(int i=1;i<centroidMap1.size();i++)
			{
				ArrayList<Double> list1=centroidMap1.get(i);
				ArrayList<Double> list2=centroidMap2.get(i);
				for(int j=0;j<list1.size();j++)
				{
					if(Math.abs(list1.get(j)-list2.get(j))>0.001) //scope for pruning
						return false;
				}
			}
			
		}
		return true;
	}
	
	public static HashMap<Integer, ArrayList<Double>> updateCentroid(HashMap<Integer, ArrayList<Double>> centroids)	{
		
		HashMap<Integer, ArrayList<Double>> localCentroid=new HashMap<Integer, ArrayList<Double>>();
		for(Map.Entry<Integer, ArrayList<Integer>> ent : clusters.entrySet())
		{
			ArrayList<Integer> list=ent.getValue();
			double sum=0;
			for(int i=0;i< geneInfoMap.get(list.get(0)).size();i++)
			{
				
				sum=0;
				for(int elem:list)
				{
					sum = sum+ geneInfoMap.get(elem).get(i);
				}
				sum=sum/list.size();
				if(!localCentroid.containsKey(ent.getKey()))
				{
					ArrayList<Double> newList=new ArrayList<Double>();
					newList.add(i,sum);
					localCentroid.put(ent.getKey(), newList);
				}
				else
					localCentroid.get(ent.getKey()).add(i,sum);
			}
		}
		
		return localCentroid;
	}

	public static boolean checkProcessed(int sampleSize, boolean[] isProcessed) {

		for (int i = 1; i < sampleSize; i++) {
			if (!isProcessed[i])
				return false;
		}
		return true;
	}
	public static ArrayList<Integer> findNeighbors(int geneId, double eps, ArrayList<Integer> neighbors) {

		for (int i = 1; i <= geneInfoMap.size(); i++) {
			if (distMatrix[geneId][i] <= eps) {
				neighbors.add(i);
			}
		}
		return neighbors;
	}
	@SuppressWarnings("rawtypes")
	public double calculateEps(ArrayList geneSet, int minPts) {
		double eps = 0;
		
		for(int i = 0; i < geneSet.size(); i++) {
			double dist[] = new double[geneSet.size() - 1];
			int count = 0;
			for(int j = 0; j < geneSet.size(); j++) {
				if(j != i) {
					dist[count++] = ((DBScanImpl) geneSet.get(i)).eucDist(geneSet.get(j));
				}
			}
			Arrays.sort(dist);
			eps+= dist[minPts - 1];
		}
		//System.out.println(eps);
		return eps/geneSet.size();
	}
	private HashMap<Integer, Integer> writeToFile(String dbScanFileName) {
		
		HashMap<Integer, Integer> tempClusterMap = new HashMap<Integer, Integer>();
		for (Map.Entry<Integer, ArrayList<Integer>> ent : dbScanClusters.entrySet()) {
			
			ArrayList<Integer> list = ent.getValue();
			for (int elem : list) {
				tempClusterMap.put(elem, ent.getKey());
				//fwr.write(elem + "," + ent.getKey());
				//.write("\n");
			}
		}
		//fwr.close();
		return tempClusterMap;
	}
	@SuppressWarnings("unchecked")
	public double eucDist(Object gene){
		double eucDist = 0;

		for(int index = 0; index < ((HashMap<Integer, Integer>) gene).size(); index++){
			double entry = ((HashMap<Integer, Integer>) gene).get(index) - ((HashMap<Integer, Integer>) gene).get(index);
			eucDist+= (entry*entry);
		}		

		return Math.sqrt(eucDist); 
	}
}
