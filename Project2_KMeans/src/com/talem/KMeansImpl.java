package com.talem;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

public class KMeansImpl {
	
	CommonMethods objCommon;
	
	static HashMap<Integer, ArrayList<Double>> geneInfoMap;
	static HashMap<Integer, ArrayList<Double>> centroidMap;
	static HashMap<Integer, ArrayList<Integer>> groundTruthMap;
	static HashMap<Integer, Integer> geneIdGroundTruthMap;
	static int k;
	static int rotations;
	
	static HashMap<Integer, ArrayList<Integer>> clusters;
	static double[][] distMatrix;
	public KMeansImpl(HashMap<Integer, ArrayList<Double>> geneInfoMap,
			HashMap<Integer, ArrayList<Double>> centroidMap, HashMap<Integer, ArrayList<Integer>> groundTruthMap,
			HashMap<Integer, Integer> geneIdGroundTruthMap, int k, int rotations){
		objCommon = new CommonMethods();
		KMeansImpl.geneInfoMap = geneInfoMap;
		KMeansImpl.centroidMap = centroidMap;
		KMeansImpl.groundTruthMap = groundTruthMap;
		KMeansImpl.geneIdGroundTruthMap = geneIdGroundTruthMap;
		KMeansImpl.k = k;
		KMeansImpl.rotations = rotations;
	}
	
	public void KMEansImplementation() {

		distMatrix = new double[k+1][geneInfoMap.size()+1];
		for(int i=1;i<=k;i++)
		{
			for(int j=1;j<geneInfoMap.size();j++)
			{
				distMatrix[i][j]= objCommon.getEuclideanDistance(centroidMap.get(i), geneInfoMap.get(j));
			}
		}
		clusters=new HashMap<Integer, ArrayList<Integer>>();
		
		HashMap<Integer, ArrayList<Double>> mapAtI = new HashMap<Integer, ArrayList<Double>>();
		HashMap<Integer, ArrayList<Double>> mapAtJ = new HashMap<Integer, ArrayList<Double>>();
		ArrayList<Integer> mergedList1=new ArrayList<Integer>();
		ArrayList<Integer> mergedList2=new ArrayList<Integer>();
		for(Map.Entry<Integer, ArrayList<Double>> ent:mapAtI.entrySet())
		{
			mergedList1.add(ent.getKey());
		}
		for(Map.Entry<Integer, ArrayList<Double>> ent:mapAtJ.entrySet())
		{
			mergedList2.add(ent.getKey());
		}
		for(Map.Entry<Integer, ArrayList<Double>> ent:mapAtJ.entrySet())
		{
			mapAtI.put(ent.getKey(), ent.getValue());
			
		}
		
		for(int i=1;i<=geneInfoMap.size();i++)
		{
			double min = Double.MAX_VALUE;
			int minRow = 0;
			for(int j=1;j<=k; j++)
			{
				if(distMatrix[j][i]<min)
				{
					min = distMatrix[j][i];
					minRow = j;
				}
			}
			if(!clusters.containsKey(minRow))
			{
				ArrayList<Integer> list=new ArrayList<Integer>();
				list.add(i);
				clusters.put(minRow, list);
			}
			else
			{
				clusters.get(minRow).add(i);
			}
		}
		for(Integer ent:mergedList1)
		{
			//mergedList1.add(ent.getKey());
			initialClusters(mergedList1, ent, mergedList2);
		}
		for(Integer ent:mergedList1)
		{
			//mergedList2.add(ent.getKey());
			initialClusters(mergedList2, ent, mergedList1);
		}
		if(rotations == 100){
			while(!checkCentroidEquality(centroidMap, updateCentroid(centroidMap)))
			{
				centroidMap = updateCentroid(centroidMap);
				updatedDistanceMatrix(centroidMap);
			}
		}
		else {
			for(int i=0; i<rotations; i++){
				if(checkCentroidEquality(centroidMap, updateCentroid(centroidMap))){
					break ;
				}
				centroidMap = updateCentroid(centroidMap);
				updatedDistanceMatrix(centroidMap);
			}
		}
		
		/*
		while(!checkCentroidEquality(centroidMap, updateCentroid(centroidMap)))
		{
			centroidMap = updateCentroid(centroidMap);
			updatedDistanceMatrix(centroidMap);
		}
		*/
		String KMeansFilename = "";
		HashMap<Integer, Integer> tempClusterMap = writeToFile(KMeansFilename);
		
		TreeMap<Integer, Integer> gtTreeMap=new TreeMap<Integer, Integer>(geneIdGroundTruthMap);
		Set<Integer> gtKeys=gtTreeMap.keySet();
		ArrayList<Integer> gtKeysList=new ArrayList<Integer>();
		Iterator<Integer> gtIter=gtKeys.iterator();
		while(gtIter.hasNext())
			gtKeysList.add(gtIter.next());
		
		System.out.println("Jaccard value is: " + objCommon.calculateJaccardCoeff(gtKeysList, tempClusterMap, 
															geneIdGroundTruthMap));
		System.out.println("Silhoutte Index is: "+ objCommon.calculateSilhoutte(geneInfoMap, clusters));
		
		try {
			String PCAFilename = "KMeansPCA_new_dataset_1.txt";
			objCommon.generatePCAFile(clusters, geneInfoMap, PCAFilename);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
	}
	
	@SuppressWarnings({ "unchecked", "rawtypes" })
	private ArrayList initialClusters(ArrayList inputs, int clusters, ArrayList<Integer> clusterIndex) {
		ArrayList initialClusters = new ArrayList();
		Map<Integer, Object> initialGenes = new HashMap<Integer,Object>();


		int count = 0;
		for (int index : clusterIndex){
			if(!initialGenes.containsKey(index)){
				initialGenes.put(index, inputs.get(count));
			}
			count++;			
		}
	
		for(int index: initialGenes.keySet()){
			initialClusters.add(initialGenes.get(index));
		}
		//System.out.println("Initial Clusters are: "+initialClusters);
		
		return initialClusters;
	}
	
	private HashMap<Integer, Integer> writeToFile(String KMeansFilename) {
		HashMap<Integer, Integer> tempClusterMap = new HashMap<Integer, Integer>();
		
		for(Map.Entry<Integer, ArrayList<Integer>> ent:clusters.entrySet())
		{
			ArrayList<Integer> list=ent.getValue();
			
			for(int elem:list)
			{
				tempClusterMap.put(elem, ent.getKey());
			}
		}
		return tempClusterMap;
	}

	private void updatedDistanceMatrix(HashMap<Integer, ArrayList<Double>> centroidMap2) {
		for(int i=1;i<centroidMap2.size();i++)
		{
			for(int j=1;j<geneInfoMap.size();j++)
			{
				distMatrix[i][j] = objCommon.getEuclideanDistance(centroidMap2.get(i), geneInfoMap.get(j));
			}
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
					sum = sum + geneInfoMap.get(elem).get(i);
				}
				sum = sum/list.size();
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
	
}
