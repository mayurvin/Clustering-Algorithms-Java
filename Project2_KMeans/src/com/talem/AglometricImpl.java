package com.talem;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

public class AglometricImpl {
	CommonMethods objCommon;
	
	static ArrayList<HashMap<Integer, ArrayList<Double>>> listOfMaps;
	static HashMap<Integer, ArrayList<Double>> geneInfoMap;
	static HashMap<Integer, Integer> geneIdGroundTruthMap;
	static HashMap<Integer, ArrayList<Integer>> aglmClusters;
	
	static int k;
	
	public AglometricImpl(ArrayList<HashMap<Integer, ArrayList<Double>>> listOfMaps, HashMap<Integer, ArrayList<Double>> geneInfoMap,
			HashMap<Integer, Integer> geneIdGroundTruthMap, HashMap<Integer, ArrayList<Integer>> aglmClusters, int k){
		objCommon = new CommonMethods();
		AglometricImpl.listOfMaps = listOfMaps;
		AglometricImpl.geneInfoMap = geneInfoMap;
		AglometricImpl.geneIdGroundTruthMap = geneIdGroundTruthMap;
		AglometricImpl.aglmClusters = aglmClusters;
		AglometricImpl.k = k;
	}
	
	
	public void AglometricImplementation(){
		int i=0,j=0;
		int minI=-1,minJ=-1;
		double minDistance=Double.MAX_VALUE;
		
		while(listOfMaps.size()>k)
		{
			i=0;j=0;
			minI=-1;minJ=-1;
			minDistance=Double.MAX_VALUE;
			for(i=0;i<listOfMaps.size()-1;i++)
			{
				for(j=i+1;j<listOfMaps.size();j++)
				{
					Set<Integer> p1Keys=listOfMaps.get(i).keySet();
					Set<Integer> p2Keys=listOfMaps.get(j).keySet();
					ArrayList<Integer> p1KeysList=new ArrayList<Integer>();
					ArrayList<Integer> p2KeysList=new ArrayList<Integer>();
					Iterator<Integer> iter1=p1Keys.iterator();
					Iterator<Integer> iter2=p2Keys.iterator();
					
					HashMap<Integer, ArrayList<Double>> mapAtI = listOfMaps.get(i);
					HashMap<Integer, ArrayList<Double>> mapAtJ = listOfMaps.get(j);
					
					while(iter1.hasNext())
						p1KeysList.add(iter1.next());
					while(iter2.hasNext())
						p2KeysList.add(iter2.next());
					
					for(int p=0;p<p1KeysList.size();p++)
					{
						for(int q=0;q<p2KeysList.size();q++)
						{
							if(objCommon.getEuclideanDistance(mapAtI.get(p1KeysList.get(p)), mapAtJ.get(p2KeysList.get(q)))<minDistance)
							{
								minDistance= objCommon.getEuclideanDistance(mapAtI.get(p1KeysList.get(p)), mapAtJ.get(p2KeysList.get(q)));
								minI=i;
								minJ=j;
							}
						}
					}
				}
			}
	
			
			HashMap<Integer, ArrayList<Double>> mapAtI=listOfMaps.get(minI);
			HashMap<Integer, ArrayList<Double>> mapAtJ=listOfMaps.get(minJ);
			ArrayList<ArrayList<Integer>> listOfLists=new ArrayList<ArrayList<Integer>>();
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
			listOfLists.add(mergedList1);
			listOfLists.add(mergedList2);
			//mergedGenes.put(minI+1, listOfLists);
			listOfMaps.remove(minJ);
		}
		
		String agglometricFileName = "Gene-ClusterId-Map-Aglm.txt";
		HashMap<Integer, Integer> tempClusterMap = writeToFile(agglometricFileName);
		//HashMap<Integer, Integer> tempClusterMap = new HashMap<Integer, Integer>();
		int clusterIndex=1;
				
		//creating map of Clusters
		aglmClusters = new HashMap<Integer, ArrayList<Integer>>();
		clusterIndex=1;
		for(HashMap<Integer, ArrayList<Double>> eachMap:listOfMaps)
		{
			ArrayList<Integer> cluserList=new ArrayList<Integer>();
			for(Map.Entry<Integer, ArrayList<Double>> ent:eachMap.entrySet())
			{
				cluserList.add(ent.getKey());
			}
			aglmClusters.put(clusterIndex, cluserList);
			clusterIndex++;
		}
		
		TreeMap<Integer, Integer> gtTreeMap = new TreeMap<Integer, Integer>(geneIdGroundTruthMap);
		Set<Integer> gtKeys=gtTreeMap.keySet();
		ArrayList<Integer> gtKeysList=new ArrayList<Integer>();
		Iterator<Integer> gtIter=gtKeys.iterator();
		while(gtIter.hasNext())
			gtKeysList.add(gtIter.next());
		
		System.out.println("Jaccard coefficient is: " + objCommon.calculateJaccardCoeff(gtKeysList, tempClusterMap, geneIdGroundTruthMap));
		System.out.println("Silhoutte Index is: " + objCommon.calculateSilhoutte(geneInfoMap, aglmClusters));
		
		try {
			String filename = "AlgmClusters_new_dataset_2.txt";
			objCommon.generatePCAFile(aglmClusters, geneInfoMap, filename);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}


	private HashMap<Integer, Integer> writeToFile(String agglometricFileName) {
		int clusterIndex = 1;
		HashMap<Integer, Integer> tempClusterMap = new HashMap<Integer, Integer>();
		for(HashMap<Integer, ArrayList<Double>> eachMap:listOfMaps)
		{
			for(Map.Entry<Integer, ArrayList<Double>> ent:eachMap.entrySet())
			{
				tempClusterMap.put(ent.getKey(), clusterIndex);
	
			}
			clusterIndex++;
		}
		return tempClusterMap;
	}
}
