package com.talem;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

public class CommonMethods {
	
	
	public double calculateJaccardCoeff(ArrayList<Integer> gtKeysList, HashMap<Integer, Integer> tempClusterMap,
			HashMap<Integer, Integer> geneIdGroundTruthMap)
	{
		int m11=0;
		int m10=0,m01=0;
		
		int i,j;
		for(i=0;i<gtKeysList.size()-1;i++)
		{
			for(j=i+1;j<gtKeysList.size();j++)
			{
				int p1=gtKeysList.get(i);
				int p2=gtKeysList.get(j);
				if(geneIdGroundTruthMap.get(p1)==geneIdGroundTruthMap.get(p2) && 
						tempClusterMap.get(p1)==tempClusterMap.get(p2)){
					m11++;
				}
				else if(geneIdGroundTruthMap.get(p1)==geneIdGroundTruthMap.get(p2) && 
						tempClusterMap.get(p1)!=tempClusterMap.get(p2)){
					m10++;
				}	
				else if(geneIdGroundTruthMap.get(p1)!=geneIdGroundTruthMap.get(p2) && 
						tempClusterMap.get(p1)==tempClusterMap.get(p2)){
					m01++;
				}	
				
			}
		}
			
		return ((double)m11/(double)(m11+m01+m10));
	}
	
	
	
	public double calculateSilhoutte(HashMap<Integer, ArrayList<Double>> geneInfoMap,
			HashMap<Integer, ArrayList<Integer>> clusters)
	{
		double[][] distMat = new double[geneInfoMap.size() + 1][geneInfoMap.size() + 1];
		for (int i = 1; i <= geneInfoMap.size(); i++) {
			for (int j = 1; j <= geneInfoMap.size(); j++) {
				distMat[i][j] = getEuclideanDistance(geneInfoMap.get(i),geneInfoMap.get(j));
			}
		}
		
		Set<Integer> clusterIdkeys=clusters.keySet();
		ArrayList<Integer> clusterIdKeysList=new ArrayList<Integer>();
		Iterator<Integer> iterator=clusterIdkeys.iterator();
		while(iterator.hasNext())
		{
			clusterIdKeysList.add(iterator.next());
		}
		
		double ai=0;
		int clusterId=clusterIdKeysList.get(0);
		ArrayList<Integer> geneList=clusters.get(clusterId);
		int chosenPoint=geneList.get(0);
		for(int point:geneList)
		{
			if(point==chosenPoint)
				;
			else
			{
				ai=ai+distMat[chosenPoint][point];
			}
		}
		ai=ai/(double)(geneList.size()-1);
		
		ArrayList<Double> biList=new ArrayList<Double>();
		for(Map.Entry<Integer, ArrayList<Integer>> ent:clusters.entrySet())
		{
			double bi=0;
			if(ent.getKey()==clusterId)
				;	
			else
			{
				ArrayList<Integer> list=ent.getValue();
				for(int elem:list)
				{
					bi=bi+distMat[chosenPoint][elem];
				}
				bi=bi/(double)list.size();
				biList.add(bi);
			}
			
			
		}
		
		Collections.sort(biList);
		double silhoutte=0;
		if(!biList.isEmpty())
			silhoutte=Math.abs((biList.get(0)-ai)/(Math.max(biList.get(0), ai)));
		return silhoutte;
		
	}
	
	public void generatePCAFile(HashMap<Integer, ArrayList<Integer>> clusters, HashMap<Integer, ArrayList<Double>> geneInfoMap,
			String filename) throws IOException
	{
		FileWriter fwr=new FileWriter(new File(filename));
		for(Map.Entry<Integer, ArrayList<Integer>> ent:clusters.entrySet())
		{
			int size=ent.getValue().size();
			ArrayList<Integer> list=ent.getValue();
			for(int i=0;i<size;i++)
			{
				fwr.write(list.get(i)+"\t");
				fwr.write(ent.getKey()+"\t");
				ArrayList<Double> expressions = geneInfoMap.get(list.get(i));
				for(double expression:expressions){
					fwr.write(expression+"\t");
				}
				fwr.write("\n");
			}
		}
		
		fwr.close();
		System.out.println("PCA file Generated!");
	}
	
	public double getEuclideanDistance(ArrayList<Double> arr1, ArrayList<Double> arr2)
	{
		double distance=0.0;
		
		for(int i=0;i<arr1.size();i++)
		{
			distance=distance+(Math.pow(arr1.get(i)-arr2.get(i), 2));
		}
		distance=Math.sqrt(distance);
		
		return distance;
	}

}
