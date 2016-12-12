package com.talem;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;
import java.util.Scanner;

public class ClustringDriver {

	public static void main(String[] args) {
		HashMap<Integer, ArrayList<Integer>> groundTruthMap=new HashMap<Integer, ArrayList<Integer>>();
		HashMap<Integer, Integer> geneIdGroundTruthMap = new HashMap<Integer,Integer>();
		HashMap<Integer, ArrayList<Double>> geneInfoMap = new HashMap<Integer, ArrayList<Double>>();
		
		ArrayList<HashMap<Integer, ArrayList<Double>>> listOfMaps=new ArrayList<HashMap<Integer,ArrayList<Double>>>();
		HashMap<Integer, ArrayList<Integer>> aglmClusters = new HashMap<Integer, ArrayList<Integer>>();

		int sampleSize = 0;
		
		Scanner sc = new Scanner(System.in);
		System.out.println("Enter file name: ");
		String inputFileName = sc.next();
		System.out.println("Enter value of K: ");
		int k = sc.nextInt();
		System.out.println("Enter number of iterations/100 for infinite: ");
		int rotation = sc.nextInt();
				
		try {
			
			FileReader fileReader = new FileReader(inputFileName);
			BufferedReader bReader = new BufferedReader(fileReader);
			String line = "";
			int initialClustId = 1;
			
			while ((line = bReader.readLine()) != null)
			{
				sampleSize++;
			}
			fileReader = new FileReader(inputFileName);
			bReader = new BufferedReader(fileReader);
			while((line = bReader.readLine()) != null){
				String lSplit[]=line.split("\t");
				if (lSplit.length == 2) {
					sampleSize = Integer.parseInt(lSplit[0]);
					continue;
				}
				int geneId=Integer.parseInt(line.split("\t")[0]);
			    int groundTruthId=Integer.parseInt(line.split("\t")[1]);
			    ArrayList<Double> list=new ArrayList<Double>();
			    int i=0;
			    for(int j=2;j<lSplit.length;j++)
			    {
			    	list.add(i, Double.parseDouble(lSplit[j]));
			    	i++;
			    }
			    if(!geneInfoMap.containsKey(geneId)){
			    	geneInfoMap.put(geneId, list);
			    }
			    if(groundTruthMap.containsKey(groundTruthId)){
			    	groundTruthMap.get(groundTruthId).add(geneId);
			    }
			    else {
			    	ArrayList<Integer> geneIdList = new ArrayList<Integer>();
			    	geneIdList.add(geneId);
			    	groundTruthMap.put(groundTruthId, geneIdList);
			    }
			    if(!geneIdGroundTruthMap.containsKey(geneId)){
			    	geneIdGroundTruthMap.put(geneId, groundTruthId);
			    }
			    
			    /////////////////Aglometric//////////////////////////
				HashMap<Integer, ArrayList<Double>> geneExpressions=new HashMap<Integer, ArrayList<Double>>();
				geneExpressions.put(geneId,list);
				listOfMaps.add(geneExpressions); 
				
				ArrayList<Integer> individualGeneList=new ArrayList<Integer>();
				individualGeneList.add(geneId);
				aglmClusters.put(initialClustId,individualGeneList);
				initialClustId++;
				/////////////////////////////////////////////////////

			}
			
			// Uncomment for K-Means Clustering algorithm
//			HashMap<Integer, ArrayList<Double>> centroidMap = generateInitialCentroids(k, geneInfoMap, sc);
//			KMeansImpl objKMeans = new KMeansImpl(geneInfoMap,centroidMap,groundTruthMap, geneIdGroundTruthMap, k, rotation);
//			objKMeans.KMEansImplementation();
			
			// Uncomment for Hierarchial Aglometric Clustering algorithm
//			AglometricImpl objAgmt = new AglometricImpl(listOfMaps, geneInfoMap, geneIdGroundTruthMap, aglmClusters, k);
//			objAgmt.AglometricImplementation();
			
			// Uncomment for Density Based Clustering algorithm 
			DBScanImpl objDBScan = new DBScanImpl(geneInfoMap, geneIdGroundTruthMap);
			objDBScan.DBScanImplementation(sampleSize);
			
			bReader.close();
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
	}
	

	public static HashMap<Integer, ArrayList<Double>> generateInitialCentroids(int k, 
									HashMap<Integer, ArrayList<Double>> geneInfoMap, Scanner sc){
		HashMap<Integer, ArrayList<Double>> centroidMap = new HashMap<Integer, ArrayList<Double>>();
		Random rand = new Random();
		//Scanner sc = new Scanner(System.in);
		int[] arr = new int[k];
		System.out.println("Enter random to generate clusters using rand(): ");
		String str = sc.next();
		//Random rand = new Random();
		if(str.equals("random")){
			for(int i=0;i<k;i++){
				arr[i] = rand.nextInt(geneInfoMap.size());
			}
		}
		else {
			System.out.println("Enter initial cluster values: ");
			for(int i=0;i<k;i++){
				arr[i] = sc.nextInt();
			}
		}
		ArrayList<Integer> geneIdSet = new ArrayList<Integer>(geneInfoMap.keySet());
//		Scanner scan = new Scanner(System.in);
//		System.out.println("Enter initial cluster values: ");
		for(int i=1;i<=k;i++){
			//int randomGeneId = geneIdSet.get(random.nextInt(geneInfoMap.size()));
			int randomGeneId = geneIdSet.get(arr[i-1]);
			if(randomGeneId>-1){
				centroidMap.put(i, geneInfoMap.get(randomGeneId));
			}
		}
		return centroidMap;
	}
	

}
