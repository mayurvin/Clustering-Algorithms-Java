public class Driver {
	public static int num_of_clusters = 5;
	public static int num_of_columns;
	public static Map<Integer, Integer> gene_cluster_kmeans;
	public static Map<Integer, Integer> gene_cluster_hierarchical;
	
	public static void main(String[] args) {
		// Generate a list of items.
		List<GeneExpression> geneSet = new ArrayList<GeneExpression>();	    
		String filename = args[0];
		if(filename.compareTo("cho.txt") == 0)
			num_of_clusters = 5;
		if(filename.compareTo("iyer.txt") == 0) {
			num_of_clusters = 10;
		}
			
		
		FileOp io = new FileOp(filename);
		geneSet = io.createInputs();
		num_of_columns = geneSet.get(0).size();
		Map<Integer, Integer> external_index = io.getExternalIndex();
		
		//Test K-Means algorithm
				
		List<Integer> clusterIndex = new ArrayList<Integer>();
		Iterator<Entry<Integer, Integer>> it = external_index.entrySet().iterator();
		
		while (it.hasNext()) {
			 Entry<Integer, Integer> entry = (Entry<Integer, Integer>) it.next();
	        clusterIndex.add(entry.getValue());
	    }
		System.out.println("Executing KMeans \n");
		KMeans(geneSet, clusterIndex, external_index);
		//End of K-Means  
		
		//Test DBScan
		System.out.println("\nExecuting DBScan \n");
		DBScan(geneSet, external_index);
		//End DBScan
		
		//Test Aggleromative clustering
		System.out.println("\nExecuting Hierarchical clustering \n");
		Hierarchical(geneSet, external_index);
		//End Aggleromative clustering
		
		System.out.println("\nExecuting Ensembling using re-labelling and voting \n");
		Map<Integer, Integer> gene_cluster_ensemble = ClusterEnsemble.cluster_ensemble(gene_cluster_kmeans, DBScanCluster.gene_cluster_dbscan, gene_cluster_hierarchical, num_of_clusters);
		
		ExternalIndexValidation externalIndexTest = new ExternalIndexValidation();
		System.out.println("External Index in Ensemble = " + externalIndexTest.validate(gene_cluster_ensemble, external_index));
		
		InternalIndexValidation internalIndexTest = new InternalIndexValidation();
		System.out.println("Internal Index in Ensemble = " + internalIndexTest.validate(gene_cluster_ensemble, geneSet));
		
	}
	
	public static void KMeans(List<GeneExpression> geneSet, List<Integer> clusterIndex,Map<Integer, Integer> external_index) {
		KMeans clusterer = new KMeans();
		
		
		Map<Integer, ArrayList<Integer>> results = clusterer.cluster(geneSet, num_of_clusters, clusterIndex);
		
		int cluster_id = 0;
		gene_cluster_kmeans = new HashMap<Integer,Integer>();
		Iterator<Entry<Integer, ArrayList<Integer>>> it = results.entrySet().iterator();
		while (it.hasNext()) {
	        Entry<Integer, ArrayList<Integer>> entry = (Entry<Integer, ArrayList<Integer>>) it.next();
	        List<Integer> gene_list = entry.getValue();
	        //System.out.println(gene_list);
	        for(int i = 0; i < gene_list.size(); i++) {
	        	gene_cluster_kmeans.put(gene_list.get(i), cluster_id);
	        }
	        cluster_id++;
	    }
		
		//System.out.println(gene_cluster_kmeans);
		
		ExternalIndexValidation externalIndexTest = new ExternalIndexValidation();
		System.out.println("External Index in K-Means = " + externalIndexTest.validate(gene_cluster_kmeans, external_index));
		
		InternalIndexValidation internalIndexTest = new InternalIndexValidation();
		System.out.println("Internal Index in K-Means = " + internalIndexTest.validate(gene_cluster_kmeans, geneSet));
	}
	
	public static void DBScan(List<GeneExpression> geneSet, Map<Integer, Integer> external_index) {
		DBScanCluster dbscanTest = new DBScanCluster();
		
		int minPts = num_of_clusters;
		double eps = dbscanTest.calculateEps(geneSet, minPts);
		//System.out.println("eps = " + eps);
		
		if(num_of_clusters == 5)
			eps = eps * 3;
		else
			eps = eps * 2.3;
		dbscanTest.DBScan(geneSet, eps, minPts);
		//System.out.println("Cluster size = " + DBScanCluster.clusterList.size());
		
		//System.out.println(DBScanCluster.gene_cluster_dbscan.size());
		//System.out.println(external_index.size());
		
		//System.out.println(DBScanCluster.gene_cluster_dbscan.containsValue(0));
		
		ExternalIndexValidation externalIndexTest = new ExternalIndexValidation();
		System.out.println("External Index in DBScan = " + externalIndexTest.validate(DBScanCluster.gene_cluster_dbscan, external_index));
		
		InternalIndexValidation internalIndexTest = new InternalIndexValidation();
		System.out.println("Internal Index in DBScan = " + internalIndexTest.validate(DBScanCluster.gene_cluster_dbscan, geneSet));
	}
	
	public static void Hierarchical(List<GeneExpression> geneSet, Map<Integer, Integer> external_index) {
		
		HierarchicalClustering test = new HierarchicalClustering();
		test.formClusters2(geneSet);
		//System.out.println(HierarchicalClustering.cluster_map.size());
		//System.out.println(HierarchicalClustering.cluster_map);
		
		int cluster_id = 0;
		gene_cluster_hierarchical = new HashMap<Integer,Integer>();
		Iterator<Entry<Integer, ArrayList<Integer>>> it = HierarchicalClustering.cluster_map.entrySet().iterator();
		
		while (it.hasNext()) {
	        Entry<Integer, ArrayList<Integer>> entry = (Entry<Integer, ArrayList<Integer>>) it.next();
	        List<Integer> gene_list = entry.getValue();
	        for(int i = 0; i < gene_list.size(); i++) {
	        	gene_cluster_hierarchical.put(gene_list.get(i), cluster_id);
	        }
	        cluster_id++;
	    }
		
		//System.out.println(gene_cluster_hierarchical.values());
		
		ExternalIndexValidation externalIndexTest = new ExternalIndexValidation();
		System.out.println("External Index in Hierarchical = " + externalIndexTest.validate(gene_cluster_hierarchical, external_index));
		
		InternalIndexValidation internalIndexTest = new InternalIndexValidation();
		System.out.println("Internal Index in Hierarchical = " + internalIndexTest.validate(gene_cluster_hierarchical, geneSet));
	}
}

public class FileOp {
	private File fileHandler;
	private Scanner scan;
	private int rows, columns;
	private Map<Integer,Integer> externalIndex;
	private List<GeneExpression> inputs;
	
	
	public FileOp(String nameOfFile){
		fileHandler = new File(nameOfFile);		
		try {
			scan = new Scanner(fileHandler);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}		
	}

	public List<GeneExpression> createInputs(){

		List<Double> geneData;		
		this.inputs = new ArrayList<GeneExpression>();		
		GeneExpression gene;
		int id;

		// ignore the first line;
		/*String line = scan.nextLine();
		String[] columns = line.split("\t");
		this.rows = Integer.parseInt(columns[0]);
		this.columns = Integer.parseInt(columns[1]);
		*/

		this.externalIndex = new TreeMap<Integer,Integer>();		

		while(scan.hasNext()){

			geneData = new ArrayList<Double>();
			String line = scan.nextLine();
			String[] splits = line.split("\t");

			int extIndex = Integer.parseInt((splits[1]));
			
			// exclude genes which are Outliers
			if(extIndex > -1){
				id = Integer.parseInt(splits[0]);

				this.externalIndex.put(id,extIndex);

				for(int index = 2; index<splits.length; index++){
					geneData.add(Double.parseDouble(splits[index]));
				}
				
				gene = new GeneExpression(geneData, id);
				this.inputs.add(gene);
				
			}
		}
		//System.out.println("input: "+inputs);
		return this.inputs;
	}

	public int getRowSize(){
		return this.rows;
	}

	public int getColumnSize(){
		return this.columns;
	}

	public Map<Integer,Integer> getExternalIndex(){
		return this.externalIndex;
	}

	/*
	public static void main(String[] args){
		FileOp f = new FileOp("/Users/saurabhtalbar/Documents/workspace/ClusteringAlgo/src/cho.txt");
		f.createInputs();
	}
	*/
}
public class GeneExpression {

	private List<Double> genes;
	private int id;

	public GeneExpression(List<Double> genes, int id){
		// size = 16 | 12
		this.genes = genes;
		this.id = id;
		//System.out.println("Gene Expression Array Initialized: "+ genes);
	}
	
	public GeneExpression average(List<GeneExpression> geneSet){
		List<Double> avgSet = new ArrayList<Double>();
		double avg;
		//System.out.println("geneSet: " +geneSet.size());

		for(int index1=0 ; index1 < Driver.num_of_columns; index1++){
			avg = 0;
			for(int index2 = 0; index2< geneSet.size(); index2++){
				avg += geneSet.get(index2).get(index1);				
			}
			avgSet.add(index1, avg/geneSet.size());
		}		
		return new GeneExpression(avgSet, -1);
	}

	public double eucDist(GeneExpression gene){
		double eucDist = 0;

		for(int index = 0; index < gene.size(); index++){
			double entry = this.genes.get(index) - gene.get(index);
			eucDist+= (entry*entry);
		}		

		return Math.sqrt(eucDist); 
	}

	public int getId(){
		return this.id;
	}

	public Double get(int index){
		return this.genes.get(index);
	}

	
	public int size(){
		return this.genes.size();
	}

	public String toString(){
		return this.genes.toString();
	}

	public void add(Double d){
		this.genes.add(d);
	}
}