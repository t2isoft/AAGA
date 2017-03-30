package algorithms;

import java.awt.Point;
import java.util.ArrayList;
import java.util.Random;

public class DefaultTeam {

	public final int K = 5;
	
	public ArrayList<ArrayList<Point>> calculKMeans(ArrayList<Point> points) {
		ArrayList<ArrayList<Point>> kmeans = new ArrayList<ArrayList<Point>>();
		
		/* ---------- Initialize the means */
		/* KMeans++ */
		Point[] centers = kmeansPlusPlus(points);
		
		/* Pure random */
		//Point[] centers = randomMeans(points, points.size());
		
		for (int iter = 0; iter < 10; iter++) {
			
			/* ---------- Clear all the means for the new iteration ---------- */
			
			kmeans = new ArrayList<ArrayList<Point>>();
			for (int i = 0; i < K; i++)
				kmeans.add(new ArrayList<Point>());
			
			/* ---------- Assignment step ---------- */
			
			/* Assign each point to the clusters (Voronoi diagram) */
			for (Point p: points) {
				double min = Double.POSITIVE_INFINITY;
				int candidate = -1;
				
				for (int i = 0; i < K; i++) {
					double squareDist = p.distance(centers[i]);
					
					if (squareDist < min) {
						candidate = i;
						min = squareDist;
					}
				}
				kmeans.get(candidate).add(p);
			}
			
			/* ---------- Update step ---------- */
			
			/* Compute the new centroids */
			for (int i = 0; i < K; i++) {
				double x = 0, y = 0;
				int card = kmeans.get(i).size();
				
				for (Point p: kmeans.get(i)) {
					x += p.getX();
					y += p.getY();
				}
				centers[i] = new Point((int)(x / card), (int)(y / card));
			}
		}
		return kmeans;
	}
	
	/* KMeans random choice of the initial means */
	public Point[] randomMeans(ArrayList<Point> points, int size) {
		Random rand = new Random();
		Point[] centers = new Point[K];
		
		
		for (int i = 0; i < K; i++)
			centers[i] = points.get(rand.nextInt(size));
		
		return centers;
	}
	
	/* KMeans++ randomly chooses the centers at the beginning
	 * Algorithm from wikipedia */
	public Point[] kmeansPlusPlus(ArrayList<Point> points) {
		int size = points.size();	
		Random rand = new Random();
		int center = rand.nextInt(size);
		Point[] centers = new Point[K];
		
		centers[0] = points.get(center);
		
		for (int k = 1; k < K; k++) {
			double sum = 0;
			double[] D = new double[size];
			
			for (int i = 0; i < size; i++) {
				Point p = points.get(i);
				double min = centers[0].distance(p);
				
				for (int j = 1; j < k; j++) {
					double d = centers[j].distance(p);
					
					if (d < min)
						min = d;
				}
				D[i] = min;
				sum += min;
			}
			
			for (int i = 0; i < size; i++)
				D[i] /= sum;
			
			centers[k] = points.get(weightedRandom(D));
		}
		return centers;
	}
	
	/* For KMeans++ algorithm */
	public int weightedRandom(double[] D) {
		Random rand = new Random();
		double uniform = rand.nextDouble();
		double cumulativeSum = 0;
		
		for (int i = 0; i < D.length; i++) {
			cumulativeSum += D[i];
			if (uniform < cumulativeSum)
				return i;
		}
		
		// Should not occur
		return rand.nextInt(D.length);
	}

	public ArrayList<ArrayList<Point>> calculKMeansBudget(ArrayList<Point> points) {
		ArrayList<Point> rouge = new ArrayList<Point>();
		ArrayList<Point> verte = new ArrayList<Point>();

		for (int i=0;i<points.size()/2;i++){
			rouge.add(points.get(i));
			verte.add(points.get(points.size()-i-1));
		}
		if (points.size()%2==1) rouge.add(points.get(points.size()/2));

		ArrayList<ArrayList<Point>> kmeans = new ArrayList<ArrayList<Point>>();
		kmeans.add(rouge);
		kmeans.add(verte);


		/* WRITE */



		/* WRITE */


		return kmeans;
	}
}
