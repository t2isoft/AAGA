package algorithms;

import java.awt.Point;
import java.util.ArrayList;
import java.util.Random;

public class DefaultTeam {

	public ArrayList<Point> calculCoupeMax(ArrayList<Point> inpoints, int edgeThreshold) {
		ArrayList<Point> points = (ArrayList<Point>)inpoints.clone();
		ArrayList<Point> result = new ArrayList<Point>();

		Random rand= new Random();
		for (int i=0;i<points.size();i++) if (rand.nextBoolean()) result.add(points.get(i));

		System.out.println("Resultat coupe: score = "+scoreCoupe(result, inpoints, edgeThreshold));
		return result;
	}
	
	public ArrayList<Point> calculEntraineur(ArrayList<Point> inpoints, int edgeThreshold) {
		ArrayList<Point> points = (ArrayList<Point>)inpoints.clone();
		ArrayList<Point> result = new ArrayList<Point>();
		int size = points.size();
		
		
		if (true) return result;
		
		
		// Arrange the points
		for (int k = 0; k < size; k++) {
			int bestIndex = -1;
			int lowestScore = Integer.MAX_VALUE;
			
			// Test all the points, we want the lowest score for the cup
			for (int i = 0; i < size; i++) {
				if (result.contains(points.get(i)))
					continue;
				
				ArrayList<Point> tmpCup = (ArrayList<Point>)result.clone();
				tmpCup.add(points.get(i));
				int tmpScore = scoreCoupe(tmpCup, points, edgeThreshold);
				if (tmpScore < lowestScore) {
					bestIndex = i;
					lowestScore = tmpScore;
				}
			}
			
			result.add(points.get(bestIndex));
			System.out.println("Lowest tmp score : " + lowestScore);
		}

		System.out.println("Resultat alignement : score = " + scoreEntraineur(result,edgeThreshold));
		System.out.println("Nombre dans alignement : " + result.size());
		if (result.size() == points.size())
			System.out.println("OK");
		else
			System.out.println("NON OK");
		System.out.println(""+inpoints.containsAll(result));
		return result;
	}

	//UTILITIES
	private ArrayList<Point> localSearch(ArrayList<Point> result, int edgeThreshold) {
		int currentScore = scoreEntraineur(result, edgeThreshold);
		int size = result.size();
		
		// Local-searching to improve the result
		for (int i = 0; i < size; i++) {
			for (int j = i + 1; j < size; j++) {
				ArrayList<Point> resultCopy = (ArrayList<Point>)result.clone();
				Point tmp = resultCopy.get(i);
				resultCopy.set(i, resultCopy.get(j));
				resultCopy.set(j, tmp);
				int tmpScore = scoreEntraineur(resultCopy, edgeThreshold);
				
				if (tmpScore < currentScore)
					return resultCopy;
			}
		}
		return null;
	}
	
	private int scoreCoupe(ArrayList<Point> coupe, ArrayList<Point> points, int edgeThreshold) {
		int s=0;
		for (Point p:coupe){
			for (Point q:points){
				if (coupe.contains(q)) continue;
				if (p.distance(q)<=edgeThreshold) s++;
			}
		}
		return s;
	}
	
	private int scoreEntraineur(ArrayList<Point> alignement, int edgeThreshold) {
		int max=0;
		for (int i=0;i<alignement.size();i++) {
			int s=scoreVi(alignement,i,edgeThreshold);
			if (max<s) max=s;
		}
		return max;
	}
	
	private int scoreVi(ArrayList<Point> alignement, int i, int edgeThreshold) {
		int s=0;
		for (int j=0;j<=i;j++){
			for (int k=i+1;k<alignement.size();k++){
				if (alignement.get(j).distance(alignement.get(k))<=edgeThreshold) s++;
			}
		}
		return s;
	}
	
	private int degree(ArrayList<Point> points, int i, int edgeThreshold){
		int d=-1;
		for (int j=0;j<points.size();j++) if (points.get(j).distance(points.get(i))<=edgeThreshold) d++;
		return d;
	}
}
