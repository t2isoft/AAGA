package algorithms;

import java.awt.Point;

public class Edge implements Comparable<Edge> {
	
	Point p1 ; 
	Point p2;
	double poids ; 

	public Edge(Point p1, Point p2, double poid){
		this.p1=p1;
		this.p2=p2;
		this.poids=poid ;
		
	}

	@Override
	public int compareTo(Edge e2) {
		return Double.compare(this.poids, e2.poids);
	}
	
	
	@Override
	public String toString() {
		return "Edge: p1: " + p1.toString() + " - p2: "+p2.toString();
	}

}
