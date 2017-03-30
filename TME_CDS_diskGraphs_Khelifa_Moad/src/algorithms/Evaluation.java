package algorithms;
import java.util.ArrayList;
import java.awt.Point;

public class Evaluation {
	
  private boolean isMember(ArrayList<Point> points, Point p){
    for (Point point:points) if (point.equals(p)) return true; return false;
  }
  
//  public boolean isValide(ArrayList<Point> origPoints, ArrayList<Point> ensDom, int edgeThreshold){
//   
//	  ArrayList<Point> points = new ArrayList<Point>(origPoints);
//	  for(Point p: ensDom){
//		  points.removeAll(neighbor(p, points, edgeThreshold));
//		  points.remove(p);
//		  if(points.isEmpty()) return true;
//	  }
//	  
//    return points.isEmpty();
//  }
  
  public boolean isValide(ArrayList<Point> domSet, ArrayList<Point> points, int edgeThreshold){
	    for (Point p:points){
	      boolean isDom=false;
	      for (Point q:domSet) 
	    	  if (isDom=p.distance(q)<edgeThreshold) 
	    		  break;
	      if (!isDom) return false;
	    }
	    return true;
	  }
  public ArrayList<Point> neighbor(Point p, ArrayList<Point> vertices, int edgeThreshold){
    ArrayList<Point> result = new ArrayList<Point>();

    for (Point point:vertices) if (point.distance(p)<edgeThreshold && !point.equals(p)) result.add((Point)point.clone());

    return result;
  }
}
