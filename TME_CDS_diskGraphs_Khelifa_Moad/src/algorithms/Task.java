package algorithms;

import java.awt.Point;
import java.util.ArrayList;
import java.util.Collections;
import java.util.concurrent.Callable;

public class Task implements Callable<ArrayList<Point>> {

	private ArrayList<Point> points;
	
	public Task(ArrayList<Point> points) {
		this.points=points;
	}
	@Override
	public ArrayList<Point> call() throws Exception {
	
		ArrayList<Point> result = new ArrayList<Point>();
		result = DefaultTeam.glouton(points);
		
		System.out.println("Score glouton= "+result.size());
		
		ArrayList<Point> resultFinal = new ArrayList<Point>(result);

		for(int i=0; i< 10; i++){
			
			Collections.shuffle(points);
			Collections.shuffle(result);
			ArrayList<Point> result2 = DefaultTeam.localSearchingWhile(points,result); 
			result = DefaultTeam.glouton(points);
			
			if(result2.size() < resultFinal.size()){
				resultFinal=result2;
				//result=resultFinal;
			}
			
			
			//System.out.println(i+" : "+result.size());
			System.out.println(i+ " : result2.size()=" +result2.size()+" resultFinal.size()="+resultFinal.size());
			
			//if(resultFinal.size()<78) break;
			
		}

		System.out.println(""+resultFinal.size());
		return resultFinal;
	}

}
