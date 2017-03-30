package algorithms;

import java.awt.Point;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Random;

public class Generateur {

	private static int nombreDePoints = 1000;
	private static int maxWidth = 1400;
	private static int maxHeight = 700;
	static byte edgeThreshold = 55;

	public static void main(String[] args) {
		Evaluation e = new Evaluation();
		try {
			for(int i = 0; i < 100; ++i) {
				PrintStream output = new PrintStream(new FileOutputStream("tests/input" + i + ".points"));
				Random generator = new Random();
				ArrayList<Point> points = new ArrayList<Point>();

				Point point = new Point(maxWidth/2, maxHeight/2);

				while(points.size()!=nombreDePoints){
					System.out.println(points.size());
					int x;
					int y;
					do {
						x = generator.nextInt(maxWidth);
						y = generator.nextInt(maxHeight);
						point = new Point(x, y);
					} while(!connecte(point,points) );
					points.add(point);	
					output.println(Integer.toString(x) + " " + Integer.toString(y));
				}
				output.close();
			}
		} catch (FileNotFoundException ex) {
			ex.printStackTrace();
		}
	}

	private static boolean connecte(Point point, ArrayList<Point> points) {
		if(points.isEmpty()) return true;
		for(Point p: points){
			if(point.distance(p)<=(double)edgeThreshold) return true;
		}
		return false;
	}

}
