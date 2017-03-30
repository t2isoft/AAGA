package algorithms;

import java.awt.Point;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;

public class Experimentation {

	public static void main(String[] args) throws FileNotFoundException {
		PrintStream output = new PrintStream(new FileOutputStream("Experimentation"));
		ArrayList<Point> points = new ArrayList<>();
		DefaultTeam dt = new DefaultTeam();
		for(int i=0; i<100;i++){
			points=DefaultTeam.readFromFile("tests/input"+i+".points");
			output.println(dt.calculConnectedDominatingSet(points, 55).size() +" "+
			dt.calculConnectedDominatingSetArticle(points, 55).size());
		}
		output.close();
	}

}
