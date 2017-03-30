package algorithms;

import java.awt.Point;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

public class Kruskal {
	
	private static int compteur = 0;

	public static ArrayList<Edge> calcul_krusk(ArrayList<Point> pts, ArrayList<Edge> ledge) {

		HashMap<Point, Integer> pet = new HashMap<Point, Integer>();

		for (Point p : pts)
			pet.put(p, -1);

		Collections.sort((List<Edge>) ledge);

		ArrayList<Edge> solutions = new ArrayList<Edge>();
		for (Edge e : ledge) {
			cycleExists(solutions, e, pet, pts);

		}
		return solutions ;
		//Point first = solutions.get(0).p1;
		//System.out.println(solutions.size());
		//return convert(first, solutions);

	}

	public Tree2D convert(Point p, ArrayList<Edge> solution) {
		if (solution.size() == 0)
			return new Tree2D(p, new ArrayList<Tree2D>());
		else {
			ArrayList<Tree2D> fils = new ArrayList<Tree2D>();
			ArrayList<Edge> edge_p = new ArrayList<Edge>();
			for (Edge e : solution) {
				if (p.equals(e.p1) || p.equals(e.p2))
					edge_p.add(e);

			}
			for (Edge e : edge_p) {
				if (e.p1.equals(p)) {

					ArrayList<Edge> copySol = new ArrayList<Edge>();
					copySol.addAll(solution);
					copySol.removeAll(edge_p);
					fils.add(convert(e.p2, copySol));
				}
				if ((e.p2.equals(p))) {
					ArrayList<Edge> copySol = new ArrayList<Edge>(solution);
					copySol.removeAll(edge_p);
					fils.add(convert(e.p1, copySol));
				}

			}

			return new Tree2D(p, fils);
		}
	}

	public static void cycleExists(List<Edge> solutions, Edge e,
			HashMap<Point, Integer> pet, ArrayList<Point> pts) {

		// int et = pet.get(solutions.get(0).p1);
		if (pet.get(e.p1) == -1 && pet.get(e.p2) == -1) {
			pet.put(e.p1, compteur);
			pet.put(e.p2, compteur);
			solutions.add(e);
			compteur++;

		} else if (pet.get(e.p1) == -1 && pet.get(e.p2) != -1) {

			pet.put(e.p1, pet.get(e.p2));
			solutions.add(e);

		} else if (pet.get(e.p1) != -1 && pet.get(e.p2) == -1) {

			pet.put(e.p2, pet.get(e.p1));
			solutions.add(e);

		} else if (pet.get(e.p1) == pet.get(e.p2)) {

		} else {
			if (pet.get(e.p1) != pet.get(e.p2)) {

				for (Point p : pts) {
					if (p.equals(e.p1))
						continue;
					if (pet.get(e.p1) == pet.get(p))
						pet.put(p, pet.get(e.p2));
				}
				pet.put(e.p1, pet.get(e.p2));

				solutions.add(e);
			}
		}

	}
}
