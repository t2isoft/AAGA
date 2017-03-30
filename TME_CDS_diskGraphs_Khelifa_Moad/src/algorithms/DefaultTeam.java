package algorithms;

import java.awt.Point;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.FutureTask;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;

import algorithms.Vertex.Color;

public class DefaultTeam {
	public static int infini = Integer.MAX_VALUE;
	public static int edgeThreshold;
	public ArrayList<Point> calculConnectedDominatingSet(ArrayList<Point> points, int edgeThreshold) {
		this.edgeThreshold = edgeThreshold;
		Kruskal kruskal = new Kruskal();

		ArrayList<Point> ds =  multiThread(points, 2);
		//ArrayList<Point> mis = maximalIndependentSet(points);		
		//ArrayList<Point> mis = MIS(points);		
		Tree2D tree2D = calculSteiner(points, edgeThreshold, ds);

		ArrayList<Point> res = edgesTosommets(kruskal(edgesTosommets(treeToEdges(tree2D))));
		//	return mis;
		//ArrayList<Point> res = lialAlgo(mis, points);
		//ArrayList<Point> res = algoritmeA(mis, points);
		
		System.out.println("isDom="+(new Evaluation()).isValide(res, points, edgeThreshold) +" res.size="+res.size());
		return res;
	}

	public ArrayList<Point> calculConnectedDominatingSetArticle(ArrayList<Point> points, int edgeThreshold) {
		this.edgeThreshold = edgeThreshold;

		//		ArrayList<Point> mis = maximalIndependentSet(points);		
		ArrayList<Point> mis = MIS(points);
		ArrayList<Point> res = algoritmeA(mis, points);
		System.out.println("isDom="+(new Evaluation()).isValide(res, points, edgeThreshold));
		return res;
	}

	private ArrayList<Point> maximalIndependentSet(ArrayList<Point> points) {
		ArrayList<Point> pointsCopie = (ArrayList<Point>)points.clone();
		ArrayList<Point> I = new ArrayList<>();
		Point u ; Evaluation e = new Evaluation();
		int i=0;
		while (!pointsCopie.isEmpty()) {
			u = pointsCopie.get(0);
			I.add(u);
			pointsCopie.removeAll(e.neighbor(u, pointsCopie, edgeThreshold));
			pointsCopie.remove(u);
		}
		return I;
	}

	private ArrayList<Point> MIS(ArrayList<Point> points) {
		ArrayList<Vertex> V = new ArrayList<>();
		Evaluation e = new Evaluation();
		ArrayList<Point> mis = new ArrayList<>();
		int i=0;

		for(Point p: points){
			Vertex v=new Vertex(p, i++);
			V.add(v);
		}

		Vertex s = V.get(0);
		s.color=Color.BLACK;
		mis.add(s);
		for(Vertex u:neighbor(s, V, edgeThreshold)){
			u.color=Color.GREY;

			ArrayList<Vertex> tmp = neighbor(u, V, edgeThreshold);
			for(Vertex p:tmp){
				p.setActif(true);
			}
		}
		while( (s=existeWhite(V))!=null ){

			Vertex actifMax=null ;
			for(Vertex u: V){
				if(u.color==Color.WHITE && u.isActif()){ 
					actifMax = u;
					break;
				}
			}
			if(actifMax==null) continue;
			for(Vertex u: V){
				if(u.actif && u.color==Color.WHITE){
					if(neighborWhite(u,V).size() > neighborWhite(actifMax,V).size()){
						actifMax=u;
					}
				}
			}
			actifMax.color=Color.BLACK;
			mis.add(actifMax);
			s=actifMax;

				for(Vertex u:neighbor(s, V, edgeThreshold)){
					u.color=Color.GREY;
	
					ArrayList<Vertex> tmp = neighbor(u, V, edgeThreshold);
					for(Vertex p:tmp){
						p.setActif(true);
					}
			}
		}
		return mis;
		//return VtoPoint(V);
	}

	private ArrayList<Vertex> neighborWhite(Vertex u, ArrayList<Vertex> v) {
		ArrayList<Vertex> neib = new ArrayList<>();
		for(Vertex r:neighbor(u, v, edgeThreshold)){
			if(r.color==Color.WHITE) neib.add(r);
		}
		return neib;
	}

	private Vertex existeWhite(ArrayList<Vertex> V) {
		for(Vertex v:V){
			if(v.color==Color.WHITE) return v;
		}
		return null;
	}

	public ArrayList<Point> algoritmeA(ArrayList<Point> MIS, ArrayList<Point> points){


		ArrayList<Vertex> vs = new ArrayList<Vertex>();
		

		for(Point p: points){
			if(!MIS.contains(p)){
				Vertex g = new Vertex(p,0);
				g.color=Color.GREY;
				vs.add(g);
			}
			else{
				Vertex b = new Vertex(p,0);
				b.color=Color.BLACK;
				vs.add(b);
			}
		}
		System.out.println(vs.size());
		
		lialAlgo(vs);
		System.out.println(vs.size());
		return VBBtoPoint(vs);

	}
	public void lialAlgo(ArrayList<Vertex> vs){
		Map<Integer, Set<Vertex>> comps = new HashMap<>();
		Integer compId = 0;
		for(compId = 0; compId < vs.size(); compId++){
			if(vs.get((int)compId).color == Color.BLACK){
				Set<Vertex> s = new HashSet<>();
				s.add(vs.get((int)compId));
				comps.put(compId, s);
			}
		}

		for(int i = 5; i > 1; i--){
			boolean existGrey;
			do{
				existGrey = false;
				for(Vertex g : vs){
					if(g.color == Color.GREY){
						Set<Integer> foundComp = new HashSet<>();
						for(Vertex b : neighbor(g, vs, edgeThreshold)){
							if(b.color == Color.BLACK){
								foundComp.add(getComp(comps, b));
							}
						}
						if(foundComp.size() >= i){
							existGrey = true;
							g.color = Color.BLUE;
							Set<Vertex> newSet = new HashSet<Vertex>();
							for(Integer c : foundComp){
								newSet.addAll(comps.get(c));
								comps.remove(c);
							}
							comps.put(++compId, newSet);
						}
					}
				}
			}while(existGrey);
		}
	}
	public int getComp(Map<Integer, Set<Vertex>> comp, Vertex v){
		for(Integer key : comp.keySet()){
			if(comp.get(key).contains(v))
				return key;
		}
		
		return -1;
	}

	private Point grayAdja(int i, ArrayList<Point> grey, ArrayList<Point> black, ArrayList<Point> blue) {
		Evaluation e = new Evaluation();
		ArrayList<Point> black_blue = new ArrayList<>(black);
		black_blue.addAll(blue);
		for(Point g: grey){	
			if(e.neighbor(g, black_blue, edgeThreshold).size()<=i){
				return g;
			}
		}

		return null;
	}

	public ArrayList<Vertex> neighbor(Vertex p, ArrayList<Vertex> vertices, int edgeThreshold){
		ArrayList<Vertex> result = new ArrayList<Vertex>();

		for (Vertex point:vertices) {
			if (point.distance(p)<edgeThreshold && !point.equals(p)) 
				result.add(point);
		}

		return result;
	}

	ArrayList<Point> DM(ArrayList<Point> points){
		ArrayList<Vertex> V = new ArrayList<>();
		Evaluation e = new Evaluation();
		int i=0;

		for(Point p: points){
			Vertex v=new Vertex(p, i++);
			v.whiteCount=e.neighbor(p, points, edgeThreshold).size();
			V.add(v);
		}

		for(Vertex v:V){
			if(v.color==Color.WHITE){
				boolean	bool=true;
				for(Vertex u: neighbor(v, V, edgeThreshold)){
					if(! (v.whiteCount>=u.whiteCount) )
						bool=false;
				}
				if(bool==true){
					v.color=Color.BLACK;
					for(Vertex u1: neighbor(v, V, edgeThreshold)){
						if(u1.color==Color.WHITE){
							u1.color=Color.GREY;
						}
					}
				}
			}
		}
		for(Vertex v:V){
			if(v.color==Color.WHITE){
				int high=0;

				int max = Integer.MIN_VALUE;
				ArrayList<Vertex> neighbor=neighbor(v, V, edgeThreshold);
				for(int j=0; j<neighbor.size();j++ ){
					Vertex u= neighbor.get(j);
					if(u.whiteCount>max){
						max = u.whiteCount;
						high=j;
					}
				}
				V.get(high).color=Color.BLACK;
				for(Vertex u:neighbor(V.get(high), V, edgeThreshold)){
					if(u.color==Color.WHITE){
						u.color=Color.GREY;
					}
				}
			}
		}


		ArrayList<Point> ds = new ArrayList<>();

		for(Vertex v : V){
			if(v.color==Color.BLACK || v.color==Color.GREY){
				ds.add(new Point(v.x, v.y));
			}
		}
		return ds;
	}

	public ArrayList<Point> VBBtoPoint(ArrayList<Vertex> V) {
		ArrayList<Point> ds = new ArrayList<>();

		for(Vertex v : V){
			if(v.color==Color.BLACK || v.color==Color.BLUE ){
				ds.add(new Point(v.x, v.y));
				//ds.add(v);
			}
				
		}
		return ds;
	}
	public ArrayList<Point> VtoPoint(ArrayList<Vertex> V) {
		ArrayList<Point> ds = new ArrayList<>();

		for(Vertex v : V){
			if(v.color==Color.BLACK ){
				ds.add(new Point(v.x, v.y));
			}
		}
		return ds;
	}
	private ArrayList<Point> tree2DtoPoint(Tree2D tree2d) {
		ArrayList<Edge> edges = treeToEdges(tree2d);
		ArrayList<Point> res = new ArrayList<>();

		for(Edge e : edges){
			res.add(e.p1);
			res.add(e.p2);
		}
		return res;
	}
	private ArrayList<Point> edgesTosommets(ArrayList<Edge> kruskal) {

		Set<Point> set = new HashSet<Point>();
		for(Edge e: kruskal){
			set.add(e.p1);
			set.add(e.p2);
		}
		System.out.println("edgesTosommets set.size="+set.size());
		return new ArrayList<Point>(set);
	}



	public Tree2D calculSteiner(ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints) {
		Kruskal krusk = new Kruskal();

		double[][] D = new double[points.size()][points.size()];
		int[][] P = new int[points.size()][points.size()];

		for (int i = 0; i < points.size(); i++)
			for (int j = 0; j < points.size(); j++) {
				Point p = points.get(i);
				Point q = points.get(j);
				double distance = p.distance(q);
				if (distance > edgeThreshold) {
					D[i][j] = Double.POSITIVE_INFINITY;
				} else {
					D[i][j] = distance;
				}
				P[i][j] = j;
			}

		for (int k = 0; k < D.length; k++) {
			for (int i = 0; i < D.length; i++) {
				for (int j = 0; j < D.length; j++) {
					if (D[i][k] + D[k][j] < D[i][j]) {
						D[i][j] = D[i][k] + D[k][j];
						P[i][j] = P[i][k];
					}
				}
			}
		}
		ArrayList<Edge> ledge = new ArrayList<Edge>();// liste des aretes des
		// points rouges
		// calcule d'arbre couvrant des points rouges
		for (Point p : hitPoints)
			for (Point q : hitPoints) {
				if (p.equals(q))
					continue;
				ledge.add(new Edge(p, q, D[points.indexOf(p)][points.indexOf(q)]));
			}
		ArrayList<Edge> t1 = Kruskal.calcul_krusk(hitPoints, ledge);
		// solutionsHS pour eliminer les doublons
		HashSet<Edge> solutionsHS = new HashSet<Edge>();
		for (Edge e : Kruskal.calcul_krusk(points, t1)) {
			solutionsHS.add(e);
		}
		ArrayList<Edge> solutions = new ArrayList<Edge>(solutionsHS);
		ArrayList<Edge> ledge_rouge_final = local_search(solutions, points);
		ArrayList<Edge> ledge_final = ajout_plus_court_chemin(ledge_rouge_final, D, P, points);
		Point first = ledge_final.get(0).p1;
		Tree2D arbre = krusk.convert(first, ledge_final);
		return arbre;

	}
	private boolean contains (ArrayList<Edge> edges, Point p, Point q){
		for(Edge e: edges){
			if(e.p1.equals(p) && e.p2.equals(q) ||
					e.p1.equals(q) && e.p2.equals(p) ) return true;
		}
		return false;
	}
	private ArrayList<Edge> kruskal(ArrayList<Point> points) {

		ArrayList<Edge> edges = new ArrayList<Edge>();
		Edge e;
		for(int i=0;i<points.size();i++){
			for(int j=0;j<points.size();j++){
				if(points.get(i).equals(points.get(j)) || contains(edges, points.get(i), points.get(j))) continue;
				e = new Edge(points.get(i), points.get(j),points.get(i).distance( points.get(j)) );
				edges.add(e);
			}
		}
		Collections.sort(edges);

		ArrayList<Edge> kruskal = new ArrayList<Edge>();
		Edge current ;

		NameTag forest = new NameTag(points);

		while (edges.size()!=0){
			current=edges.remove(0);
			if(forest.tag(current.p1)!=forest.tag(current.p2)){
				kruskal.add(current);
				forest.reTag(forest.tag(current.p1), forest.tag(current.p2));
			}
		}

		return kruskal;
	}
	// local searching
	public ArrayList<Edge> local_search(ArrayList<Edge> ledge, ArrayList<Point> pts) {
		Point first = ledge.get(0).p1;
		Kruskal krusk = new Kruskal();
		Tree2D arbre = krusk.convert(first, ledge);
		double oldscore = score(arbre);
		double newscore = score(arbre) - 1;
		int i = 0;
		while (oldscore > newscore || i < 50) {
			i++;
			oldscore = newscore;
			ledge = local_opt(ledge, pts);
			ledge = Kruskal.calcul_krusk(pts, ledge);
			first = ledge.get(0).p1;
			arbre = krusk.convert(first, ledge);
			newscore = score(arbre);
		}
		return ledge;
	}
	private ArrayList<Edge> local_opt(ArrayList<Edge> ledge, ArrayList<Point> pts) {
		ArrayList<Point> tag = extrairepoints(ledge);
		for (Point p : tag) {
			ArrayList<Point> res = new ArrayList<Point>();
			res = voisins(p, ledge);
			if (res.size() == 1)
				continue;
			Point f = Fermat(p, res.get(0), res.get(1));
			ArrayList<Point> pointsproches = plusproche(f, pts);
			Point g = plus_optimal(pointsproches, p, res);

			if (g.distance(p) + g.distance(res.get(0)) + g.distance(res.get(1)) < p.distance(res.get(0))
					+ p.distance(res.get(1))) {
				// recuperer les arrete a supprimer
				ArrayList<Edge> supp = asupprimer(ledge, p, res.get(0), res.get(1));
				for (int i = 0; i < supp.size(); i++)
					ledge.remove(ledge.indexOf(supp.get(i)));
				ledge.add(new Edge(g, p, g.distance(p)));
				ledge.add(new Edge(g, res.get(0), g.distance(res.get(0))));
				ledge.add(new Edge(g, res.get(1), g.distance(res.get(1))));

			}
		}
		return ledge;
	}
	// calcul point-fermat de 3 points111
	public Point Fermat(Point A, Point B, Point C) {
		double AB = A.distance(B);
		double BC = C.distance(B);
		double CA = A.distance(C);
		double angleA = angle(A, B, C);
		double angleB = angle(B, A, C);
		double angleC = angle(C, B, A);
		if (((angleA * 180 / Math.PI) < 120) && ((angleB * 180 / Math.PI) < 120) && ((angleC * 180 / Math.PI) < 120)) {
			Point F = new Point(
					(int) ((BC * A.x / Math.sin(angleA + Math.PI / 3) + CA * B.x / Math.sin(angleB + Math.PI / 3)
							+ AB * C.x / Math.sin(angleC + Math.PI / 3))
							/ (BC / Math.sin(angleA + Math.PI / 3) + CA / Math.sin(angleB + Math.PI / 3)
									+ AB / Math.sin(angleC + Math.PI / 3))),
									(int) ((BC * A.y / Math.sin(angleA + Math.PI / 3) + CA * B.y / Math.sin(angleB + Math.PI / 3)
											+ AB * C.y / Math.sin(angleC + Math.PI / 3))
											/ (BC / Math.sin(angleA + Math.PI / 3) + CA / Math.sin(angleB + Math.PI / 3)
													+ AB / Math.sin(angleC + Math.PI / 3))));
			return F;
		} else if ((angleA * 180 / Math.PI) >= 120)
			return A;
		else if ((angleB * 180 / Math.PI) >= 120)
			return B;
		else
			return C;
	}
	public double angle(Point C, Point B, Point A) {
		double a = C.distance(B);
		double b = C.distance(A);
		double c = B.distance(A);
		return Math.acos((Math.pow(a, 2) + Math.pow(b, 2) - Math.pow(c, 2)) / (2 * a * b));
	}
	// renvoie la liste des point qui sont connectes direct avec le point p
	private ArrayList<Point> voisins(Point p, ArrayList<Edge> ledge) {
		ArrayList<Point> res = new ArrayList<Point>();
		for (Edge e : ledge) {
			if (e.p1.equals(p)) {
				if (!res.contains(e.p2))
					res.add(e.p2);
			} else if (e.p2.equals(p)) {
				if (!res.contains(e.p1))
					res.add(e.p1);
			}
		}
		return res;
	}
	// 5 points les plus proches
	public ArrayList<Point> plusproche(Point b, ArrayList<Point> pts) {
		ArrayList<Point> res = new ArrayList<Point>();
		ArrayList<Point> pts_aux = new ArrayList<Point>(pts);
		while (res.size() < 50) {
			double dist_min = Double.MAX_VALUE;
			Point pres = null;
			for (Point p : pts_aux)
				if (p.distance(b) < dist_min) {
					dist_min = p.distance(b);
					pres = p;
				}
			res.add(pres);
			pts_aux.remove(pres);
		}
		return res;
	}

	// calcule du chemin en prenant en compte les points bleu
	private ArrayList<Edge> ajout_plus_court_chemin(ArrayList<Edge> t0, double[][] D, int[][] P,
			ArrayList<Point> points) {

		ArrayList<Edge> t1 = new ArrayList<Edge>();
		for (Edge a : t0) {
			if (points.get(P[points.indexOf(a.p1)][points.indexOf(a.p2)]).equals(a.p2)) {

				t1.add(a);
			} else {
				int i = P[points.indexOf(a.p1)][points.indexOf(a.p2)];
				Edge e = new Edge(a.p1, points.get(i), D[points.indexOf(a.p1)][i]);
				t1.add(e);

				int j;
				while (i != points.indexOf(a.p2)) {
					j = P[points.indexOf(points.get(i))][points.indexOf(a.p2)];

					Edge e2 = new Edge(points.get(i), points.get(j), D[points.indexOf(a.p1)][j]);

					t1.add(e2);
					i = j;

				}

			}

		}
		return t1;
	}
	private Point plus_optimal(ArrayList<Point> pointsproches, Point p, ArrayList<Point> res) {
		Point g = null;
		double min = Double.MAX_VALUE;
		for (Point proche : pointsproches)
			if (proche.distance(p) + proche.distance(res.get(0)) + proche.distance(res.get(1)) < min) {
				min = proche.distance(p) + proche.distance(res.get(0)) + proche.distance(res.get(1));
				g = proche;
			}
		return g;
	}
	// extraire tout les points de la liste d'arretes ledge
	private ArrayList<Point> extrairepoints(ArrayList<Edge> ledge) {
		ArrayList<Point> res = new ArrayList<Point>();
		for (Edge e : ledge) {
			if (!res.contains(e.p1))
				res.add(e.p1);
			if (!res.contains(e.p2))
				res.add(e.p2);
		}
		return res;
	}
	// renvoie la liste d'arrete composee par ces points et qui existe deja

	private ArrayList<Edge> asupprimer(ArrayList<Edge> ledge, Point p, Point b, Point c) {
		ArrayList<Edge> res = new ArrayList<Edge>();
		for (Edge e : ledge) {
			if (e.p1.equals(p) && e.p2.equals(b))
				res.add(e);
			if (e.p1.equals(p) && e.p2.equals(c))
				res.add(e);
			if (e.p1.equals(b) && e.p2.equals(p))
				res.add(e);
			if (e.p1.equals(c) && e.p2.equals(p))
				res.add(e);
		}
		return res;
	}

	public double score(ArrayList<Edge> e) {
		double res = 0;
		for (Edge s : e)
			res += s.poids;
		return res;
	}
	private double score(Tree2D tree) {
		double res = 0;
		for (Tree2D t : tree.getSubTrees())
			res += t.getRoot().distance(tree.getRoot()) + score(t);
		return res;
	}


	private static ArrayList<Edge> treeToEdges(Tree2D tree) {
		ArrayList<Edge> edges = new ArrayList<Edge>();
		for (Tree2D subTree : tree.getSubTrees()) {
			edges.add(new Edge(tree.getRoot(),subTree.getRoot(), tree.getRoot().distance(subTree.getRoot())));
			edges.addAll(treeToEdges(subTree));
		}
		return edges;
	}
	private static Tree2D edgesToTree(ArrayList<Edge> edges, Point root) {
		ArrayList<Edge> remainder = new ArrayList<Edge>();
		ArrayList<Point> subTreeRoots = new ArrayList<Point>();
		Edge current;
		while (edges.size()!=0) {
			current = edges.remove(0);
			if (current.p1.equals(root)) {
				subTreeRoots.add(current.p2);
			} else {
				if (current.p2.equals(root)) {
					subTreeRoots.add(current.p1);
				} else {
					remainder.add(current);
				}
			}
		}

		ArrayList<Tree2D> subTrees = new ArrayList<Tree2D>();
		for (Point subTreeRoot: subTreeRoots) subTrees.add(edgesToTree((ArrayList<Edge>)remainder.clone(),subTreeRoot));

		return new Tree2D(root, subTrees);
	}
	public ArrayList<Point> multiThread(ArrayList<Point> points, int n) {

		ArrayList<FutureTask<ArrayList<Point>>> ft = new ArrayList<>();
		ArrayList<Point> res ;
		for(int i=0; i<n; i++){
			ft.add(new FutureTask<>(new Task((ArrayList<Point>)points.clone())));
		}

		ExecutorService executor = Executors.newFixedThreadPool(n);

		for(int i=0; i<n; i++){
			executor.execute(ft.get(i));
		}

		boolean fin = true;

		while(true){
			fin = true;
			for(int i=0; i<n; i++){
				if(!ft.get(i).isDone()){
					fin = false;
				}
			}

			if(fin){
				try {
					res = ft.get(0).get();
					for(int i=1; i<n; i++){
						if(ft.get(i).get().size() < res.size()){
							res = ft.get(i).get();
						}
					}
					return res;
				} catch (InterruptedException | ExecutionException e) {

					e.printStackTrace();
				}
			}

		}

	}

	public static ArrayList<Point> glouton(ArrayList<Point> points) {
		Evaluation e = new Evaluation();
		ArrayList<Point> pointsCopie = (ArrayList<Point>) points.clone();
		ArrayList<Point> result = new ArrayList<Point>();
		Point  u;

		while(!pointsCopie.isEmpty()){
			//			Collections.shuffle(pointsCopie);
			u= getMaxVoisin(pointsCopie);
			result.add(u);

			ArrayList<Point> voisins= e.neighbor(u, pointsCopie, edgeThreshold);

			pointsCopie.remove(u);
			pointsCopie.removeAll(voisins);
		}

		return result;
	}




	private static Point getMaxVoisin(ArrayList<Point> points) {
		Evaluation e = new Evaluation();
		ArrayList<Point> voisins = e.neighbor(points.get(0), points,edgeThreshold);

		int nbVoisin = voisins.size();
		//Point MaxVoisin = points.get((new Random()).nextInt(points.size()));
		Point MaxVoisin = points.get(0);
		for(int i = 1 ; i < points.size(); i++){
			voisins=e.neighbor(points.get(i), points,edgeThreshold);
			if(nbVoisin < voisins.size()){
				nbVoisin = voisins.size();
				MaxVoisin = points.get(i);
			}
		}

		return MaxVoisin;
	}
	public static ArrayList<Point> localSearchingWhile(ArrayList<Point> points, ArrayList<Point> result) {
		ArrayList<Point> resultMin=(ArrayList<Point>)result.clone();

		resultMin = ameliorer(points, result);

		int i=0;

		do{
			i++;
			result = resultMin;
			resultMin = ameliorerRuturn(points, result);
			//			System.out.println(i+" : "+resultMin.size() );

		}while( resultMin.size() < result.size());

		return resultMin;
	}

	private static ArrayList<Point> ameliorer(ArrayList<Point> points,ArrayList<Point> result) {
		ArrayList<Point> pointsCopie = (ArrayList<Point>)points.clone();
		ArrayList<Point> ensDom =  (ArrayList<Point>)result.clone();
		Evaluation e = new Evaluation();


		for(int i=0; i < ensDom.size(); i++){
			for(int j=i+1; j < ensDom.size(); j++){

				Point p = ensDom.get(i);
				Point q = ensDom.get(j);

				if (p.distance(q) > 3.5 * edgeThreshold) continue;

				//				pointsCopie.removeAll(ensDom);

				for(Point z: pointsCopie){

					if( !(p.distance(z) <= 2.5 * edgeThreshold && q.distance(z) <= 2.5 * edgeThreshold)) continue;

					ensDom.add(z);
					ensDom.remove(p) ;	
					ensDom.remove(q);

					if(e.isValide(ensDom,points,edgeThreshold)) {
						break;
					}
					else{
						ensDom.add(p);
						ensDom.add(q);
						ensDom.remove(z);
					}	
				}
				//				pointsCopie = new ArrayList<Point>(points);
			}
		}

		return ensDom;
	}

	private static ArrayList<Point> ameliorerRuturn(ArrayList<Point> points, ArrayList<Point> result) {

		ArrayList<Point> pointsCopie = (ArrayList<Point>)points.clone();
		ArrayList<Point> ensDom =  (ArrayList<Point>)result.clone();
		Evaluation e = new Evaluation();


		for(int i=0; i < ensDom.size(); i++){
			for(int j=i+1; j < ensDom.size(); j++){

				Point p = ensDom.get(i);
				Point q = ensDom.get(j);

				if (p.distance(q) > 3.5 * edgeThreshold) continue;

				//				pointsCopie.removeAll(ensDom);

				for(Point z: pointsCopie){

					if( !(p.distance(z) <= 2.5 * edgeThreshold && q.distance(z) <= 2.5 * edgeThreshold)) continue;

					ensDom.add(z);
					ensDom.remove(p) ;	
					ensDom.remove(q);

					if(e.isValide(ensDom,points,edgeThreshold)) {
						return ensDom;
					}
					else{
						ensDom.add(p);
						ensDom.add(q);
						ensDom.remove(z);
					}	
				}
				//				pointsCopie = new ArrayList<Point>(points);
			}
		}

		return ensDom;
	}

	//FILE PRINTER
	private void saveToFile(String filename,ArrayList<Point> result){
		int index=0;
		try {
			while(true){
				BufferedReader input = new BufferedReader(new InputStreamReader(new FileInputStream(filename+Integer.toString(index)+".points")));
				try {
					input.close();
				} catch (IOException e) {
					System.err.println("I/O exception: unable to close "+filename+Integer.toString(index)+".points");
				}
				index++;
			}
		} catch (FileNotFoundException e) {
			printToFile(filename+Integer.toString(index)+".points",result);
		}
	}
	private void printToFile(String filename,ArrayList<Point> points){
		try {
			PrintStream output = new PrintStream(new FileOutputStream(filename));
			int x,y;
			for (Point p:points) output.println(Integer.toString((int)p.getX())+" "+Integer.toString((int)p.getY()));
			output.close();
		} catch (FileNotFoundException e) {
			System.err.println("I/O exception: unable to create "+filename);
		}
	}

	//FILE LOADER
	public static ArrayList<Point> readFromFile(String filename) {
		String line;
		String[] coordinates;
		ArrayList<Point> points=new ArrayList<Point>();
		try {
			BufferedReader input = new BufferedReader(
					new InputStreamReader(new FileInputStream(filename))
					);
			try {
				while ((line=input.readLine())!=null) {
					coordinates=line.split("\\s+");
					points.add(new Point(Integer.parseInt(coordinates[0]),
							Integer.parseInt(coordinates[1])));
				}
			} catch (IOException e) {
				System.err.println("Exception: interrupted I/O.");
			} finally {
				try {
					input.close();
				} catch (IOException e) {
					System.err.println("I/O exception: unable to close "+filename);
				}
			}
		} catch (FileNotFoundException e) {
			System.err.println("Input file not found.");
		}
		return points;
	}
}
