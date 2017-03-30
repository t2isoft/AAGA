package algorithms;

import java.awt.Point;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.Exchanger;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.FutureTask;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.lang.reflect.Array;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;

public class DefaultTeam {
	public static  int edgeThreshold;
	public static double epsilone = 0.01;

	public static double taux =  epsilone +1;

	public ArrayList<Point> calculDominatingSet(ArrayList<Point> points, int edgeThreshold) {
		this.edgeThreshold = edgeThreshold;


		//		ArrayList<Point> ensDomMin = localSearching2(points,decoupeMult(points, 5));
		//		saveToFile("EnsDom_"+ensDomMin.size(), ensDomMin);
		//		return ensDomMin;
		//		

		//return localSearching2(points,decoupeMult(points, 5));

		
		//return decoupe(points);

		//return glouton(points);

		//return localSearching2(points, decoupeArticle(points));

		//return decoupeArticle(points);


		//return localSearchingWhile(points, result);
		
		
		
/**
		ArrayList<Point> result = new ArrayList<Point>();
		result = glouton(points);
		
		System.out.println("Score glouton= "+result.size());
		
		ArrayList<Point> resultFinal = new ArrayList<Point>(result);

		for(int i=0; i< 100; i++){
			
			Collections.shuffle(points);
			Collections.shuffle(result);
			ArrayList<Point> result2 = localSearchingWhile(points,result); 
			result = glouton(points);
			
			if(result2.size() < resultFinal.size()){
				resultFinal=result2;
				//result=resultFinal;
			}
			
			
			//System.out.println(i+" : "+result.size());
			System.out.println(i+ " : result2.size()=" +result2.size()+" resultFinal.size()="+resultFinal.size());

			
		}

		System.out.println(""+resultFinal.size());
		return resultFinal;
**/
		
		return multiThread(points, 8);
		
		//return result;

		//return localSearching(points);
		//return localSearching3(points,localSearching(points));

		//return readFromFile("scoreTMEEnsDom4.points");

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

	private ArrayList<Point> decoupeArticle(ArrayList<Point> points) {

		ArrayList<Point> pointsCopie = new ArrayList<Point>(points);
		Point v;
		int r=0, k=1;

		//HashMap<String, Set<Point>> Ni = new HashMap<String, Set<Point>>();
		//HashMap<String, ArrayList<Point>> Ni = new HashMap<String, ArrayList<Point>>();
		HashMap<String, ArrayList<Point>> Nk = new HashMap<String, ArrayList<Point>>();
		//Set<Point> Nr,Nr2;
		ArrayList<Point> Nr,Nr2;
		ArrayList<Point> D1, D2;
		do{
			v = pointsCopie.get(0);
			do{
				//System.out.println("AVANT Nr");
				Nr = voisinageEti(r, v, pointsCopie);
				//System.out.println("Nr="+Nr);

				Nr2 = voisinageEti(r+2, v, pointsCopie);

				//System.out.println("Nr2="+Nr2);


				//				Ni.put("N"+r, Nr);
				//				Ni.put("N"+r+2, Nr2);

				//				D1=glouton(new ArrayList<Point>(Nr));
				//				D2=glouton(new ArrayList<Point>(Nr2));

				//								D1=localSearching(new ArrayList<Point>(Nr));
				//								D2=localSearching(new ArrayList<Point>(Nr2));

				D1=localSearching2(new ArrayList<Point>(Nr),localSearching(new ArrayList<Point>(Nr)));
				D2=localSearching2(new ArrayList<Point>(Nr2),localSearching(new ArrayList<Point>(Nr2)));

				//				D1=ensDomMin(new ArrayList<Point>(Nr));
				//				D2=ensDomMin(new ArrayList<Point>(Nr2));

				//				
				r++;
				System.out.println(pointsCopie.size() +" r=" +r);
			}while(D2.size() > taux * D1.size());

			Nk.put("N"+k, D2);
			pointsCopie.removeAll(Nr2);

			k++;
		}while( ! pointsCopie.isEmpty());

		ArrayList<Point> ensDom = new ArrayList<Point>();

		for(String key : Nk.keySet())
			ensDom.addAll(Nk.get(key));

		return ensDom;
	}

	private ArrayList<Point> ensDomMin(ArrayList<Point> arrayList) {



		return null;
	}

	private Set<Point> voisinage1(int r, Point v, ArrayList<Point> pointsCopie,
			HashSet<Point> hashSet) {
		hashSet.add(v);
		if(r==0) return hashSet;

		for(Point p: pointsCopie){
			if(v.distance(p)<=r*edgeThreshold)
				hashSet.add(p);
		}
		return hashSet;
	}

	private ArrayList<Point> decoupeMult(ArrayList<Point> points, int nbParts){

		int xMax=Integer.MIN_VALUE, yMax=Integer.MIN_VALUE, xMin=Integer.MAX_VALUE, yMin=Integer.MAX_VALUE;

		for(Point p: points){
			if(p.getX() > xMax) xMax=(int)p.getX();
			if(p.getX() < xMin) xMin=(int)p.getX();
			if(p.getY() > yMax) yMax=(int)p.getY();
			if(p.getY() < yMin) yMin=(int)p.getY();
		}

		int pasx  = (int) ((xMax-xMin)/Math.sqrt(nbParts)); 
		int pasy = (int) ((yMax-yMin)/Math.sqrt(nbParts));

		Map<String, ArrayList<Point>> partitions =new  HashMap<String,ArrayList<Point>>();


		for(int i=xMin; i<xMax; i+=pasx){
			for(int j=yMin; j<yMax; j+=pasy){
				ArrayList<Point> set = partitions.get(i+" "+j);
				for(Point p: points){
					if(i<=p.getX() && p.getX()< i+pasx && j<=p.getY() && p.getY()<j+pasy){
						if(set!=null && !set.isEmpty()){
							set.add(p);
						}else{
							set=new ArrayList<Point>();
							set.add(p);
						}		
					}
				}
				if(set!=null && !set.isEmpty())
					partitions.put(i+" "+j, localSearching(set));
			}	
		}

		//System.out.println(partitions);
		//	xMax=0;
		//	for(String k:partitions.keySet()){
		//		
		//		if(partitions.get(k) != null)
		//			xMax+=partitions.get(k).size();
		//		else
		//			partitions.remove(k);
		//		System.out.println(k + " "+ partitions.get(k));
		//		
		//	}
		//	System.out.println(xMax);

		Set<Point> ensDom = new HashSet<Point>();
		for(String k:partitions.keySet()){
			if(partitions.get(k) == null) continue;
			ensDom.addAll(localSearching(partitions.get(k)));
			//ensDom.addAll(partitions.get(k));
		}

		return new ArrayList<Point>(ensDom);
	}
	private ArrayList<Point> decoupe(ArrayList<Point> points) {

		ArrayList<Point> s1 = new ArrayList<>();
		ArrayList<Point> s2 = new ArrayList<>();
		ArrayList<Point> s3 = new ArrayList<>();
		ArrayList<Point> s4 = new ArrayList<>();


		ArrayList<Point> res1 = new ArrayList<>();
		ArrayList<Point> res2 = new ArrayList<>();
		ArrayList<Point> res3 = new ArrayList<>();
		ArrayList<Point> res4 = new ArrayList<>();

		Set<Point> result = new HashSet<Point>();

		double xmoyen = 100;
		double ymoyen = 100;

		for(Point p: points){
			xmoyen+=p.getX();
			ymoyen+=p.getY();
		}
		xmoyen=xmoyen/points.size();
		ymoyen=ymoyen/points.size();

		for(Point p: points ){
			if(p.getX()<=xmoyen && p.getY()>=ymoyen) s1.add(p);
			else if(p.getX()>xmoyen && p.getY()>=ymoyen) s2.add(p);
			else if(p.getX()<=xmoyen && p.getY()<ymoyen) s3.add(p);
			else if(p.getX()>xmoyen && p.getY()<ymoyen) s4.add(p);
		}


		res1=localSearchingWhile(s1,glouton(s1));
		res2=localSearchingWhile(s2,glouton(s2));
		res3=localSearchingWhile(s3,glouton(s3));
		res4=localSearchingWhile(s4,glouton(s4));

		result.addAll(localSearching2(s1, res1));
		result.addAll(localSearching2(s2,res2));
		result.addAll(localSearching2(s3,res3));
		result.addAll(localSearching2(s4,res4));

		return localSearching2(points, new ArrayList<Point>(result));
		//return ameliorer(points, new ArrayList<>(result));
	}


	private ArrayList<Point> localSearching2(ArrayList<Point> points, ArrayList<Point> ensDom){
		ArrayList<Point> result = new ArrayList<Point>();

		System.out.println("localSearching2");

		result = ensDom;//glouton(points);
		//result = glouton(points);
		int currentScore = result.size();
		//System.out.println("Score glouton = "+currentScore);
		int i=1;


		result = ameliorer(points,result);

		int newScore = result.size();
		//System.out.println("Debut score = "+result.size());

		while( newScore < currentScore){
			result = ameliorer(points, result);
			currentScore = newScore;
			newScore = result.size();

			System.out.println(i+" : "+result.size());
			//saveToFile("scoreTMEEnsDomDecoupe", result);
			i++;
		}

		return result;

	}

	private ArrayList<Point> localSearching(ArrayList<Point> points) {

		ArrayList<Point> result = new ArrayList<Point>();


		result = glouton(points);
		int currentScore = result.size();
		System.out.println("Score glouton = "+currentScore);
		int i=1;


		result = ameliorer(points,result);
		int newScore = result.size();

		System.out.println("Score debut = "+result.size());


		while( newScore < currentScore){

			result = ameliorer(points, result);

			currentScore = newScore;
			newScore = result.size();

			System.out.println(i+" : "+result.size());

			i++;			

		}

		return result;
	}


	private ArrayList<Point> localSearching3(ArrayList<Point> points, ArrayList<Point> ensDom) {

		ArrayList<Point> result = new ArrayList<Point>();


		result = ensDom;
		int currentScore = result.size();
		System.out.println("Score glouton = "+currentScore);
		int i=1;


		result = ameliorer3(points,result);
		int newScore = result.size();

		System.out.println("Score debut = "+result.size());


		while( newScore < currentScore){

			result = ameliorer3(points, result);

			currentScore = newScore;
			newScore = result.size();

			System.out.println(i+" : "+result.size());

			i++;			

		}

		return result;
	}


	private ArrayList<Point> ameliorerOld(ArrayList<Point> points,ArrayList<Point> enDom) {
		Evaluation e = new Evaluation();
		ArrayList<Point> ensDomCopie = new ArrayList<Point>(enDom);
		ArrayList<Point> ensDomCopie1 = new ArrayList<Point>(enDom);

		for(Point p : ensDomCopie1){
			for(Point q : ensDomCopie1){
				if(p.equals(q)) continue;

				if (p.distance(q)>3*edgeThreshold)
					continue;

				for(Point z: points){

					if( !(p.distance(z) <= edgeThreshold && q.distance(z) <= edgeThreshold))
						continue;

					//Point z = points.get((new Random()).nextInt(points.size()));
					if(p.equals(z) || q.equals(z)) continue;

					if( ensDomCopie.remove(p) && ensDomCopie.remove(q))
						ensDomCopie.add(z);

					if(e.isValide(ensDomCopie,points,edgeThreshold)) {
						if( enDom.remove(p) && enDom.remove(q))
							enDom.add(z);	
						return enDom;
					}
					else{
						ensDomCopie.add(p);
						ensDomCopie.add(q);
						ensDomCopie.remove(z);
					}
					ensDomCopie = new ArrayList<Point>(enDom);
					ensDomCopie1= new ArrayList<Point>(enDom);
					//System.out.println(fvs.size());
				}
			}
		}

		return enDom;
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



	private ArrayList<Point> ameliorer3(ArrayList<Point> points,ArrayList<Point> enDom) {
		//ArrayList<Point> pointsCopie = new ArrayList<Point>(points);
		Evaluation e = new Evaluation();
		ArrayList<Point> ensDomCopie = new ArrayList<Point>(enDom);
		ArrayList<Point> ensDomCopie1 = new ArrayList<Point>(enDom);
		//ArrayList<Point> enDom = new ArrayList<Point>(enDomGlouton);

		//	Collections.shuffle(ensDomCopie1);

		for(int i=0; i < ensDomCopie1.size(); i++){
			for(int j=i+1; j < ensDomCopie1.size(); j++){
				for(int k=j+1; k < ensDomCopie1.size(); k++){
					//		
					//		for(Point p:ensDomCopie1){
					//			for(Point q:ensDomCopie1){

					Point p = ensDomCopie1.get(i);
					Point q = ensDomCopie1.get(j);
					Point r = ensDomCopie1.get(k);
					//if(p.equals(q)) continue;

					if (p.distance(q) > 3*edgeThreshold) continue;

					//pointsCopie.removeAll(enDom);

					for(Point z: points){
						for(Point x: points){
							if( !(p.distance(z) <= 4*edgeThreshold && q.distance(z) <= 4*edgeThreshold)) continue;

							//Point z = points.get((new Random()).nextInt(points.size()));
							//if(p.equals(z) || q.equals(z)) continue;

							//if( ensDomCopie.remove(p) && ensDomCopie.remove(q))
							ensDomCopie.add(z);
							ensDomCopie.add(x);
							ensDomCopie.remove(p) ;	
							ensDomCopie.remove(q);
							ensDomCopie.remove(r);


							if(e.isValide(ensDomCopie,points,edgeThreshold)) {
								//if( enDom.remove(p) && enDom.remove(q)){
								enDom.remove(p) ;
								enDom.remove(q);
								enDom.remove(r);
								enDom.add(z);
								enDom.add(x);


								//System.out.println("remplacer "+enDom.size());
								return enDom;
								//break;
								//}
							}
							else{
								ensDomCopie.add(p);
								ensDomCopie.add(q);
								ensDomCopie.add(r);
								ensDomCopie.remove(z);
								ensDomCopie.remove(x);
							}
							ensDomCopie = new ArrayList<Point>(enDom);
							ensDomCopie1= new ArrayList<Point>(enDom);	
						}
					}
				}
			}
		}

		return enDom;
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
	private ArrayList<Point> readFromFile(String filename) {
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

	private Set<Point> voisinage(int r,Point p, ArrayList<Point >points, Set<Point> voisinageCurr){
		Evaluation e = new Evaluation();
		//System.out.println("voisinage r="+r);
		if(r==0){
			voisinageCurr.add(p);
			return voisinageCurr;
		}
		else{
			ArrayList<Point> vois = e.neighbor(p, points, edgeThreshold);
			//if(vois.isEmpty()){voisinageCurr.add(p); return voisinageCurr;}

			for(Point u: vois){
				voisinageCurr.addAll(voisinage(r-1, u, points, voisinageCurr));
			}
			voisinageCurr.add(p);
			return voisinageCurr;
		}
	}

	private ArrayList<Point> voisinageEti(int r,Point p, ArrayList<Point >points){
		Evaluation e = new Evaluation();
		ArrayList<Point> voisins = new ArrayList<Point>();
		int indx=1, indx_1=0, nbVoisin=0;


		voisins.add(p);

		while(r>0){

			nbVoisin=0;

			for(int i=indx_1; i<indx; i++ ){

				for(Point u: e.neighbor(voisins.get(i), points, edgeThreshold)){
					if(!voisins.contains(u)){
						voisins.add(u);
						nbVoisin++;
					}
				}	
			}
			indx_1=indx;
			indx+=nbVoisin;
			r--;
		}




		return voisins;
	}
}
