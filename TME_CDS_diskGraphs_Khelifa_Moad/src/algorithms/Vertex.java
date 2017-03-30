package algorithms;

import java.awt.Point;

public class Vertex extends Point {
	
	enum Color{WHITE, GREY, BLACK, BLUE};
	
	Color color;
	int  id;
	boolean actif=false;
	int whiteCount;
	
	public Vertex(Point arg0, int id) {
	super(arg0);
	color=Color.WHITE;
	this.id=id;
	whiteCount=0;
	
	}
	
	public void setActif(boolean b){
		actif=b;
	}

	public boolean isActif() {
		return actif;
	}
	
	
}
