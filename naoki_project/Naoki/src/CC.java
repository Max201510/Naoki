import java.util.ArrayList;

public class CC {
	int cindex;
	ArrayList<Edge> edges;
	long length;

	public CC(ArrayList<Edge> edges){
		this.cindex=edges.get(0).cindex;
		this.edges=edges;

		for(int i=0;i<edges.size();i++){
			this.length+=edges.get(i).length;
		}
	}
}
// a CC is an edge arraylist