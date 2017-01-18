import java.io.File;
import java.util.ArrayList;

public class Naoki {
	static int boundTheshold = 0;
	static String folders = "/home/max/Matlab/New/Naoki/cc";// connected component folder

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		File folder = new File(folders);
		File[] listOfFiles = folder.listFiles();

		for (int i = 0; i < listOfFiles.length; i++) {
			System.out.println(listOfFiles[i]);

			ArrayList<CC> ccs = new ArrayList<CC>();// to store all connected components
			
			In in = new In(listOfFiles[i]);
			Out points = new Out("./Points/"+ listOfFiles[i].getName().substring(0,listOfFiles[i].getName().length() - 3) + ".point"+ boundTheshold);
			Out edges = new Out("./Edges/"+listOfFiles[i].getName().substring(0,listOfFiles[i].getName().length() - 3) + ".edge"+ boundTheshold);
			Out boundariesSet = new Out("./Boundaries/"	+ listOfFiles[i].getName().substring(0,listOfFiles[i].getName().length() - 3) + ".bounSet"+ boundTheshold);

			ccs = SingleImageProcess(in, points, edges);
			System.out.println("There are "+ccs.size()+" boundaries in "+ listOfFiles[i].getName().substring(0,listOfFiles[i].getName().length() - 3));
			BoundarySets(ccs, boundariesSet);
		}
	}

	/**
	 * @param in: image name
	 * @param pointsThes: image threshold points output file name
	 * @param edgesThes: edges output file name
	 * 
	 * @return cc: all boundaries of an image, a boundary is an ArrayList of Edge
	 */
	public static ArrayList<CC> SingleImageProcess(In in, Out pointsThes,Out edgesThes) {
		ArrayList<CC> ccs = new ArrayList<CC>();
		int numberOfCC = Integer.parseInt(in.readLine());
		int rowOfBW = Integer.parseInt(in.readLine());
		int colOfBW = Integer.parseInt(in.readLine());
		String temp = null;
		int line = 3;

		while ((temp = in.readLine()) != null) {// every line is a connected component
			ArrayList<Edge> edges = new ArrayList<Edge>();
			line++;

			int length = temp.split("\t").length;
			ArrayList<PointCor> pc = new ArrayList<PointCor>();
			ArrayList<PointCor> pcTemp = new ArrayList<PointCor>();
			for (int j = 0; j < length; j++) {
				long v = Long.parseLong(temp.split("\t")[j]);
				long y = v / rowOfBW + 1;// column
				long x = v % rowOfBW;// line
				pc.add(new PointCor(v, x, y));
				pcTemp.add(new PointCor(v, x, y));
			}

			int eindex = 0;
			int window = 0;
			int edirection = 0;

			while (pc.size() > 1) {// while process a line
				if (pc.get(0).x == pc.get(1).x)	edirection = 0;// vertical
				if (pc.get(0).y == pc.get(1).y)	edirection = 1;// horizontal
				if (Math.abs(pc.get(0).y - pc.get(1).y) == 1 && Math.abs(pc.get(0).x - pc.get(1).x) == 1) edirection = 2;// inclined

				if (edirection == 0) {
					for (window = 1; window < pc.size(); window++) {
						if (pc.get(0).x != pc.get(window).x) break;// now, window equals to the first location of a new edge
					}
				}
				if (edirection == 1) {
					for (window = 1; window < pc.size(); window++) {
						if (pc.get(0).y != pc.get(window).y) break;
					}
				}
				if (edirection == 2) {
					for (window = 1; window < pc.size(); window++) {
						if (Math.abs(pc.get(window - 1).y - pc.get(window).y) != 1	|| Math.abs(pc.get(window - 1).x - pc.get(window).x) != 1)	break;
					}
				}

				// add this processed edge to variable edges
				eindex++;
				if (edirection == 2) {
					edges.add(new Edge(line - 3, eindex, edirection, window * 1.414));// each element in edges is an edge of a component
				} else {
					edges.add(new Edge(line - 3, eindex, edirection, window));// each element in edges in an edge of a component
				}

				// remove this processed edge
				for (int j = 0; j < window; j++) {
					pc.remove(0);
				}
			}
			if (pc.size() == 1) {
				edges.add(new Edge(line - 3, eindex, 2, 1));// each element in edges is an edge f a component
				pc.remove(0);
			}

			CC tempCCThes = new CC(edges);
			if (tempCCThes.length >= boundTheshold) {
				for (int j = 0; j < pcTemp.size(); j++) {
					pointsThes.println(tempCCThes.cindex + "\t"	+ pcTemp.get(j).value + " " + "\t"+ pcTemp.get(j).x + "\t" + pcTemp.get(j).y);
				}
				for (int j = 0; j < edges.size(); j++) {
					edgesThes.println(edges.get(j).cindex + " "	+ edges.get(j).index + " " + edges.get(j).direction	+ " " + edges.get(j).length);
				}
				ccs.add(tempCCThes);// each element in cc is a component in an image
			}
		}
		return ccs;
	}
	
	/**
	 * 
	 * @param ccs: all boundaries of an image, a boundary is an ArrayList of Edge
	 * @param boundariesSetThes: output file name
	 * @param BS: boundaries series
	 * 
	 * function: change cc into boundaries series and write them into file
	 */
	public static void BoundarySets(ArrayList<CC> ccs, Out boundariesSetThes) {
		ArrayList<Boundary> bs = new ArrayList<Boundary> ();
		for (int i = 0; i < ccs.size(); i++) {
			ArrayList<Double> bound = new ArrayList<Double>();
			for (int j = 0; j < ccs.get(i).edges.size(); j++) {
				bound.add(ccs.get(i).edges.get(j).length - ccs.get(i).edges.get((j + 1) % ccs.get(i).edges.size()).length);
			}
			Boundary b = new Boundary(bound);
			bs.add(b);
		}
		// System.out.println("bs size: "+BS.bs.size());
		for (int i = 0; i < bs.size(); i++) {
			boundariesSetThes.print(i + 1 + "\t");
			for (int j = 0; j < bs.get(i).bound.size(); j++) {
				boundariesSetThes.print(bs.get(i).bound.get(j) + "\t");
			}
			boundariesSetThes.println();
		}
	}

}
