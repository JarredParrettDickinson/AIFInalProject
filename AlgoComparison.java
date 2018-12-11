
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.time.Duration;
import java.time.Instant;
import java.util.*;
import java.awt.*;

public class AlgoComparison {

	public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
		PrintWriter writer = new PrintWriter("PATest.txt", "UTF-8");
		writer.println("outBFS,outBi,Dij");
		String location = "/Users/jarredparrett/Desktop/USA-road-d.NY.gr";
		Scanner sc = new Scanner(new File(location));
		String line0 = sc.nextLine();// c 9th DIMACS Implementation Challenge: Shortest Paths
		String line1 = sc.nextLine();// c 9th DIMACS Implementation Challenge: Shortest Paths
		String line2 = sc.nextLine();// c TIGER/Line graph USA-road-d.NY
		String line3 = sc.nextLine();// c
		String line4 = sc.nextLine();// p sp 264346 733846
		String line5 = sc.nextLine();// c graph contains 264346 nodes and 733846 arcs
		String line6 = sc.nextLine();// c
		String[] line4Split = line4.split(" ");

		WeightedGraph.Graph searchSpace = new WeightedGraph.Graph(Integer.parseInt(line4Split[2]) + 1);
		int ii = 0;
		while (sc.hasNextLine()) {
			if (ii++ % 10000 == 0) {
				// System.out.println(ii);
			}
			String line = sc.nextLine();
			String[] line_S = line.split(" ");
			int source = Integer.parseInt(line_S[1]);
			int dest = Integer.parseInt(line_S[2]);
			int val = Integer.parseInt(line_S[3]);
			searchSpace.addEdge(source, dest, val);
		}
		String location_lat = "/Users/jarredparrett/Desktop/USA-road-d.NY.co";
		Scanner sc2 = new Scanner(new File(location_lat));
		String line_0 = sc2.nextLine();// c 9th DIMACS Implementation Challenge: Shortest Paths
		String line_1 = sc2.nextLine();// c 9th DIMACS Implementation Challenge: Shortest Paths
		String line_2 = sc2.nextLine();// c TIGER/Line nodes coords for graph USA-road-d.NY
		String line_3 = sc2.nextLine();// c
		String line_4 = sc2.nextLine();// aux sp co 264346
		String line_5 = sc2.nextLine();// c graph contains 264346 nodes
		String line_6 = sc2.nextLine();// c
		String[] line_4Split = line_4.split(" ");
		int arrLen = Integer.parseInt(line4Split[2]);
		Point[] pointLocations = new Point[arrLen + 1];
		while (sc2.hasNextLine()) {

			if (ii++ % 10000 == 0) {
				// System.out.println(ii);
			}
			String cur = sc2.nextLine();
			String[] curSplit = cur.split(" ");
			int point = Integer.parseInt(curSplit[1]);
			int x = Integer.parseInt(curSplit[2]);
			int y = Integer.parseInt(curSplit[3]);
			pointLocations[point] = new Point(x, y);
		}
		sc.close();
		sc2.close();

		Random rand = new Random();
		int start = rand.nextInt((arrLen) + 1) + 0; // select start
		System.out.println("Start " + start);
		for (int i = 0; i < 50; i++) {
			int end = rand.nextInt((arrLen) + 1) + 0;
			double[] dij = dijkstraGraphSearch(start, end, searchSpace);
			System.out.println("Dij Path Len: " + dij[0] + " Dij Time: " + dij[1]);
			// System.out.println("End " + end);
			double[] outBFS = breathFirstSearch(start, end, searchSpace);
			System.out.println("BFS Path Len: " + outBFS[0] + " BFS Time: " + outBFS[1]);
			double[] outBi = biDirectionalSearch(start, end, searchSpace);
			System.out.println("Bi-Di Path Len: " + outBi[0] + " Bi-Di Time: " + outBi[1]);
			// + outBi[0]
			// + " Out Bi 1 " + outBi[1]);
			// long[] outA = aStarGraphSearch(start, end, searchSpace, pointLocations);
			// System.out.println("Astar 0 " + outA[0] + " Astar 1 " + outA[1] + " Astar 2 "
			// + outA[2]);
			// long[] aStarBi = aStarGraphSearchBi(start, end, searchSpace, pointLocations);
			// System.out.println("AstarBi 0 " + aStarBi[0] + " AstarBi 1 " + aStarBi[1]);
			// long[] aStarBiOther = aStarGraphSearchBiApproach(start,end, searchSpace,
			// pointLocations);
			// System.out.println("AstarBiOther 0 " + aStarBiOther[0] + " AstarBiOther 1 " +
			// aStarBiOther[1]);

		}

		writer.close();
	}

	/**
	 * Breath First Graph Search
	 * 
	 * @param start
	 * @param destination
	 * @param searchSpace
	 * @return
	 */
	static double[] breathFirstSearch(int start, int destination, WeightedGraph.Graph searchSpace) {
		Instant begin = Instant.now();
		long[] visitedNum = new long[searchSpace.verticies + 1];
		for (int i = 0; i < visitedNum.length; i++) {
			visitedNum[i] = Integer.MAX_VALUE;
		}
		// HashSet<Integer> visitedCheck = new HashSet<Integer>();
		Vertex star = new Vertex(start, 0, null);
		// long expanded = 0;
		// long explored = 1;
		double[] output = new double[10];
		// boolean visited[] = new boolean[searchSpace.verticies + 1];
		Queue<Vertex> q = new LinkedList<Vertex>();
		// visited[start] = true;
		q.add(star);
		visitedNum[start] = 0;
		while (q.size() != 0) {
			star = q.poll();
			// expanded++;
			Iterator<WeightedGraph.Edge> i = searchSpace.adjList[star.node].listIterator();
			while (i.hasNext()) {
				WeightedGraph.Edge next = i.next();
				if (visitedNum[next.destination] != Integer.MAX_VALUE
						|| (next.weight + star.distance) > visitedNum[next.destination]) {
					continue;
				}
				// explored++;
				int n = next.destination;
				visitedNum[n] = next.weight + star.distance;
				Vertex nVert = new Vertex(n, (star.distance + next.weight), star);
				q.add(nVert);
				// visitedCheck.add(nVert);
				// visitedNum[n] = nVert.distance;
				if (visitedNum[destination] != Integer.MAX_VALUE) {
					output[0] = visitedNum[destination];
					// output[1] = explored;
					Instant finish = Instant.now();
					long timeElapsed = Duration.between(begin, finish).toMillis();
					output[1] = timeElapsed;
					return output;
				}

			}
		}
		return output;
	}

	/**
	 * Bi-directional Graph Search
	 * 
	 * @param start
	 * @param destination
	 * @param searchSpace
	 * @return
	 */
	static double[] biDirectionalSearch(int start, int destination, WeightedGraph.Graph searchSpace) {
		Instant begin = Instant.now();
		// HashSet<Integer> visitedStartCheck = new HashSet<Integer>();
		// HashSet<Integer> visitedDestCheck = new HashSet<Integer>();
		// long expanded = 0;
		// long explored = 1;
		long[] visitedNumStart = new long[searchSpace.verticies + 1];
		for (int i = 0; i < visitedNumStart.length; i++) {
			visitedNumStart[i] = Integer.MAX_VALUE;
		}
		long[] visitedNumDest = new long[searchSpace.verticies + 1];
		for (int i = 0; i < visitedNumDest.length; i++) {
			visitedNumDest[i] = Integer.MAX_VALUE;
		}
		visitedNumStart[start] = 0;
		visitedNumDest[destination] = 0;
		double[] output = new double[10];
		// boolean visitedStart[] = new boolean[searchSpace.verticies + 1];
		// boolean visitedDest[] = new boolean[searchSpace.verticies + 1];
		Queue<Vertex> startQ = new LinkedList<Vertex>();
		Queue<Vertex> destQ = new LinkedList<Vertex>();
		Vertex star = new Vertex(start, 0, null);
		Vertex dest = new Vertex(destination, 0, null);

		startQ.add(star);
		destQ.add(dest);
		// visitedStart[start] = true;
		// visitedDest[destination] = true;

		while (startQ.size() != 0 && destQ.size() != 0) {
			// start BFS
			Vertex orig = startQ.poll();
			// expanded++;
			Iterator<WeightedGraph.Edge> i = searchSpace.adjList[orig.node].listIterator();
			while (i.hasNext()) {
				WeightedGraph.Edge next = i.next();
				if (visitedNumStart[next.destination] != Integer.MAX_VALUE
						|| (next.weight + star.distance) > visitedNumStart[next.destination]) {
					continue;
				}
				// explored++;
				int n = next.destination;
				visitedNumStart[n] = next.weight + star.distance;
				Vertex nVert = new Vertex(n, (star.distance + next.weight), star);
				startQ.add(nVert);
			}
			// dest BFS
			Vertex dest1 = destQ.poll();
			//expanded++;
			Iterator<WeightedGraph.Edge> j = searchSpace.adjList[dest1.node].listIterator();
			while (j.hasNext()) {
				WeightedGraph.Edge next = j.next();
				if (visitedNumDest[next.destination] != Integer.MAX_VALUE
						|| (next.weight + dest1.distance) > visitedNumDest[next.destination]) {
					continue;
				}
				//explored++;
				int n = next.destination;
				visitedNumDest[n] = next.weight + dest1.distance;
				Vertex nVert = new Vertex(n, (dest1.distance + next.weight), dest1);
				destQ.add(nVert);
			}
			
			
			for (int k = 0; k < visitedNumDest.length; k++) {
				if (visitedNumStart[k] < Integer.MAX_VALUE && visitedNumDest[k] < Integer.MAX_VALUE) {
					Instant finish = Instant.now();
					long timeElapsed = Duration.between(begin, finish).toMillis();
					output[0] = (visitedNumStart[k] + visitedNumDest[k]);
					output[1] = timeElapsed;
					return output;
				}

			}

		}
		return null;
	}

	/**
	 * dijkstraGraphSearch
	 * 
	 * @param start
	 * @param destination
	 * @param searchSpace
	 * @return
	 */
	static double[] dijkstraGraphSearch(int start, int destination, WeightedGraph.Graph searchSpace) {
		Instant begin = Instant.now();
		HashSet<Integer> visitedCheck = new HashSet<Integer>();
		double[] output = new double[10];
		long minEdge[] = new long[searchSpace.verticies];
		for (int i = 0; i < minEdge.length; i++) {
			minEdge[i] = Integer.MAX_VALUE;
		}
		PriorityQueue<Vertex> pq = new PriorityQueue<Vertex>();
		minEdge[start] = 0;
		Vertex vertStart = new Vertex(start, 0, null);
		pq.add(vertStart);
		while (pq.size() != 0) {
			Vertex cur = pq.poll();
			if (visitedCheck.contains(cur.node)) {
				continue;
			}
			visitedCheck.add(cur.node);
			minEdge[cur.node] = cur.distance;

			if (cur.node == destination) {
				Instant finish = Instant.now();
				long timeElapsed = Duration.between(begin, finish).toMillis();
				output[0] = cur.distance;
				output[1] = timeElapsed;
				return output;
			}
			for (WeightedGraph.Edge curEdge : searchSpace.adjList[cur.node]) {

				long pathDist = curEdge.weight + minEdge[curEdge.source];

				pq.add(new Vertex(curEdge.destination, pathDist, cur));
			}
		}
		return output;
	}

	/**
	 * aStarGraphSearch
	 */
	static long[] aStarGraphSearch(int start, int destination, WeightedGraph.Graph searchSpace,
			Point[] pointLocations) {
		HashSet<Integer> visited = new HashSet<Integer>();
		long expanded = 0;
		long explored = 1;
		long minEdge[] = new long[searchSpace.verticies + 1];
		for (int i = 0; i < minEdge.length; i++) {
			minEdge[i] = Integer.MAX_VALUE;
		}
		minEdge[start] = 0;
		long[] output = new long[10];
		PriorityQueue<VertexA> pq = new PriorityQueue<VertexA>();
		double distance = longLatDist(pointLocations[start], pointLocations[destination], "N");
		VertexA sVert = new VertexA(start, distance, null);
		pq.add(sVert);
		while (pq.size() != 0) {
			VertexA cur = pq.poll();
			expanded++;
			Iterator<WeightedGraph.Edge> i = searchSpace.adjList[cur.node].listIterator();
			while (i.hasNext()) {
				WeightedGraph.Edge curEdge = i.next();
				if (visited.contains(curEdge.destination)) {
					continue;
				}
				explored++;
				double newDist = longLatDist(pointLocations[curEdge.destination], pointLocations[destination], "N");
				if (newDist == Double.NaN) {
					newDist = 10;// resolve this
				}
				long pathDist = curEdge.weight + minEdge[curEdge.source];
				minEdge[curEdge.destination] = pathDist;
				double value = newDist + pathDist;
				pq.add(new VertexA(curEdge.destination, value, cur));
				visited.add(curEdge.destination);

			}
			if (minEdge[destination] != Integer.MAX_VALUE) {
				output[0] = expanded;
				output[1] = explored;
				output[2] = minEdge[destination];
				return output;
			}
		}
		return null;
	}

	/**
	 * aStarGraphSearch
	 */
	static long[] aStarGraphSearchBi(int start, int destination, WeightedGraph.Graph searchSpace,
			Point[] pointLocations) {
		HashSet<Integer> visitedStartCheck = new HashSet<Integer>();
		HashSet<Integer> visitedDestCheck = new HashSet<Integer>();

		long expanded = 0;
		long explored = 1;
		long minEdgeFromStart[] = new long[searchSpace.verticies + 1];
		for (int i = 0; i < minEdgeFromStart.length; i++) {
			minEdgeFromStart[i] = Integer.MAX_VALUE;
		}
		long minEdgeFromDest[] = new long[searchSpace.verticies + 1];
		for (int i = 0; i < minEdgeFromDest.length; i++) {
			minEdgeFromDest[i] = Integer.MAX_VALUE;
		}
		boolean visitedStart[] = new boolean[searchSpace.verticies + 1];
		boolean visitedDest[] = new boolean[searchSpace.verticies + 1];
		PriorityQueue<VertexA> pqStart = new PriorityQueue<VertexA>();
		PriorityQueue<VertexA> pqDest = new PriorityQueue<VertexA>();
		long[] output = new long[10];
		double distance = longLatDist(pointLocations[start], pointLocations[destination], "N");

		VertexA sVert = new VertexA(start, distance, null);
		VertexA dVert = new VertexA(destination, distance, null);
		visitedStartCheck.add(start);
		visitedDestCheck.add(destination);

		pqStart.add(sVert);
		pqDest.add(dVert);

		visitedStart[start] = true;
		visitedDest[destination] = true;

		while (pqStart.size() != 0 && pqDest.size() != 0) {
			// start aStarGraph
			expanded++;
			VertexA curStart = pqStart.poll();
			Iterator<WeightedGraph.Edge> i = searchSpace.adjList[curStart.node].listIterator();
			while (i.hasNext()) {
				WeightedGraph.Edge curEdge = i.next();
				if (visitedStartCheck.contains(curEdge.destination)) {
					continue;
				}
				explored++;
				double newDist = longLatDist(pointLocations[curEdge.destination], pointLocations[destination], "N");
				if (newDist == Double.NaN) {
					newDist = 100;
				}
				long pathDist = curEdge.weight + minEdgeFromStart[curEdge.source];
				minEdgeFromStart[curEdge.destination] = pathDist;
				double value = newDist + pathDist;
				if (!visitedStart[curEdge.destination]) {
					visitedStart[curEdge.destination] = true;
					pqStart.add(new VertexA(curEdge.destination, value, curStart));
					visitedStartCheck.add(curEdge.destination);
				}
			}

			VertexA curDest = pqDest.poll();
			expanded++;
			Iterator<WeightedGraph.Edge> i2 = searchSpace.adjList[curDest.node].listIterator();
			while (i2.hasNext()) {
				WeightedGraph.Edge curEdge = i2.next();
				if (visitedDestCheck.contains(curEdge.destination)) {
					continue;
				}
				explored++;
				double newDist = longLatDist(pointLocations[curEdge.destination], pointLocations[start], "N");
				if (newDist == Double.NaN) {
					newDist = 100;
				}
				long pathDist = curEdge.weight + minEdgeFromDest[curEdge.source];
				minEdgeFromDest[curEdge.destination] = pathDist;
				double value = newDist + pathDist;
				if (!visitedDest[curEdge.destination]) {
					visitedDest[curEdge.destination] = true;
					pqDest.add(new VertexA(curEdge.destination, value, curDest));
					visitedDestCheck.add(curEdge.destination);
				}
			}
			for (int k = 0; k < visitedStart.length; k++) {
				if (visitedStart[k] == true && visitedDest[k] == true) {
					output[0] = expanded;
					output[1] = explored;
					return output;
				}
			}
		}
		return null;
	}

	/**
	 * aStarGraphSearch
	 */
	static long[] aStarGraphSearchBiApproach(int start, int destination, WeightedGraph.Graph searchSpace,
			Point[] pointLocations) {
		HashSet<Integer> visitedStartCheck = new HashSet<Integer>();
		HashSet<Integer> visitedDestCheck = new HashSet<Integer>();

		long expanded = 0;
		long explored = 1;
		long minEdgeFromStart[] = new long[searchSpace.verticies + 1];
		for (int i = 0; i < minEdgeFromStart.length; i++) {
			minEdgeFromStart[i] = Integer.MAX_VALUE;
		}
		long minEdgeFromDest[] = new long[searchSpace.verticies + 1];
		for (int i = 0; i < minEdgeFromDest.length; i++) {
			minEdgeFromDest[i] = Integer.MAX_VALUE;
		}
		boolean visitedStart[] = new boolean[searchSpace.verticies + 1];
		boolean visitedDest[] = new boolean[searchSpace.verticies + 1];
		PriorityQueue<VertexA> pqStart = new PriorityQueue<VertexA>();
		PriorityQueue<VertexA> pqDest = new PriorityQueue<VertexA>();
		long[] output = new long[10];
		double distance = longLatDist(pointLocations[start], pointLocations[destination], "N");

		VertexA sVert = new VertexA(start, distance, null);
		VertexA dVert = new VertexA(destination, distance, null);
		visitedStartCheck.add(start);
		visitedDestCheck.add(destination);

		pqStart.add(sVert);
		pqDest.add(dVert);

		visitedStart[start] = true;
		visitedDest[destination] = true;

		while (pqStart.size() != 0 && pqDest.size() != 0) {
			// start aStarGraph
			expanded++;
			VertexA curStart = pqStart.poll();
			Iterator<WeightedGraph.Edge> i = searchSpace.adjList[curStart.node].listIterator();
			while (i.hasNext()) {
				WeightedGraph.Edge curEdge = i.next();
				if (visitedStartCheck.contains(curEdge.destination)) {
					continue;
				}
				explored++;
				VertexA reach = pqDest.peek();
				double newDist = longLatDist(pointLocations[curEdge.destination], pointLocations[reach.node], "N");
				if (newDist == Double.NaN) {
					newDist = 100;
				}
				long pathDist = curEdge.weight + minEdgeFromStart[curEdge.source];
				minEdgeFromStart[curEdge.destination] = pathDist;
				double value = newDist + pathDist;
				if (!visitedStart[curEdge.destination]) {
					visitedStart[curEdge.destination] = true;
					pqStart.add(new VertexA(curEdge.destination, value, curStart));
					visitedStartCheck.add(curEdge.destination);
				}
			}

			VertexA curDest = pqDest.poll();
			expanded++;
			Iterator<WeightedGraph.Edge> i2 = searchSpace.adjList[curDest.node].listIterator();
			while (i2.hasNext()) {
				WeightedGraph.Edge curEdge = i2.next();
				if (visitedDestCheck.contains(curEdge.destination)) {
					continue;
				}
				explored++;
				VertexA reach = pqStart.peek();
				double newDist = longLatDist(pointLocations[curEdge.destination], pointLocations[reach.node], "N");
				if (newDist == Double.NaN) {
					newDist = 100;
				}
				long pathDist = curEdge.weight + minEdgeFromDest[curEdge.source];
				minEdgeFromDest[curEdge.destination] = pathDist;
				double value = newDist + pathDist;
				if (!visitedDest[curEdge.destination]) {
					visitedDest[curEdge.destination] = true;
					pqDest.add(new VertexA(curEdge.destination, value, curDest));
					visitedDestCheck.add(curEdge.destination);
				}
			}
			for (int k = 0; k < visitedStart.length; k++) {
				if (visitedStart[k] == true && visitedDest[k] == true) {
					output[0] = expanded;
					output[1] = explored;
					return output;
				}
			}
		}
		return null;
	}

	/**
	 * Coordinate Distance
	 */
	static double longLatDist(Point p1, Point p2, String unit) {
		double theta = p1.x - p2.x;
		double dist = Math.sin(deg2rad(p1.y)) * Math.sin(deg2rad(p1.y))
				+ Math.cos(deg2rad(p1.x)) * Math.cos(deg2rad(p1.x)) * Math.cos(deg2rad(theta));
		dist = Math.acos(dist);
		dist = rad2deg(dist);
		dist = dist * 60 * 1.1515;
		if (unit == "K") {
			dist = dist * 1.609344;
		} else if (unit == "N") {
			dist = dist * 0.8684;
		}
		return (dist);
	}

	/**
	 * Degree to Rad(ical dude)
	 * 
	 * @param deg
	 * @return
	 */
	private static double deg2rad(double deg) {
		return (deg * Math.PI / 180.0);
	}

	/**
	 * Rad(ical dude) to Degree
	 * 
	 * @param rad
	 * @return
	 */
	private static double rad2deg(double rad) {
		return (rad * 180.0 / Math.PI);
	}

	/**
	 * Graph Class
	 * 
	 * @author jarredparrett
	 *
	 */
	public static class WeightedGraph {

		static class Edge {
			int source;
			int destination;
			int weight;

			public Edge(int source, int destination, int weight) {
				this.source = source;
				this.destination = destination;
				this.weight = weight;
			}
		}

		static class Graph {
			int verticies;
			LinkedList<Edge>[] adjList;

			Graph(int verticies) {
				this.verticies = verticies;
				adjList = new LinkedList[verticies];
				for (int i = 0; i < verticies; i++) {
					adjList[i] = new LinkedList<>();
				}
			}

			public void addEdge(int source, int destination, int weight) {
				Edge edge = new Edge(source, destination, weight);
				adjList[source].addFirst(edge);
			}
		}
	}

	/**
	 * Vertex class for Dijsktra
	 * 
	 * @author jarredparrett
	 *
	 */
	static class Vertex implements Comparable<Vertex> {
		int node;
		long distance;
		Vertex parent;

		public Vertex(int start, long distance, Vertex parent) {
			this.node = start;
			this.distance = distance;
			this.parent = parent;
		}

		public int compareTo(Vertex o) {
			if (o.distance < distance) {
				return 1;
			} else if (this.distance > o.distance) {
				return 0;
			} else {
				return -1;
			}
		}

		public long getLengthOfPath() {
			long outVal = this.distance;
			Vertex curParent = this.parent;
			while (curParent != null) {
				outVal = +parent.distance;
				curParent = curParent.parent;
			}
			return outVal;
		}
	}

	/**
	 * Vertex class for Dijsktra
	 * 
	 * @author jarredparrett
	 *
	 */
	static class VertexA implements Comparable<VertexA> {
		int node;
		double distance;
		VertexA parent;

		public VertexA(int start, double distance, VertexA parent) {
			this.node = start;
			this.distance = distance;
			this.parent = parent;
		}

		public int compareTo(VertexA o) {
			if (o.distance < distance) {
				return 1;
			} else if (this.distance > o.distance) {
				return 0;
			} else {
				return -1;
			}
		}

		public double getLengthOfPath() {
			double outVal = this.distance;
			VertexA curParent = this.parent;
			while (curParent != null) {
				outVal = +parent.distance;
				curParent = curParent.parent;
			}
			return outVal;
		}
	}
}
