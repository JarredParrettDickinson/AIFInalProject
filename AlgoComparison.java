import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.time.Duration;
import java.time.Instant;
import java.util.*;

//current - working
public class AlgoComparison {

	public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
//		File fileright = new File("out6.txt");
//		PrintWriter writer = new PrintWriter(fileright);
//		writer.println(
//				"BFS,BFSTIME,BFSBI,BFSBITIME,Dij,DijTIME,A*,A*TIME,START,START LOCATION,END,END LOCATION,EUCLIDEAN DISTANCE,LONGLAT DISTANCE");
		String location = "/Users/jarredparrett/Desktop/USA-road-d.NE.gr";
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
			String line = sc.nextLine();
			String[] line_S = line.split(" ");
			int source = Integer.parseInt(line_S[1]);
			int dest = Integer.parseInt(line_S[2]);
			int val = Integer.parseInt(line_S[3]);
			searchSpace.addEdge(source, dest, val);
		}
		String location_lat = "/Users/jarredparrett/Desktop/USA-road-d.NE.co";
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
			String cur = sc2.nextLine();
			String[] curSplit = cur.split(" ");
			int point = Integer.parseInt(curSplit[1]);
			double x = Double.parseDouble(curSplit[3]);
			double y = Double.parseDouble(curSplit[2]);
			pointLocations[point] = new Point(x, y);
		}
		sc.close();
		sc2.close();

		Random rand = new Random();
		for (int i = 0; i < 10; i++) {
			int start = 414206;//rand.nextInt((arrLen) + 1) + 0;
			int end = 626524;//rand.nextInt((arrLen) + 1) + 0;

			double[] dij = dijkstraGraphSearch(start, end, searchSpace);
			double[] outBFS = breathFirstSearch(start, end, searchSpace);
			double[] outBi = BFSBI2(start, end, searchSpace);
			double[] outA = aStarGraphSearch(start, end, searchSpace, pointLocations, 0.0);
			double outBFSPathDist = outBFS[0];
			double outBFSTime = outBFS[1];
			double outBiPathDist = outBi[0];
			double outBiTime = outBi[1];
			double dijPathDist = dij[0];
			double dijTime = dij[1];
			double outAStarPathDist = outA[0];
			double outAStarTime = outA[1];
			System.out.println("Round Dij " + i + dijTime);
			System.out.println("Round A* " + i +outAStarTime);
			
			
//			writer.print(outBFSPathDist);
//			writer.print(",");
//			writer.print(outBFSTime);
//			writer.print(",");
//			writer.print(outBiPathDist);
//			writer.print(",");
//			writer.print(outBiTime);
//			writer.print(",");
//			writer.print(dijPathDist);
//			writer.print(",");
//			writer.print(dijTime);
//			writer.print(",");
//			writer.print(outAStarPathDist);
//			writer.print(",");
//			writer.print(outAStarTime);
//			writer.print(",");
//			writer.print(start);
//			writer.print(",");
//			writer.print((pointLocations[start].x + " " + pointLocations[start].y));
//			writer.print(",");
//			writer.print(end);
//			writer.print(",");
//			writer.print((pointLocations[end].x + " " + pointLocations[end].y));
//			writer.print(",");
//			writer.print(longLatDist(pointLocations[start], pointLocations[end], "M"));
//			writer.print(",");
//			writer.print(0.621371 * (Haversine.distance(pointLocations[start], pointLocations[end])));
//			writer.println();

		}

//		writer.close();
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
		Instant startTime = Instant.now();
		long[] visitedNum = new long[searchSpace.verticies + 1];
		for (int i = 0; i < visitedNum.length; i++) {
			visitedNum[i] = Integer.MAX_VALUE;
		}
		VertexA startVertex = new VertexA(start, 0, null, 0);
		double[] output = new double[10];
		Queue<VertexA> q = new LinkedList<VertexA>();
		q.add(startVertex);
		visitedNum[start] = 0;
		while (q.size() != 0) {
			startVertex = q.poll();
			Iterator<WeightedGraph.Edge> i = searchSpace.adjList[startVertex.node].listIterator();
			while (i.hasNext()) {
				WeightedGraph.Edge next = i.next();
				if (visitedNum[next.destination] != Integer.MAX_VALUE
						|| (next.weight + startVertex.weightForQueue) > visitedNum[next.destination]) {
					continue;
				}
				int n = next.destination;
				visitedNum[n] = (long) (next.weight + startVertex.weightForQueue);
				VertexA nVert = new VertexA(n, (startVertex.weightForQueue + next.weight), startVertex, 0);
				q.add(nVert);
				if (visitedNum[destination] != Integer.MAX_VALUE) {
					output[0] = visitedNum[destination];
					Instant endTime = Instant.now();
					long timeElapsed = Duration.between(startTime, endTime).toMillis();
					output[1] = timeElapsed;
					return output;
				}

			}
		}
		return output;
	}

	/**
	 * Bidirectional BFS Search
	 * 
	 * @param start
	 * @param destination
	 * @param searchSpace
	 * @return
	 */
	public static double[] BFSBI2(int start, int destination, WeightedGraph.Graph searchSpace) {
		Instant startTime = Instant.now();
		long[] visitedNumStart = new long[searchSpace.verticies + 1];
		for (int i = 0; i < visitedNumStart.length; i++) {
			visitedNumStart[i] = Integer.MAX_VALUE;
		}
		long[] visitedNumDest = new long[searchSpace.verticies + 1];
		for (int i = 0; i < visitedNumDest.length; i++) {
			visitedNumDest[i] = Integer.MAX_VALUE;
		}
		VertexA startVertex = new VertexA(start, 0, null, 0);
		VertexA destVertex = new VertexA(destination, 0, null, 0);
		double[] output = new double[10];
		Queue<VertexA> qStart = new LinkedList<VertexA>();
		Queue<VertexA> qDest = new LinkedList<VertexA>();
		qStart.add(startVertex);
		qDest.add(destVertex);
		visitedNumStart[start] = 0;
		visitedNumDest[destination] = 0;
		while (qStart.size() != 0 || qDest.size() != 0) {
			startVertex = qStart.poll();
			destVertex = qDest.poll();
			Iterator<WeightedGraph.Edge> i = searchSpace.adjList[startVertex.node].listIterator();
			while (i.hasNext()) {
				WeightedGraph.Edge next = i.next();
				if (visitedNumStart[next.destination] != Integer.MAX_VALUE
						|| (next.weight + startVertex.weightForQueue) > visitedNumStart[next.destination]) {
					continue;
				}
				int n = next.destination;
				visitedNumStart[n] = (long) (next.weight + startVertex.weightForQueue);
				VertexA nVert = new VertexA(n, (startVertex.weightForQueue + next.weight), startVertex, 0);
				if (visitedNumDest[n] != Integer.MAX_VALUE) {
					Instant endTime = Instant.now();
					output[0] = visitedNumDest[n] + visitedNumStart[n];
					long timeElapsed = Duration.between(startTime, endTime).toMillis();
					output[1] = timeElapsed;
					return output;
				}
				qStart.add(nVert);
			}
			Iterator<WeightedGraph.Edge> j = searchSpace.adjList[destVertex.node].listIterator();
			while (j.hasNext()) {
				WeightedGraph.Edge next = j.next();
				int n = next.destination;
				if (visitedNumDest[next.destination] != Integer.MAX_VALUE
						|| (next.weight + destVertex.weightForQueue) > visitedNumDest[next.destination]) {
					continue;
				}
				visitedNumDest[n] = (long) (next.weight + destVertex.weightForQueue);
				VertexA nVert = new VertexA(n, (destVertex.weightForQueue + next.weight), destVertex, 0);
				if (visitedNumStart[n] != Integer.MAX_VALUE) {
					Instant endTime = Instant.now();
					output[0] = visitedNumDest[n] + visitedNumStart[n];
					long timeElapsed = Duration.between(startTime, endTime).toMillis();
					output[1] = timeElapsed;
					return output;
				}
				qDest.add(nVert);
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
		PriorityQueue<VertexA> pq = new PriorityQueue<VertexA>();
		minEdge[start] = 0;
		VertexA vertStart = new VertexA(start, 0, null, 0);
		pq.add(vertStart);
		int total = 0;
		while (pq.size() != 0) {
			total++;
			VertexA cur = pq.poll();
			if (visitedCheck.contains(cur.node)) {
				continue;
			}
			visitedCheck.add(cur.node);
			minEdge[cur.node] = (long) cur.weightForQueue;

			if (cur.node == destination) {
				Instant finish = Instant.now();
				long timeElapsed = Duration.between(begin, finish).toMillis();
				output[0] = cur.weightForQueue;
				output[1] = timeElapsed;
				return output;
			}
			for (WeightedGraph.Edge curEdge : searchSpace.adjList[cur.node]) {

				long pathDist = curEdge.weight + minEdge[curEdge.source];

				if (visitedCheck.contains(curEdge.destination)) {
					continue;
				}
				pq.add(new VertexA(curEdge.destination, pathDist, cur, 0));
			}
		}
		return output;
	}

	/**
	 * aStarGraphSearch
	 */
	static double[] aStarGraphSearch(int start, int destination, WeightedGraph.Graph searchSpace,
			Point[] pointLocations, double factor) {
		Instant begin = Instant.now();
		HashMap<Integer, Double> visited = new HashMap<Integer, Double>();

		double[] output = new double[10];
		PriorityQueue<VertexA> pq = new PriorityQueue<VertexA>();
		double distanceToDestination = longLatDist(pointLocations[start], pointLocations[destination], "M");
		VertexA startVertex = new VertexA(start, distanceToDestination, null, 0);
		pq.add(startVertex);
		while (pq.size() != 0) {
			VertexA currentVertex = pq.poll();
			Iterator<WeightedGraph.Edge> i = searchSpace.adjList[currentVertex.node].listIterator();
			if (visited.containsKey(currentVertex.node)) {
				continue;
			}
			visited.put(currentVertex.node, currentVertex.phyDist);
			if (visited.containsKey(destination)) {
				Instant finish = Instant.now();

				long timeElapsed = Duration.between(begin, finish).toMillis();
				output[0] = visited.get(destination);
				output[1] = timeElapsed;
				return output;
			}
			while (i.hasNext()) {
				WeightedGraph.Edge curEdge = i.next();
				double pathFromStartToThisNode = curEdge.weight + visited.get(currentVertex.node);

				if (visited.containsKey(curEdge.destination)) {
					continue;
				}
				double straightLineDistToDestination = longLatDist(pointLocations[curEdge.destination],
						pointLocations[destination], "M");

				double weightForQueue = straightLineDistToDestination + pathFromStartToThisNode;

				pq.add(new VertexA(curEdge.destination, weightForQueue, currentVertex, pathFromStartToThisNode));

			}

		}
		return null;
	}

	/**
	 * Coordinate Distance
	 */
	static double longLatDist(Point p1, Point p2, String unit) {
		double xVal = Math.pow(p2.x - p1.x, 2);
		double yVal = Math.pow(p2.y - p1.y, 2);
		double dVal = Math.sqrt(xVal + yVal);
		return dVal;
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
	 * Vertex class
	 * 
	 * @author jarredparrett
	 *
	 */
	static class VertexA implements Comparable<VertexA> {
		int node;
		double weightForQueue;
		VertexA parent;
		double phyDist;

		public VertexA(int start, double weightForQueue, VertexA parent, double phyDist) {
			this.node = start;
			this.weightForQueue = weightForQueue;
			this.parent = parent;
			this.phyDist = phyDist;
		}

		public int compareTo(VertexA o) {
			if (o.weightForQueue < weightForQueue) {
				return 1;
			} else if (this.weightForQueue < o.weightForQueue) {
				return -1;
			} else {
				return 0;
			}
		}

		public double getLengthOfPath() {
			double outVal = this.weightForQueue;
			VertexA curParent = this.parent;
			while (curParent != null) {
				outVal = +parent.weightForQueue;
				curParent = curParent.parent;
			}
			return outVal;
		}
	}

	static class Point {
		double x;
		double y;

		public Point(double x, double y) {
			this.x = x;
			this.y = y;
		}
	}

	/**
	 * https://github.com/jasonwinn/haversine/blob/master/Haversine.java
	 * 
	 * @author Jason Winn
	 *
	 */
	public static class Haversine {
		private static final int EARTH_RADIUS = 6371;

		public static double distance(Point p1, Point p2) {

			double startLat = p1.x / 10000;
			double startLong = p1.y / 10000;
			double endLat = p2.x / 10000;
			double endLong = p2.y / 10000;
			double dLat = Math.toRadians((endLat - startLat));
			double dLong = Math.toRadians((endLong - startLong));

			startLat = Math.toRadians(startLat);
			endLat = Math.toRadians(endLat);

			double a = haversin(dLat) + Math.cos(startLat) * Math.cos(endLat) * haversin(dLong);
			double c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
			System.out.println((EARTH_RADIUS * c) * 621371);
			return EARTH_RADIUS * c; // <-- d
		}

		public static double haversin(double val) {
			return Math.pow(Math.sin(val / 2), 2);
		}
	}
}
