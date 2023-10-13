/*
 * Search.java
 *
 */

import java.io.*;
import java.lang.Math;
import java.util.*;

/**
 * This program finds path from source city to destination city using various
 * search algorithms such as BFS, DFS and A*.
 *
 * @author Anurag Kallurwar ak6491@rit.edu
 */


class Vertex implements Comparable<Vertex> {
    // This class holds the Vertex of the graph(map)

    public final String city;   // CityName
    public final String state;  // StateName
    public final double latitude;  // Latitude
    public final double longitude;  // Longitude
    public Double heuristics;   // Variable to hold heuristics estimate: h(v)
    public Double totalDistance; // Variable to hold total distance estimate:
    // g(v)
    public Double estimate; // Variable to hold total estimate: f(v)

    /**
     * Constructor
     * @param data Split line from city.dat file
     */
    public Vertex(String[] data) {
        city = data[0];
        state = data[1];
        latitude = Double.parseDouble(data[2]);
        longitude = Double.parseDouble(data[3]);
        heuristics = 0.0;
        totalDistance = Double.MAX_VALUE;
        estimate = Double.MAX_VALUE;
    }

    @Override
    public String toString() {
        return city;
    }

    @Override
    public int compareTo(Vertex o) {
        return city.toLowerCase().compareTo(o.city.toLowerCase());
    }
} // Vertex

class Edge implements Comparable<Edge> {
    // This class holds Edge for a given Vertex in the graph(map)

    public final Vertex toVertex; // Edge to Vertex
    public final double distance; // Distance (Weight)

    /**
     * Constructor
     * @param toVertex Destination Vertex (Edge to vertex)
     * @param distance Weight of edge calculated using latitudes and
     *                 longitudes of source and destination vertex
     */
    public Edge(Vertex toVertex, double distance) {
        this.toVertex = toVertex;
        this.distance = distance;
    }

    @Override
    public String toString() {
        return "[" + toVertex + ": " + distance + "]";
    }

    @Override
    public int compareTo(Edge o) {
        return toVertex.city.toLowerCase().compareTo(o.toVertex.city.toLowerCase());
    }
} // Edge

public class Search {
    // This class implements different search algorithms on the input map

    // HashMap to hold the graph(map)
    private final Map<Vertex, List<Edge>> graph = new HashMap<>();

    /**
     * Calculates Edge Weight
     * @param source Source vertex
     * @param destination destination vertex
     * @return the edge weight
     */
    public static double calculateEdgeDistance(Vertex source,
                                               Vertex destination) {
        return Math.sqrt(Math.pow(source.latitude - destination.latitude
                , 2) + Math.pow(source.longitude - destination.longitude, 2))
                * 100;
    }

    /**
     * Function to add Vertex to graph
     * @param input Input city and details
     */
    private void addVertex(String input) {
        String[] data = input.split("\\s+");
        // Removing null or empty values
        Arrays.asList(data).removeIf(str -> str == null || "".equals(str));
        Vertex source = new Vertex(data);
        graph.put(source, new LinkedList<>());
    }

    /**
     * Function to find vertex
     * @param city City Vertex to be searched
     * @return Vertex with given city
     */
    private Vertex getVertex(String city) {
        for (Vertex vertex : graph.keySet()) {
            if (city.equals(vertex.city))
                return vertex;
        }
        return null;
    }

    /**
     * Function to add Edge to a Vertex of the graph
     * @param input Input edge details
     */
    private void addEdge(String input) throws NullPointerException {
        String[] data = input.split("\\s+");
        Arrays.asList(data).removeIf(str -> str == null || "".equals(str));

        Vertex source = getVertex(data[0]);
        Vertex destination = getVertex(data[1]);
        if (source == null || destination == null) // If vertices not found
            throw new NullPointerException(data[0] + " or " + data[1] +
                    " Vertex was not present in city.dat");

        // Calculate weight
        double distance = calculateEdgeDistance(source, destination);
        // Add edge for both vertices
        graph.get(source).add(new Edge(destination, distance));
        graph.get(destination).add(new Edge(source, distance));
    }

    /**
     * Function to find edge between two vertices
     * @param fromVertex source vertex
     * @param toVertex destination vertex
     * @return Edge
     */
    private Edge getEdge(Vertex fromVertex, Vertex toVertex) {
        for (Edge edge: graph.get(fromVertex))
            if (edge.toVertex.compareTo(toVertex) == 0)
                return edge;
        return null;
    }

    /**
     * Sort Edges for every vertex in the graph
     * @param desc whether descending order = true/false
     */
    private void sortEdges(boolean desc) {
        if (!desc)
            for (Vertex key: graph.keySet()) {
                 Collections.sort(graph.get(key));
            }
        else
            for (Vertex key: graph.keySet()) {
                Collections.reverse(graph.get(key));
            }
    }

    /**
     * Count of Vertices
     * @return count of vertices
     */
    private int countOfVertices() {
        return graph.keySet().size();
    }

    /**
     * Count of Edges
     * @return count of edges
     */
    private int countOfEdges() {
        int count = 0;
        for (Vertex key: graph.keySet()) {
            count += graph.get(key).size();
        }
        return count;
    }

    /**
     * Pretty prints the graph
     */
    private void printGraph() {
        for (Vertex key: graph.keySet()) {
            System.out.print(key + " : [ ");
            for (Edge edge: graph.get(key))
                System.out.print(edge);
            System.out.println(" ]");
        }
    }

    /**
     * Function implements BFS algorithm
     * @param start Start vertex
     * @param destination destination vertex
     * @return Predecessors map visited vertices
     */
    private Map<Vertex, Vertex> bfs(Vertex start, Vertex destination) {
        // Predecessors Map / Visited / Closed List
        Map<Vertex, Vertex> predecessors = new HashMap<>();
        // Queue / Open List
        Queue<Vertex> queue = new LinkedList<>();
        // Initializing
        predecessors.put(start, null);
        queue.add(start);

        // Till Queue is empty
        while (!queue.isEmpty()) {
            // Dequeue from queue
            Vertex vertex = queue.remove();
            // If vertex is destination (Late Goal Test)
            if (vertex.compareTo(destination) == 0)
                break;
            // Find neighbors
            for (Edge edge: graph.get(vertex))
                if (!predecessors.containsKey(edge.toVertex)) {
                    // If neighbor not visited
                    predecessors.put(edge.toVertex, vertex);
                    // Enqueue neighbor to queue
                    queue.add(edge.toVertex);
                    // If neighbor is destination (Early Goal Test)
                    if (edge.toVertex.compareTo(destination) == 0) {
                        queue = new LinkedList<>();
                        break;
                    }
                }
        }
        return predecessors;
    }

    /**
     * Function implements DFS algorithm
     * @param start Start vertex
     * @param destination destination vertex
     * @return Predecessors map visited vertices
     */
    private Map<Vertex, Vertex> dfs(Vertex start, Vertex destination) {
        // Predecessors Map / Visited / Closed List
        Map<Vertex, Vertex> predecessors = new HashMap<>();
        // Stack / Open List
        Stack<Vertex> stack = new Stack<>();
        // Initializing
        predecessors.put(start, null);
        stack.push(start);

        // Till Stack is Empty
        while (!stack.isEmpty()) {
            // Pop from stack
            Vertex vertex = stack.pop();
            // If vertex is destination (Late Goal Test)
            if (vertex.compareTo(destination) == 0)
                break;
            // Finding neighbors
            for (Edge edge: graph.get(vertex)) {
                if (!predecessors.containsKey(edge.toVertex)) {
                    // If neighbor not visited
                    predecessors.put(edge.toVertex, vertex);
                    // Push neighbor to stack
                    stack.push(edge.toVertex);
                    // If neighbor is destination (Early Goal Test)
                    if(edge.toVertex.compareTo(destination) == 0) {
                        stack = new Stack<>();
                        break;
                    }
                }
            }
        }
        return predecessors;
    }

    /**
     * Calculates Straight-Line Distance (SLD) heuristics with respect to
     * destination vertex
     * @param destination destination vertex
     */
    public void calculateHeuristics(Vertex destination) {
        for (Vertex vertex : graph.keySet()) {
            vertex.heuristics = calculateEdgeDistance(vertex, destination);
        }
    }

    /**
     * Function implements A* algorithm
     * @param start Start vertex
     * @param destination destination vertex
     * @return Predecessors map visited vertices
     */
    private Map<Vertex, Vertex> aStar(Vertex start, Vertex destination) {
        // Predecessors Map / Visited / Closed List
        Map<Vertex, Vertex> predecessors = new HashMap<>();
        // Priority Queue / Open List
        PriorityQueue<Vertex> openQueue =
                new PriorityQueue<>((a, b) -> a.estimate.compareTo(b.estimate));
        // Initializing
        start.totalDistance = 0.0;
        start.estimate = start.totalDistance + start.heuristics;
        predecessors.put(start, null);
        openQueue.add(start);

        // Till priority queue is empty
        while (!openQueue.isEmpty()) {
            // Dequeue from priority queue
            Vertex vertex = openQueue.poll();
            // If vertex is destination (Late Goal Test)
            if(vertex.compareTo(destination) == 0)
                break;
            // Finding neighbors
            for (Edge edge: graph.get(vertex)) {
                // Calculating estimate for neighbor
                Vertex neighbor = edge.toVertex;
                double currentDistance = edge.distance + vertex.totalDistance;
                // If new estimate is less than old estimate
                if (currentDistance < neighbor.totalDistance) {
                    openQueue.remove(neighbor);
                    neighbor.totalDistance = currentDistance;
                    neighbor.estimate =
                            neighbor.totalDistance + neighbor.heuristics;
                    // Update neighbor to priority queue and update parent
                    predecessors.put(neighbor, vertex);
                    openQueue.add(neighbor);
                }
            }
        }
        return predecessors;
    }

    /**
     * Calculates total path distance from source to destination vertex
     * @param order path from source to destination
     * @return integer distance in miles
     */
    private int calculateTotalDistance(List<Vertex> order) {
        double distance = 0;
        for (int i = 0; i < (order.size() - 1); i++) {
            Edge edge;
            // Finding edge and its distance
            if ((edge = getEdge(order.get(i), order.get(i + 1))) != null)
                distance += edge.distance;
            else {
                distance = 0;
                break;
            }
        }
        return (int) Math.round(distance);
    }

    /**
     * Creates the pretty output format
     * @param predecessors Predecessors map
     * @param start Start vertex
     * @param destination Destination vertex
     * @param type Type of Algorithm
     *             1    : BFS
     *             0    : A*
     *             -1   : DFS
     * @return pretty output
     */
    public String createOutput(Map<Vertex, Vertex> predecessors, Vertex start,
                               Vertex destination, int type) {
        // Retracing the path
        LinkedList<Vertex> cities = new LinkedList<>();
        Vertex current = destination;
        cities.add(destination);
        while (current.compareTo(start) != 0) {
            if (predecessors.containsKey(current)) {
                current = predecessors.get(current);
                cities.addFirst(current);
            }
            else {
                cities = null;
                break;
            }
        }

        // Pretty output
        String output = "\n";
        if (cities != null) {
            if (type > 0)
                output += "Breadth-First Search Results: \n";
            else if (type < 0)
                output += "Depth-First Search Results: \n";
            else
                output += "A* Search Results: \n";
            for (Vertex city : cities)
                output += city + "\n";
            output += "That took " + (cities.size() - 1) + " hops to find.\n";
            output += "Total distance = " + calculateTotalDistance(cities) +
                    " miles.\n";
        }
        return (output + "\n");
    }

    /**
     * Print output to STDOut or a given output file
     * @param output pretty output
     * @param printOutTo fileName of output file
     * @throws Exception if error while writing to file
     */
    private void printOutput(String output, String printOutTo)
            throws Exception {
        if (printOutTo.isEmpty()) // printOutTo empty => STDOut
            System.out.print(output);
        else { // Print to output file
            try (BufferedWriter out =
                         new BufferedWriter(new FileWriter(printOutTo))) {
                out.write(output);
            }
            catch (IOException ex) {
                throw new Exception(ex);
            }
        }
    }

    /**
     * Find Path from source to destination using various paths
     * @param start source city
     * @param destination destination city
     * @param printOutTo fileName of output file
     * @throws Exception if any error occurs
     */
    public void findPath(String start, String destination,
                         String printOutTo) throws Exception {
        Vertex startV = getVertex(start);
        Vertex endV = getVertex(destination);
        // If vertices not found
        if (startV == null)
            throw new NullPointerException("No such city: " + start);
        else if (endV == null)
            throw new NullPointerException("No such city: " + destination);

        // Map to hold predecessors map
        Map<Vertex, Vertex> predecessors;
        String output = "";

        // BFS
        sortEdges(false); // Edges in ascending
        predecessors = bfs(startV, endV);
        output += createOutput(predecessors, startV, endV, 1);

        // DFS
        sortEdges(true); // Edges in descending
        predecessors = dfs(startV, endV);
        output += createOutput(predecessors, startV, endV, -1);

        // A*
        calculateHeuristics(endV); // Calculating heuristics
        predecessors = aStar(startV, endV);
        output += createOutput(predecessors, startV, endV, 0);

        // Print output
        printOutput(output, printOutTo);
    }

    /**
     * Read file according to filename
     * @param fileName name of file
     * @return file content
     * @throws NullPointerException if missing file content
     * @throws IOException if issue while reading file or file not exists
     */
    private String[] readFile(String fileName) throws NullPointerException,
            IOException {
        try (BufferedReader input =
                     new BufferedReader(new FileReader((fileName)))) {
            String line;
            if (fileName.equals("city.dat")) {
                while ((line = input.readLine()) != null)
                    addVertex(line);
            } else if (fileName.equals("edge.dat")) {
                while ((line = input.readLine()) != null)
                    addEdge(line);
            } else {
                String[] inputCities = new String[2];
                int count = 0;
                while ((line = input.readLine()) != null && count < 2)
                    inputCities[count++] = line.trim();
                return inputCities;
            }
        } catch (IOException ex) {
            throw new IOException(fileName);
        }
        return null;
    }

    /**
     * Main method
     * @param args Command Line Arguments
     */
    public static void main(String[] args) {
        Search search = new Search();
        String printOutput = "";
        try {
            String[] inputCities = new String[2];
            if (args.length == 2) { // checking number of arguments
                if (args[0].trim().equals("-")) {
                    // STD Input
                    Scanner scan = new Scanner(System.in);
                    inputCities[0] = scan.nextLine().trim();
                    inputCities[1] = scan.nextLine().trim();
                }
                else {
                    // File Input
                    inputCities = search.readFile(args[0].trim());
                }
                if (!args[1].trim().equals("-")) {
                    // File Output
                    printOutput = args[1].trim();
                }
                if (inputCities == null) { // If input values missing
                    throw new NullPointerException("Cities missing " +
                            "in input");
                }
                for(int i = 0; i < 2; i++) { // If input values missing
                    if(inputCities[i] == null) {
                        throw new NullPointerException("Cities missing " +
                                "in input");
                    }
                }
            }
            else // if missing arguments
                throw new Exception("Usage: java Search inputFile outputFile");
            search.readFile("city.dat");
            search.readFile("edge.dat");
            // Find Paths
            search.findPath(inputCities[0], inputCities[1], printOutput);
            /*
            search.printGraph();
            System.out.println("Vertices: " + search.countOfVertices());
            System.out.println("Directed Edges: " + search.countOfEdges());
            search.findShortestPath("Olympia", "SaltLakeCity");
            search.findShortestPath("Olympia", "Olympia");
            search.findShortestPath("Denver", "Boston");
            */
        } catch (IOException e) {
            System.err.println("File not found: " + e.getMessage());
        } catch (Exception e) {
            System.err.println(e.getMessage());
        }
    }
} // Search
