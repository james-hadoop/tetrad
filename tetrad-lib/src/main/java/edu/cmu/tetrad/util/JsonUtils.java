package edu.cmu.tetrad.util;

import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.graph.EdgeTypeProbability.EdgeType;
import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Dec 9, 2016 5:43:47 PM
 *
 * @author Chirayu (Kong) Wongchokprasitti, PhD
 */
public class JsonUtils {

    public static Graph parseJSONObjectToTetradGraph(String jsonResponse) {
        return JsonUtils.parseJSONObjectToTetradGraph(new JSONObject(jsonResponse));
    }

    public static Graph parseJSONObjectToTetradGraph(JSONObject jObj) {
        if (!jObj.isNull("graph")) {
            return JsonUtils.parseJSONObjectToTetradGraph(jObj.getJSONObject("graph"));
        }

        // Node
        List<Node> nodes = JsonUtils.parseJSONArrayToTetradNodes(jObj.getJSONArray("nodes"));
        EdgeListGraph graph = new EdgeListGraph(nodes);

        // Edge
        Set<Edge> edges = JsonUtils.parseJSONArrayToTetradEdges(graph, jObj.getJSONArray("edgesSet"));
        for (Edge edge : edges) {
            graph.addEdge(edge);
        }

        // ambiguousTriples
        Set<Triple> ambiguousTriples = JsonUtils.parseJSONArrayToTetradTriples(jObj.getJSONArray("ambiguousTriples"));
        for (Triple triple : ambiguousTriples) {
            graph.addAmbiguousTriple(triple.getX(), triple.getY(), triple.getZ());
        }

        // underLineTriples
        Set<Triple> underLineTriples = JsonUtils.parseJSONArrayToTetradTriples(jObj.getJSONArray("underLineTriples"));
        for (Triple triple : underLineTriples) {
            graph.addUnderlineTriple(triple.getX(), triple.getY(), triple.getZ());
        }

        // dottedUnderLineTriples
        Set<Triple> dottedUnderLineTriples = JsonUtils.parseJSONArrayToTetradTriples(jObj.getJSONArray("dottedUnderLineTriples"));
        for (Triple triple : dottedUnderLineTriples) {
            graph.addDottedUnderlineTriple(triple.getX(), triple.getY(), triple.getZ());
        }

        // stuffRemovedSinceLastTripleAccess
        boolean stuffRemovedSinceLastTripleAccess = jObj.getBoolean("stuffRemovedSinceLastTripleAccess");
        graph.setStuffRemovedSinceLastTripleAccess(stuffRemovedSinceLastTripleAccess);

        // highlightedEdges
        Set<Edge> highlightedEdges = JsonUtils.parseJSONArrayToTetradEdges(graph, jObj.getJSONArray("highlightedEdges"));
        for (Edge edge : highlightedEdges) {
            graph.setHighlighted(edge, true);
        }

        return graph;
    }

    public static Set<Triple> parseJSONArrayToTetradTriples(JSONArray jArray) {
        Set<Triple> triples = new HashSet<>();

        for (int i = 0; i < jArray.length(); i++) {
            Triple triple = JsonUtils.parseJSONArrayToTetradTriple(jArray.getJSONObject(i));
            triples.add(triple);
        }

        return triples;
    }

    public static Triple parseJSONArrayToTetradTriple(JSONObject jObj) {
        Node x = JsonUtils.parseJSONObjectToTetradNode(jObj.getJSONObject("x"));
        Node y = JsonUtils.parseJSONObjectToTetradNode(jObj.getJSONObject("y"));
        Node z = JsonUtils.parseJSONObjectToTetradNode(jObj.getJSONObject("z"));

        return new Triple(x, y, z);
    }

    public static Set<Edge> parseJSONArrayToTetradEdges(Graph graph, JSONArray jArray) {
        Set<Edge> edges = new HashSet<>();

        for (int i = 0; i < jArray.length(); i++) {
            Edge edge = JsonUtils.parseJSONObjectToTetradEdge(graph, jArray.getJSONObject(i));
            edges.add(edge);
        }

        return edges;
    }

    public static Edge parseJSONObjectToTetradEdge(Graph graph, JSONObject jObj) {
        Node node1 = graph.getNode(jObj.getJSONObject("node1").getString("name"));
        Node node2 = graph.getNode(jObj.getJSONObject("node2").getString("name"));
        Endpoint endpoint1 = Endpoint.TYPES[jObj.getJSONObject("endpoint1").getInt("ordinal")];
        Endpoint endpoint2 = Endpoint.TYPES[jObj.getJSONObject("endpoint2").getInt("ordinal")];
        Edge edge = new Edge(node1, node2, endpoint1, endpoint2);

        try {
            // properties
            JSONArray jArray = jObj.getJSONArray("properties");
            if (jArray != null) {
                for (int i = 0; i < jArray.length(); i++) {
                    edge.addProperty(JsonUtils.parseJSONObjectToEdgeProperty(jArray.getString(i)));
                }
            }
        } catch (JSONException e) {
            // TODO Auto-generated catch block
            //e.printStackTrace();
        }

        try {
            // properties
            JSONArray jArray = jObj.getJSONArray("edgeTypeProbabilities");
            if (jArray != null) {
                for (int i = 0; i < jArray.length(); i++) {
                    edge.addEdgeTypeProbability(JsonUtils.parseJSONObjectToEdgeTypeProperty(jArray.getJSONObject(i)));
                }
            }
        } catch (JSONException e) {
            // TODO Auto-generated catch block
            //e.printStackTrace();
        }

        return edge;
    }

    public static EdgeTypeProbability parseJSONObjectToEdgeTypeProperty(JSONObject jObj) {
        String _edgeType = jObj.getString("edgeType");
        EdgeType edgeType = EdgeType.nil;
        switch (_edgeType) {
            case "ta":
                edgeType = EdgeType.ta;
                break;
            case "at":
                edgeType = EdgeType.at;
                break;
            case "ca":
                edgeType = EdgeType.ca;
                break;
            case "ac":
                edgeType = EdgeType.ac;
                break;
            case "cc":
                edgeType = EdgeType.cc;
                break;
            case "aa":
                edgeType = EdgeType.aa;
                break;
            case "tt":
                edgeType = EdgeType.tt;
                break;
        }
        double probability = jObj.getDouble("probability");
        EdgeTypeProbability edgeTypeProbability = new EdgeTypeProbability(edgeType, probability);

        try {
            // properties
            JSONArray jArray = jObj.getJSONArray("properties");
            if (jArray != null) {
                for (int i = 0; i < jArray.length(); i++) {
                    edgeTypeProbability.addProperty(JsonUtils.parseJSONObjectToEdgeProperty(jArray.getString(i)));
                }
            }
        } catch (JSONException e) {
            // TODO Auto-generated catch block
            //e.printStackTrace();
        }
        return edgeTypeProbability;
    }

    public static Edge.Property parseJSONObjectToEdgeProperty(String prop) {
        if (prop.equalsIgnoreCase("dd")) {
            return Edge.Property.dd;
        }
        if (prop.equalsIgnoreCase("nl")) {
            return Edge.Property.nl;
        }
        if (prop.equalsIgnoreCase("pd")) {
            return Edge.Property.pd;
        }
        if (prop.equalsIgnoreCase("pl")) {
            return Edge.Property.pl;
        }
        return null;
    }

    public static List<Node> parseJSONArrayToTetradNodes(JSONArray jArray) {
        List<Node> nodes = new ArrayList<>();

        for (int i = 0; i < jArray.length(); i++) {
            Node node = JsonUtils.parseJSONObjectToTetradNode(jArray.getJSONObject(i));
            nodes.add(node);
        }

        return nodes;
    }

    public static Node parseJSONObjectToTetradNode(JSONObject jObj) {
        JSONObject nodeType = jObj.getJSONObject("nodeType");
        int ordinal = nodeType.getInt("ordinal");
        int centerX = jObj.getInt("centerX");
        int centerY = jObj.getInt("centerY");
        String name = jObj.getString("name");

        GraphNode graphNode = new GraphNode(name);
        graphNode.setNodeType(NodeType.TYPES[ordinal]);
        graphNode.setCenter(centerX, centerY);

        return graphNode;
    }

}
