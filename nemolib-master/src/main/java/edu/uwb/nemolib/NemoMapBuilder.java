package edu.uwb.nemolib;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * @Author: Zican Li
 * @Date：10/27/2020 5:36 PM
 */
public class NemoMapBuilder {
    private List<Integer> mappedHNodes;
    private BufferedWriter bufferedWriter;
    private HashMap<Integer, String> mapIndexToString;
    public NemoMapBuilder(String filename) {
        try {
            bufferedWriter = new BufferedWriter(new FileWriter(filename));
        }catch (Exception e){
            System.out.println(e.getMessage());
        }

    }

    /**
     * Method to find the most constrained neighboring node of mapped nodes in the query graph.
     *
     * @param partialMap:List[int] - the current partial mapping of query graph to target graph
     * @param queryGraph:Graph     - the query graph
     * @return int - number corresponding to most constrained node
     */
    private int getMostConstrainedNeighbor(List<Integer> partialMap, Graph queryGraph) {
        // get all neoghbors of already mapped nodes
        // NO Duplicates
        List<Integer> neighborList = new ArrayList<>();
        neighborList = chooseNeighborsOfRange(partialMap, queryGraph, neighborList);

        if (neighborList == null || neighborList.size() == 0) {
            return -1;
        } else if (neighborList.size() == 1) {
            return neighborList.get(0);
        }
        // 2D list to create pairs
        List<List<Integer>> constrainRank = new ArrayList<>();
        for (int i = 0; i < neighborList.size(); i++) {
            constrainRank.add(new ArrayList<>());
            constrainRank.get(i).add(0);
            constrainRank.get(i).add(neighborList.get(i));
        }
        for (int i = 0; i < neighborList.size(); i++) {
            for (Integer vertex : getNeighors(queryGraph, constrainRank.get(i).get(1))) {
                if (partialMap.contains(vertex)) {
                    constrainRank.get(i).set(0, constrainRank.get(i).get(0) + 1);
                }
            }
        }
        // Rank neighbor nodes with most already-mapped neighbors
        constrainRank.sort((a, b) -> {
            if (!a.get(0).equals(b.get(0))) {
                return b.get(0) - a.get(0);
            } else {
                return b.get(1) - a.get(1);
            }
        });

        int highestNeighborMapped = constrainRank.get(0).get(0);
        int count = neighborList.size();
        for (int i = 1; i < neighborList.size(); i++) {
            if (constrainRank.get(i).get(0) < highestNeighborMapped) {
                if (i == 1) {
                    return constrainRank.get(0).get(1);
                }
                count = i;
                break;
            }
        }
        // Rank neighbor nodes with highest degree
        for (int i = 0; i < count; i++) {
            constrainRank.get(i).set(0, queryGraph.getAdjacencyList(constrainRank.get(i).get(1)).size());
        }
        constrainRank = constrainRank.subList(0, count);
        constrainRank.sort((a, b) -> {
            if (!a.get(0).equals(b.get(0))) {
                return b.get(0) - a.get(0);
            } else {
                return b.get(1) - a.get(1);
            }
        });

        int highestDegree = constrainRank.get(0).get(0);
        for (int i = 1; i < count; i++) {
            if (constrainRank.get(i).get(0) < highestDegree) {
                if (i == 1) {
                    return constrainRank.get(0).get(1);
                }
                count = i;
                break;
            }
        }
        // rank neighbor nodes wth largest neighbor degree sequence
        for (int i = 0; i < count; i++) {
            int temp = 0;
            for (Integer neighOfPotential : getNeighors(queryGraph, constrainRank.get(i).get(1))) {
                temp += queryGraph.getAdjacencyList(neighOfPotential).size();
            }
            constrainRank.get(i).set(0, temp);
        }
        constrainRank.sort((a, b) -> {
            if (!a.get(0).equals(b.get(0))) {
                return b.get(0) - a.get(0);
            } else {
                return b.get(1) - a.get(1);
            }
        });
        return constrainRank.get(0).get(1);
    }

    /**
     * @param targetNodes:List[int]  - the IDs of the target set of nodes
     * @param inputGraph:Graph       - the graph to be searched for motif
     * @param neighborList:List[int[ - the reference to the return list of neighbors
     * @return List[int] - modified neighborList
     */
    private List<Integer> chooseNeighborsOfRange(List<Integer> targetNodes, Graph inputGraph, List<Integer> neighborList) {
        for (Integer node : targetNodes) {
            for (Integer neighbor : getNeighors(inputGraph, node)) {
                if (!targetNodes.contains(neighbor)) {
                    neighborList.add(neighbor);
                }
            }
        }
        neighborList = new ArrayList<Integer>(new HashSet<Integer>(neighborList));
        Collections.sort(neighborList);
        return neighborList;
    }

    /**
     * Method to check if a neighbor node n of the target graph could be mapped to a node m of the query graph
     *
     * @param inputGraph:Graph       - target graph
     * @param n:int                  - ID of the node n in the target graph
     * @param partialMap:Dict[int,   int] - the current partial mapping from query graph to target graph #dit
     * @param neighborsOfM:List[int] - the list of neighbors of node m to the query graph
     * @return boolean - True if node n can be mapped to node m, otherwise false
     */
    private boolean isNeighborIncompatible(Graph inputGraph, int n, HashMap<Integer, Integer> partialMap, List<Integer> neighborsOfM) {
        for (Integer node : partialMap.keySet()) {
//            AdjacencyList neighborsOfNode = inputGraph.getAdjacencyList(partialMap.get(node));
            List<Integer> neighborsOfNode = getNeighors(inputGraph, partialMap.get(node));
            if (neighborsOfM.contains(node)) {
                if (!neighborsOfNode.contains(n)) {
                    return true;
                }
            } else {
                if (neighborsOfNode.contains(n)) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * Method to check if a mapping from node m of query graph to node n of target graph satisfy the symmetry-breaking conditions
     *
     * @param fixed:int              - the representative node from each equivalence class
     * @param nodesToCheck:List[int] - the symmetry-breaking conditions
     * @param partialMap:Dict[int,   int] - the current partial mapping from query graph to target graph
     * @param m:int                  - ID number of node m of query graph
     * @param n:int                  - ID number of node n of target graph
     * @return boolean - True if the symmetry-breaking condition is satisfied and the mapping is okay, False == mapping not okay
     */
    private boolean checkSymmetryBreak(int fixed, List<Integer> nodesToCheck, HashMap<Integer, Integer> partialMap, int m, int n) {
        if (!nodesToCheck.contains(m) || m != fixed) {
            return true;
        }
        int fixedLabel;
        if (m == fixed) {
            fixedLabel = n;
        } else {
            fixedLabel = partialMap.get(fixed);
        }

        if (m == fixed) {
            for (Integer node : nodesToCheck) {
                if (partialMap.containsKey(node)) {
                    if (partialMap.get(node) < fixedLabel) {
                        return false;
                    }
                }
            }
            return true;
        } else {
            return n >= fixedLabel;
        }
    }

    /**
     * Helper function to check if the list of keys of obj1 (D) is equal to obj2 (H)
     * Equal if all elements of object 1's keys are present in object 2,
     * and the elements don't have to be in the same order between objects
     *
     * @param obj1:List[int] - vectorList of queryGraph
     * @param obj2:List[int] - list of keys
     * @return boolean - isEqual
     */
    private boolean equalDtoH(List<Integer> obj1, List<Integer> obj2) {
        Collections.sort(obj1);
        Collections.sort(obj2);
        return obj1.equals(obj2);
    }

    private Graph createIsomorphicGraphs(Graph inputGraph, List<Integer> vertexList) {
        Graph newGraph = new Graph(true);
        Map<String, Integer> nameToIndex = new HashMap<>();
        for(Integer source: vertexList){
            CompactHashSet.Iter iter = inputGraph.getAdjacencyList(source).iterator();
            while (iter.hasNext()) {
                int target = iter.next();
                if(target < 0 && vertexList.contains(Math.abs(target))){
//                    newGraph.addEdge(Math.abs(target), source);
                    int fromIndex = newGraph.getOrCreateIndex(String.valueOf(Math.abs(target)), nameToIndex);
                    int toIndex = newGraph.getOrCreateIndex(String.valueOf(source), nameToIndex);
                    newGraph.getAdjacencyList(fromIndex).add(toIndex);
                    newGraph.getAdjacencyList(toIndex).add((-1) * fromIndex);
                }
            }
        }
        newGraph.setNameToIndexMap(nameToIndex);
        return newGraph;
    }

    private boolean isIsomorphicGraphSimilar(List<Integer> vertexList, Graph inputGraph, Graph queryGraph) {
        Graph newGraph = createIsomorphicGraphs(inputGraph, vertexList);
        List<List<Integer>> a = getFromToCount(newGraph);
        List<List<Integer>> b = getFromToCount(queryGraph);
        return a.equals(b);
    }

    /**
     * Method to find the symmetry-breaking conditions by Grochow-Kellis.
     * *****NOTE*****: should combine this with Algorithm2_Modified_For_Equivalence_Class()
     *
     * @param theMappings:List[List[int]] -
     * @param condition:Dict[int,         List[int]] - the symmetry break condition tha
     *                                    t we are testing in this iteration
     * @param equivalenceClass:           set[int] -
     * @return Dict[int, List[int]] - returns the symmetry break condition that was found
     */
    private HashMap<Integer, List<List<Integer>>> findCondition(List<List<Integer>> theMappings, HashMap<Integer, List<List<Integer>>> condition, HashSet<Integer> equivalenceClass) {
        if (theMappings.size() == 1) {
            return condition;
        }
        HashMap<Integer, HashSet<Integer>> equivalenceFilter = new LinkedHashMap<>();
        for (List<Integer> maps : theMappings) {
            for (int i = 0; i < maps.size(); i++) {
                equivalenceFilter.putIfAbsent(i, new HashSet<>());
                equivalenceFilter.get(i).add(maps.get(i));
            }
        }
        Integer filterKey = equivalenceFilter.entrySet().iterator().next().getKey();
        int maxSize = equivalenceFilter.get(filterKey).size();

        HashSet<Integer> temp = new HashSet<>();
        if (equivalenceClass.size() == 0) {
            temp = new HashSet<>(equivalenceFilter.get(filterKey));
        } else {
            // TODO question?
//            Integer classItem = equivalenceClass.iterator().next();
//            temp = new HashSet<>(equivalenceClass.get(classItem));
        }

        for (Integer entry : equivalenceFilter.keySet()) {
            if (equivalenceFilter.get(entry).size() > 1) {
                if (equivalenceFilter.get(entry).size() > maxSize) {
                    maxSize = equivalenceFilter.get(entry).size();
                    temp = new HashSet<>(equivalenceFilter.get(entry));
                }
            }
        }
        equivalenceClass.retainAll(temp);

        List<Integer> sortedTemp = new ArrayList<>(temp);
        Collections.sort(sortedTemp);
        Integer fixedNode = sortedTemp.get(0);
        if (!condition.containsKey(fixedNode)) {
            condition.putIfAbsent(fixedNode, new ArrayList<>());
        }
        condition.get(fixedNode).add(sortedTemp);

        List<List<Integer>> newMappings = new ArrayList<>();
        for (List<Integer> maps : theMappings) {
            for (int i = 0; i < maps.size(); i++) {
                if (maps.get(i).equals(fixedNode) && maps.get(i).equals(mappedHNodes.get(i))) {
                    newMappings.add(new ArrayList<>(maps));
                }
            }
        }
        findCondition(newMappings, condition, equivalenceClass);
        return condition;
    }

    /**
     * Method to count all of the isomorphic extensions (no duplicates) of a partial map between the query graph and the target graph
     * @param partialMap:Dict[int, int] - the current partial mapping from query graph to target graph #is a dictionary
     * @param queryGraph:Graph - reference to the query graph
     * @param inputGraph:Graph - reference to the target graph
     * @param symBreakCondition:Dict[int, List[List[int]] - set of symmetry-breaking conditions
     * @return int - representing the count of all the isomorphic extensions
     */
    private int isomorphicExtention(HashMap<Integer, Integer> partialMap, Graph queryGraph, Graph inputGraph, HashMap<Integer, List<List<Integer>>> symBreakCondition) {
        int listOfIsomorphisms = 0; // # tracks number of isomorphisms
        List<Integer> partialMapValuesG = new ArrayList<>();
        List<Integer> partialMapKeysH = new ArrayList<>();

        // extract list of keys and list values from paritalMap
        for (Integer maps : partialMap.keySet()) {
            partialMapValuesG.add(partialMap.get(maps));
            partialMapKeysH.add(maps);
        }

        List<Integer> mapValueOriginal = new ArrayList<>(partialMapValuesG);
        List<Integer> mapKeyOriginal = new ArrayList<>(partialMapKeysH);

        Collections.sort(partialMapValuesG);
        Collections.sort(partialMapKeysH);

        if (equalDtoH(getVertexList(queryGraph), partialMapKeysH)) {
            return 1;
        }

        int m = getMostConstrainedNeighbor(partialMapKeysH, queryGraph);
        if (m < 0) {
            return 0;
        }

        List<Integer> neighborsofM = getNeighors(queryGraph, m);
        int bestMappedNeighborOfM = 0;
        for (Integer neighbor : neighborsofM) {
            if (partialMap.containsKey(neighbor)) {
                bestMappedNeighborOfM = neighbor;
                break;
            }
        }

        List<Integer> possibleMappingNodes = new ArrayList<>();
        for (Integer node : getNeighors(inputGraph, partialMap.get(bestMappedNeighborOfM))) {
            if (!partialMapValuesG.contains(node)) {
                possibleMappingNodes.add(node);
            }
        }

        int partialMapKeysHSize = partialMapKeysH.size();
        for (int i = 0; i < partialMapKeysHSize; i++) {
            List<Integer> neighborsOfMappedGNode = getNeighors(inputGraph, mapValueOriginal.get(i));
            List<Integer> temp = new ArrayList<>();
            if (neighborsofM.contains(mapKeyOriginal.get(i))) {
                for (Integer node : possibleMappingNodes) {
                    if (neighborsOfMappedGNode.contains(node)) {
                        temp.add(node);
                    }
                }
                possibleMappingNodes = new ArrayList<>(temp);
            } else {
                for(Integer node: possibleMappingNodes){
                    if(!neighborsOfMappedGNode.contains(node)){
                        temp.add(node);
                    }
                }
                possibleMappingNodes = new ArrayList<>(temp);
            }
        }

        for(Integer n: possibleMappingNodes){
            if(!isNeighborIncompatible(inputGraph, n, partialMap, neighborsofM)){
                boolean skip = false;
                for(Integer condition: symBreakCondition.keySet()){
                    if(!checkSymmetryBreak(symBreakCondition.get(condition).get(0).get(0), symBreakCondition.get(condition).get(0), partialMap, m, n)){
                        skip = true;
                        break;
                    }
                }
                if(skip){
                    continue;
                }
                HashMap<Integer, Integer> newPartialMap = new HashMap<>(partialMap);
                newPartialMap.put(m, n);

                int subList = isomorphicExtention(newPartialMap, queryGraph, inputGraph, symBreakCondition);
                if(subList == 1){
                    if(newPartialMap.values().size() == queryGraph.getSize()){
                        if(inputGraph.getDir()){

                            if(isIsomorphicGraphSimilar(new ArrayList<>(newPartialMap.values()), inputGraph, queryGraph)){
                                // write the directed graph
                                List<String> graphlet = new ArrayList<>();
                                for(Integer node: new ArrayList<>(newPartialMap.values())){
                                    graphlet.add(mapIndexToString.get(node));
                                }
                                try {
                                    bufferedWriter.write(graphlet.toString() + "\n");
                                    bufferedWriter.flush();
                                } catch (IOException e){
                                    System.out.println(e.getMessage());
                                }

                            }else{
                                subList = 0;
                            }
                        }else{
                            // write the undirected Graph
                            List<String> graphlet = new ArrayList<>();
                            for(Integer node: new ArrayList<>(newPartialMap.values())){
                                graphlet.add(mapIndexToString.get(node));
                            }
                            try {
                                bufferedWriter.write(graphlet.toString() + "\n");
                                bufferedWriter.flush();
                            } catch (IOException e){
                                System.out.println(e.getMessage());
                            }
                        }
                    }
                }
                listOfIsomorphisms += subList;
            }
        }
        return listOfIsomorphisms;

    }

    /**
     * Helper method to find all of the isomorphic extensions of a partial map between the query graph and itself
     * @param partialMap:dict[int, int] - a partial map
     * @param queryGraph:Graph
     * @param inputGraph:Graph - same as query graph
     * @return List[List[int]] -
     */
    private List<List<Integer>> isomorphicExtensionForEquivalenceClass(HashMap<Integer, Integer> partialMap, Graph queryGraph, Graph inputGraph) {
        List<List<Integer>> result = new ArrayList<>();
        List<List<Integer>> listOfIsomorphisms = new ArrayList<>();
        List<Integer> partialMapValuesG = new ArrayList<>();
        List<Integer> partialMapKeysH = new ArrayList<>();

        // extract list of keys and list of values from partialMap
        for(Integer maps: partialMap.keySet()){
            partialMapValuesG.add(partialMap.get(maps));
            partialMapKeysH.add(maps);
        }

        List<Integer> mapValueOriginal = new ArrayList<>(partialMapValuesG);
        List<Integer> mapKeyOriginal = new ArrayList<>(partialMapKeysH);

        Collections.sort(partialMapValuesG);
        Collections.sort(partialMapKeysH);

        if(equalDtoH(getVertexList(queryGraph), partialMapKeysH)){
            mappedHNodes = new ArrayList<>(mapKeyOriginal);
            result.add(mapValueOriginal);
            return result;
        }

        int m = getMostConstrainedNeighbor(partialMapKeysH, queryGraph);
        if(m < 0){
            return listOfIsomorphisms;
        }

        List<Integer> neighbourRange = new ArrayList<>();
        neighbourRange = chooseNeighborsOfRange(partialMapValuesG, inputGraph, neighbourRange);

        List<Integer> neighborsOfM = getNeighors(queryGraph, m);
        for(Integer n: neighbourRange){
            if(!isNeighborIncompatible(inputGraph, n, partialMap, neighborsOfM)){
                HashMap<Integer, Integer> newPartialMap = new HashMap<>(partialMap);
                newPartialMap.put(m, n);
                List<List<Integer>> subList = isomorphicExtensionForEquivalenceClass(newPartialMap, queryGraph, inputGraph);

                listOfIsomorphisms.addAll(subList);
            }
        }
        return listOfIsomorphisms;
    }

    /**
     * Method to find the symmetry-breaking conditions by Grochow-Kellis. It starts by choosing one node to be the anchor point and create conditions from
     * @param queryGraph:Graph - reference to query graph
     * @param inputGraph:Graph - reference to input graph
     * @param fixedNode:int - the node we choose to be fixed as the anchor for symmetry (might not be needed??)
     * @return Dict[int, List[List[int]]] - a set of symmetry-breaking conditions for each represented node from each equivalance class
     */
    private HashMap<Integer, List<List<Integer>>> algorithm2_modified_for_equivalance_class(Graph queryGraph, Graph inputGraph, int fixedNode) {
        List<Integer> vertexList = getVertexList(queryGraph);
        Integer h = vertexList.iterator().next();

        List<Integer> inputGraphDegSeq = getNodesSortedByDegree(inputGraph, queryGraph.getAdjacencyList(h).size());
        List<List<Integer>> theMappings = new ArrayList<>();
        mappedHNodes = new ArrayList<>();

        for(Integer item: inputGraphDegSeq){
            HashMap<Integer, Integer> f = new HashMap<>();
            f.put(h, item);
            List<List<Integer>> mappings = isomorphicExtensionForEquivalenceClass(f, queryGraph, queryGraph);
            theMappings.addAll(mappings);
        }
        System.out.println("The Mappings: " + theMappings);
        HashMap<Integer, List<List<Integer>>> condition = new HashMap<>();
        HashSet<Integer> equivalenceClass = new HashSet<>();
        return findCondition(theMappings, condition, equivalenceClass);
    }

    /**
     * Method to use NemoMap algorithm (i.e. Algorithm 5 from the NemoMap paper)
     * ***Modified from Grochow-Kelis algorithm***
     * Implemented in C++ by Tien Huynh
     * For more information please see the research paper of NemoMap and/or Grochow-Kellis' paper
     * "Network Motif Discovery Using Subgraph Enumeration and Symmetry-Breaking"
     * @param queryGraph:Graph - reference to query graph
     * @param inputGraph:Graph - reference to input graph
     * @param h:int - the starting node h of query graph -
     *        (should be the most constrained node of H -> first rank by out-degree; second rank by neighbor degree sequence)
     * @return int - The count of all of possible mappings of the query graph to the target graph
     */
    public int algorithm2_modified(Graph queryGraph, Graph inputGraph, int h) {
        mapIndexToString = getMapIndexToString(inputGraph);
        HashMap<Integer, List<List<Integer>>> condition = algorithm2_modified_for_equivalance_class(queryGraph, queryGraph, h);
        List<Integer> inputGraphDegSeq = getNodesSortedByDegree(inputGraph, queryGraph.getAdjacencyList(h).size());
        int mappingCount = 0;
        HashMap<Integer, Integer> f = new HashMap<>();

        for(Integer value: inputGraphDegSeq){
            f.put(h, value);
            mappingCount += isomorphicExtention(f, queryGraph, inputGraph, condition);
        }
        try {
            bufferedWriter.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return mappingCount;
    }

    private HashMap<Integer, String> getMapIndexToString(Graph inputGraph){
        HashMap<Integer, String> indexToNameMap = new HashMap<>();
        for(Map.Entry<String, Integer> entry : ((HashMap<String, Integer>)inputGraph.getNameToIndexMap()).entrySet()){
            indexToNameMap.put(entry.getValue(), entry.getKey());
        }
        return indexToNameMap;
    }

    private List<List<Integer>> getFromToCount(Graph graph){
        List<List<Integer>> graphFromToCount = new ArrayList<>();
        for(int node = 1; node <= graph.getSize(); node++){
            graphFromToCount.add(new ArrayList<>());
            CompactHashSet.Iter iter = graph.getAdjacencyList(node).iterator();
            int fromCount = 0;
            int toCount = 0;
            while (iter.hasNext()) {
                int target = iter.next();
                if(target > 0){
                    fromCount++;
                }else{
                    toCount++;
                }
            }
            graphFromToCount.get(node-1).add(fromCount);
            graphFromToCount.get(node-1).add(toCount);
        }
        graphFromToCount.sort((a, b) -> {
            if (!a.get(0).equals(b.get(0))) {
                return a.get(0).compareTo(b.get(0));
            } else {
                return a.get(1).compareTo(b.get(1));
            }
        });
        return graphFromToCount;
    }

    public List<Integer> getNeighors(Graph graph, Integer index) {
        List<Integer> neighborList = new ArrayList<>();
        CompactHashSet.Iter iter = graph.getAdjacencyList(index).iterator();
        while (iter.hasNext()) {
            neighborList.add(Math.abs(iter.next()));
        }
        return neighborList;
    }

    public List<Integer> getVertexList(Graph graph){
        List<Integer> vertexList = new ArrayList<>();
        for (int i = 1; i <= graph.getSize(); i++) {
            vertexList.add(i);
        }
        return vertexList;
    }

    public List<Integer> getNodesSortedByDegree(Graph graph, Integer degreeCutOff){
        List<Integer> nodesSortedByDegree = new ArrayList<>();
        List<Pair> temp = new ArrayList<>();
        for(Integer vertex : getVertexList(graph)){
            int degree = graph.getAdjacencyList(vertex).size();
            if(degree >= degreeCutOff){
                temp.add(new Pair(vertex, degree));
            }
        }
        temp.sort(Comparator.comparingInt(Pair::getDegree));
        for(Pair pair: temp){
            nodesSortedByDegree.add(pair.getVertex());
        }
        return nodesSortedByDegree;
    }

    private static class Pair{
        private int vertex;
        private int degree;
        Pair(int vertex, int degree){
            this.vertex = vertex;
            this.degree = degree;
        }

        public int getVertex() {
            return vertex;
        }

        public int getDegree() {
            return degree;
        }
    }
}
