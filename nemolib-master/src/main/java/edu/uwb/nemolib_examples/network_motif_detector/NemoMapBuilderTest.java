package edu.uwb.nemolib_examples.network_motif_detector;

import edu.uwb.nemolib.Graph;
import edu.uwb.nemolib.GraphParser;
import edu.uwb.nemolib.NemoMapBuilder;

import java.io.IOException;
import java.util.List;

/**
 * @Author: Zican Li
 * @Date：10/29/2020 5:32 PM
 */
public class NemoMapBuilderTest {

    public static void main(String[] args) {
        Graph inputGraph = null;
        Graph queryGraph = null;

        String inputGraphName = args[0];
        String queryGraphName = args[1];
        System.out.println("inputGraphName = " + args[0]);
        System.out.println("queryGraphName = " + args[1]);

        try{
            inputGraph = GraphParser.parse(inputGraphName, false);
            queryGraph = GraphParser.parse(queryGraphName, false);
        }catch (IOException e){
            System.err.println("Could not process " + inputGraphName);
            System.err.println(e);
            System.exit(-1);
        }

        System.out.println("Print vertex to index map =" + inputGraph.getNameToIndexMap());
        System.out.println("Print vertex to index map =" + queryGraph.getNameToIndexMap());

        NemoMapBuilder nemoMapBuilder = new NemoMapBuilder("NemoCollection.txt");
        List<Integer> h = nemoMapBuilder.getNodesSortedByDegree(queryGraph, 0);
        int h1 = h.get(h.size() - 1);
        int total_mappings = nemoMapBuilder.algorithm2_modified(queryGraph, inputGraph, h1);
        System.out.println("Number of Mappings: " + total_mappings);
    }
}
