package edu.uwb.nemolib_examples.network_motif_detector;

import edu.uwb.nemolib.*;


import java.io.*;
import java.util.*;

/**
 * @Author: Zican Li
 * @Date：11/12/2020 11:33 PM
 */
public class DirectMethodBuilderTest {
//    private static final String generateAdjMatrixCommand = "nemolib-master\\src\\main\\resources\\showg.exe -a subgraphs.txt";

    public static int adjMatrixToNumber(List<String> adj){
        StringBuilder num_code = new StringBuilder();
        for(String line: adj){
            num_code.append(line);
        }
        return Integer.parseInt(num_code.toString(), 2);
    }

    public static void main(String[] args) throws IOException{

        Graph inputGraph = null;
        String inputGraphName = args[0];
        int motifSize = Integer.parseInt(args[1]);
        try{
            inputGraph = GraphParser.parse(inputGraphName, false);
        }catch (IOException e){
            System.err.println("Could not process " + inputGraphName);
            System.err.println(e.getMessage());
            System.exit(-1);
        }

        List<String> g6s = new ArrayList<>();
        List<Integer> motifs = new ArrayList<>();
        Map<String, Integer> g6_to_number = new HashMap<>();

        SubgraphCollection subgraphCount = new SubgraphCollection();
        SubgraphEnumerator targetGraphESU = new ESU();
        TargetGraphAnalyzer targetGraphAnalyzer = new TargetGraphAnalyzer(targetGraphESU, subgraphCount);
        Map<String, Double> targetLabelToRelativeFrequency =
                (HashMap<String, Double>) targetGraphAnalyzer.analyze(inputGraph, motifSize);

        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter("subgraphs.txt"));
            for (String key : targetLabelToRelativeFrequency.keySet()) {
                writer.write(key);
                g6s.add(key);
                writer.write('\n');
            }
            writer.close();
        }catch (IOException e){
            System.err.println(e.getMessage());
        }


        Runtime rt = Runtime.getRuntime();
        Process proc = rt.exec("nemolib-master\\src\\main\\resources\\showg.exe -a subgraphs.txt");

        BufferedReader stdInput = new BufferedReader(new
                InputStreamReader(proc.getInputStream()));

        // Read the output from the command
        List<String> adjMatrixs = new ArrayList<>();
        String s = null;
        while ((s = stdInput.readLine()) != null) {
            if(!s.equals("")){
                adjMatrixs.add(s);
            }
        }

        int k = 0;
        for(String g6: g6s){
            List<String> adj = new ArrayList<>();
            for(int i = 0; i <= motifSize; i++){
                if(adjMatrixs.get(k * (motifSize + 1) + i).startsWith("Graph")){
                    continue;
                }
                adj.add(adjMatrixs.get(k * (motifSize + 1) + i));
            }
            k++;
            int number_code = adjMatrixToNumber(adj);
            g6_to_number.put(g6, number_code);
            motifs.add(number_code);
        }

        long[] degree_r = DirectMethodBuilder.getDegrees(inputGraph);
        long[] degree_c = DirectMethodBuilder.getDegrees(inputGraph);

        int num_samples = 100000;

        if(inputGraph.getDir()){
            DirectMethodBuilder directMethodBuilder = new DirectMethodBuilder(inputGraph.getSize(), degree_r, degree_c, new Random());

        }else{
            DirectMethodBuilder directMethodBuilder = new DirectMethodBuilder(inputGraph.getSize(), degree_r, new Random());
            double[] res = new double[motifs.size()];
            for(int i = 0; i != motifs.size(); i++){
                System.out.println(   "Motif " + (i + 1) + " of " + motifs.size());
                res[i] = directMethodBuilder.motif_sample_log(motifs.get(i), (short) motifSize,  num_samples);
            }
//            for(int i = 0; i != motifs.size(); i++){
//                System.out.println("res[i]: " + res[i]);
//            }

            double max = -50000000.0;
            for (int i = 0; i!= motifs.size(); ++i) {
                if (res[i] > max && !(res[i] != res[i]))  //test for NaN!
                { max = res[i];	}
            }
            double sum = 0.0;
            for (int i = 0; i!= motifs.size(); ++i) {
                if (!(res[i] != res[i]))  //test for NaN!
                {
                    //cout << res[i];
                    res[i] = Math.exp(res[i] - max);
                    sum += res[i];
                }
            }

            for (int i = 0; i!= motifs.size(); ++i) {
                System.out.print("   Motif " + motifs.get(i) + " :");
                if (!(res[i] != res[i]))
                    System.out.println(res[i] / sum);
                else
                    System.out.println("Not found");
            }
        }
    }
}
