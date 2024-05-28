import java.io.*;
import java.util.*;
import java.util.HashMap;
import java.util.stream.Collectors;

public class A2_G13_t2 {
    static double eps;
    static int mu;
    static int dim;
    static Map<String, Point> database = new HashMap<>();

    private static class Point {
        public ArrayList<Double> coord;
        public Point(ArrayList<Double> coord) {
            this.coord = coord;
        }
        // The Euclidean Distance between a point
        public double dist(Point p) {
            double ret=0.0;
            for(int i=0;i<dim;++i) {
                // get() has O(1) complexity for ArrayList (but O(n) for LinkedList)
                // https://www.baeldung.com/java-collections-complexity#arraylist
                ret += coord.get(i) * p.coord.get(i);
            }
            return ret;
        }

        @Override
        public String toString() {
            return "Point{" +
                    "coord=" + coord +
                    '}';
        }
    }

    public static void manual(String msg) {
        System.out.println("Usage: java <bytecode file> <csv file> [mu] [epsilon]");
        System.out.println("Message: $"+msg);
    }
    public static void main(String[] args) {
//        System.out.println(args.length);
//        for(String s: args) System.out.println(s); // ./resources/artd-31.csv
        String path = args[0], line;
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            line = br.readLine();
            int idx = line.indexOf(',');
            String title=line.substring(0, idx);
            ArrayList<Double> vals = Arrays.stream(line.substring(idx+1).split(",")).map(Double::parseDouble).collect(Collectors.toCollection(ArrayList::new));
            database.put(title, new Point(vals));
            while ((line = br.readLine()) != null) {
                idx = line.indexOf(',');
                title=line.substring(0, idx);
                vals = Arrays.stream(line.substring(idx+1).split(",")).map(Double::parseDouble).collect(Collectors.toCollection(ArrayList::new));
                database.put(title, new Point(vals));
            }
        } catch (IOException e) {
            e.printStackTrace();
            return;
        }

//        int i=0;
//        for(Map.Entry<String, Point> entry: database.entrySet()) {
//            System.out.println(entry.getKey() + " - " + entry.getValue());
//            if(i++>=5) break;
//        }
        
    }
}
