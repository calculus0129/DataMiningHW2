import java.io.*;
import java.util.*;
import java.util.HashMap;
import java.util.stream.Collectors;

public class A2_G13_t2 {
    static double eps=0.5;
    static int mu=5;
    static int dim;
    static Map<String, Point> database = new HashMap<>();

    private static class Point {
        public ArrayList<Double> coord;
        public Point(ArrayList<Double> coord) {
            this.coord = coord;
        }
        // The Euclidean Distance between two points
        public double dist(Point p) {
            double ret=0.0;
            for(int i=0;i<dim;++i) {
                // get() has O(1) complexity for ArrayList (but O(n) for LinkedList)
                // https://www.baeldung.com/java-collections-complexity#arraylist
                ret += (coord.get(i) - p.coord.get(i))*(coord.get(i) - p.coord.get(i));
            }
            return Math.sqrt(ret);
        }

        public static double dist(Point p, Point q) {
            double ret=0.0;
            for(int i=0;i<dim;++i) {
                ret += (p.coord.get(i) - q.coord.get(i))*(p.coord.get(i) - q.coord.get(i));
            }
            return Math.sqrt(ret);
        }

        @Override
        public String toString() {
            return "Point{" +
                    "coord=" + coord +
                    '}';
        }
    }

    private static ArrayList<Set<String>> dbscan(Map<String, Point> db, double eps, int mu) {
        ArrayList<Set<String>> ret = new ArrayList<>();
        Map<String, Set<String>> graph = makeGraph(db, eps);
        Set<String> corePoints = getCores(graph, mu);
        Set<String> visited = new HashSet<>(), discovered = new HashSet<>();
        LinkedList<String> stack = new LinkedList<>();
        for(String s: corePoints) if(!visited.contains(s)) {
            stack.add(s);
            visited.add(s);
            discovered.add(s);
            Set<String> cluster = new HashSet<>();
            while (!stack.isEmpty()) {
                String node = stack.removeLast();
                cluster.add(node);
                if(corePoints.contains(node)) {
                    visited.add(node);
                    for (String k : graph.get(node)) if(!discovered.contains(k)) {
                        discovered.add(k);
                        if (corePoints.contains(k)) {
                            stack.add(k);
                        }
                    }
                }
            }
            ret.add(cluster);
        }

        return ret;
    }

    private static Map<String, Set<String>> makeGraph(Map<String, Point> db, double eps) {
        Map<String, Set<String>> ret = new HashMap<>();
        for(Map.Entry<String, Point> p: db.entrySet()) for(Map.Entry<String, Point> q: db.entrySet()) if(Point.dist(p.getValue(), q.getValue())<=eps && !p.getKey().equals(q.getKey())) {
            ret.computeIfAbsent(p.getKey(), k->new HashSet<>()).add(q.getKey());
//            ret.computeIfAbsent(q.getKey(), k->new HashSet<>()).add(p.getKey());
        }
        return ret;
    }

    private static Set<String> getCores(Map<String, Set<String>> graph, int mu) {
        Set<String> ret = new HashSet<>();
        for(String s: graph.keySet()) if(graph.get(s).size()>=mu) ret.add(s);
        return ret;
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
            dim = vals.size();
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

        ArrayList<Set<String>> clusters = dbscan(database, eps, mu);

        System.out.printf("Number of clusters : %d\n", clusters.size());
        System.out.printf("Number of noise: %d%n",
                clusters.stream()
                        .mapToInt(Set::size)
                        .sum());
        for(int i=0,e=clusters.size();i<e;++i) {
            System.out.printf("Cluster #%d\t=>\t", i+1);
            for(String s: clusters.get(i)) System.out.print(s+String.valueOf(' '));
            System.out.println();
        }

//        int i=0;
//        for(Map.Entry<String, Point> entry: database.entrySet()) {
//            System.out.println(entry.getKey() + " - " + entry.getValue());
//            if(i++>=5) break;
//        }

    }
}
