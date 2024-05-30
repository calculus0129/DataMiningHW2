import java.io.*;
import java.util.*;
import java.util.HashMap;
import java.util.stream.Collectors;
import static java.util.Collections.binarySearch;

public class A2_G13_t2 {
    static double eps;
    static int mu;
    static int dim;
    // key representing the 'actual' cluster, value representing the set of points.
    static Map<String, Set<Point>> database = new HashMap<>();

    private static class Point {
        public ArrayList<Double> coord;
        public Point(ArrayList<Double> coord) {
            this.coord = coord;
        }
        // The Euclidean Distance between two points
        final public double dist(final Point p) {
            double ret=0.0;
            for(int i=0;i<dim;++i) {
                // get() has O(1) complexity for ArrayList (but O(n) for LinkedList)
                // https://www.baeldung.com/java-collections-complexity#arraylist
                ret += (coord.get(i) - p.coord.get(i))*(coord.get(i) - p.coord.get(i));
            }
            return Math.sqrt(ret);
        }

        public static double dist(final Point p, final Point q) {
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

        @Override
        public boolean equals(Object obj) {
            var o = (Point) obj;
            if(coord.size() != o.coord.size()) return false;
            for(int i=0,e=coord.size();i<e;++i) if(coord.get(i) != o.coord.get(i)) return false;
            return true;
        }

        // necessary in practice to use the HashMap
        @Override
        public int hashCode() {
            return Objects.hash(coord);
        }
    }

    // The plain dbscan. This returns the ArrayList of cluster sets of points. Noise points are excluded from their union.
    private static ArrayList<Set<Point>> dbscan(Collection<Point> db, double eps, int mu) {
        ArrayList<Set<Point>> ret = new ArrayList<>();
        Map<Point, Set<Point>> graph = makeGraph(db, eps);
        Set<Point> corePoints = getCores(graph, mu);
        Set<Point> discovered = new HashSet<>();

        for(Point s: corePoints) if(!discovered.contains(s)) {
            discovered.add(s);
            ret.add(findCluster(graph, s, corePoints, discovered));
        }

        return ret;
    }

    // Creating an edge-list graph of points of the dataset in which two points are incident iff their distance <= eps.
    private static Map<Point, Set<Point>> makeGraph(Collection<Point> db, double eps) {
        Map<Point, Set<Point>> ret = new HashMap<>();
        for(Point p: db) for(Point q: db) if(Point.dist(p, q)<=eps && !p.equals(q)) {
            ret.computeIfAbsent(p, k->new HashSet<>()).add(q);
//            ret.computeIfAbsent(q.getKey(), k->new HashSet<>()).add(p.getKey());
        }
        return ret;
    }

    // Finding all core points based on the edge-list graph.
    private static Set<Point> getCores(Map<Point, Set<Point>> graph, int mu) {
        Set<Point> ret = new HashSet<>();
        for(Point s: graph.keySet()) if(graph.get(s).size()>=mu-1) ret.add(s);
        return ret;
    }

    // Finding cluster of a core point.
    private static Set<Point> findCluster(Map<Point, Set<Point>> graph, Point s, Set<Point> corePoints, Set<Point> discovered) {
        LinkedList<Point> stack = new LinkedList<>();
        stack.push(s);
        Set<Point> cluster = new HashSet<>();
        while (!stack.isEmpty()) {
            Point node = stack.pop();
            cluster.add(node);
            if(corePoints.contains(node)) {
                for (Point nb : graph.get(node)) if(!discovered.contains(nb)) {
                    discovered.add(nb);
                    stack.push(nb);
                }
            }
        }
        return cluster;
    }

    // Simple manual printing function
    private static void manual(String msg) {
        System.out.println("Usage: java <bytecode file> <csv file> [mu] [epsilon]");
        System.out.println("Message: $"+msg);
    }

    // returns the optimal epsilon estimated, given data points and mu value.
    private static double epsEstimate(final Collection<Point> db, int mu) {
        int n = db.size(), opt=0;
        ArrayList<Double> k_dist = getKdists(db, mu); // This is equivalent to (mu-1)-dist list.
        double avgDiff = (k_dist.get(0)-k_dist.get(n-1))/(n-1), dist = Double.MAX_VALUE, diff;
        for(int i=0,e=n-1;i<e;++i) {
            diff = Math.abs(k_dist.get(i)-k_dist.get(i+1) - avgDiff);
            if(diff<dist) {
                dist = diff;
                opt = i;
            }
        }
        return k_dist.get(opt);
    }

    // returns the k-dist values sorted in monotonically decreasing order.
    private static ArrayList<Double> getKdists(final Collection<Point> db, int k) {
        ArrayList<Double> ret = new ArrayList<>();
        double kd;
        int idx;
        for(Point p: db) {
            kd = kDist(db, p, k);
            idx = binarySearch(ret, kd, Comparator.reverseOrder());
            ret.add(idx<0?~idx:idx, kd);
        }
        return ret;
    }

    // returns the k-dist value of a point p. Here, we consider the points including p itself as possible neighbors, so this is, 'technically', 'k-1'-dist.
    private static double kDist(final Collection<Point> db, Point p, int k) {
        ArrayList<Double> dists = new ArrayList<>();
        int cnt=0, idx;
        double d;
        for(Point q: db) {
            d = Point.dist(p, q);
            idx = binarySearch(dists, d);
            dists.add(idx<0?~idx:idx, d);
            if(cnt++>=k) {
                dists.remove(k);
            }
        }
        return dists.get(k-1);
    }

    public static void main(String[] args) {
//        System.out.println(args.length);
//        for(String s: args) System.out.println(s); // ./resources/artd-31.csv
        String path = args[0];
        Map<Point, String> names = new HashMap<>();
        // This is the boolean of whether it is originated from the given example datasets.
        boolean isExample = args[0].substring(args[0].lastIndexOf('/')+1).startsWith("art") && args[0].endsWith(".csv");
        // Input
        if(isExample) {
            if(getExamples(path, names)) {

                System.out.println("Data read successfully!");
//                int i=0;
//                for(Map.Entry<Point, String> entry: names.entrySet()) {
//                    System.out.println(entry.getKey() + " - " + entry.getValue());
//                    if(i++>=5) break;
//                }

            }
            else {
                // TODO
                return;
            }
        }

        // If mu or something is given.
        if(args.length==2) {
            try {
                mu = Integer.parseInt(args[1]);
            }
            catch(NumberFormatException e) {
                eps = Double.parseDouble(args[1]);
            }
        }
        if(args.length >= 3) {
            try {
                mu = Integer.parseInt(args[1]);
                eps = Double.parseDouble(args[2]);
            }
            catch(NumberFormatException e) {
                eps = Double.parseDouble(args[1]);
                mu = Integer.parseInt(args[2]);
                return;
            }
        }
        Set<Point> data = database.values().stream().reduce(new HashSet<>(), (a, b) -> {
            a.addAll(b);
            return a;
        });
        if(mu == 0) {
            mu = dim<<1;
            System.out.println("Estimated MinPts: "+mu);
        }
        if(eps == 0.0) {
            eps = epsEstimate(data, mu);
            System.out.println("Estimated eps: "+eps);
        }
//        System.out.println("size of data: "+data.size());

        ArrayList<Set<Point>> clusters = dbscan(data, eps, mu);

        System.out.printf("Number of clusters : %d\n", clusters.size());
        System.out.printf("Number of noise: %d%n",
                clusters.stream()
                        .mapToInt(Set::size)
                        .sum());
        if(isExample) {
            for(int i=0,e=clusters.size();i<e;++i) {
                System.out.printf("Cluster #%d\t=>\t", i+1);
                for(Point s: clusters.get(i)) System.out.print(names.get(s)+String.valueOf(' '));
                System.out.println();
            }
        }
    }

    // Called for data from the examples, etc.
    private static boolean getExamples(String path, Map<Point, String> names) {
        String line;
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            line = br.readLine();
            int idx = line.indexOf(','), lastidx = line.lastIndexOf(',');
            String name = line.substring(0, idx);
            ArrayList<Double> vals = Arrays.stream(line.substring(idx + 1, lastidx).split(",")).map(Double::parseDouble).collect(Collectors.toCollection(ArrayList::new));
            dim = vals.size();
            String id = line.substring(lastidx+1);
            database.computeIfAbsent(id, k->new HashSet<>()).add(new Point(vals));
            names.put(new Point(vals), name);
            while ((line = br.readLine()) != null) {
                idx = line.indexOf(',');
                lastidx = line.lastIndexOf(',');
                name = line.substring(0, idx);
                id = line.substring(lastidx+1);
                vals = Arrays.stream(line.substring(idx + 1, lastidx).split(",")).map(Double::parseDouble).collect(Collectors.toCollection(ArrayList::new));
                database.computeIfAbsent(id, k->new HashSet<>()).add(new Point(vals));
                names.put(new Point(vals), name);
            }
        } catch (IOException e) {
            e.printStackTrace();
            return false;
        }
        return true;
    }
}
