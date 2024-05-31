import java.io.*;
import java.util.*;
import java.util.HashMap;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static java.util.Collections.binarySearch;

public class A2_G13_t2 {
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
    private static Map<Point, Set<Point>> makeGraph(final Collection<Point> db, double eps) {
        Map<Point, Set<Point>> ret = new HashMap<>();
        for(Point p: db) for(Point q: db) if(Point.dist(p, q)<=eps && !p.equals(q)) {
            ret.computeIfAbsent(p, k->new HashSet<>()).add(q);
        }
        return ret;
    }

    // Finding all core points based on the edge-list graph.
    private static Set<Point> getCores(final Map<Point, Set<Point>> graph, int mu) {
        Set<Point> ret = new HashSet<>();
        for(Point s: graph.keySet()) if(graph.get(s).size()>=mu-1) ret.add(s);
        return ret;
    }

    // Finding cluster of a core point.
    private static Set<Point> findCluster(final Map<Point, Set<Point>> graph, Point s, Set<Point> corePoints, Set<Point> discovered) {
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
    private static List<Double> epsEstimates(final Collection<Point> db, int mu) {
        int n = db.size(), opt=0;
        ArrayList<Double> k_dist = getKdists(db, mu); // This is equivalent to (mu-1)-dist list.
        List<Integer> opts = List.of(cand1(k_dist), cand2(k_dist), cand3(k_dist), cand4(k_dist));
        return opts.stream().map(k_dist::get).collect(Collectors.toList());
    }

    // Candidate 1: Find the index i of point pi whose 'instantaneous rate of change', pi+1.y-pi.y is most similar to the average rate of change.
    private static int cand1(ArrayList<Double> k_dist) {
        int ret=0, n = k_dist.size();
        double avgDiff = (k_dist.get(0)-k_dist.get(n-1))/(n-1), dist = Double.MAX_VALUE, diff;
        for(int i=0,e=n-1;i<e;++i) {
            diff = Math.abs(k_dist.get(i)-k_dist.get(i+1) - avgDiff);
            if(diff<dist) {
                dist = diff;
                ret = i;
            }
        }
        return ret;
    }
    // Candidate 2: Find the index i of point pi whose ratio slope(p0pi)/slope(pipn-1) is maximized!
    private static int cand2(ArrayList<Double> k_dist) {
        int ret=0, n = k_dist.size();
        double maxRatio = 0.0, ratio;
        for(int i=1,e=n-1;i<e;++i) {
            ratio = (k_dist.get(0)-k_dist.get(i))*(n-1-i)/((k_dist.get(i)-k_dist.get(n-1))*i);
            if(maxRatio<ratio) {
                maxRatio = ratio;
                ret = i;
            }
        }
        return ret;
    }
    // Candidate 3 : Find the index i of point pi whose angle pi-1pipi+1 is maximized!
    private static int cand3(ArrayList<Double> k_dist) {
        int ret=0, n = k_dist.size();
        double maxTan = 0.0, tan, a, b;
        for(int i=1,e=n-1;i<e;++i) {
            // Use the tangent subtraction law: a-b/(1+ab)
            a = k_dist.get(i-1)-k_dist.get(i);
            b = k_dist.get(i)-k_dist.get(i+1);
            tan = (a-b) / (1+a*b);
            if(maxTan<tan) {
                maxTan = tan;
                ret = i;
            }
        }
        return ret;
    }
    // Candidate 4 : Find the index i of point pi whose slope ratio slope(pi-1pi)/slope(pipi+1) is maximized!
    private static int cand4(ArrayList<Double> k_dist) {
        int ret=0, n = k_dist.size();
        double maxRatio = 0.0, ratio;
        for(int i=1,e=n-1;i<e;++i) {
            ratio = (k_dist.get(i-1)-k_dist.get(i))/(k_dist.get(i)-k_dist.get(i+1));
            if(maxRatio<ratio) {
                maxRatio = ratio;
                ret = i;
            }
        }
        return ret;
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

    // Conducting the experiment!
    public static void multiRun(String[] args) throws Exception {
        File file;
        for(int fidx=1,e=args.length;fidx<e;++fidx) {
            file = new File(args[fidx]);
            if(file.isDirectory()) {
                System.out.printf("%s is a directory!\n", file.toPath());
                List<List<Double>> accMatrix = new ArrayList<>(), precMatrix = new ArrayList<>();
                for(String fileName: file.list((f, fname)->fname.endsWith(".csv"))) {
                    String path = file.toPath().toString()+'/'+fileName;
                    boolean isExample = path.substring(path.lastIndexOf('/')+1).startsWith("art") && path.endsWith(".csv");

                    readFile(path, isExample);
                    int mu = dim<<1;
                    Set<Point> data = database.values().stream().reduce(new HashSet<>(), (a, b) -> {
                        a.addAll(b);
                        return a;
                    });

                    List<Double> epsList = epsEstimates(data, mu);
                    List<List<Integer>> confMatrix = epsList.stream()
                            .map(eps->getConfusion(database.values(), dbscan(data, eps, mu)))
                            .collect(Collectors.toList());

                    List<Double> accs = confMatrix.stream().map(confusion->{
                        // Evaluation
                        int tp = confusion.get(0), fp = confusion.get(1), tn = confusion.get(2), fn = confusion.get(3);
                        return (double)(tp+tn)/(tp+fp+tn+fn);
                    }).collect(Collectors.toList());
                    accMatrix.add(accs);

                    List<Double> precs = confMatrix.stream().map(confusion->{
                        // Evaluation
                        int tp = confusion.get(0), fp = confusion.get(1), tn = confusion.get(2), fn = confusion.get(3);
                        return (double)(tp)/(tp+fp);
                    }).collect(Collectors.toList());
                    precMatrix.add(precs);

                    double minEps = epsList.stream().min(Double::compareTo).get();
                    List<Integer> epsMinIndices = IntStream.range(0, epsList.size())
                            .filter(idx -> epsList.get(idx).equals(minEps))
                            .boxed()  // Convert from int to Integer
                            .collect(Collectors.toList());  // Correct way to collect into a list
                    List<Integer> accMaxIndices = IntStream.range(0, accs.size())
                            .filter(idx -> accs.get(idx).equals(accs.stream().max(Double::compareTo).get())) // Fixed missing parenthesis
                            .boxed()
                            .collect(Collectors.toList());
                    List<Integer> precMaxIndices = IntStream.range(0, precs.size())
                            .filter(idx -> precs.get(idx).equals(precs.stream().max(Double::compareTo).get())) // Fixed missing parenthesis
                            .boxed()
                            .collect(Collectors.toList());

//                    epsList.forEach(d -> System.out.printf("% 5.3f\t", d));
                    System.out.printf("\t%.3f\t%.3f\t%.3f\t%.3f\t%s\t%s\t\t" +
                                    "\t%.3f\t%.3f\t%.3f\t%.3f\t%s\t%s\t\t" +
                                    "%s\n",
                            accs.get(0), accs.get(1), accs.get(2), accs.get(3), accMaxIndices.toString(), accMaxIndices.equals(epsMinIndices),
                            precs.get(0), precs.get(1), precs.get(2), precs.get(3), precMaxIndices.toString(), precMaxIndices.equals(epsMinIndices),
                            epsMinIndices);
                }
                System.out.println("Accuracy:");
                printStats(accMatrix);
                System.out.println("Precision:");
                printStats(precMatrix);
            }
            else {
                System.out.printf("%s is NOT a directory!", file.toPath());
            }
        }
    }
    private static void printStats(List<List<Double>> matrix) {

        // Calculate and print the statistics for each candidate
        for (int i = 0; i < 4; i++) { // Assuming four candidates
            final int index = i;
            double mean = matrix.stream()
                    .mapToDouble(list -> list.get(index))
                    .average()
                    .orElse(0.0);
            double variance = matrix.stream()
                    .mapToDouble(list -> Math.pow(list.get(index) - mean, 2))
                    .average()
                    .orElse(0.0);
            double stdev = Math.sqrt(variance);

            System.out.printf("Candidate %d -> Mean: %.3f, Stdev: %.3f%n", i+1, mean, stdev);
        }
    }

    public static void main(String[] args) throws Exception {
        boolean isTest = Objects.equals(args[0], "-test");
        if(isTest) {
            multiRun(args);
        } else {
            singleRun(args);
        }

    }


    private static void singleRun(String[] args) throws Exception {
        String path = args[0];
        boolean isExample = path.substring(path.lastIndexOf('/')+1).startsWith("art") && path.endsWith(".csv");
        Map<Point, String> names = readFile(path, isExample);
        // This is the boolean of whether it is originated from the given example datasets.

        int mu = 0;
        double eps = 0.0;
        // If mu or something is given.
        if(args.length==2) {
//            if(Objects.equals(args[1], "test")) {
//                test();
//            }
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
            }
        }
        Set<Point> data = database.values().stream().reduce(new HashSet<>(), (a, b) -> {
            a.addAll(b);
            return a;
        });
        if(mu == 0) {
            mu = dim<<1;
            System.out.println("Estimated MinPts : "+mu);
        }
        if(eps == 0.0) {
            // Single epsEstimate value? Let's first take the minimum of the four candidates.
            eps = epsEstimates(data, mu).stream().min(Double::compareTo).get();
            System.out.println("Estimated eps: "+eps);
        }
        runAndEvaluate(data, mu, eps);
    }

    private static Map<Point, String> readFile(String path, boolean isExample) throws Exception {
        // This is the boolean of whether it is originated from the given example datasets.
        Map<Point, String> names = new HashMap<>();
        database.clear();
        // Input
        if(isExample) {
            if(getExamples(path, names)) {

//                System.out.println("Data read successfully!");
//                int i=0;
//                for(Map.Entry<Point, String> entry: names.entrySet()) {
//                    System.out.println(entry.getKey() + " - " + entry.getValue());
//                    if(i++>=5) break;
//                }

            }
            else {
                // TODO
                System.out.println("Data NOT read!");
                throw new Exception("Data NOT read!");
            }
        } else { // datasets from https://www.kaggle.com/datasets/joonasyoon/clustering-exercises
            if(getKaggle(path)) {
//                System.out.println("Data read successfully!");
            }
            else {
                // TODO
                System.out.println("Data NOT read!");
                throw new Exception("Data NOT read!");
            }
        }
        return names;
    }

    private static void runAndEvaluate(Set<Point> data, int mu, double eps) {

//        System.out.println("size of data: "+data.size());

        ArrayList<Set<Point>> clusters = dbscan(data, eps, mu);

        System.out.printf("Number of clusters : %d\n", clusters.size());
        System.out.printf("Number of noise : %d%n",
                clusters.stream()
                        .mapToInt(Set::size)
                        .sum());
//        if(isExample) {
//            ArrayList<List<String>> clusterLabels = clusters.stream()
//                    .map(cluster -> cluster.stream()
//                            .map(p -> names.get(p)) // Replace each identifier with its corresponding name
//                            .collect(Collectors.toList())) // Collect names into a List to form a transformed cluster
//                    .collect(Collectors.toCollection(ArrayList::new)); // Collect transformed clusters into an ArrayList
//
//            sortCollectionOfLists(clusterLabels);
//
//            // Display sorted clusters
//            for (int i = 0, e = clusterLabels.size(); i < e; ++i) {
//                System.out.printf("Cluster #%d\t=>\t", i + 1);
//                clusterLabels.get(i).forEach(s -> System.out.print(s + " "));
//                System.out.println();
//            }
//        }

        // Evaluation

        List<Integer> confusion = getConfusion(database.values(), clusters);
        int tp = confusion.get(0), fp = confusion.get(1), tn = confusion.get(2), fn = confusion.get(3);
        System.out.println("TP: " + confusion.get(0) + ", FP: " + confusion.get(1) + ", TN: " + confusion.get(2) + ", FN: " + confusion.get(3));
        System.out.printf("Accuracy: %d/%d (%.3f%%) \n", tp+tn, tp+fp+tn+fn, 100*(double)(tp+tn)/(tp+fp+tn+fn));
    }

    // Called for data from the given examples with additional label for each data (which is not necessary).
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


    // Called for data from the given examples with additional label for each data (which is not necessary).
    private static boolean getKaggle(String path) {
        String line;
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            line = br.readLine();
            line = br.readLine();
            int lastidx = line.lastIndexOf(',');
            ArrayList<Double> vals = Arrays.stream(line.substring(0, lastidx).split(",")).map(Double::parseDouble).collect(Collectors.toCollection(ArrayList::new));
            dim = vals.size();
            String id = line.substring(lastidx+1);
            database.computeIfAbsent(id, k->new HashSet<>()).add(new Point(vals));
            while ((line = br.readLine()) != null) {
                lastidx = line.lastIndexOf(',');
                id = line.substring(lastidx+1);
                vals = Arrays.stream(line.substring(0, lastidx).split(",")).map(Double::parseDouble).collect(Collectors.toCollection(ArrayList::new));
                database.computeIfAbsent(id, k->new HashSet<>()).add(new Point(vals));
            }
        } catch (IOException e) {
            e.printStackTrace();
            return false;
        }
        return true;
    }

    // The TP, FP, TN, FN result
    private static List<Integer> getConfusion(final Collection<Set<Point>> ans, final Collection<Set<Point>> clusters) {
        Map<Point, Integer> ansMap = buildMap(ans);
        Map<Point, Integer> clusterMap = buildMap(clusters);

        int tp = 0, fp = 0, tn = 0, fn = 0;

        List<Point> points = new ArrayList<>();
        for (Set<Point> set : ans) {
            points.addAll(set);
        }

        for (int i = 0, j, e=points.size(); i < e; i++) {
            for (j = i + 1; j < e; j++) {
                Point p = points.get(i);
                Point q = points.get(j);

                boolean sameAns = Objects.equals(ansMap.get(p), ansMap.get(q));
                boolean sameCluster = Objects.equals(clusterMap.get(p), clusterMap.get(q)); //&& clusterMap.get(p)!=null;

                if (sameCluster && sameAns) tp++;
                if (sameCluster && !sameAns) fp++;
                if (!sameCluster && !sameAns) tn++;
                if (!sameCluster && sameAns) fn++;
            }
        }

        return List.of(tp, fp, tn, fn);
    }

    // Literally, build a map that maps points to its arbitrarily given corresponding cluster id.
    private static Map<Point, Integer> buildMap(Collection<Set<Point>> collection) {
        Map<Point, Integer> pointToSetMap = new HashMap<>();
        int setId = 0;
        for (Set<Point> set : collection) {
            for (Point p : set) {
                pointToSetMap.put(p, setId);
            }
            setId++;
        }
        return pointToSetMap;
    }

    // Sort the output of the results of the given example datasets.
    private static void sortCollectionOfLists(ArrayList<List<String>> collection) {
        // Define a custom comparator that compares based on the numeric part of the strings
        Comparator<String> numericComparator = (s1, s2) -> {
            int num1 = Integer.parseInt(s1.substring(1)); // Extract number part, skip 'p'
            int num2 = Integer.parseInt(s2.substring(1)); // Extract number part, skip 'p'
            return Integer.compare(num1, num2); // Compare the numbers
        };

        // Sort each list in the main collection
        for (List<String> list : collection) {
            list.sort(numericComparator); // Sort each individual list
        }

        // Sort the collection of lists based on the lexicographic order of the first element
        collection.sort((list1, list2) -> {
            String first1 = list1.get(0); // Get first element of list1
            String first2 = list2.get(0); // Get first element of list2
            return numericComparator.compare(first1, first2); // Compare using numericComparator
        });
    }
}
