import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class A1_G13_t1 {
    /* DATA STRUCTURES */
    // point class for representing 2d points
    public static class point {
        public int label = 0;
        public double x;
        public double y;
        // constructor
        public point(double x, double y) {
            this.x = x;
            this.y = y;
        }
        public point(point other) {
            this.x = other.x;
            this.y = other.y;
        }
        // accessor; returns squared distance with other point instance
        public double dist_sq(point other) {
            double x_diff = other.x - this.x;
            double y_diff = other.y - this.y;
            return x_diff*x_diff + y_diff*y_diff;
        }
        // accessor; compare x, y value with other point instance
        public boolean compare(point other) {
            double x_diff = this.x - other.x;
            if (x_diff < 0) x_diff = -x_diff;
            double y_diff = this.y - other.y;
            if (y_diff < 0) y_diff = -y_diff;
            // returns true if difference is low enough
            return x_diff < 0.001 && y_diff < 0.001;
        }
        // modifier; change x, y value same as other point instance
        public void change(point other) {
            this.x = other.x;
            this.y = other.y;
        }
    }
    // min class
    public static class min {
        public double value;
        public int index;
        public min (double v, int i) {
            this.value = v;
            this.index = i;
        }
    }
    // cluster class for final output
    public static class cluster {
        private List<String> pointIds;
        public cluster() {
            pointIds = new ArrayList<>();      
        }
        public void add(String p) {
            pointIds.add(p);
        }
        public void print() {
            int n = pointIds.size();
            for (int i = 0; i < n; i++) {
                System.out.print(pointIds.get(i) + " ");
            }
            System.out.println();
        }
    }
    // compare two point lists
    public static boolean compare (point[] p1, point[] p2) {
        for (int i = 0; i < p1.length; i++) {
            if (!p1[i].compare(p2[i])) return false;
        }
        return true;
    }
    // copy point list to destination point list
    public static void copylist (point[] dst, point[] src) {
        for (int i = 0; i < dst.length; i++) {
            dst[i].change(src[i]);
        }
    }
    /* ALGORITHM */
    // SHORTEST-DISTANCE
    public static min SHORTEST_DISTANCE (point p, point[] clusters, int size) {
        double min = p.dist_sq(clusters[0]);
        int idx = 0;
        for (int i = 1; i < size; i++) {
            double dist = p.dist_sq(clusters[i]);
            if (dist < min) {
                min = dist;
                idx = i;
            }
        }
        return new min(min, idx);
    }
    // INITIALIZE-K-CENTERS
    public static point[] INITIALIZE_K_CENTERS (List<point> points, int k, Random rand) {
        point[] clusters = new point[k];
        int n = points.size();
        // initialize first center of cluster as randomly chosen data point
        clusters[0] = new point(points.get(rand.nextInt(n)));
        // initialize the other centers 
        for (int i = 1; i < k; i++) {
            // cumulative distribution
            double[] P = new double[n];
            double sum = 0;
            for (int j = 0; j < n; j++) {
                sum += SHORTEST_DISTANCE(points.get(j), clusters, i).value;
                P[j] = sum;
            }
            // choose random number between 0 and sum
            double random = rand.nextDouble()*sum;
            for (int j = 0; j < n; j++) {
                if (P[j] >= random) {
                    clusters[i] = new point(points.get(j));
                    break;
                }
            }
        }
        return clusters;
    }
    // K-MEANS++
    public static point[] K_MEANS_PlusPlus (List<point> points, int k, Random rand) {
        int n = points.size();
        point[] clusters = INITIALIZE_K_CENTERS(points, k, rand);
        point[] prev = new point[k];
        for (int i = 0; i < prev.length; i++) {
            prev[i] = new point(clusters[i]);
        }
        do {
            // save cluster(centers) state
            copylist(prev, clusters);
            // for calculating C.O.M(center of mass)
            int[] num = new int[k];
            double[] x_sum = new double[k];
            double[] y_sum = new double[k];
            // assign each points to cluster
            for (int i = 0; i < n; i++) {
                point cur = points.get(i);
                int label = SHORTEST_DISTANCE(points.get(i), clusters, k).index;
                cur.label = label;
                x_sum[label] += cur.x;
                y_sum[label] += cur.y;
                num[label] += 1;
            }
            // update cluster k centers with C.O.M of each cluster
            for (int i = 0; i < k; i++) {
                clusters[i].x = x_sum[i]/num[i];
                clusters[i].y = y_sum[i]/num[i];
            }
        } while (compare(prev, clusters));
        return clusters;
    }

    public static void main(String[] args) throws IOException {
        int k = 1;
        boolean estimateK = false;
        if (args.length == 2) {
            k = Integer.parseInt(args[1]);
        } else estimateK = true;
        BufferedReader file = new BufferedReader(new FileReader(args[0]));
        // n: number of items in the file
        int n = 0;
        String temp;
        // answer labels
        List<String> ids = new ArrayList<>();
        List<Integer> labels = new ArrayList<>();
        List<point> points = new ArrayList<>();
        // initialize dataset of n points
        while ((temp = file.readLine()) != null) {
            n++;
            String[] line = temp.split(",");
            ids.add(line[0]);
            labels.add(Integer.parseInt(line[3]));
            points.add(new point(Double.parseDouble(line[1]), Double.parseDouble(line[2])));
        }
        file.close();
        // Random variable
        Random rand = new Random(42);
        // estimate k
        if (estimateK) {
            // k = 
            System.out.println("estimated k: " + k);
        }
        // clustering
        K_MEANS_PlusPlus(points, k, rand);

        // gather clusters
        cluster[] clusters = new cluster[k];
        for (int i = 0; i < k; i++) {
            clusters[i] = new cluster();
        }
        for (int i = 0; i < n; i++) {
            clusters[points.get(i).label].add(ids.get(i)); 
        }
        // print out clusters
        for (int i = 0; i < k; i++) {
            System.out.print("Cluster #" + (i+1) + " => ");
            clusters[i].print();
            System.out.println();
        }
        // save results to csv
        BufferedWriter writer = new BufferedWriter(new FileWriter("res.csv"));
        writer.write("x,y,label\n");
        for (int i = 0; i < n; i++) {
            point cur = points.get(i);
            writer.write(cur.x + "," + cur.y + "," + cur.label + "\n");
        }
        writer.close();
        // calculate accuracy 0.9860236 0.9561234
        // number of TP, TN / FP, FN
        int[] acc = {0,0};
        for (int i = 0; i < n-1; i++) {
            int exp_label = points.get(i).label;
            int ans_label = labels.get(i);
            for (int j = i+1; j < n; j++) {
                int exp_label_c = points.get(j).label;
                int ans_label_c = labels.get(j);
                boolean exp_comp = exp_label == exp_label_c;
                boolean ans_comp = ans_label == ans_label_c;
                if (exp_comp == ans_comp) acc[0] += 1; // TP or TN
                else acc[1] += 1; // FP or FN
            }
        }
        float accuracy = (float) acc[0] / (acc[0]+acc[1]);
        System.out.println("accuracy: " + accuracy);
    }
}