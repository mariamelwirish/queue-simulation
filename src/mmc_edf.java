import java.util.*;
import java.io.*;

/**
 * M/M/c queue simulator with deadlines and non-preemptive EDF service:
 * - Exponential interarrival times (rate lambda)
 * - Exponential service times (rate mu)
 * - Exponential deadlines/patience (rate DEADLINE_RATE)
 * - Infinite buffer (no blocking)
 * - Non-preemptive EDF discipline on the waiting queue:
 *      * Customers that find all servers busy join a priority queue ordered
 *        by earliest absolute deadlineTime.
 *      * When a server becomes idle, it takes the waiting customer with the
 *        smallest deadlineTime (Earliest Deadline First).
 *
 * For each rho in [RHO_START, RHO_END] with step RHO_STEP, it:
 *  1. Sets lambda = rho * c * mu
 *  2. Runs a single simulation until RUN_COMPLETED customers depart
 *  3. Estimates Ns, Nq, Ds, Dq, utilization
 *  4. Estimates pExpired = P(served customer already expired at service start)
 *  5. Writes one summary row to CSV and prints a console summary
 *
 * No closed-form theory is computed here (non-FIFO discipline).
 */
public class mmc_edf {

    // CONTROL PARAMETERS

    // Number of completed customers per rho value
    static final int RUN_COMPLETED = 100000;

    // rho sweep range
    static final double RHO_START = 0.05;
    static final double RHO_END   = 0.95;
    static final double RHO_STEP  = 0.05;

    // Service rate and number of servers
    static final double MU = 1;
    static final int C = 1;

    // Deadline (patience) rate
    static final double DEADLINE_RATE = 0.02;

    // Output CSV path
    static final String CSV_FILE = "results/csv/mmc_edf.csv";

    // RANDOM NUMBER GENERATOR

    static final Random RNG = new Random(12345);

    // EVENT TYPES

    private static enum EventType {
        ARRIVAL,
        DEPARTURE
    }

    // EVENT CLASS

    /**
     * Generic event:
     *  - time: when the event occurs
     *  - type: ARRIVAL or DEPARTURE
     *  - serverId: which server (for departures), or -1 for arrivals
     */
    private static class Event implements Comparable<Event> {
        double time;
        EventType type;
        int serverId;

        Event(double time, EventType type, int serverId) {
            this.time = time;
            this.type = type;
            this.serverId = serverId;
        }

        @Override
        public int compareTo(Event other) {
            return Double.compare(this.time, other.time);
        }
    }

    // SERVICE DEPARTURE EVENT

    /**
     * DEPARTURE event that remembers arrivalTime of the served customer
     * so we can compute system and queueing delays.
     */
    private static class ServiceEvent extends Event {
        double arrivalTime;

        ServiceEvent(double time, int serverId, double arrivalTime) {
            super(time, EventType.DEPARTURE, serverId);
            this.arrivalTime = arrivalTime;
        }
    }

    // SERVER CLASS

    /**
     * Represents a single server:
     *   - busy = true if currently serving a customer
     */
    private static class Server {
        boolean busy = false;
    }

    // CUSTOMER CLASS

    /**
     * Customer waiting in the EDF queue:
     *  - arrivalTime: arrival instant
     *  - deadlineTime: absolute expiry time (arrival + Exp(DEADLINE_RATE))
     */
    private static class Customer {
        double arrivalTime;
        double deadlineTime;

        Customer(double arrivalTime, double deadlineTime) {
            this.arrivalTime  = arrivalTime;
            this.deadlineTime = deadlineTime;
        }
    }

    // SCHEDULING COMPARATOR FOR EDF

    /**
     * Orders customers by earliest deadlineTime (EDF).
     * Ties broken by earlier arrivalTime.
     */
    private static final Comparator<Customer> EDF_COMPARATOR =
            Comparator.<Customer>comparingDouble(c -> c.deadlineTime)
                    .thenComparingDouble(c -> c.arrivalTime);

    // SUMMARY CLASS

    /**
     * Holds per-rho statistics (simulation only):
     *  - avgNs, avgNq, avgDs, avgDq
     *  - utilization
     *  - completed
     *  - pExpired: probability that a served customer was expired at service start
     */
    private static class Summary {
        double avgNs, avgNq, avgDs, avgDq;
        double utilization;
        int    completed;
        double pExpired;
    }

    // SUMMARY ROW FOR CSV

    /**
     * Row written to CSV for each rho.
     */
    private static class SummaryRow {
        double rho, lambda, mu;
        int    c;

        double Ns_sim;
        double Nq_sim;
        double Ds_sim;
        double Dq_sim;

        double utilization_sim;

        double pExpired_sim;
    }

    // EXPONENTIAL SAMPLING

    /**
     * Exponential sample with rate = rate using inverse CDF.
     */
    private static double expSample(double rate) {
        double u = RNG.nextDouble();
        if (u == 0.0) u = 1e-16;
        return -Math.log(1.0 - u) / rate;
    }

    // FIND IDLE SERVER

    private static int findIdleServer(Server[] servers) {
        for (int i = 0; i < servers.length; i++) {
            if (!servers[i].busy) return i;
        }
        return -1;
    }

    // SINGLE RUN SIMULATOR

    /**
     * Runs one simulation with EDF scheduling in the waiting queue:
     *  - Input: lambda, mu, c, deadlineRate
     *  - Uses event-driven simulation
     *  - Returns Summary with Ns, Nq, Ds, Dq, utilization, pExpired
     */
    private static Summary simulate(double lambda,
                                    double mu,
                                    int c,
                                    double deadlineRate) {

        PriorityQueue<Event> fel = new PriorityQueue<>();
        Queue<Customer> edfQ = new PriorityQueue<>(EDF_COMPARATOR);
        Server[] servers = new Server[c];
        for (int i = 0; i < c; i++) servers[i] = new Server();

        double t = 0.0;
        double lastT = 0.0;

        int Ns = 0;
        int Nq = 0;

        double areaNs = 0.0;
        double areaNq = 0.0;
        double[] serverBusyArea = new double[c];

        double totalSystemDelay = 0.0;
        double totalQueueDelay  = 0.0;
        int    completed        = 0;

        int expiredServed = 0;

        // First arrival at t = 0
        fel.add(new Event(0.0, EventType.ARRIVAL, -1));

        while (completed < RUN_COMPLETED) {

            Event e = fel.poll();
            double oldT = t;
            t = e.time;

            double dt = t - oldT;
            areaNs += Ns * dt;
            areaNq += Nq * dt;
            for (int i = 0; i < c; i++) {
                if (servers[i].busy) serverBusyArea[i] += dt;
            }

            if (e.type == EventType.ARRIVAL) {
                // New arrival
                Ns++;

                // Next arrival
                double inter = expSample(lambda);
                fel.add(new Event(t + inter, EventType.ARRIVAL, -1));

                // Deadline
                double deadlineTime = t + expSample(deadlineRate);
                Customer cust = new Customer(t, deadlineTime);

                int idle = findIdleServer(servers);
                if (idle >= 0) {
                    // Immediate service
                    servers[idle].busy = true;

                    if (t >= cust.deadlineTime) {
                        expiredServed++;
                    }

                    double serviceTime = 1;
                    double departTime  = t + serviceTime;
                    fel.add(new ServiceEvent(departTime, idle, cust.arrivalTime));
                } else {
                    // Join EDF queue
                    edfQ.add(cust);
                    Nq++;
                }

            } else { // DEPARTURE
                ServiceEvent se = (ServiceEvent) e;
                int sid = se.serverId;

                completed++;
                Ns--;
                servers[sid].busy = false;

                double Ds_i = t - se.arrivalTime;
                totalSystemDelay += Ds_i;

                if (!edfQ.isEmpty()) {
                    Customer next = edfQ.poll();
                    Nq--;

                    double serviceStart = t;

                    if (serviceStart >= next.deadlineTime) {
                        expiredServed++;
                    }

                    double Dq_i = serviceStart - next.arrivalTime;
                    totalQueueDelay += Dq_i;

                    servers[sid].busy = true;
                    double serviceTime = 1;
                    double departTime  = serviceStart + serviceTime;
                    fel.add(new ServiceEvent(departTime, sid, next.arrivalTime));
                }
            }
        }

        Summary S = new Summary();
        S.completed = completed;

        double totalTime = t;
        S.avgNs = areaNs / totalTime;
        S.avgNq = areaNq / totalTime;

        double busySum = 0.0;
        for (int i = 0; i < c; i++) busySum += serverBusyArea[i];
        S.utilization = busySum / (c * totalTime);

        S.avgDs = totalSystemDelay / completed;
        S.avgDq = totalQueueDelay  / completed;

        S.pExpired = expiredServed / (double) completed;

        return S;
    }

    // CSV WRITER

    private static void writeCSV(List<SummaryRow> rows) {
        try {
            File f = new File(CSV_FILE);
            f.getParentFile().mkdirs();

            PrintWriter pw = new PrintWriter(new FileWriter(f));

            pw.println("rho,lambda,mu,c,"
                    + "Ns_sim,"
                    + "Nq_sim,"
                    + "Ds_sim,"
                    + "Dq_sim,"
                    + "utilization_sim,"
                    + "pExpired_sim");

            for (SummaryRow r : rows) {
                pw.printf(Locale.US,
                        "%.4f,%.6f,%.6f,%d,"
                                + "%.6f,"
                                + "%.6f,"
                                + "%.6f,"
                                + "%.6f,"
                                + "%.6f,"
                                + "%.6f%n",
                        r.rho, r.lambda, r.mu, r.c,
                        r.Ns_sim,
                        r.Nq_sim,
                        r.Ds_sim,
                        r.Dq_sim,
                        r.utilization_sim,
                        r.pExpired_sim);
            }

            pw.close();
        } catch (IOException e) {
            System.err.println("Error writing CSV: " + e.getMessage());
        }
    }

    // SUMMARY PRINTER

    private static void printSummary(List<SummaryRow> rows) {
        System.out.println("==============================================================");
        System.out.println(" c | rho |   Ns   |   Nq   |   Ds   |   Dq   | util | pExpired");
        System.out.println("==============================================================");
        for (SummaryRow r : rows) {
            System.out.printf(Locale.US,
                    "%2d | %.2f | %6.4f | %6.4f | %6.4f | %6.4f | %5.3f | %7.4f%n",
                    r.c, r.rho,
                    r.Ns_sim,
                    r.Nq_sim,
                    r.Ds_sim,
                    r.Dq_sim,
                    r.utilization_sim,
                    r.pExpired_sim);
        }
        System.out.println("==============================================================");
    }

    // MAIN DRIVER

    public static void main(String[] args) {

        List<SummaryRow> rows = new ArrayList<>();

        for (double rho = RHO_START; rho <= RHO_END + 1e-9; rho += RHO_STEP) {
            double lambda = rho * C * MU;

            Summary S = simulate(lambda, MU, C, DEADLINE_RATE);

            SummaryRow R = new SummaryRow();
            R.rho = rho;
            R.lambda = lambda;
            R.mu = MU;
            R.c = C;

            R.Ns_sim = S.avgNs;
            R.Nq_sim = S.avgNq;
            R.Ds_sim = S.avgDs;
            R.Dq_sim = S.avgDq;

            R.utilization_sim = S.utilization;
            R.pExpired_sim    = S.pExpired;

            rows.add(R);

            System.out.printf(Locale.US,
                    "EDF: rho=%.2f -> Ns=%.4f, pExpired=%.6f%n",
                    rho, S.avgNs, S.pExpired);
        }

        writeCSV(rows);
        printSummary(rows);
    }
}
