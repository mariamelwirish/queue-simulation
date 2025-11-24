import java.util.*;
import java.io.*;

/**
 * M/M/c queue simulator with deadlines and FIFO service:
 * - Exponential interarrival times (rate lambda)
 * - Exponential service times (rate mu)
 * - Exponential deadlines/patience (rate DEADLINE_RATE)
 * - Infinite buffer (no blocking)
 * - FIFO service
 *
 * For each rho in [RHO_START, RHO_END] with step RHO_STEP, it:
 *  1. Sets lambda = rho * c * mu
 *  2. Runs a single simulation until RUN_COMPLETED customers depart
 *  3. Estimates Ns, Nq, Ds, Dq
 *  4. Computes theoretical M/M/c (Erlang-C) values (no deadlines)
 *  5. Estimates pExpired = P(served customer already expired at service start)
 *  6. Writes one summary row to CSV
 *  7. Prints a brief console summary
 */
public class mmc_fifo {

    // CONTROL PARAMETERS

    // Number of completed customers per rho value
    static final int RUN_COMPLETED = 100000;

    // rho sweep range: 0.05, 0.10, ..., 0.95
    static final double RHO_START = 0.05;
    static final double RHO_END   = 0.95;
    static final double RHO_STEP  = 0.05;

    // Service rate and number of servers
    static final double MU = 1;
    static final int C = 1;

    // Deadline (patience) rate
    static final double DEADLINE_RATE = 0.02;

    // Output CSV path
    static final String CSV_FILE = "results/csv/mmc_fifo.csv";

    // RANDOM NUMBER GENERATOR

    static final Random RNG = new Random(12345);

    // EVENT TYPES

    private static enum EventType {
        ARRIVAL,
        DEPARTURE
    }

    // EVENT CLASS

    /**
     * Represents a generic event in the simulation:
     *  - time: when the event happens
     *  - type: ARRIVAL or DEPARTURE
     *  - serverId: for departures, which server finishes service (-1 for arrivals)
     *
     * Events are stored in a PriorityQueue ordered by time.
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
     * Specialized DEPARTURE event that remembers the arrivalTime
     * of the customer being served, so we can compute:
     *   - system delay = departure_time - arrival_time
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
     * Represents a single customer in the FIFO queue:
     *  - arrivalTime: when they arrived
     *  - deadlineTime: absolute time they expire
     *                  (arrivalTime + Exp(DEADLINE_RATE))
     */
    private static class Customer {
        double arrivalTime;
        double deadlineTime;

        Customer(double arrivalTime, double deadlineTime) {
            this.arrivalTime  = arrivalTime;
            this.deadlineTime = deadlineTime;
        }
    }

    // SUMMARY CLASS

    /**
     * Holds statistics for one (lambda, mu, c):
     *  - avgNs, avgNq, avgDs, avgDq
     *  - utilization
     *  - completed
     *  - Ns_th, Nq_th, Ds_th, Dq_th (theory)
     *  - pExpired: probability a served customer was expired at service start
     */
    private static class Summary {
        double avgNs, avgNq, avgDs, avgDq;
        double utilization;
        int    completed;

        double Ns_th, Nq_th, Ds_th, Dq_th;

        double pExpired;
    }

    // SUMMARY ROW FOR CSV

    /**
     * Row written to CSV for each rho:
     *  - rho, lambda, mu, c
     *  - Ns_sim, Ns_th
     *  - Nq_sim, Nq_th
     *  - Ds_sim, Ds_th
     *  - Dq_sim, Dq_th
     *  - utilization_sim
     *  - pExpired_sim
     */
    private static class SummaryRow {
        double rho, lambda, mu;
        int    c;

        double Ns_sim, Ns_th;
        double Nq_sim, Nq_th;
        double Ds_sim, Ds_th;
        double Dq_sim, Dq_th;

        double utilization_sim;

        double pExpired_sim;
    }

    // EXPONENTIAL SAMPLING

    /**
     * Draws a sample from an exponential distribution with rate = rate
     * via the inverse CDF method:
     *   X = -ln(1-U)/rate,  U ~ Uniform(0,1).
     */
    private static double expSample(double rate) {
        double u = RNG.nextDouble();
        if (u == 0.0) u = 1e-16;       // avoid log(0)
        return -Math.log(1.0 - u) / rate;
    }

    // FIND IDLE SERVER

    /**
     * Finds the index of the first idle server.
     * Returns:
     *  - index in [0, c-1] if a free server exists
     *  - -1 if all servers are busy
     */
    private static int findIdleServer(Server[] servers) {
        for (int i = 0; i < servers.length; i++) {
            if (!servers[i].busy) return i;
        }
        return -1;
    }

    // SINGLE RUN SIMULATOR

    /**
     * Runs a single M/M/c simulation with deadlines and FIFO:
     *  - Input: lambda, mu, c, deadlineRate
     *  - Uses event-driven simulation (arrival/departure events)
     *  - Tracks Ns, Nq, Ds, Dq, utilization, pExpired
     *  - Stops after RUN_COMPLETED customers have completed service
     *  - Returns a Summary object with all statistics
     */
    private static Summary simulate(double lambda,
                                    double mu,
                                    int c,
                                    double deadlineRate) {

        PriorityQueue<Event> fel = new PriorityQueue<>();
        Queue<Customer> fifoQ = new LinkedList<>();
        Server[] servers = new Server[c];
        for (int i = 0; i < c; i++) servers[i] = new Server();

        double t = 0.0;          // current time
        double lastT = 0.0;      // time of last event

        int Ns = 0;              // number in system
        int Nq = 0;              // number in queue

        double areaNs = 0.0;     // integral of Ns(t) over time
        double areaNq = 0.0;     // integral of Nq(t) over time
        double[] serverBusyArea = new double[c];

        double totalSystemDelay = 0.0;
        double totalQueueDelay  = 0.0;
        int    completed        = 0;

        int expiredServed = 0;   // how many served customers were expired at service start

        // First arrival at t = 0
        fel.add(new Event(0.0, EventType.ARRIVAL, -1));

        while (completed < RUN_COMPLETED) {

            Event e = fel.poll();
            double oldT = t;
            t = e.time;

            // Update integrals for Ns, Nq, and server busy times
            double dt = t - oldT;
            areaNs += Ns * dt;
            areaNq += Nq * dt;
            for (int i = 0; i < c; i++) {
                if (servers[i].busy) serverBusyArea[i] += dt;
            }

            if (e.type == EventType.ARRIVAL) {
                // New arrival
                Ns++;

                // Schedule next arrival
                double inter = expSample(lambda);
                fel.add(new Event(t + inter, EventType.ARRIVAL, -1));

                // Customer's deadline
                double deadlineTime = t + expSample(deadlineRate);
                Customer cust = new Customer(t, deadlineTime);

                // Check for idle server
                int idle = findIdleServer(servers);
                if (idle >= 0) {
                    // Immediate service
                    servers[idle].busy = true;

                    // Check expiry at service start
                    if (t >= cust.deadlineTime) {
                        expiredServed++;
                    }

                    double serviceTime = 1;
                    double departTime  = t + serviceTime;
                    fel.add(new ServiceEvent(departTime, idle, cust.arrivalTime));

                    // Queue delay = 0
                } else {
                    // Join FIFO queue
                    fifoQ.add(cust);
                    Nq++;
                }

            } else { // DEPARTURE
                ServiceEvent se = (ServiceEvent) e;
                int sid = se.serverId;

                // One customer departs system
                completed++;
                Ns--;
                servers[sid].busy = false;

                // System delay = departure - arrival
                double Ds_i = t - se.arrivalTime;
                totalSystemDelay += Ds_i;

                // If queue non empty, start service for next in FIFO
                if (!fifoQ.isEmpty()) {
                    Customer next = fifoQ.poll();
                    Nq--;

                    double serviceStart = t;

                    // Check expiry at service start
                    if (serviceStart >= next.deadlineTime) {
                        expiredServed++;
                    }

                    // Queue delay = service start - arrival
                    double Dq_i = serviceStart - next.arrivalTime;
                    totalQueueDelay += Dq_i;

                    // Start service
                    servers[sid].busy = true;
                    double serviceTime = 1;
                    double departTime  = serviceStart + serviceTime;
                    fel.add(new ServiceEvent(departTime, sid, next.arrivalTime));
                }
            }
        }

        // Build summary
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

    // FACTORIAL (FOR ERLANG-C)

    /**
     * Computes n! for small n using a simple loop.
     */
    private static long factorial(int n) {
        long f = 1L;
        for (int i = 2; i <= n; i++) f *= i;
        return f;
    }

    // THEORETICAL M/M/c (ERLANG-C)

    /**
     * Fills in S.Ns_th, S.Nq_th, S.Ds_th, S.Dq_th using standard
     * M/M/c (Erlang-C) formulas, ignoring deadlines:
     *  - Input: lambda, mu, c
     *  - Uses offered load a = lambda / mu and rho = lambda / (c * mu)
     *  - If rho >= 1, sets theory values to NaN
     */
    private static void computeTheory(Summary S, double lambda, double mu, int c) {

        double a   = lambda / mu;
        double rho = lambda / (c * mu);

        if (rho >= 1.0) {
            S.Nq_th = Double.NaN;
            S.Ns_th = Double.NaN;
            S.Dq_th = Double.NaN;
            S.Ds_th = Double.NaN;
            return;
        }

        // P0
        double sum = 0.0;
        for (int n = 0; n <= c - 1; n++) {
            sum += Math.pow(a, n) / factorial(n);
        }
        double cTerm = Math.pow(a, c) / (factorial(c) * (1 - rho));
        double P0    = 1.0 / (sum + cTerm);

        // Probability of waiting
        double Pwait = cTerm * P0;

        // Nq and Ns
        S.Nq_th = (Pwait * rho) / (1 - rho);
        S.Ns_th = S.Nq_th + a;

        // Delays
        S.Dq_th = S.Nq_th / lambda;
        S.Ds_th = S.Dq_th + 1.0 / mu;
    }

    // CSV WRITER

    /**
     * Writes all rows to CSV with columns:
     *  rho, lambda, mu, c,
     *  Ns_sim, Ns_th,
     *  Nq_sim, Nq_th,
     *  Ds_sim, Ds_th,
     *  Dq_sim, Dq_th,
     *  utilization_sim,
     *  pExpired_sim
     */
    private static void writeCSV(List<SummaryRow> rows) {
        try {
            File f = new File(CSV_FILE);
            f.getParentFile().mkdirs();

            PrintWriter pw = new PrintWriter(new FileWriter(f));

            pw.println("rho,lambda,mu,c,"
                    + "Ns_sim,Ns_th,"
                    + "Nq_sim,Nq_th,"
                    + "Ds_sim,Ds_th,"
                    + "Dq_sim,Dq_th,"
                    + "utilization_sim,"
                    + "pExpired_sim");

            for (SummaryRow r : rows) {
                pw.printf(Locale.US,
                        "%.4f,%.6f,%.6f,%d,"
                                + "%.6f,%.6f,"      // Ns_sim, Ns_th
                                + "%.6f,%.6f,"      // Nq_sim, Nq_th
                                + "%.6f,%.6f,"      // Ds_sim, Ds_th
                                + "%.6f,%.6f,"      // Dq_sim, Dq_th
                                + "%.6f,%.6f%n",    // utilization_sim, pExpired_sim
                        r.rho, r.lambda, r.mu, r.c,
                        r.Ns_sim, r.Ns_th,
                        r.Nq_sim, r.Nq_th,
                        r.Ds_sim, r.Ds_th,
                        r.Dq_sim, r.Dq_th,
                        r.utilization_sim, r.pExpired_sim);

            }

            pw.close();
        } catch (IOException e) {
            System.err.println("Error writing CSV: " + e.getMessage());
        }
    }

    // SUMMARY PRINTER

    /**
     * Prints a compact summary table of key statistics to the console.
     */
    private static void printSummary(List<SummaryRow> rows) {
        System.out.println("==============================================================");
        System.out.println(" c | rho |   Ns(sim)   Ns(th) |   Nq(sim)   Nq(th) | pExpired");
        System.out.println("==============================================================");
        for (SummaryRow r : rows) {
            System.out.printf(Locale.US,
                    "%2d | %.2f | %8.4f %8.4f | %8.4f %8.4f | %8.4f%n",
                    r.c, r.rho,
                    r.Ns_sim, r.Ns_th,
                    r.Nq_sim, r.Nq_th,
                    r.pExpired_sim);
        }
        System.out.println("==============================================================");
    }

    // MAIN DRIVER

    /**
     * Steps:
     *  1. For rho from RHO_START to RHO_END:
     *      a) Set lambda = rho * C * MU
     *      b) Run the simulation with FIFO + deadlines
     *      c) Compute Erlang-C theory (no deadlines)
     *      d) Build a SummaryRow and add to list
     *  2. Write CSV
     *  3. Print console summary
     */
    public static void main(String[] args) {

        List<SummaryRow> rows = new ArrayList<>();

        for (double rho = RHO_START; rho <= RHO_END + 1e-9; rho += RHO_STEP) {
            double lambda = rho * C * MU;

            Summary S = simulate(lambda, MU, C, DEADLINE_RATE);
            computeTheory(S, lambda, MU, C);

            SummaryRow R = new SummaryRow();
            R.rho = rho;
            R.lambda = lambda;
            R.mu = MU;
            R.c = C;

            R.Ns_sim = S.avgNs;
            R.Ns_th  = S.Ns_th;
            R.Nq_sim = S.avgNq;
            R.Nq_th  = S.Nq_th;
            R.Ds_sim = S.avgDs;
            R.Ds_th  = S.Ds_th;
            R.Dq_sim = S.avgDq;
            R.Dq_th  = S.Dq_th;

            R.utilization_sim = S.utilization;
            R.pExpired_sim    = S.pExpired;

            rows.add(R);

            System.out.printf(Locale.US,
                    "FIFO: rho=%.2f -> Ns_sim=%.4f, Ns_th=%.4f, pExpired=%.6f%n",
                    rho, S.avgNs, S.Ns_th, S.pExpired);
        }

        writeCSV(rows);
        printSummary(rows);
    }
}
