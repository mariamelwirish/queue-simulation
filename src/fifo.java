import java.util.*;
import java.io.*;

/**
 * M/M/c queue simulator with deadlines and FIFO service.
 *
 * 
 * BASELINE MODEL (when all flags disabled):
 * 
 * - Exponential interarrival times (rate lambda)
 * - Exponential service times (rate mu)
 * - Exponential deadlines/patience (rate DEADLINE_RATE)
 * - Infinite buffer (no blocking)
 * - FIFO service
 *
 * 
 * DETERMINISTIC MODE (DDC) — enabled via flags:
 * 
 * - Deterministic interarrival times (1 / lambda)
 * - Deterministic service times (1 / mu)
 * - Deadlines remain exponential
 *
 * 
 * DROP MODE — enabled via flag:
 * 
 * - Jobs whose deadlines expire before service are DROPPED
 * - Dropped jobs do not consume service
 *
 * 
 * UNIFIED DEADLINE METRIC:
 * 
 * - pExpired = P(job expires before service begins)
 *   • If DROP_EXPIRED = false → expired jobs are still served
 *   • If DROP_EXPIRED = true  → expired jobs are dropped
 *
 * This file serves as a controlled extension of the baseline simulator.
 * All changes are guarded by flags so the original behavior is preserved
 * and can be used as a reference for comparison.
 *
 * 
 * For each rho in [RHO_START, RHO_END] with step RHO_STEP:
 * 
 *  1. Sets lambda = rho * c * mu
 *  2. Runs the simulation until RUN_COMPLETED jobs complete service
 *  3. Estimates Ns, Nq, Ds, Dq
 *  4. Computes Erlang-C M/M/c theory (baseline reference)
 *  5. Estimates pExpired
 *  6. Writes one summary row to CSV
 *
 * NOTE:
 * - Erlang-C theory applies strictly to the baseline M/M/c model
 * - In deterministic or drop modes, theory values are kept only
 *   as a reference baseline
 */
public class fifo {

    
    // CONTROL PARAMETERS
    

    // Number of completed jobs per rho value
    static final int RUN_COMPLETED = 1000000;

    // rho sweep range
    static final double RHO_START = 0.05;
    static final double RHO_END   = 0.95;
    static final double RHO_STEP  = 0.05;

    // Service rate and number of servers
    static final double MU = 1.0;
    static final int C = 1;

    // Deadline (patience) rate
    static final double DEADLINE_RATE = 0.1;

    
    // MODE FLAGS (baseline preserved)
    

    // Deterministic arrival process
    static final boolean DETERMINISTIC_ARRIVAL = false;

    // Deterministic service times
    static final boolean DETERMINISTIC_SERVICE = false;

    // Drop jobs whose deadlines expire before service
    static final boolean DROP_EXPIRED = false;

    
    // CSV FILE SELECTION
    

    static final String CSV_FILE;
    static {
        // Build a mode-aware filename to avoid overwriting.
        if (DROP_EXPIRED) {
            CSV_FILE = "../results/csv/drop_fifo.csv";
        } else if (DETERMINISTIC_ARRIVAL || DETERMINISTIC_SERVICE) {
            CSV_FILE = "../results/csv/ddc_fifo.csv";
        } else {
            CSV_FILE = "../results/csv/mmc_fifo.csv";
        }
    }

    // Random number generator (used only in exponential mode)
    static final Random RNG = new Random(12345);

    
    // EVENT TYPES
    

    private static enum EventType {
        ARRIVAL,
        DEPARTURE
    }

    
    // EVENT CLASS
    

    /**
     * Generic event stored in the future event list (FEL).
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
     * Departure event that remembers arrival time
     * for delay calculations.
     */
    private static class ServiceEvent extends Event {
        double arrivalTime;

        ServiceEvent(double time, int serverId, double arrivalTime) {
            super(time, EventType.DEPARTURE, serverId);
            this.arrivalTime = arrivalTime;
        }
    }

    // SERVER & CUSTOMER CLASSES

    private static class Server {
        boolean busy = false;
    }

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
     * Holds statistics for a single (lambda, mu, c) run.
     */
    private static class Summary {
        double avgNs, avgNq, avgDs, avgDq;
        double utilization;
        int completed;

        double Ns_th, Nq_th, Ds_th, Dq_th;

        double pExpired;
    }

    /**
     * Row written to CSV for each rho.
     */
    private static class SummaryRow {
        double rho, lambda, mu;
        int c;

        double Ns_sim, Ns_th;
        double Nq_sim, Nq_th;
        double Ds_sim, Ds_th;
        double Dq_sim, Dq_th;

        double utilization_sim;
        double pExpired_sim;
    }

    // EXPONENTIAL SAMPLING

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
     * Runs a single FIFO simulation with deadlines.
     *
     * INVARIANTS:
     *  - Ns counts ONLY jobs that are either in service or waiting in queue
     *  - Each job contributes:
     *        exactly one Ns++ (when admitted)
     *        exactly one Ns-- (when leaving system)
     *
     * DEADLINE HANDLING:
     *  - If DROP_EXPIRED = false:
     *        expired jobs are still served
     *  - If DROP_EXPIRED = true:
     *        expired jobs are removed immediately and never served
     *
     * pExpired:
     *  - Probability a job expires before service begins
     *  - Counted regardless of whether the job is later served or dropped
     */
    private static Summary simulate(double lambda, double mu, int c, double deadlineRate) {

        PriorityQueue<Event> fel = new PriorityQueue<>();
        Queue<Customer> fifoQ = new LinkedList<>();

        Server[] servers = new Server[c];
        for (int i = 0; i < c; i++) servers[i] = new Server();

        double t = 0.0;
        int Ns = 0;        // number in system
        int Nq = 0;        // number in queue

        double areaNs = 0.0;
        double areaNq = 0.0;
        double[] serverBusyArea = new double[c];

        double totalSystemDelay = 0.0;
        double totalQueueDelay  = 0.0;

        int arrivals  = 0;
        int completed = 0;
        int expired   = 0;

        // Schedule first arrival
        fel.add(new Event(0.0, EventType.ARRIVAL, -1));

        while (arrivals < RUN_COMPLETED || Ns > 0) {
        
            Event e = fel.poll();
            double oldT = t;
            t = e.time;

            // ---- time integrals ----
            double dt = t - oldT;
            areaNs += Ns * dt;
            areaNq += Nq * dt;
            for (int i = 0; i < c; i++)
                if (servers[i].busy)
                    serverBusyArea[i] += dt;

            
            // ARRIVAL EVENT
            
            if (e.type == EventType.ARRIVAL) {
                if (arrivals >= RUN_COMPLETED)
                    continue;
                arrivals++;

                // Schedule next arrival
                double interArrival =
                        DETERMINISTIC_ARRIVAL ? 1.0 / lambda : expSample(lambda);
                fel.add(new Event(t + interArrival, EventType.ARRIVAL, -1));

                // Create customer
                Customer cust = new Customer(t, t + expSample(deadlineRate));



                // -------- ADMISSION POINT --------
                // Job is now in the system
                Ns++;

                int idle = findIdleServer(servers);
                if (idle >= 0) {
                    // Start service immediately
                    servers[idle].busy = true;

                    double serviceTime =
                            DETERMINISTIC_SERVICE ? 1.0 / mu : expSample(mu);

                    fel.add(new ServiceEvent(t + serviceTime, idle, cust.arrivalTime));
                } else {
                    // Join FIFO queue
                    fifoQ.add(cust);
                    Nq++;
                }
            }

            
            // DEPARTURE EVENT
            
            else {

                ServiceEvent se = (ServiceEvent) e;
                int sid = se.serverId;

                // Job leaves system after service
                completed++;
                Ns--;
                servers[sid].busy = false;

                totalSystemDelay += t - se.arrivalTime;

                // Try to start service for next job in queue
                while (!fifoQ.isEmpty()) {

                    Customer next = fifoQ.poll();
                    Nq--;

                    // Deadline check at service start
                    if (t >= next.deadlineTime) {
                        expired++;

                        if (DROP_EXPIRED) {
                            // Job leaves system permanently
                            Ns--;
                            continue;
                        }
                        // else: serve expired job
                    }

                    totalQueueDelay += t - next.arrivalTime;

                    servers[sid].busy = true;

                    double serviceTime =
                            DETERMINISTIC_SERVICE ? 1.0 / mu : expSample(mu);

                    fel.add(new ServiceEvent(t + serviceTime, sid, next.arrivalTime));
                    break;
                }
            }

            // ---- SAFETY CHECK (debug only; can remove later) ----
            if (Ns < 0) {
                throw new RuntimeException("Ns < 0 at time " + t);
            }
        }

        
        // BUILD SUMMARY
        

        Summary S = new Summary();
        double totalTime = t;

        S.avgNs = areaNs / totalTime;
        S.avgNq = areaNq / totalTime;
        S.avgDs = totalSystemDelay / arrivals;
        S.avgDq = totalQueueDelay / arrivals;

        double busySum = 0.0;
        for (double b : serverBusyArea) busySum += b;
        S.utilization = busySum / (c * totalTime);

        S.pExpired = expired / (double) arrivals;
        System.out.print(expired + " ");
        S.completed = completed;
        System.out.println(completed + " " + arrivals);
        return S;
    }

    
    // ERLANG-C THEORY (BASELINE REFERENCE)

    private static long factorial(int n) {
        long f = 1L;
        for (int i = 2; i <= n; i++) f *= i;
        return f;
    }

    private static void computeTheory(Summary S, double lambda, double mu, int c) {

        double a   = lambda / mu;
        double rho = lambda / (c * mu);

        if (rho >= 1.0) {
            S.Nq_th = S.Ns_th = S.Dq_th = S.Ds_th = Double.NaN;
            return;
        }

        double sum = 0.0;
        for (int n = 0; n <= c - 1; n++)
            sum += Math.pow(a, n) / factorial(n);

        double cTerm = Math.pow(a, c) / (factorial(c) * (1 - rho));
        double P0 = 1.0 / (sum + cTerm);
        double Pwait = cTerm * P0;

        S.Nq_th = (Pwait * rho) / (1 - rho);
        S.Ns_th = S.Nq_th + a;
        S.Dq_th = S.Nq_th / lambda;
        S.Ds_th = S.Dq_th + 1.0 / mu;
    }

    // CSV WRITER & MAIN DRIVER

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
            R.pExpired_sim = S.pExpired;

            rows.add(R);
        }

        writeCSV(rows);
    }

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
                        + "pExpired_sim");


            for (SummaryRow r : rows) {
               pw.printf(Locale.US,
                        "%.4f,%.6f,%.6f,%d,"
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
                        r.pExpired_sim);

            }

            pw.close();
        } catch (IOException e) {
            System.err.println("Error writing CSV: " + e.getMessage());
        }
    }
}
