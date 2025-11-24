import java.util.*;
import java.io.*;

/**
 * M/M/c queue simulator with:
 * - FIFO Service
 * - Exponential interarrival and service times
 * - Infinite buffer (no blocking)
 *
 * For each rho in [RHO_START, RHO_END] with step RHO_STEP, it:
 *  1. Sets lambda = rho * c * mu
 *  2. Runs a single simulation until RUN_COMPLETED customers depart
 *  3. Estimates Ns, Nq, Ds, Dq, utilization
 *  4. Computes theoretical M/M/c (Erlang-C) values
 *  5. Writes one summary row to CSV
 *  6. Stores results and prints a grand summary table with % errors
 */
public class mmc_infinite {

    // CONTROL PARAMETERS
    // Number of completed customers per rho value
    static final int RUN_COMPLETED = 1000000;
    // rho sweep range: 0.05, 0.10, ..., 0.95
    static final double RHO_START = 0.05;
    static final double RHO_END   = 0.95;
    static final double RHO_STEP  = 0.05;
    // Output CSV path
    static final String CSV_FILE = "results/csv/mmc_infinite.csv";

    // EVENT TYPES
    private static enum EventType {
        ARRIVAL,
        DEPARTURE
    }

    //EVENT CLASS
    /**
     * Represents a generic event in the simulation:
     *  - time: when the event happens
     *  - type: ARRIVAL or DEPARTURE
     *  - serverId: for departures, which server is finishing service (-1 for arrivals)
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

    // SERVICE DEPARTURE EVENT with arrivalTime
    /**
     * Specialized event for DEPARTURE that also remembers the arrivalTime
     * of the customer being served.
     *
     * This lets us compute:
     *   - system delay = departure_time - arrival_time
     *   - queue delay = service_start_time - arrival_time
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

    // SUMMARY CLASS
    /**
     * Holds per-rho results:
     *  - avgNs: time-average number in system
     *  - avgNq: time-average number in queue
     *  - avgDs: average system delay (arrival to departure)
     *  - avgDq: average queue delay (arrival to service start)
     *  - utilization: average fraction of time servers are busy
     *  - completed: how many departures were simulated
     *
     * Also holds the theoretical M/M/c values (Erlang-C):
     *  - Ns_th, Nq_th, Ds_th, Dq_th
     */
    private static class Summary {
        double avgNs, avgNq, avgDs, avgDq;
        double utilization;
        int completed;

        // theory
        double Ns_th, Nq_th, Ds_th, Dq_th;
    }

    // For summary table printed at the end
    /**
     * A flattened record used to print the summary table:
     *  - rho, lambda, mu, c
     *  - sim and theory for Ns, Nq, Ds, Dq
     *  - percentage error for each metric
     */
    private static class SummaryRow {
        double rho, lambda, mu;
        int c;

        double Ns_sim, Ns_th, Ns_errp;
        double Nq_sim, Nq_th, Nq_errp;
        double Ds_sim, Ds_th, Ds_errp;
        double Dq_sim, Dq_th, Dq_errp;
    }

    // EXPONENTIAL SAMPLING
    /**
     * Draws a sample from an exponential distribution with rate = rate.
     * Uses the inverse CDF method:
     *   X = -ln(1-U)/rate, where U ~ Uniform(0,1).
     *
     * Used for both:
     *   - interarrival times (rate = lambda)
     *   - service times (rate = mu)
     */
    private static double expSample(double rate) {
        return -Math.log(1 - Math.random()) / rate;
    }

    /**
     * Finds the index of the first free server in the array.
     * Returns:
     *   - index in [0, c-1] if a free server exists
     *   - -1 if all servers are busy
     */
    private static int findFreeServer(Server[] servers) {
        for (int i = 0; i < servers.length; i++) {
            if (!servers[i].busy) return i;
        }
        return -1;
    }

    // THEORETICAL VALUES (Erlang C)
    /**
     * Computes n! as double via a simple loop.
     * n is small here (up to c), so this is safe and fast.
     */
    private static double factorial(int n) {
        double r = 1.0;
        for (int i = 2; i <= n; i++) r *= i;
        return r;
    }
    /**
     * Given lambda, mu, c:
     *  - Computes M/M/c theoretical metrics using the Erlang-C formula:
     *      a   = lambda / mu
     *      rho = lambda / (c * mu)
     *      P0  = normalizing constant
     *      Pwait = probability that an arrival waits
     *      Nq_th  = (Pwait * rho) / (1 - rho)
     *      Ns_th  = Nq_th + a
     *      Dq_th  = Nq_th / lambda
     *      Ds_th  = Dq_th + 1/mu
     *
     *  - Stores them in the Summary S.
     *
     * If rho >= 1, the system is unstable and theory is set to NaN.
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

        // Compute P0
        double sum = 0.0;
        for (int n = 0; n <= c - 1; n++) {
            sum += Math.pow(a, n) / factorial(n);
        }

        double cTerm = Math.pow(a, c) / (factorial(c) * (1 - rho));
        double P0 = 1.0 / (sum + cTerm);

        // Probability of waiting
        double Pwait = cTerm * P0;

        // Nq and Ns
        S.Nq_th = (Pwait * rho) / (1 - rho);
        S.Ns_th = S.Nq_th + a;

        // Delays
        S.Dq_th = S.Nq_th / lambda;
        S.Ds_th = S.Dq_th + 1.0 / mu;
    }


    // MAIN DRIVER
    /**
     * Steps:
     *  1. Create/open the CSV file and write header.
     *  2. Initialize a list to store per-rho summaries for the final table.
     *  3. For rho from RHO_START to RHO_END (step RHO_STEP):
     *      a) Set c, mu, lambda = rho * c * mu
     *      b) Run the simulation for RUN_COMPLETED completions
     *      c) Compute theory for this (lambda, mu, c)
     *      d) Append summary row to CSV
     *      e) Store stats in allSummaries (includes % error)
     *      f) Print quick 'finished rho' message
     *  4. After the loop, print a summary table:
     *      - For each rho: show sim vs theory, and % error for Ns, Nq, Ds, Dq
     *  5. Print a final "Done" message.
     */
    public static void main(String[] args) throws Exception {

        // 1) Create CSV and write header
        PrintWriter pw = new PrintWriter(new FileWriter(CSV_FILE));
        pw.println("rho,lambda,mu,c,avg_Ns,avg_Nq,avg_Ds,avg_Dq,utilization,completed");
        pw.close();

        // 2) List to hold all per-rho summaries
        List<SummaryRow> allSummaries = new ArrayList<>();

        // 3) Sweep over rho values
        for (double rho = RHO_START; rho <= RHO_END + 1e-12; rho += RHO_STEP) {
            // a) Set c, mu, lambda = rho * c * mu
            int c = 3;
            double mu = 0.75;
            double lambda = rho * c * mu;

            // b) Run one simulation for this rho
            Summary S = runSingleSimulation(lambda, mu, c);

            // c) Compute theory for this (lambda, mu, c)
            computeTheory(S, lambda, mu, c);

            // d) Append summary to CSV
            appendSummaryToCSV(S, lambda, mu, c);

            // e) Save data for final grand summary table
            SummaryRow R = new SummaryRow();
            R.rho = rho;
            R.lambda = lambda;
            R.mu = mu;
            R.c = c;

            R.Ns_sim = S.avgNs; R.Ns_th = S.Ns_th;
            R.Nq_sim = S.avgNq; R.Nq_th = S.Nq_th;
            R.Ds_sim = S.avgDs; R.Ds_th = S.Ds_th;
            R.Dq_sim = S.avgDq; R.Dq_th = S.Dq_th;

            R.Ns_errp = 100 * Math.abs(R.Ns_sim - R.Ns_th) / S.Ns_th;
            R.Nq_errp = 100 * Math.abs(R.Nq_sim - R.Nq_th) / S.Nq_th;
            R.Ds_errp = 100 * Math.abs(R.Ds_sim - R.Ds_th) / S.Ds_th;
            R.Dq_errp = 100 * Math.abs(R.Dq_sim - R.Dq_th) / S.Dq_th;

            allSummaries.add(R);

            // f) Quick progress info
            System.out.printf("Finished rho=%.2f\n", rho);
        }

        // 4) Print grand summary table
        System.out.println("\n=====================================================");
        System.out.println(" SUMMARY TABLE (Simulation vs Theory, % Error)");
        System.out.println("=====================================================");

        System.out.printf(
                "%5s | %10s %10s %10s | %10s %10s %10s | %10s %10s %10s | %10s %10s %10s\n",
                "rho",
                "Ns_sim","Ns_th","Ns_err%",
                "Nq_sim","Nq_th","Nq_err%",
                "Ds_sim","Ds_th","Ds_err%",
                "Dq_sim","Dq_th","Dq_err%"
        );

        for (SummaryRow R : allSummaries) {
            System.out.printf(
                    "%.2f  | %10.5f %10.5f %10.3f | %10.5f %10.5f %10.3f | %10.5f %10.5f %10.3f | %10.5f %10.5f %10.3f\n",
                    R.rho,
                    R.Ns_sim, R.Ns_th, R.Ns_errp,
                    R.Nq_sim, R.Nq_th, R.Nq_errp,
                    R.Ds_sim, R.Ds_th, R.Ds_errp,
                    R.Dq_sim, R.Dq_th, R.Dq_errp
            );
        }

        System.out.println("=====================================================");
        System.out.println("Done. CSV written: " + CSV_FILE);
    }

    // SINGLE SIMULATION
    /**
     * Runs a single M/M/c FIFO simulation for given (lambda, mu, c).
     *
     * Steps:
     *  1. Initialize:
     *      - Future Event List (fel) as a min-heap of events
     *      - FIFO queue for waiting customers (stores arrival times)
     *      - c servers, all idle
     *      - Clock variables t, lastT
     *      - State Ns (# in system), Nq (# in queue)
     *      - Time-integral accumulators areaNs, areaNq
     *      - Per-server busy time integrals
     *      - Counters: totalSystemDelay, totalQueueDelay, completed
     *  2. Schedule the first ARRIVAL at t=0.
     *  3. While completed < RUN_COMPLETED:
     *      a) Pop the earliest event from fel.
     *      b) Advance time to event time; integrate Ns, Nq, and server busy times.
     *      c) If ARRIVAL:
     *           - Ns++
     *           - Schedule next arrival at t + Exp(lambda)
     *           - If a server is free:
     *                 start service immediately, schedule DEPARTURE
     *             else:
     *                 push arrival time into FIFO queue, Nq++
     *      d) If DEPARTURE:
     *           - completed++
     *           - Ns--
     *           - Use arrivalTime to accumulate system delay
     *           - Free the server
     *           - If queue not empty:
     *                 pop next arrivalTime from FIFO
     *                 Nq--
     *                 start service on this server
     *                 accumulate queue delay
     *                 schedule a new DEPARTURE for that customer
     *  4. After the loop:
     *      - Compute avgNs = areaNs / totalTime
     *      - Compute avgNq = areaNq / totalTime
     *      - Compute avgDs = totalSystemDelay / completed
     *      - Compute avgDq = totalQueueDelay / completed
     *      - Compute utilization from serverBusyArea
     *      - Return Summary S
     */
    private static Summary runSingleSimulation(double lambda, double mu, int c) {

        PriorityQueue<Event> fel = new PriorityQueue<>();
        Queue<Double> fifoQ = new LinkedList<>();
        Server[] servers = new Server[c];
        for (int i = 0; i < c; i++) servers[i] = new Server();

        double t = 0.0, lastT = 0.0;
        int Ns = 0, Nq = 0;
        double areaNs = 0.0, areaNq = 0.0;
        double[] serverBusyArea = new double[c];

        double totalSystemDelay = 0.0;
        double totalQueueDelay = 0.0;
        int completed = 0;

        // First arrival at t = 0
        fel.add(new Event(0.0, EventType.ARRIVAL, -1));

        while (completed < RUN_COMPLETED) {

            Event e = fel.poll();
            double tNew = e.time;
            double dt = tNew - lastT;

            // Integrate Ns and Nq over [lastT, tNew)
            areaNs += Ns * dt;
            areaNq += Nq * dt;

            // Integrate per-server busy time
            for (int i = 0; i < c; i++)
                if (servers[i].busy) serverBusyArea[i] += dt;

            t = tNew;
            lastT = tNew;

            // HANDLE ARRIVAL
            if (e.type == EventType.ARRIVAL) {

                Ns++;

                // Schedule next arrival
                fel.add(new Event(t + expSample(lambda), EventType.ARRIVAL, -1));

                int free = findFreeServer(servers);
                if (free != -1) {
                    // Serve immediately
                    servers[free].busy = true;
                    double st = expSample(mu);
                    fel.add(new ServiceEvent(t + st, free, t));
                } else {
                    // Join FIFO queue
                    fifoQ.add(t);
                    Nq++;
                }

                // HANDLE DEPARTURE
            } else {

                ServiceEvent de = (ServiceEvent) e;
                completed++;
                Ns--;

                // System delay = departure time - arrival time
                totalSystemDelay += (t - de.arrivalTime);

                int sid = de.serverId;
                servers[sid].busy = false;

                // Serve next customer from queue if any
                if (!fifoQ.isEmpty()) {

                    double arrTime = fifoQ.poll();
                    Nq--;

                    servers[sid].busy = true;

                    // Queue delay = service start - arrival time
                    totalQueueDelay += (t - arrTime);

                    double st = expSample(mu);
                    fel.add(new ServiceEvent(t + st, sid, arrTime));

                    // IMPORTANT: do NOT Ns++ here.
                    // Customer was already counted in Ns upon arrival.
                }
            }
        }

        // Build summary
        Summary S = new Summary();
        double totalTime = lastT;

        S.avgNs = areaNs / totalTime;
        S.avgNq = areaNq / totalTime;
        S.avgDs = totalSystemDelay / completed;
        S.avgDq = totalQueueDelay / completed;

        double busySum = 0.0;
        for (int i = 0; i < c; i++) busySum += serverBusyArea[i];
        S.utilization = busySum / (c * totalTime);

        S.completed = completed;

        return S;
    }

    // CSV APPEND
    /**
     * Appends a single summary row to the CSV file for a given (lambda, mu, c).
     *
     * Columns:
     *   rho, lambda, mu, c,
     *   avg_Ns, avg_Nq, avg_Ds, avg_Dq,
     *   utilization, completed
     *
     * This is what the MATLAB script reads to produce Ns, Nq, Ds, Dq vs rho plots.
     */
    private static void appendSummaryToCSV(Summary S, double lambda, double mu, int c) {
        try {
            FileWriter fw = new FileWriter(CSV_FILE, true);
            PrintWriter pw = new PrintWriter(fw);

            double rho = lambda / (c * mu);

            pw.printf(Locale.US,
                    "%.4f,%.6f,%.6f,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%d\n",
                    rho, lambda, mu, c,
                    S.avgNs, S.avgNq, S.avgDs, S.avgDq,
                    S.utilization, S.completed);

            pw.close();
        }
        catch (IOException e) {
            System.out.println("CSV write error: " + e.getMessage());
        }
    }

}
