import java.util.*;
import java.io.*;

/**
 * Queue simulator with deadlines and non-preemptive SJF service.
 *
 * BASELINE MODEL:
 * - Interarrival: Exp(lambda) unless DETERMINISTIC_ARRIVALS
 * - Deadlines: Exp(DEADLINE_RATE) (always exponential)
 * - Job size: L ~ Geometric(p) in bits
 * - Scheduling: non-preemptive SJF on waiting queue (smallest L first)
 *
 * SERVICE MODEL (new):
 * - If LENGTH_BASED_SERVICE = true: serviceTime = L / R  (R in bits/time)
 *   => per-job effective service rate = R / L
 * - If LENGTH_BASED_SERVICE = false: service ~ Exp(mu) unless DETERMINISTIC_SERVICE
 *
 * FIFO BASELINE CONVENTIONS APPLIED:
 * - Arrivals capped strictly at RUN_COMPLETED (no more scheduled after cap)
 * - Main loop drains: while (arrivals < RUN_COMPLETED || Ns > 0)
 * - No expiry upon arrival (deadline checks occur only at service start / while waiting)
 * - pExpired definition:
 *      pExpired = (# jobs whose deadline expired BEFORE service start) / (# arrivals)
 * - DROP_EXPIRED flag controls whether expired jobs are dropped or still served
 *   WITHOUT changing pExpired definition.
 * - Ns accounting:
 *      Ns++ once on admission
 *      Ns-- on service completion, and additionally on drop (if DROP_EXPIRED)
 *
 * NEW (OPTIONAL): PREEMPTIVE SJF (SRPT-like)
 * - If PREEMPTIVE_SJF = true:
 *      A newly arriving job can preempt a running job if its (remaining) service time
 *      is smaller than the running job's remaining service time.
 *      The preempted job is put back into the SJF queue with remaining service time
 *      (so it does NOT restart from the beginning).
 * - To handle stale DEPARTURE events after preemption, we use a per-server "version" token.
 */
public class sjf {

    // CONTROL FLAGS (DDC MODE + FIFO CONVENTIONS)

    /**
     * If true: interarrival time = 1/lambda (deterministic).
     * If false: interarrival ~ Exp(lambda).
     */
    static final boolean DETERMINISTIC_ARRIVALS = false;

    /**
     * If true (and LENGTH_BASED_SERVICE=false): service time = 1/mu (deterministic).
     * If false (and LENGTH_BASED_SERVICE=false): service ~ Exp(mu).
     *
     * Note: when LENGTH_BASED_SERVICE=true, service time is L/R (deterministic given L),
     * so this flag is ignored by design.
     */
    static final boolean DETERMINISTIC_SERVICE  = false;

    /**
     * FIFO baseline behavior:
     * If true: drop jobs that are already expired at service start.
     * If false: still serve them (but they still count as expired for pExpired).
     */
    static final boolean DROP_EXPIRED = false;

    /**
     * New service model:
     * If true: serviceTime = L / R (size-based service).
     * If false: revert to Exp(mu) / deterministic 1/mu based on flag.
     */
    static final boolean LENGTH_BASED_SERVICE = true;

    /**
     * NEW: Preemptive SJF mode.
     * If true: enable preemption based on smallest remaining service time.
     * If false: baseline non-preemptive SJF.
     */
    static final boolean PREEMPTIVE_SJF = false;

    // CONTROL PARAMETERS

    // Number of arrivals per rho value (strict cap)
    static final int RUN_COMPLETED = 1000000;

    // rho sweep range
    static final double RHO_START = 0.05;
    static final double RHO_END   = 0.95;
    static final double RHO_STEP  = 0.05;

    // Service parameter (kept for consistency + to derive R)
    static final double MU = 1.0;
    static final int C = 1;

    // Deadline (patience) rate
    static final double DEADLINE_RATE = 0.1;

    static final double MEAN_LENGTH_BITS = 100;
    static final double GEO_P            = 1.0 / MEAN_LENGTH_BITS;

    /**
     * Transmission/processing rate R (bits/time).
     * Chosen to preserve mean service time:
     *   E[S] = E[L]/R  and  E[S] (old) = 1/MU  =>  R = MU * E[L]
     */
    static final double R = MU * MEAN_LENGTH_BITS;

    // OUTPUT CSV (FILENAME DEPENDS ON FLAGS)

    static final String CSV_FILE;
    static {
        // Build a mode-aware filename to avoid overwriting. (Same style as EDF naming)
        if (PREEMPTIVE_SJF) {
            CSV_FILE = "../results/csv/preemptive_sjf.csv";
        }
        else if (DROP_EXPIRED) {
            CSV_FILE = "../results/csv/drop_sjf.csv";
        }
        else if (DETERMINISTIC_ARRIVALS || DETERMINISTIC_SERVICE) {
            CSV_FILE = "../results/csv/ddc_sjf.csv";
        }
        else {
            CSV_FILE = "../results/csv/mmc_sjf.csv";
        }
    }

    // RANDOM NUMBER GENERATOR

    static final Random RNG = new Random(12345);

    // EVENT TYPES

    private static enum EventType {
        ARRIVAL,
        DEPARTURE
    }

    // EVENT CLASS

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

    private static class ServiceEvent extends Event {
        double arrivalTime;

        // NEW: version token to ignore stale departure events after preemption.
        int version;

        ServiceEvent(double time, int serverId, double arrivalTime, int version) {
            super(time, EventType.DEPARTURE, serverId);
            this.arrivalTime = arrivalTime;
            this.version = version;
        }
    }

    // SERVER CLASS

    private static class Server {
        boolean busy = false;

        // NEW (needed for preemptive mode)
        Customer current = null;       // currently running job
        double serviceStartTime = 0.0; // start time of current service segment
        int version = 0;              // invalidates stale departures
    }

    // CUSTOMER CLASS

    /**
     * Customer for SJF queue:
     *  - arrivalTime: arrival instant
     *  - deadlineTime: absolute expiry time
     *  - lengthBits: packet/job size L in bits (priority key + service time if LR mode)
     *
     * NEW (for preemption):
     *  - remainingService: remaining service time so preempted jobs resume (no restart)
     *  - lastQueueEnterTime: for correct queue delay accumulation across multiple queue entries
     */
    private static class Customer {
        double arrivalTime;
        double deadlineTime;
        int lengthBits;

        // NEW
        double remainingService = Double.NaN;
        double lastQueueEnterTime = Double.NaN;

        Customer(double arrivalTime, double deadlineTime, int lengthBits) {
            this.arrivalTime  = arrivalTime;
            this.deadlineTime = deadlineTime;
            this.lengthBits   = lengthBits;
        }
    }

    // SJF COMPARATOR

    /**
     * Orders customers by shortest lengthBits; ties broken by earlier arrival.
     * NOTE: In preemptive mode we use remainingService at runtime for preemption decisions,
     * but the waiting queue still uses lengthBits as the SJF priority key (baseline behavior).
     */
    private static final Comparator<Customer> SJF_COMPARATOR =
            Comparator.<Customer>comparingDouble(c -> c.remainingService)
                    .thenComparingDouble(c -> c.arrivalTime);

    // SUMMARY CLASS

    private static class Summary {
        double avgNs, avgNq, avgDs, avgDq;
        double utilization;
        int    servedCompleted;   // number of actual service completions
        int    arrivals;          // number of admitted arrivals (== RUN_COMPLETED at end)
        double pExpired;          // per FIFO definition (expired before service start / arrivals)
    }

    // SUMMARY ROW

    private static class SummaryRow {
        double rho, lambda, mu;
        int    c;

        double Ns_sim;
        double Nq_sim;
        double Ds_sim;
        double Dq_sim;

        double pExpired_sim;
    }

    // SAMPLING HELPERS

    private static double expSample(double rate) {
        double u = RNG.nextDouble();
        if (u == 0.0) u = 1e-16;
        return -Math.log(1.0 - u) / rate;
    }

    /**
     * Geometric(p) with support {1,2,...}.
     */
    private static int geometricSample(double p) {
        double u = RNG.nextDouble();
        if (u == 0.0) u = 1e-16;
        double oneMinusP = 1.0 - p;
        int k = (int) Math.ceil(Math.log(1.0 - u) / Math.log(oneMinusP));
        if (k < 1) k = 1;
        return k;
    }

    /**
     * Interarrival sampling:
     * - Baseline: Exp(lambda)
     * - DDC: 1/lambda
     */
    private static double sampleInterarrival(double lambda) {
        if (DETERMINISTIC_ARRIVALS) return 1.0 / lambda;
        return expSample(lambda);
    }

    /**
     * Service time sampling:
     * - If LENGTH_BASED_SERVICE: L / R
     * - Else:
     *    - DDC: 1/mu
     *    - Baseline: Exp(mu)
     */
    private static double serviceTimeFor(Customer cust, double mu) {
        if (LENGTH_BASED_SERVICE) {
            return cust.lengthBits / R;
        }
        if (DETERMINISTIC_SERVICE) return 1.0 / mu;
        return expSample(mu);
    }

    // FIND IDLE SERVER

    private static int findIdleServer(Server[] servers) {
        for (int i = 0; i < servers.length; i++) {
            if (!servers[i].busy) return i;
        }
        return -1;
    }

    // HELPER: current remaining service on a server at time now

    private static double currentRemainingOn(Server s, double now) {
        if (s == null || s.current == null) return Double.POSITIVE_INFINITY;
        double elapsed = now - s.serviceStartTime;
        if (elapsed < 0) elapsed = 0;
        double rem = s.current.remainingService - elapsed;
        if (rem < 0) rem = 0;
        return rem;
    }

    // HELPER: start service on server sid for job cust at time now
    // - Applies expiry-at-service-start rule
    // - Drops expired (if DROP_EXPIRED) and returns false (server remains idle)
    // - Otherwise schedules departure using remainingService and increments version
    // - Accumulates queue delay using lastQueueEnterTime (for correct Dq under preemption)

    private static boolean startService(
            int sid,
            Customer cust,
            double now,
            double mu,
            PriorityQueue<Event> fel,
            Server[] servers,
            MutableCounters ctr
    ) {
        // Expiry check at service start (FIFO definition; no expiry upon arrival)
        if (now >= cust.deadlineTime) {
            ctr.expiredBeforeServiceStart++;

            if (DROP_EXPIRED) {
                // drop immediately at service start
                ctr.Ns--;
                if (ctr.Ns < 0) throw new IllegalStateException("Ns < 0 after drop");
                return false;
            }
        }

        // Queue delay accumulation (works even if job was queued multiple times)
        if (!Double.isNaN(cust.lastQueueEnterTime)) {
            ctr.totalQueueDelay += now - cust.lastQueueEnterTime;
            cust.lastQueueEnterTime = Double.NaN;
        }

        // Sample service once per job (needed to resume after preemption)
        if (Double.isNaN(cust.remainingService)) {
            cust.remainingService = serviceTimeFor(cust, mu);
        }

        Server s = servers[sid];
        s.busy = true;
        s.current = cust;
        s.serviceStartTime = now;

        // Version bump invalidates any old departure for this server (preemption-safe)
        s.version++;
        int v = s.version;

        // Schedule departure based on remaining service
        fel.add(new ServiceEvent(now + cust.remainingService, sid, cust.arrivalTime, v));
        return true;
    }

    // HELPER: choose a preemption victim server for arriving job
    // - Preempt if arriving remaining service < some running remaining service
    // - Choose the running job with the LARGEST remaining service (worst) to preempt

    private static int choosePreemptionVictim(Server[] servers, Customer arriving, double now) {
        int victim = -1;
        double worstRem = Double.NEGATIVE_INFINITY;

        for (int i = 0; i < servers.length; i++) {
            Server s = servers[i];
            if (!s.busy || s.current == null) continue;

            double rem = currentRemainingOn(s, now);

            // Preempt only if arriving is shorter (smaller remaining service)
            if (arriving.remainingService < rem) {
                if (rem > worstRem) {
                    worstRem = rem;
                    victim = i;
                }
            }
        }
        return victim;
    }

    // HELPER: preempt server sid at time now and enqueue the running job back to sjfQ
    // - Subtract elapsed service from remainingService
    // - Put job back to queue with lastQueueEnterTime = now
    // - Do NOT change Ns (job stays in system)
    // - Increment Nq

    private static void preemptToQueue(
            int sid,
            double now,
            PriorityQueue<Customer> sjfQ,
            Server[] servers,
            MutableCounters ctr
    ) {
        Server s = servers[sid];
        if (s.current == null) return;

        Customer running = s.current;

        // Update remaining service
        double rem = currentRemainingOn(s, now);
        running.remainingService = rem;

        // Enqueue back
        running.lastQueueEnterTime = now;
        sjfQ.add(running);
        ctr.Nq++;

        // Make server idle; old departure becomes stale due to version change when new service starts
        s.busy = false;
        s.current = null;
    }

    // MUTABLE COUNTERS (to keep code changes minimal)

    private static class MutableCounters {
        int Ns;
        int Nq;

        double totalSystemDelay;
        double totalQueueDelay;

        int arrivals;
        int servedCompleted;

        int expiredBeforeServiceStart;
    }

    // SINGLE RUN SIMULATOR (FIFO CONVENTIONS APPLIED)

    private static Summary simulate(double lambda,
                                    double mu,
                                    int c,
                                    double deadlineRate) {

        PriorityQueue<Event> fel = new PriorityQueue<>();
        PriorityQueue<Customer> sjfQ = new PriorityQueue<>(SJF_COMPARATOR);

        Server[] servers = new Server[c];
        for (int i = 0; i < c; i++) servers[i] = new Server();

        double t = 0.0;

        double areaNs = 0.0;
        double areaNq = 0.0;
        double[] serverBusyArea = new double[c];

        MutableCounters ctr = new MutableCounters();
        ctr.Ns = 0;
        ctr.Nq = 0;
        ctr.arrivals = 0;
        ctr.servedCompleted = 0;
        ctr.totalSystemDelay = 0.0;
        ctr.totalQueueDelay  = 0.0;
        ctr.expiredBeforeServiceStart = 0;

        // First arrival at t = 0
        fel.add(new Event(0.0, EventType.ARRIVAL, -1));

        // FIFO baseline drain condition
        while (ctr.arrivals < RUN_COMPLETED || ctr.Ns > 0) {

            Event e = fel.poll();
            double oldT = t;
            t = e.time;

            // ---- time integrals ----
            double dt = t - oldT;
            areaNs += ctr.Ns * dt;
            areaNq += ctr.Nq * dt;
            for (int i = 0; i < c; i++) {
                if (servers[i].busy) serverBusyArea[i] += dt;
            }

            if (e.type == EventType.ARRIVAL) {

                // Strict arrivals cap (FIFO convention)
                if (ctr.arrivals >= RUN_COMPLETED) {
                    continue;
                }
                ctr.arrivals++;

                // Schedule next arrival only if still under cap
                if (ctr.arrivals < RUN_COMPLETED) {
                    double inter = sampleInterarrival(lambda);
                    fel.add(new Event(t + inter, EventType.ARRIVAL, -1));
                }

                // Create customer (no expiry upon arrival assumption preserved)
                double deadlineTime = t + expSample(deadlineRate);
                int lengthBits = geometricSample(GEO_P);
                Customer cust = new Customer(t, deadlineTime, lengthBits);

                // Sample service once per job (needed for preemption comparisons)
                // (This is safe even in non-preemptive mode; it does not change the model.)
                cust.remainingService = serviceTimeFor(cust, mu);

                // Admit into system
                ctr.Ns++;

                int idle = findIdleServer(servers);
                if (idle >= 0) {
                    // Immediate service start at time t
                    Customer candidate = cust;

                    // If candidate is dropped (DROP_EXPIRED), try to fill idle server from queue
                    while (candidate != null) {
                        boolean served = startService(idle, candidate, t, mu, fel, servers, ctr);
                        if (served) break;

                        // Candidate dropped => try another from queue
                        if (!sjfQ.isEmpty()) {
                            candidate = sjfQ.poll();
                            ctr.Nq--;
                        } else {
                            candidate = null;
                        }
                    }

                } else {
                    // All servers busy:
                    // - If PREEMPTIVE_SJF enabled, potentially preempt based on remaining service time
                    // - Otherwise, join SJF queue
                    if (PREEMPTIVE_SJF) {

                        // Choose victim if arriving is shorter than some running job
                        int victim = choosePreemptionVictim(servers, cust, t);

                        if (victim >= 0) {
                            // Preempt the victim job back into queue with remaining service
                            preemptToQueue(victim, t, sjfQ, servers, ctr);

                            // Start serving arriving job on the victim server.
                            Customer candidate = cust;

                            while (candidate != null) {
                                boolean served = startService(victim, candidate, t, mu, fel, servers, ctr);
                                if (served) break;

                                // Candidate dropped => fill from queue (SJF order)
                                if (!sjfQ.isEmpty()) {
                                    candidate = sjfQ.poll();
                                    ctr.Nq--;
                                } else {
                                    candidate = null;
                                }
                            }

                        } else {
                            // No preemption => enqueue arriving job
                            cust.lastQueueEnterTime = t;
                            sjfQ.add(cust);
                            ctr.Nq++;
                        }

                    } else {
                        // Baseline: join SJF queue
                        cust.lastQueueEnterTime = t;
                        sjfQ.add(cust);
                        ctr.Nq++;
                    }
                }

            } else { // DEPARTURE

                ServiceEvent se = (ServiceEvent) e;
                int sid = se.serverId;

                // Ignore stale departure events (after preemption/reschedule)
                if (se.version != servers[sid].version) {
                    continue;
                }

                // Service completion (of the current job on this server)
                ctr.servedCompleted++;
                ctr.Ns--;
                if (ctr.Ns < 0) throw new IllegalStateException("Ns < 0 after completion");

                servers[sid].busy = false;

                // System delay for served job
                double Ds_i = t - se.arrivalTime;
                ctr.totalSystemDelay += Ds_i;

                // Clear server current
                servers[sid].current = null;

                // Now try to start service for next queued jobs (may drop expired ones)
                while (!sjfQ.isEmpty()) {
                    Customer next = sjfQ.poll();
                    ctr.Nq--;

                    boolean served = startService(sid, next, t, mu, fel, servers, ctr);
                    if (served) {
                        break;
                    }
                    // If dropped, keep looping to see if another job can start
                }
            }
        }

        Summary S = new Summary();
        S.arrivals = ctr.arrivals;
        S.servedCompleted = ctr.servedCompleted;

        double totalTime = t;
        S.avgNs = areaNs / totalTime;
        S.avgNq = areaNq / totalTime;

        double busySum = 0.0;
        for (int i = 0; i < c; i++) busySum += serverBusyArea[i];
        S.utilization = busySum / (c * totalTime);

        // NOTE: delays are averages over served completions (like FIFO baseline)
        S.avgDs = ctr.totalSystemDelay / ctr.arrivals;
        S.avgDq = ctr.totalQueueDelay / ctr.arrivals;

        // FIFO definition: expired BEFORE service start / arrivals
        S.pExpired = (ctr.arrivals > 0) ? (ctr.expiredBeforeServiceStart / (double) ctr.arrivals) : 0.0;

        System.out.println(ctr.arrivals + " " + ctr.servedCompleted);

        return S;
    }

    // CSV WRITER

    private static void writeCSV(List<SummaryRow> rows) {
        try {
            File f = new File(CSV_FILE);
            f.getParentFile().mkdirs();

            PrintWriter pw = new PrintWriter(new FileWriter(f));

            // IMPORTANT: keep the SAME CSV columns you specified
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

    // SUMMARY PRINTER

    private static void printSummary(List<SummaryRow> rows) {
        System.out.println("==============================================================");
        System.out.println(" c | rho |   Ns   |   Nq   |   Ds   |   Dq   | pExpired");
        System.out.println("==============================================================");
        for (SummaryRow r : rows) {
            System.out.printf(Locale.US,
                    "%2d | %.2f | %6.4f | %6.4f | %6.4f | %6.4f | %7.4f%n",
                    r.c, r.rho,
                    r.Ns_sim,
                    r.Nq_sim,
                    r.Ds_sim,
                    r.Dq_sim,
                    r.pExpired_sim);
        }
        System.out.println("==============================================================");
        System.out.println("CSV written to: " + CSV_FILE);
        System.out.println("Deterministic arrivals: " + DETERMINISTIC_ARRIVALS
                + " | Deterministic service: " + DETERMINISTIC_SERVICE
                + " | Length-based service (L/R): " + LENGTH_BASED_SERVICE
                + " | DROP_EXPIRED: " + DROP_EXPIRED
                + " | PREEMPTIVE_SJF: " + PREEMPTIVE_SJF);
        if (LENGTH_BASED_SERVICE) {
            System.out.println("R (bits/time) = " + R + "  (derived as MU * MEAN_LENGTH_BITS)");
        }
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

            R.pExpired_sim = S.pExpired;

            rows.add(R);

            // System.out.printf(Locale.US, "SJF: rho=%.2f -> Ns=%.4f, pExpired=%.6f%n", rho, S.avgNs, S.pExpired);
        }

        writeCSV(rows);
        printSummary(rows);
    }
}
