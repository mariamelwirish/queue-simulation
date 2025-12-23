import java.util.*;
import java.io.*;

/**
 * M/M/c queue simulator with deadlines and non-preemptive EDF service.
 *
 * BASELINE MODEL:
 * - Exponential interarrival times (rate lambda)
 * - Exponential service times (rate mu)
 * - Exponential deadlines (single-class)
 *
 * EXTENSION: TWO DEADLINE CLASSES
 * - Short–tight deadlines (Exp(DEADLINE_RATE_TIGHT))
 * - Long–loose deadlines (Exp(DEADLINE_RATE_LOOSE))
 * - Class chosen at arrival with probability PROB_TIGHT
 *
 * EDF always schedules by earliest absolute deadline time.
 *
 * FIFO-INHERITED RULES:
 * - No expiration upon arrival
 * - Arrival capping at RUN_COMPLETED
 * - Drain phase after arrivals stop
 * - pExpired = expired-before-service / total arrivals
 *
 * ADDITION (DROP MODE + SERVED PROBABILITY):
 * - DROP_EXPIRED controls whether expired jobs are dropped at service start
 * - pServed = completed (served) / total arrivals
 *
 * NEW: PREEMPTIVE MODE (optional)
 * - PREEMPTIVE_EDF controls whether a newly arriving job can preempt service
 *   on a busy server if it has an earlier absolute deadline than a running job.
 * - When preempting, the interrupted job is returned to the EDF queue with its
 *   remaining service time (so it does NOT restart from the beginning).
 * - Stale DEPARTURE events are ignored using a per-server "version" token.
 */
public class edf {

   
    // CONTROL PARAMETERS
   

    static final int RUN_COMPLETED = 1000000;

    static final double RHO_START = 0.05;
    static final double RHO_END   = 0.95;
    static final double RHO_STEP  = 0.05;

    static final double MU = 1;
    static final int C = 1;

   
    // DEADLINE PARAMETERS
   

    static final double DEADLINE_RATE = 0.1;

    static final boolean TWO_DEADLINE_CLASSES = true;

    static final double DEADLINE_RATE_TIGHT = 0.20;
    static final double DEADLINE_RATE_LOOSE = 0.01;

    static final double PROB_TIGHT = 0.5;

   
    // DETERMINISTIC MODE FLAGS
   

    static final boolean DETERMINISTIC_ARRIVAL = false;
    static final boolean DETERMINISTIC_SERVICE = false;

   
    // DROP MODE FLAG (NEW)
   
    // If true: jobs that are expired at service start are dropped (not served).
    // If false: expired jobs are still served (baseline behavior).
    static final boolean DROP_EXPIRED = false;

   
    // PREEMPTIVE MODE FLAG (NEW)
   
    // If true: preemptive EDF is enabled (a higher-priority arrival can preempt).
    // If false: baseline non-preemptive EDF behavior.
    static final boolean PREEMPTIVE_EDF = false;

   
    // OUTPUT CSV
   

    static final String CSV_FILE;
    static {
        // Build a mode-aware filename to avoid overwriting.
        if (PREEMPTIVE_EDF) {
            CSV_FILE = "../results/csv/preemptive_edf.csv";
        }
        else if (DROP_EXPIRED) {
            CSV_FILE = "../results/csv/drop_edf.csv";
        } else if (DETERMINISTIC_ARRIVAL || DETERMINISTIC_SERVICE) {
            CSV_FILE = "../results/csv/ddc_edf.csv";
        } else if(TWO_DEADLINE_CLASSES) {
            CSV_FILE = "../results/csv/2dc_edf.csv";
        } else {
            CSV_FILE = "../results/csv/mmc_edf.csv";
        }
    }

   
    // RNG
   

    static final Random RNG = new Random(12345);

   
    // EVENT TYPES
   

    private static enum EventType { ARRIVAL, DEPARTURE }

   
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

        public int compareTo(Event other) {
            return Double.compare(this.time, other.time);
        }
    }

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

   
    // SERVER
   

    private static class Server {
        boolean busy = false;

        // NEW: currently running job on this server (needed for preemption)
        Customer current = null;

        // NEW: service segment start time (needed to compute elapsed on preemption)
        double serviceStartTime = 0.0;

        // NEW: increments whenever this server starts a (new) service segment
        //      (used to invalidate old ServiceEvent departure events)
        int version = 0;
    }

   
    // CUSTOMER
   

    private static long NEXT_ID = 1;

    private static class Customer {
        double arrivalTime;
        double deadlineTime;
        boolean isTight;

        // NEW: remaining service time so preempted jobs resume without restarting
        double remainingService = Double.NaN; // NaN => not yet sampled

        // NEW: for correct queue delay accumulation with multiple queue entries
        double lastQueueEnterTime = Double.NaN;

        // Optional: unique id for debugging / sanity checks
        long id;

        Customer(double arrivalTime, double deadlineTime, boolean isTight) {
            this.arrivalTime  = arrivalTime;
            this.deadlineTime = deadlineTime;
            this.isTight = isTight;
            this.id = NEXT_ID++;
        }
    }

   
    // EDF COMPARATOR
   

    private static final Comparator<Customer> EDF_COMPARATOR =
            Comparator.<Customer>comparingDouble(c -> c.deadlineTime)
                      .thenComparingDouble(c -> c.arrivalTime);

   
    // SAMPLING
   

    private static double expSample(double rate) {
        double u = RNG.nextDouble();
        if (u == 0.0) u = 1e-16;
        return -Math.log(1.0 - u) / rate;
    }

    private static double sampleInterarrival(double lambda) {
        if (DETERMINISTIC_ARRIVAL) return 1.0 / lambda;
        return expSample(lambda);
    }

    private static double sampleService(double mu) {
        if (DETERMINISTIC_SERVICE) return 1.0 / mu;
        return expSample(mu);
    }

    private static int findIdleServer(Server[] servers) {
        for (int i = 0; i < servers.length; i++)
            if (!servers[i].busy) return i;
        return -1;
    }

   
    // HELPER: start service on a given server for a given job at time now.
    // - Applies deadline-at-service-start rule (no expiration upon arrival).
    // - If DROP_EXPIRED is enabled, may drop and return false (server remains idle).
    // - Otherwise schedules a departure using remaining service time.
    // - Handles queue delay using lastQueueEnterTime (for jobs taken from queue).
   

    private static boolean startServiceOnServer(
            int serverId,
            Customer job,
            double now,
            double mu,
            PriorityQueue<Event> fel,
            Server[] servers,
            MutableCounters ctr
    ) {
        // Deadline check at service start (FIFO-inherited: no check upon arrival)
        if (now >= job.deadlineTime) {
            ctr.expiredBeforeService++;

            if (DROP_EXPIRED) {
                // Drop at service start
                ctr.dropped++;
                ctr.Ns--;               // Ns decrements on drop (FIFO convention)
                return false;           // server stays idle; caller can try next job
            }
        }

        // If job is coming from queue (or preempted back into queue), accumulate queue delay
        // only when it actually starts/resumes service.
        if (!Double.isNaN(job.lastQueueEnterTime)) {
            ctr.totalQueueDelay += now - job.lastQueueEnterTime;
            job.lastQueueEnterTime = Double.NaN; // reset: it's now in service
        }

        // Sample service time once per job (if not previously sampled).
        if (Double.isNaN(job.remainingService)) {
            job.remainingService = sampleService(mu);
        }

        // Mark server busy and record state for potential future preemption.
        servers[serverId].busy = true;
        servers[serverId].current = job;
        servers[serverId].serviceStartTime = now;

        // Bump version to invalidate any previous departure for this server.
        servers[serverId].version++;
        int v = servers[serverId].version;

        // Schedule departure based on remaining service (resume-able).
        double departTime = now + job.remainingService;
        fel.add(new ServiceEvent(departTime, serverId, job.arrivalTime, v));

        ctr.servedStarted++;
        return true;
    }

   
    // HELPER: pick victim server for preemption.
    // Choose the currently running job with the LATEST deadline (worst priority),
    // among those that have deadline > arriving.deadline (i.e., can be preempted).
    // Returns victim server index, or -1 if no preemption should happen.
   

    private static int choosePreemptionVictim(Server[] servers, Customer arriving) {
        int victim = -1;
        double worstDeadline = Double.NEGATIVE_INFINITY;

        for (int i = 0; i < servers.length; i++) {
            if (!servers[i].busy || servers[i].current == null) continue;

            Customer running = servers[i].current;

            // Only preempt if arriving is higher priority (earlier deadline)
            if (arriving.deadlineTime < running.deadlineTime) {
                // pick the "worst" running job (largest deadline) to preempt
                if (running.deadlineTime > worstDeadline) {
                    worstDeadline = running.deadlineTime;
                    victim = i;
                }
            }
        }
        return victim;
    }

   
    // HELPER: preempt current job on a server at time now, enqueue it back into edfQ.
    // Adjust remaining service time so it resumes where it left off.
    // Does NOT touch Ns (preemption keeps the job in the system).
    // Increments Nq and sets lastQueueEnterTime for queue delay tracking.
   

    private static void preemptServerToQueue(
            int serverId,
            double now,
            PriorityQueue<Customer> edfQ,
            Server[] servers,
            MutableCounters ctr
    ) {
        Customer running = servers[serverId].current;
        if (running == null) return;

        // Compute elapsed service in this segment
        double elapsed = now - servers[serverId].serviceStartTime;
        if (elapsed < 0) elapsed = 0; // safety

        // Reduce remaining service (do not restart from beginning)
        running.remainingService -= elapsed;
        if (running.remainingService < 0) running.remainingService = 0;

        // Put back into EDF queue
        running.lastQueueEnterTime = now;
        edfQ.add(running);
        ctr.Nq++;

        // Make server idle (current job removed)
        servers[serverId].busy = false;
        servers[serverId].current = null;

        // Note: we do NOT decrement Ns here (job still in system).
        // Note: old departure event becomes stale due to version bump when new service starts.
    }

   
    // MUTABLE COUNTERS (to avoid changing too much of your local variables)
   

    private static class MutableCounters {
        int Ns;
        int Nq;

        int expiredBeforeService;
        int dropped;
        int servedStarted;

        double totalSystemDelay;
        double totalQueueDelay;

        int completed;
    }

   
    // SIMULATION
   

    private static Summary simulate(double lambda, double mu, int c) {

        PriorityQueue<Event> fel = new PriorityQueue<>();
        PriorityQueue<Customer> edfQ = new PriorityQueue<>(EDF_COMPARATOR);

        Server[] servers = new Server[c];
        for (int i = 0; i < c; i++) servers[i] = new Server();

        double t = 0.0;

        int arrivals = 0;

        MutableCounters ctr = new MutableCounters();
        ctr.Ns = 0;
        ctr.Nq = 0;

        double areaNs = 0.0, areaNq = 0.0;
        double[] serverBusyArea = new double[c];

        // First arrival
        fel.add(new Event(0.0, EventType.ARRIVAL, -1));

        while (arrivals < RUN_COMPLETED || ctr.Ns > 0) {

            Event e = fel.poll();
            if (e == null) break; // safety

            double oldT = t;
            t = e.time;

            double dt = t - oldT;
            areaNs += ctr.Ns * dt;
            areaNq += ctr.Nq * dt;
            for (int i = 0; i < c; i++)
                if (servers[i].busy)
                    serverBusyArea[i] += dt;

            if (e.type == EventType.ARRIVAL) {

                if (arrivals >= RUN_COMPLETED)
                    continue;

                arrivals++;
                ctr.Ns++;

                // Schedule next arrival ONLY if we still need more arrivals
                if (arrivals < RUN_COMPLETED)
                    fel.add(new Event(t + sampleInterarrival(lambda),
                                      EventType.ARRIVAL, -1));

                // Choose deadline class (if enabled)
                boolean isTight = false;
                double rate;

                if (TWO_DEADLINE_CLASSES) {
                    isTight = RNG.nextDouble() < PROB_TIGHT;
                    rate = isTight ? DEADLINE_RATE_TIGHT : DEADLINE_RATE_LOOSE;
                } else {
                    rate = DEADLINE_RATE;
                }

                // Assign deadline (still no expiration check at arrival)
                double deadlineTime = t + expSample(rate);
                Customer cust = new Customer(t, deadlineTime, isTight);

                int idle = findIdleServer(servers);
                if (idle >= 0) {
                    // Start service immediately on idle server.
                    // If dropped, try pulling from queue (EDF) to fill the idle server.
                    Customer candidate = cust;

                    while (candidate != null) {
                        boolean served = startServiceOnServer(
                                idle, candidate, t, mu, fel, servers, ctr
                        );

                        if (served) break;

                        // Candidate was dropped (DROP_EXPIRED=true). Try another from queue.
                        if (!edfQ.isEmpty()) {
                            candidate = edfQ.poll();
                            ctr.Nq--;
                        } else {
                            candidate = null;
                        }
                    }

                } else {
                    // All servers busy: either preempt (if enabled) or enqueue.
                    if (PREEMPTIVE_EDF) {
                        int victim = choosePreemptionVictim(servers, cust);

                        if (victim >= 0) {
                            // Preempt victim, enqueue its remaining work, then start cust on that server.
                            preemptServerToQueue(victim, t, edfQ, servers, ctr);

                            // Now attempt to start service with the arriving job on victim server.
                            // If arriving gets dropped at service start (DROP_EXPIRED), then try to fill
                            // the server from the EDF queue.
                            Customer candidate = cust;

                            while (candidate != null) {
                                boolean served = startServiceOnServer(
                                        victim, candidate, t, mu, fel, servers, ctr
                                );

                                if (served) break;

                                // Candidate dropped -> fill from queue
                                if (!edfQ.isEmpty()) {
                                    candidate = edfQ.poll();
                                    ctr.Nq--;
                                } else {
                                    candidate = null;
                                }
                            }

                        } else {
                            // No valid preemption => enqueue arriving job
                            cust.lastQueueEnterTime = t;
                            edfQ.add(cust);
                            ctr.Nq++;
                        }

                    } else {
                        // Baseline: join EDF queue (non-preemptive)
                        cust.lastQueueEnterTime = t;
                        edfQ.add(cust);
                        ctr.Nq++;
                    }
                }

            } else { // DEPARTURE

                ServiceEvent se = (ServiceEvent) e;

                // Ignore stale departure events (after preemption/reschedule)
                if (se.version != servers[se.serverId].version) {
                    continue;
                }

                // Complete the current job on this server
                Customer done = servers[se.serverId].current;

                // Safety: if something went wrong, skip
                if (done == null) {
                    servers[se.serverId].busy = false;
                    continue;
                }

                ctr.completed++;
                ctr.Ns--;

                // System delay for served job
                ctr.totalSystemDelay += t - done.arrivalTime;

                // Clear server
                servers[se.serverId].busy = false;
                servers[se.serverId].current = null;

                // Start next service(s) from queue, skipping dropped expired jobs if enabled
                while (!edfQ.isEmpty()) {

                    Customer next = edfQ.poll();
                    ctr.Nq--;

                    boolean served = startServiceOnServer(
                            se.serverId, next, t, mu, fel, servers, ctr
                    );

                    if (served) {
                        break; // server is busy again
                    }
                    // If dropped, continue to try another queued job
                }
            }
        }

        Summary S = new Summary();
        S.avgNs = areaNs / t;
        S.avgNq = areaNq / t;
        S.avgDs = ctr.totalSystemDelay / arrivals;
        S.avgDq = ctr.totalQueueDelay / arrivals;


        double busy = 0.0;
        for (double b : serverBusyArea) busy += b;
        S.utilization = busy / (c * t);

        S.arrivals = arrivals;
        S.completed = ctr.completed;
        S.dropped = ctr.dropped;

        // FIFO-inherited definition
        S.pExpired = ctr.expiredBeforeService / (double) arrivals;

        // NEW: probability of being served (completed / arrivals)
        S.pServed = ctr.completed / (double) arrivals;

        System.out.println(arrivals + " " + ctr.completed);

        return S;
    }

   
    // SUMMARY STRUCT
   

    private static class Summary {
        double avgNs, avgNq, avgDs, avgDq;
        double utilization;

        int arrivals;
        int completed;
        int dropped;

        double pExpired;  // expired-before-service / arrivals
        double pServed;   // served (completed) / arrivals
    }

   
    // MAIN + CSV WRITER (kept minimal; adapt to your existing pattern if needed)
   

    public static void main(String[] args) throws Exception {

        // Write header
        try (PrintWriter pw = new PrintWriter(new FileWriter(CSV_FILE))) {
           pw.println("rho,lambda,mu,c,"
                    + "Ns_sim,"
                    + "Nq_sim,"
                    + "Ds_sim,"
                    + "Dq_sim,"
                    + "pExpired_sim");
        }

        for (double rho = RHO_START; rho <= RHO_END + 1e-12; rho += RHO_STEP) {

            double lambda = rho * C * MU;

            Summary S = simulate(lambda, MU, C);

            try (PrintWriter out = new PrintWriter(new FileWriter(CSV_FILE, true))) {
                out.printf(Locale.US,
                                "%.4f,%.6f,%.6f,%d,"
                                + "%.6f,"
                                + "%.6f,"
                                + "%.6f,"
                                + "%.6f,"
                                + "%.6f%n",
                    rho, lambda, MU, C,
                    S.avgNs, S.avgNq,
                    S.avgDs, S.avgDq,
                    S.pExpired
                );

            }
        }

        System.out.println("Done. Wrote: " + CSV_FILE);
    }
}
