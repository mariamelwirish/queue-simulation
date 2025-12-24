import java.util.*;
import java.io.*;

public class sjf {

    // CONTROL VARIABLES
    // Number of arrivals.

    static final int RUN_COMPLETED = 1000000; 

    // rho sweep range.
    static final double RHO_START = 0.05;
    static final double RHO_END   = 0.95;
    static final double RHO_STEP  = 0.05;

    // Service rate and number of servers.
    static final double MU = 1.0;
    static final int C = 5;

    // Probability of small-tight packets.
    static final double PROB_TIGHT = 0.5;

    // Deadline ranges classes.
    static final double DEADLINE_RATE_TIGHT = 0.20;
    static final double DEADLINE_RATE_LOOSE = 0.01;

    // Packet size classes.
    static final double GEOM_P_TIGHT = 0.40; 
    static final double GEOM_P_LOOSE = 0.10; 

    // Mean length and rate for length-based service.
    // E[L] for Geometric(p) on {1,2,...} is 1/p
    static final double MEAN_L =
            PROB_TIGHT * (1.0 / GEOM_P_TIGHT) + (1.0 - PROB_TIGHT) * (1.0 / GEOM_P_LOOSE);

    // Preserve mean service time: serviceTime = L / R  => R = MU * E[L]
    static final double R = MU * MEAN_L;

    // Flags to control simulation modes.
    static final boolean DETERMINISTIC_ARRIVAL = false; // Deterministic arrival process
    static final boolean DETERMINISTIC_SERVICE = false; // Deterministic service times
    static final boolean LENGTH_BASED_SERVICE = true;     // Size-based service time: serviceTime = L / R
    static final boolean DROP_EXPIRED = true; // Drop expired customers at service start
    static final boolean PREEMPTIVE_SJF = false; // SJF preemptive mode

    // CSV output file (mode-aware).
    static final String CSV_FILE;
    static {
        // Build a mode-aware filename to avoid overwriting.
        if (PREEMPTIVE_SJF) {
            CSV_FILE = "../results/csv/preemptive_sjf.csv";
        } else if (DROP_EXPIRED) {
            CSV_FILE = "../results/csv/drop_sjf.csv";
        } else if (DETERMINISTIC_ARRIVAL || DETERMINISTIC_SERVICE) {
            CSV_FILE = "../results/csv/ddc_sjf.csv";
        } else {
            CSV_FILE = "../results/csv/mmc_sjf.csv";
        }
    }

    /*******************************************************************************************/
    // SIMULATION STRUCTURE HELPERS

    // Random number generators
    static final Random RNG_PKT = new Random(12345); // for arrivals, deadlines, sizes
    static final Random RNG_SVC = new Random(54321); // for service times

    // Event types.

    private static enum EventType { ARRIVAL, DEPARTURE }

    // Event class.

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

    // Service event class.
    private static class ServiceEvent extends Event {
        int version;

        ServiceEvent(double time, int serverId, int version) {
            super(time, EventType.DEPARTURE, serverId);
            this.version = version;
        }
    }

    // Server class.

    private static class Server {
        boolean busy = false;
        Customer current = null;
        double serviceStartTime = 0.0;
        int version = 0;
    }

    // Customer (packet) class.

    private static long NEXT_ID = 1;

    private static class Customer {
        double arrivalTime;
        double deadlineTime;

        boolean isTight;
        int lengthBits;

        double remainingService = Double.NaN;
        double lastQueueEnterTime = Double.NaN;

        long id;

        Customer(double arrivalTime, double deadlineTime,
                 boolean isTight, int lengthBits) {
            this.arrivalTime = arrivalTime;
            this.deadlineTime = deadlineTime;
            this.isTight = isTight;
            this.lengthBits = lengthBits;
            this.id = NEXT_ID++;
        }
    }

    // SJF comparator.
    // If true M/M/c: compare by lengthBits label.
    // Else: compare by remainingService (assigned at arrival / when queued).

    private static final Comparator<Customer> SJF_COMPARATOR = (a, b) -> {
        boolean MMC_MODE = (!DETERMINISTIC_SERVICE && !LENGTH_BASED_SERVICE);
        if (MMC_MODE) {
            int cmp = Integer.compare(a.lengthBits, b.lengthBits);
            if (cmp != 0) return cmp;
            return Double.compare(a.arrivalTime, b.arrivalTime);
        } else {
            int cmp = Double.compare(a.remainingService, b.remainingService);
            if (cmp != 0) return cmp;
            return Double.compare(a.arrivalTime, b.arrivalTime);
        }
    };

    // exponential(rate) sample using inverse-CDF.

    private static double expSample(double rate, Random rng) {
        double u = rng.nextDouble();
        if (u == 0.0) u = 1e-16;
        return -Math.log(1.0 - u) / rate;
    }

    // Geometric(p) sample using inverse-CDF.
    private static int geometricSample(double p, Random rng) {
        // X = ceil( log(1-U) / log(1-p) )
        double u = rng.nextDouble();
        if (u == 0.0) u = 1e-16;
        return (int) Math.ceil(Math.log(1.0 - u) / Math.log(1.0 - p));
    }

    private static double sampleInterarrival(double lambda) {
        return DETERMINISTIC_ARRIVAL ? 1.0 / lambda : expSample(lambda, RNG_PKT);
    }

    private static double sampleService(double mu) {
        return DETERMINISTIC_SERVICE ? 1.0 / mu : expSample(mu, RNG_SVC);
    }

    private static double serviceTimeFor(Customer job, double mu) {
        if (LENGTH_BASED_SERVICE) {
            return job.lengthBits / R;
        }
        return sampleService(mu);
    }

    // Find an idle server index, or -1 if none.
    private static int findIdleServer(Server[] servers) {
        for (int i = 0; i < servers.length; i++)
            if (!servers[i].busy) return i;
        return -1;
    }

    // Preemption helpers.
    // SRPT-like in non-MMC; in MMC we preempt based on lengthBits label.
    // (Kept same structure as EDF victim selection.)

    private static double currentRemainingOn(Server s, double now) {
        if (s == null || s.current == null) return Double.POSITIVE_INFINITY;
        double elapsed = now - s.serviceStartTime;
        if (elapsed < 0) elapsed = 0;
        double rem = s.current.remainingService - elapsed;
        if (rem < 0) rem = 0;
        return rem;
    }

    private static int choosePreemptionVictim(Server[] servers, Customer arriving, double now) {
        boolean MMC_MODE = (!DETERMINISTIC_SERVICE && !LENGTH_BASED_SERVICE);

        int victim = -1;

        if (MMC_MODE) {
            // Preempt a running job if arriving has smaller lengthBits label
            int worstLen = Integer.MIN_VALUE;

            for (int i = 0; i < servers.length; i++) {
                if (!servers[i].busy || servers[i].current == null) continue;

                int runLen = servers[i].current.lengthBits;
                if (arriving.lengthBits < runLen) {
                    if (runLen > worstLen) {
                        worstLen = runLen;
                        victim = i;
                    }
                }
            }
            return victim;
        } else {
            // Non-MMC: SRPT-like compare remaining service
            double worstRem = Double.NEGATIVE_INFINITY;

            for (int i = 0; i < servers.length; i++) {
                if (!servers[i].busy || servers[i].current == null) continue;

                double rem = currentRemainingOn(servers[i], now);

                if (arriving.remainingService < rem) {
                    if (rem > worstRem) {
                        worstRem = rem;
                        victim = i;
                    }
                }
            }
            return victim;
        }
    }

    private static void preemptServerToQueue(
            int serverId,
            double now,
            PriorityQueue<Customer> sjfQ,
            Server[] servers,
            MutableCounters ctr
    ) {
        Customer running = servers[serverId].current;
        if (running == null) return;

        boolean MMC_MODE = (!DETERMINISTIC_SERVICE && !LENGTH_BASED_SERVICE);

        if (MMC_MODE) {
            // In true M/M/c: discard remaining service (memoryless restart)
            running.remainingService = Double.NaN;
        } else {
            // Non-MMC: keep remaining work
            running.remainingService = currentRemainingOn(servers[serverId], now);
        }

        running.lastQueueEnterTime = now;
        sjfQ.add(running);
        ctr.Nq++;

        servers[serverId].busy = false;
        servers[serverId].current = null;
    }

    // Counters.

    private static class MutableCounters {
        int Ns, Nq;
        int expiredBeforeService, dropped;
        int completed;
        double totalSystemDelay, totalQueueDelay;
    }

    // Start service on server.

    private static boolean startServiceOnServer(
            int serverId,
            Customer job,
            double now,
            double mu,
            PriorityQueue<Event> fel,
            Server[] servers,
            MutableCounters ctr
    ) {
        // deadline check at service start (no expiry upon arrival)
        if (now >= job.deadlineTime) {
            ctr.expiredBeforeService++;
            if (DROP_EXPIRED) {
                ctr.dropped++;
                ctr.Ns--;
                return false;
            }
        }

        if (!Double.isNaN(job.lastQueueEnterTime)) {
            ctr.totalQueueDelay += now - job.lastQueueEnterTime;
            job.lastQueueEnterTime = Double.NaN;
        }

        if (Double.isNaN(job.remainingService)) {
            // For MMC: Exp(mu). For non-MMC: L/R or deterministic/Exp depending on flags.
            job.remainingService = serviceTimeFor(job, mu);
        }

        servers[serverId].busy = true;
        servers[serverId].current = job;
        servers[serverId].serviceStartTime = now;
        servers[serverId].version++;

        fel.add(new ServiceEvent(
                now + job.remainingService,
                serverId,
                servers[serverId].version
        ));

        return true;
    }

    /*******************************************************************************************/
    // SIMULATION CORE      

    private static Summary simulate(double lambda, double mu, int c) {

        PriorityQueue<Event> fel = new PriorityQueue<>();
        PriorityQueue<Customer> sjfQ = new PriorityQueue<>(SJF_COMPARATOR);

        Server[] servers = new Server[c];
        for (int i = 0; i < c; i++) servers[i] = new Server();

        double t = 0.0;
        int arrivals = 0;

        MutableCounters ctr = new MutableCounters();

        double areaNs = 0.0, areaNq = 0.0;
        double[] serverBusyArea = new double[c];

        // schedule first arrival
        fel.add(new Event(0.0, EventType.ARRIVAL, -1));

        while (arrivals < RUN_COMPLETED || ctr.Ns > 0) {

            Event e = fel.poll();
            double oldT = t;
            t = e.time;

            // time integrals update
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

                // schedule next arrival (packet RNG only).
                if (arrivals < RUN_COMPLETED)
                    fel.add(new Event(t + sampleInterarrival(lambda),
                            EventType.ARRIVAL, -1));

                // new packet properties (same in EDF, SJF).
                boolean isTight = RNG_PKT.nextDouble() < PROB_TIGHT;
                double rate = isTight ? DEADLINE_RATE_TIGHT : DEADLINE_RATE_LOOSE;
                double deadlineTime = t + expSample(rate, RNG_PKT);

                int lengthBits = isTight
                        ? geometricSample(GEOM_P_TIGHT, RNG_PKT)
                        : geometricSample(GEOM_P_LOOSE, RNG_PKT);

                Customer cust = new Customer(t, deadlineTime, isTight, lengthBits);

                // initialize remainingService for non-MMC mode.
                boolean MMC_MODE = (!DETERMINISTIC_SERVICE && !LENGTH_BASED_SERVICE);
                if (!MMC_MODE) {
                    cust.remainingService = serviceTimeFor(cust, mu);
                }

                // admit to system.
                ctr.Ns++;

                // check for idle server, if any.
                int idle = findIdleServer(servers);
                if (idle >= 0) { // FOUND = start service immediately
                    startServiceOnServer(idle, cust, t, mu, fel, servers, ctr);
                } else { // join FIFO queue / SJF queue
                    if (PREEMPTIVE_SJF) {
                        int victim = choosePreemptionVictim(servers, cust, t);
                        if (victim >= 0) {
                            preemptServerToQueue(victim, t, sjfQ, servers, ctr);
                            startServiceOnServer(victim, cust, t, mu, fel, servers, ctr);
                        } else {
                            cust.lastQueueEnterTime = t;
                            sjfQ.add(cust);
                            ctr.Nq++;
                        }
                    } else {
                        cust.lastQueueEnterTime = t;
                        sjfQ.add(cust);
                        ctr.Nq++;
                    }
                }

            } else {

                // DEPARTURE EVENT
                ServiceEvent se = (ServiceEvent) e;
                if (se.version != servers[se.serverId].version) continue;

                Customer done = servers[se.serverId].current;
                servers[se.serverId].busy = false;
                servers[se.serverId].current = null;

                ctr.completed++;
                ctr.Ns--;
                ctr.totalSystemDelay += t - done.arrivalTime;

                // start next from the queue, if any.
                while (!sjfQ.isEmpty()) {
                    Customer next = sjfQ.poll();
                    ctr.Nq--;
                    if (startServiceOnServer(se.serverId, next, t, mu, fel, servers, ctr))
                        break;
                }
            }
        }

        // finalize
        Summary S = new Summary();
        S.avgNs = areaNs / t;
        S.avgNq = areaNq / t;
        S.avgDs = ctr.totalSystemDelay / arrivals;
        S.avgDq = ctr.totalQueueDelay / arrivals;
        S.pExpired = ctr.expiredBeforeService / (double) arrivals;
        S.pServed = ctr.completed / (double) arrivals;

        System.out.println(ctr.completed + " " + arrivals);

        return S;
    }

    // Summary statistics class.

    private static class Summary {
        double avgNs, avgNq, avgDs, avgDq;
        double utilization;
        double pExpired, pServed;
    }

    /*******************************************************************************************/
    // DRIVER METHOD.

    public static void main(String[] args) throws Exception {

        try (PrintWriter pw = new PrintWriter(new FileWriter(CSV_FILE))) {
            pw.println("rho,lambda,mu,c,Ns_sim,Nq_sim,Ds_sim,Dq_sim,pExpired_sim, pServed");
        }

        for (double rho = RHO_START; rho <= RHO_END + 1e-12; rho += RHO_STEP) {

            double lambda = rho * C * MU;
            Summary S = simulate(lambda, MU, C);

            try (PrintWriter out = new PrintWriter(new FileWriter(CSV_FILE, true))) {
                out.printf(Locale.US,
                        "%.4f,%.6f,%.6f,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f%n",
                        rho, lambda, MU, C,
                        S.avgNs, S.avgNq,
                        S.avgDs, S.avgDq,
                        S.pExpired, S.pServed);
            }
        }
        System.out.println("Simulation completed. Results saved to " + CSV_FILE);

    }
}
