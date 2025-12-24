import java.util.*;
import java.io.*;

public class edf {

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

    // Flags to control simulation modes.

    static final boolean DETERMINISTIC_ARRIVAL = false; // Deterministic arrival process
    static final boolean DETERMINISTIC_SERVICE = false; // Deterministic service times
    static final boolean DROP_EXPIRED = true; // Drop expired customers at service start
    static final boolean PREEMPTIVE_EDF = false; // EDF preemptive mode

    // CSV output file (mode-aware).
    static final String CSV_FILE;
    static {
        // Build a mode-aware filename to avoid overwriting.
        if (PREEMPTIVE_EDF) {
            CSV_FILE = "../results/csv/preemptive_edf.csv";
        } else if (DROP_EXPIRED) {
            CSV_FILE = "../results/csv/drop_edf.csv";
        } else if (DETERMINISTIC_ARRIVAL || DETERMINISTIC_SERVICE) {
            CSV_FILE = "../results/csv/ddc_edf.csv";
        } else {
            CSV_FILE = "../results/csv/mmc_edf.csv";
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

    // Service event class (for preemption).
    // This class is used to check if an already scheduled departure event is still valid.
    // If the version of the server has changed, the event should be ignored.
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

    // EDF comparator.
    private static final Comparator<Customer> EDF_COMPARATOR =
        Comparator.<Customer>comparingDouble(c -> c.deadlineTime).thenComparingDouble(c -> c.arrivalTime);

    // exponential(rate) sample using inverse-CDF.
    private static double expSample(double rate, Random rng) {
        double u = rng.nextDouble();
        if (u == 0.0) u = 1e-16;
        return -Math.log(1.0 - u) / rate;
    }

    // Geometric(p) sample using inverse-CDF.
    private static int geometricSample(double p, Random rng) {
        double u = rng.nextDouble();
        if (u == 0.0)
            u = 1e-16;
        return (int) Math.ceil(Math.log(1.0 - u) / Math.log(1.0 - p));
    }

    // Helper to decide interarrival times.
    private static double sampleInterarrival(double lambda) {
        return DETERMINISTIC_ARRIVAL ? 1.0 / lambda : expSample(lambda, RNG_PKT);
    }
    // Helper to decide service times.
    private static double sampleService(double mu) {
        return DETERMINISTIC_SERVICE ? 1.0 / mu : expSample(mu, RNG_SVC);
    }


    /*******************************************************************************************/
    // SIMULATION CORE 

    // Find an idle server index, or -1 if none.
    private static int findIdleServer(Server[] servers) {
        for (int i = 0; i < servers.length; i++)
            if (!servers[i].busy) return i;
        return -1;
    }

    // Preemption helpers: choose victim server index, or -1 if none.
    private static int choosePreemptionVictim(Server[] servers, Customer arriving) {
        int victim = -1;
        double worstDeadline = Double.NEGATIVE_INFINITY;

        for (int i = 0; i < servers.length; i++) {
            // if (!servers[i].busy || servers[i].current == null) continue;
            Customer running = servers[i].current;
            if (arriving.deadlineTime < running.deadlineTime) {
                if (running.deadlineTime > worstDeadline) {
                    worstDeadline = running.deadlineTime;
                    victim = i;
                }
            }
        }
        return victim;
    }

    // Preempt running packet to queue and update its data.
    private static void preemptServerToQueue(int serverId, double now, PriorityQueue<Customer> edfQ, Server[] servers, MutableCounters ctr) {
        Customer running = servers[serverId].current;
        if (running == null) return;

        if (!DETERMINISTIC_SERVICE) {
            running.remainingService = Double.NaN;
        } else {
            double elapsed = now - servers[serverId].serviceStartTime;
            running.remainingService -= elapsed;
            if (running.remainingService < 0) running.remainingService = 0;
        }

        running.lastQueueEnterTime = now;
        edfQ.add(running);
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
            job.remainingService = sampleService(mu);
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


    private static Summary simulate(double lambda, double mu, int c) {

        PriorityQueue<Event> fel = new PriorityQueue<>();
        PriorityQueue<Customer> edfQ = new PriorityQueue<>(EDF_COMPARATOR);

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

                // new packet properties (same in EDF, SJF).
                boolean isTight = RNG_PKT.nextDouble() < PROB_TIGHT;
                double rate = isTight ? DEADLINE_RATE_TIGHT : DEADLINE_RATE_LOOSE;
                double deadlineTime = t + expSample(rate, RNG_PKT);
                int lengthBits = isTight
                        ? geometricSample(GEOM_P_TIGHT, RNG_PKT)
                        : geometricSample(GEOM_P_LOOSE, RNG_PKT);
                Customer cust = new Customer(t, deadlineTime, isTight, lengthBits);

                // admit to system.
                ctr.Ns++;

                // schedule next arrival (packet RNG only).
                if (arrivals < RUN_COMPLETED)
                    fel.add(new Event(t + sampleInterarrival(lambda),
                            EventType.ARRIVAL, -1));

                // check for idle server, if any.
                int idle = findIdleServer(servers);
                if (idle >= 0) { // FOUND = start service immediately
                    startServiceOnServer(idle, cust, t, mu, fel, servers, ctr);
                } else { // join FIFO queue / EDF queue
                    if (PREEMPTIVE_EDF) {
                        int victim = choosePreemptionVictim(servers, cust);
                        if (victim >= 0) {
                            preemptServerToQueue(victim, t, edfQ, servers, ctr);
                            startServiceOnServer(victim, cust, t, mu, fel, servers, ctr);
                        } else {
                            cust.lastQueueEnterTime = t;
                            edfQ.add(cust);
                            ctr.Nq++;
                        }
                    } else {
                        cust.lastQueueEnterTime = t;
                        edfQ.add(cust);
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
                while (!edfQ.isEmpty()) {
                    Customer next = edfQ.poll();
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

        double busy = 0;
        for (double b : serverBusyArea) busy += b;
        S.utilization = busy / (c * t);

        S.pExpired = ctr.expiredBeforeService / (double) arrivals;
        S.pServed = ctr.completed / (double) arrivals;
    
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
            pw.println("rho,lambda,mu,c,Ns_sim,Nq_sim,Ds_sim,Dq_sim,pExpired_sim,pServed");
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
