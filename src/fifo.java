import java.util.*;
import java.io.*;

public class fifo {

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
    static final double DEADLINE_RATE_TIGHT = 2.0;
    static final double DEADLINE_RATE_LOOSE = 0.5;
    
    // Packet size classes.
    static final double SMALL_PCKT_SIZE = 2.5; 
    static final double LARGE_PCKT_SIZE = 10; 
    static final double GEOM_P_TIGHT = 1 / SMALL_PCKT_SIZE; // for geometric distribution
    static final double GEOM_P_LOOSE = 1 / LARGE_PCKT_SIZE; // for geometric distribution

    // Flags to control simulation modes.
    static final boolean DETERMINISTIC_ARRIVAL = false; // Deterministic arrival process
    static final boolean DETERMINISTIC_SERVICE = false; // Deterministic service times
    static final boolean DROP_EXPIRED = true; // Drop expired customers at service start

    // CSV output file (mode-aware).
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

    /*******************************************************************************************/
    // SIMULATION STRUCTURE HELPERS

    // Random number generators
    static final Random RNG_PKT = new Random(12345); // for arrivals, deadlines, sizes
    static final Random RNG_SVC = new Random(54321); // for service times

    // Event types.
    private enum EventType {
        ARRIVAL,
        DEPARTURE
    }

    // Event class.
    private static class Event implements Comparable<Event> {
        double time;
        EventType type;
        int serverIndex;
        Customer customer;

        Event(double time, EventType type, int serverIndex, Customer customer) {
            this.time = time;
            this.type = type;
            this.serverIndex = serverIndex;
            this.customer = customer;
        }

        @Override
        public int compareTo(Event other) {
            return Double.compare(this.time, other.time);
        }
    }

    // Server class.
    private static class Server {
        boolean busy;
        Customer customer;

        Server() {
            busy = false;
            customer = null;
        }
    }

    // Customer (packet) class.
    static class Customer {
        double arrivalTime;
        double deadlineTime;
        boolean isTight;
        int lengthBits;

        Customer(double arrivalTime, double deadlineTime, boolean isTight, int lengthBits) {
            this.arrivalTime = arrivalTime;
            this.deadlineTime = deadlineTime;
            this.isTight = isTight;
            this.lengthBits = lengthBits;
        }
    }

    // Summary statistics class.
    private static class Summary {
        double avgNs, avgNq, avgDs, avgDq;
        double utilization;
        int completed;
        double pExpired_sim;
        double pServed;
    }

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


    /*******************************************************************************************/
    // SIMULATION CORE      

    // Find an idle server index, or -1 if none.
    private static int findIdleServer(Server[] servers) {
        for (int i = 0; i < servers.length; i++) {
            if (!servers[i].busy) return i;
        }
        return -1;
    }

    // Main simulation method.
    private static Summary simulate(double lambda, double mu, int c) {
        Summary S = new Summary();
        PriorityQueue<Event> fel = new PriorityQueue<>();
        Queue<Customer> queue = new ArrayDeque<>();
        Server[] servers = new Server[c];
        for (int i = 0; i < c; i++)
            servers[i] = new Server();

        // State variables
        double t = 0.0;
        int arrivals = 0;
        int Ns = 0;
        int Nq = 0;

        // Integrals
        double areaNs = 0.0;
        double areaNq = 0.0;

        // Trackers.
        int expiredBeforeService = 0;
        int completed = 0;

        // schedule first arrival
        fel.add(new Event(0.0, EventType.ARRIVAL, -1, null));

        while (arrivals < RUN_COMPLETED || Ns > 0) {

            Event e = fel.poll();
            double oldT = t;
            t = e.time;

            // time integrals update
            double dt = t - oldT;
            areaNs += Ns * dt;
            areaNq += Nq * dt;

            // ARRIVAL EVENT
            if (e.type == EventType.ARRIVAL) {
                if (arrivals >= RUN_COMPLETED)
                    continue;

                arrivals++;

                // new packet properties (same in EDF, SJF).
                boolean isTight = RNG_PKT.nextDouble() < PROB_TIGHT;
                double deadlineRate = isTight ? DEADLINE_RATE_TIGHT : DEADLINE_RATE_LOOSE;
                double absDeadline = t + expSample(deadlineRate, RNG_PKT);
                int lengthBits = isTight ? geometricSample(GEOM_P_TIGHT, RNG_PKT)
                        : geometricSample(GEOM_P_LOOSE, RNG_PKT);
                Customer cust = new Customer(t, absDeadline, isTight, lengthBits);

                // admit to system.
                Ns++;

                // check for idle server, if any.
                int idle = findIdleServer(servers);
                if (idle >= 0) { // FOUND = start service immediately
                    servers[idle].busy = true;
                    servers[idle].customer = cust;
                    double serviceTime = DETERMINISTIC_SERVICE ? 1.0 / mu : expSample(mu, RNG_SVC);
                    fel.add(new Event(t + serviceTime, EventType.DEPARTURE, idle, cust));
                } else { // join FIFO queue    
                    queue.add(cust);
                    Nq++;
                }

                // schedule next arrival (packet RNG only).
                double interArrival = DETERMINISTIC_ARRIVAL ? 1.0 / lambda : expSample(lambda, RNG_PKT);
                fel.add(new Event(t + interArrival, EventType.ARRIVAL, -1, null));
            }

            // DEPARTURE EVENT
            else {

                // job completes
                int s = e.serverIndex;
                servers[s].busy = false;
                servers[s].customer = null;
                Ns--;
                completed++;

                // start next from the queue, if any.
                while (!queue.isEmpty() && !servers[s].busy) {
                    Customer next = queue.poll();
                    Nq--;

                    // deadline check at service start (no expiry upon arrival)
                    if (t > next.deadlineTime) {
                        expiredBeforeService++;
                        if (DROP_EXPIRED) { // drop immediately (counts out of system)  
                            Ns--;
                        } else {
                            // still serve it
                            servers[s].busy = true;
                            servers[s].customer = next;

                            double serviceTime = DETERMINISTIC_SERVICE ? 1.0 / mu : expSample(mu, RNG_SVC);

                            fel.add(new Event(t + serviceTime, EventType.DEPARTURE, s, next));
                        }
                    } else {

                        servers[s].busy = true;
                        servers[s].customer = next;

                        double serviceTime = DETERMINISTIC_SERVICE ? 1.0 / mu : expSample(mu, RNG_SVC);

                        fel.add(new Event(t + serviceTime, EventType.DEPARTURE, s, next));
                    }
                }
            }
        }

        // finalize
        S.avgNs = areaNs / t;
        S.avgNq = areaNq / t;
        S.avgDs = S.avgNs / lambda;
        S.avgDq = S.avgNq / lambda;
        S.completed = completed;
        S.pExpired_sim = (double) expiredBeforeService / arrivals;
        S.pServed = (double) completed / arrivals;

        System.out.println(completed + " " + arrivals);
        return S;
    }

    /*******************************************************************************************/
    // DRIVER METHOD.
    public static void main(String[] args) throws Exception {
        PrintWriter out = new PrintWriter(new FileWriter(CSV_FILE));
        out.println("rho,lambda,mu,c,Ns_sim,Nq_sim,Ds_sim,Dq_sim,pExpired_sim,pServed");
        for (double rho = RHO_START; rho <= RHO_END + 1e-12; rho += RHO_STEP) {
            double lambda = rho * C * MU;
            Summary S = simulate(lambda, MU, C);
            out.printf(Locale.US,
                    "%.5f,%.5f,%.5f,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f%n",
                    rho, lambda, MU, C,
                    S.avgNs, S.avgNq, S.avgDs, S.avgDq,
                    S.pExpired_sim,
                    S.pServed);
        }
        out.close();
        System.out.println("Simulation completed. Results saved to " + CSV_FILE);

    }
}
