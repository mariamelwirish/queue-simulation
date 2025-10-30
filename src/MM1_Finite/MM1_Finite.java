package MM1_Finite;

import java.util.*;
import java.io.*;
import java.nio.file.*;

public class MM1_Finite {
    // Parameters
    private final double lambda;         // arrival rate
    private final double mu;             // service rate
    private final int    K;              // buffer capacity (max waiting)
    private final long   MAX_COMPLETED;  // stop after this many departures

    // Simulation state
    private double MC = 0.0;             // master clock (time)
    private double nextArrival;          // next external arrival time
    private double nextDeparture = Double.POSITIVE_INFINITY;
    private int n = 0;                   // # in system (in service + waiting)
    private boolean busy = false;        // server busy?

    // Track the customer currently in service
    private Double inServiceArrive = null;
    private Double inServiceStart = null;

    // FIFO queue of arrival timestamps for waiting customers
    private final Deque<Double> Q = new ArrayDeque<>();

    // RNG (Exponential sampling)
    private final Random rng;
    private double exp(double rate) {
        return -Math.log(1.0 - rng.nextDouble()) / rate;
    }

    // Stats
    private double areaNs = 0.0;         // area under curve for Ns
    private double areaNq = 0.0;         // area under curve for Nq
    private double lastT = 0.0;          // last time we added areas

    private double sumDs = 0.0;          // sum of system delays
    private double sumDq = 0.0;          // sum of queueing delays
    private long completed = 0;          // # departures so far
    private long blocked = 0;            // # blocked (dropped) arrivals

    // For Pn estimation via time-in-state
    private final ArrayList<Double> timeInState = new ArrayList<>();

    // Logging controls
    private static final boolean LOG_CONSOLE = false;
    private long iterations = 0;

    // File writers
    private BufferedWriter tsWriter = null;

    // Constructor
    public MM1_Finite(double lambda, double mu, int K, long maxCompleted) {
        this.lambda = lambda;
        this.mu = mu;
        this.K = K;
        this.MAX_COMPLETED = maxCompleted;
        this.rng = new Random();
    }

    // Ensure timeInState supports index n
    private void ensureStateSize(int idx) {
        while (timeInState.size() <= idx) timeInState.add(0.0);
    }

    // Accumulate areas/time in state for the interval [from, to)
    private void accrue(double from, double to) {
        if (to <= from) return;
        double dt = to - from;

        // Ns is current n
        areaNs += dt * n;

        // Nq is number waiting (not in service)
        int nq = Math.max(n - (busy ? 1 : 0), 0);
        areaNq += dt * nq;

        ensureStateSize(n);
        timeInState.set(n, timeInState.get(n) + dt);
    }

    // Event handlers
    private void onArrival() {
        double tA = MC;

        // Check if system is full (max capacity = K+1: K waiting + 1 in service)
        if (n >= K + 1) {
            // BLOCK/DROP the packet
            blocked++;
            // Schedule next arrival anyway
            nextArrival = MC + exp(lambda);
            return;
        }

        // Accept the arrival
        n += 1;

        if (!busy) {
            // start service immediately
            busy = true;
            inServiceArrive = tA;
            inServiceStart = MC;
            nextDeparture = MC + exp(mu);
        } else {
            // wait in FIFO queue
            Q.addLast(tA);
        }

        // schedule the next Poisson arrival
        nextArrival = MC + exp(lambda);
    }

    private void onDeparture() {
        // Customer in service departs
        double tArr = inServiceArrive;
        double tSrv = inServiceStart;
        double Ds = MC - tArr;     // system delay
        double Dq = tSrv - tArr;   // queueing delay

        sumDs += Ds;
        sumDq += Dq;
        completed += 1;


        n -= 1;

        // Start next, if any
        if (!Q.isEmpty()) {
            inServiceArrive = Q.removeFirst();
            inServiceStart = MC;
            nextDeparture = MC + exp(mu);
            busy = true;
        } else {
            inServiceArrive = null;
            inServiceStart = null;
            nextDeparture = Double.POSITIVE_INFINITY;
            busy = false;
        }
    }

    // Run loop
    public void run() {
        // Prepare results directory and CSV writers
        try {
            Path resultsDir = Paths.get("src/MM1_Finite/results");
            Files.createDirectories(resultsDir);

            Path tsPath = resultsDir.resolve("mm1_finite_timeseries.csv");
            tsWriter = Files.newBufferedWriter(tsPath);
            tsWriter.write("event_idx,MC,completed,blocked,n,q_len,busy,Ns_running,Nq_running,Ds_running,Dq_running,Pblock_running\n");

        } catch (IOException e) {
            throw new UncheckedIOException("Failed to open results CSV files", e);
        }

        try {
            // Initialize first arrival
            nextArrival = MC + exp(lambda);
            lastT = MC;

            if (LOG_CONSOLE) {
                System.out.printf("%-8s %-12s %-12s %-11s %-6s %-8s %-6s %-6s%n",
                        "Iter", "Δt", "MC", "Event", "n", "blocked", "completed", "busy");
            }

            // Main loop
            while (iterations < MAX_COMPLETED) {
                double tNext = Math.min(nextArrival, nextDeparture);
                boolean isDeparture = (nextDeparture <= nextArrival);

                accrue(lastT, tNext);
                MC = tNext;
                lastT = MC;

                if (isDeparture) onDeparture();
                else onArrival();

                iterations++;

                // Running averages
                double Ns_running = (MC > 0) ? (areaNs / MC) : 0.0;
                double Nq_running = (MC > 0) ? (areaNq / MC) : 0.0;
                double Ds_running = (completed > 0) ? (sumDs / completed) : 0.0;
                double Dq_running = (completed > 0) ? (sumDq / completed) : 0.0;

                // Blocking probability: blocked / total arrivals attempted
                long totalArrivals = completed + blocked;
                double Pblock_running = (totalArrivals > 0) ? ((double)blocked / totalArrivals) : 0.0;

                int qLen = Math.max(n - (busy ? 1 : 0), 0);

                tsWriter.write(
                        iterations + "," + MC + "," + completed + "," + blocked + "," + n + "," + qLen + "," + (busy ? 1 : 0) + ","
                                + Ns_running + "," + Nq_running + "," + Ds_running + "," + Dq_running + "," + Pblock_running + "\n"
                );

                if (LOG_CONSOLE) {
                    System.out.printf("%-8d %-12.6f %-12.6f %-11s %-6d %-8d %-6d %-6b%n",
                            iterations, tNext - lastT, MC, (isDeparture ? "DEPARTURE" : "ARRIVAL"),
                            n, blocked, completed, busy);
                }
            }
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        } finally {
            try {
                if (tsWriter != null) tsWriter.close();
            } catch (IOException ignore) {}
        }

        report();
    }

    // Reporting
    private void report() {
        double measuredTime = MC;
        double Ns = areaNs / measuredTime;
        double Nq = areaNq / measuredTime;

        double Ds = sumDs / Math.max(1, completed);
        double Dq = sumDq / Math.max(1, completed);

        double rho = lambda / mu;

        // Observed throughput (only completed packets)
        double lambdaHat = completed / measuredTime;

        // Blocking probability
        long totalArrivals = completed + blocked;
        double Pblock = (totalArrivals > 0) ? ((double)blocked / totalArrivals) : 0.0;

        // Theory (M/M/1/K) - only valid for rho != 1
        double P0_th, Ns_th, Nq_th, Pblock_th;

        if (Math.abs(rho - 1.0) == 0) {
            // Special case: rho = 1
            P0_th = 1.0 / (K + 2);
            Pblock_th = P0_th;
            Ns_th = (K + 1) / 2.0;
            Nq_th = Ns_th - (1 - P0_th);
        } else {
            // General case: rho != 1
            P0_th = (1 - rho) / (1 - Math.pow(rho, K + 2));
            Pblock_th = P0_th * Math.pow(rho, K + 1);

            // Ns calculation (assignment)
            double numerator = (K + 1) * Math.pow(rho, K + 2) - (K + 2) * Math.pow(rho, K + 1) + 1;
            double denominator = Math.pow(1 - rho, 2);
            Ns_th = (rho * (1 - rho) / (1 - Math.pow(rho, K + 2))) * (numerator / denominator);

            // Nq = Ns - (1 - P0) = Ns - probability server is busy
            Nq_th = Ns_th - (1 - P0_th);
        }

        // Effective arrival rate (accounting for blocking)
        double lambdaEff_th = lambda * (1 - Pblock_th);
        double Ds_th = (lambdaEff_th > 0) ? (Ns_th / lambdaEff_th) : 0;
        double Dq_th = (lambdaEff_th > 0) ? (Nq_th / lambdaEff_th) : 0;

        System.out.println("==== M/M/1/K Simulation Report (Finite Buffer) ====");
        System.out.printf("lambda=%.6f, mu=%.6f, K=%d, rho=%.6f%n", lambda, mu, K, rho);
        System.out.printf("completed=%d, blocked=%d, total_arrivals=%d%n", completed, blocked, totalArrivals);
        System.out.printf("measuredTime=%.6f, lambda_hat=%.6f%n", measuredTime, lambdaHat);
        System.out.println();

        // State probabilities
        System.out.println("State probabilities (Pn) from time-in-state:");
//        for (int i = 0; i < Math.min(timeInState.size(), K + 5); i++) {
//            double pn = (measuredTime > 0) ? (timeInState.get(i) / measuredTime) : 0.0;
//            if (pn > 1e-9 || i <= Math.min(10, K + 1)) {
//                System.out.printf("  P[%d] = %.9f%n", i, pn);
//            }
//        }
        System.out.println();

        // Simulation results
        System.out.printf("Ns (avg in system)       = %.6f%n", Ns);
        System.out.printf("Nq (avg in queue)        = %.6f%n", Nq);
        System.out.printf("Ds (avg system delay)    = %.6f%n", Ds);
        System.out.printf("Dq (avg queue delay)     = %.6f%n", Dq);
        System.out.printf("Pblock (blocking prob)   = %.6f (%.2f%%)%n", Pblock, Pblock * 100);
        System.out.println();

        // Little's Theorem checks (using observed lambda_hat)
        System.out.println("Little's Theorem checks (using lambda_hat):");
        System.out.printf("  Ns = λ_hat * Ds: %.6f = %.6f%n", Ns, lambdaHat * Ds);
        System.out.printf("  Nq = λ_hat * Dq: %.6f = %.6f%n", Nq, lambdaHat * Dq);
        System.out.println();

        // Theory vs Simulation
        System.out.println("M/M/1/K Theory vs Simulation:");
        System.out.printf("  Ns_theory = %.6f,  Ns_sim = %.6f,  Δ=%.6f%n", Ns_th, Ns, Math.abs(Ns_th - Ns));
        System.out.printf("  Nq_theory = %.6f,  Nq_sim = %.6f,  Δ=%.6f%n", Nq_th, Nq, Math.abs(Nq_th - Nq));
        System.out.printf("  Ds_theory = %.6f,  Ds_sim = %.6f,  Δ=%.6f%n", Ds_th, Ds, Math.abs(Ds_th - Ds));
        System.out.printf("  Dq_theory = %.6f,  Dq_sim = %.6f,  Δ=%.6f%n", Dq_th, Dq, Math.abs(Dq_th - Dq));
        System.out.printf("  Pblock_theory = %.6f,  Pblock_sim = %.6f,  Δ=%.6f%n", Pblock_th, Pblock, Math.abs(Pblock_th - Pblock));

    }

    // Driver
    public static void main(String[] args) {
        // Parameters
        double lambda = 0.18;     // arrival rate
        double mu     = 0.25;     // service rate
        int    K      = 100;      // buffer capacity (max waiting)
        double rho    = lambda / mu;

        System.out.printf("Running M/M/1/K with lambda=%.2f, mu=%.2f, K=%d, rho=%.2f%n",
                lambda, mu, K, rho);

        long maxCompleted = 900000;

        new MM1_Finite(lambda, mu, K, maxCompleted).run();
    }
}