package MM1_Infinite;

import java.util.*;
import java.io.*;
import java.nio.file.*;

public class MM1_Infinite {
    // Parameters
    private final double lambda;         // arrival rate
    private final double mu;             // service rate (server capacity)
    private final long   MAX_COMPLETED;  // stop after this many departures

    // Simulation state
    private double MC = 0.0;             // master clock (time)
    private double nextArrival;          // next external arrival time
    private double nextDeparture = Double.POSITIVE_INFINITY; // next service completion time
    private int n = 0;                   // # in system (in service + waiting)
    private boolean busy = false;        // server busy?

    // Track the customer currently in service
    private Double inServiceArrive = null;   // its arrival time
    private Double inServiceStart = null;    // when service started

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

    // For Pn estimation
    private final ArrayList<Double> timeInState = new ArrayList<>();

    private static final boolean LOG_CONSOLE = false; // console spam off by default
    private long iterations = 0;

    // File writers
    private BufferedWriter tsWriter = null;   // time series per master-clock event

    // Constructor
    public MM1_Infinite(double lambda, double mu, long maxCompleted) {
        this.lambda = lambda;
        this.mu = mu;
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

        areaNs += dt * n;
        int nq = Math.max(n - (busy ? 1 : 0), 0);
        areaNq += dt * nq;

        ensureStateSize(n);
        timeInState.set(n, timeInState.get(n) + dt);
    }

    // Event handlers
    private void onArrival() {
        n += 1;
        double tA = MC;

        if (!busy) {
            // start service immediately
            busy = true;
            inServiceArrive = tA;
            inServiceStart = MC;
            nextDeparture = MC + exp(mu);
        } else {
            // wait in FIFO queue (store arrival time for delay accounting)
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

        // Count every completed customer
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
        try {
            Path resultsDir = Paths.get("src/MM1_Infinite/results");
            Files.createDirectories(resultsDir);

            Path tsPath = resultsDir.resolve("mm1_timeseries.csv");
            tsWriter = Files.newBufferedWriter(tsPath);
            tsWriter.write("event_idx,MC,completed,n,q_len,busy,Ns_running,Nq_running,Ds_running,Dq_running\n");

        } catch (IOException e) {
            throw new UncheckedIOException("Failed to open results/* CSV files", e);
        }

        try {
            // Initialize first arrival
            nextArrival = MC + exp(lambda);
            lastT = MC;

            if (LOG_CONSOLE) {
                System.out.printf("%-8s %-12s %-12s %-11s %-6s %-6s %-6s %-12s %-12s %-11s%n",
                        "Iter", "Δt", "MC", "Event", "n", "|Q|", "busy", "nextArrival", "nextDeparture", "completed");
            }

            // Main loop
            while (completed < MAX_COMPLETED) {
                // Pick next event: min(nextArrival, nextDeparture)
                double tNext = Math.min(nextArrival, nextDeparture);
                boolean isDeparture = (nextDeparture <= nextArrival);

                // Accrue areas/time-in-state up to tNext
                double dt = tNext - lastT;
                accrue(lastT, tNext);
                MC = tNext;
                lastT = MC;

                // Dispatch event (break ties consistently: DEPARTURE first)
                if (isDeparture) onDeparture();
                else onArrival();

                // Master clock event index
                iterations++;

                // Running averages after handling the event
                double Ns_running = (MC > 0) ? (areaNs / MC) : 0.0;
                double Nq_running = (MC > 0) ? (areaNq / MC) : 0.0;
                double Ds_running = (completed > 0) ? (sumDs / completed) : 0.0;
                double Dq_running = (completed > 0) ? (sumDq / completed) : 0.0;
                int qLen = Math.max(n - (busy ? 1 : 0), 0);

                tsWriter.write(
                        iterations + "," + MC + "," + completed + "," + n + "," + qLen + "," + (busy ? 1 : 0) + ","
                                + Ns_running + "," + Nq_running + "," + Ds_running + "," + Dq_running + "\n"
                );

                if (LOG_CONSOLE) {
                    System.out.printf("%-8d %-12.6f %-12.6f %-11s %-6d %-6d %-6b %-12.6f %-12.6f %-11d%n",
                            iterations, dt, MC, (isDeparture ? "DEPARTURE" : "ARRIVAL"),
                            n, Q.size(), busy, nextArrival, nextDeparture, completed);
                }
            }
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        } finally {
            // Close writers
            try {
                if (tsWriter != null) tsWriter.close();
            } catch (IOException ignore) {}
        }

        report();
    }

    // Reporting
    private void report() {
        double measuredTime = MC;               // total simulated time
        double Ns = areaNs / measuredTime;
        double Nq = areaNq / measuredTime;

        double Ds = sumDs / Math.max(1, completed);
        double Dq = sumDq / Math.max(1, completed);

        double lambdaHat = completed / measuredTime; // observed throughput
        double rho = lambda / mu;

        // Theory (M/M/1) for rho < 1
        double NsTh = rho / (1.0 - rho);
        double NqTh = (rho * rho) / (1.0 - rho);
        double DsTh = 1.0 / (mu - lambda);
        double DqTh = NqTh / lambda; // = rho/(mu - lambda)/mu

        System.out.println("==== M/M/1 Simulation Report (Fixed Completions) ====");
        System.out.printf("lambda=%.6f, mu=%.6f, rho=%.6f %n", lambda, mu, rho);
        System.out.printf("completed=%d, measuredTime=%.6f, lambda_hat=%.6f%n", completed, measuredTime, lambdaHat);
        System.out.println();

        // Pn (state probabilities) — print a few non-zero states
        System.out.println("State probabilities (Pn) from time-in-state:");
        for (int i = 0; i < timeInState.size() && i <= 20; i++) {
            double pn = (measuredTime > 0) ? (timeInState.get(i) / measuredTime) : 0.0;
            if (pn > 1e-9 || i <= 10) {
                System.out.printf("  P[%d] = %.9f%n", i, pn);
            }
        }
        System.out.println();

        // Time averages & delays
        System.out.printf("Ns (avg in system)     = %.6f%n", Ns);
        System.out.printf("Nq (avg in queue)      = %.6f%n", Nq);
        System.out.printf("Ds (avg system delay)  = %.6f%n", Ds);
        System.out.printf("Dq (avg queue delay)   = %.6f%n", Dq);
        System.out.println();

        // Little's Theorem checks with observed throughput
        System.out.println("Little's Theorem checks (using lambda_hat):");
        System.out.printf("  Ns = λ_hat * Ds: %.6f ≈ %.6f%n", Ns, lambdaHat * Ds);
        System.out.printf("  Nq = λ_hat * Dq: %.6f ≈ %.6f%n", Nq, lambdaHat * Dq);
        System.out.println();

        // Theory vs Simulation
        System.out.println("M/M/1 Theory vs Simulation:");
        System.out.printf("  Ns_theory = %.6f,  Ns_sim = %.6f,  Δ=%.6f%n", NsTh, Ns, Math.abs(NsTh - Ns));
        System.out.printf("  Nq_theory = %.6f,  Nq_sim = %.6f,  Δ=%.6f%n", NqTh, Nq, Math.abs(NqTh - Nq));
        System.out.printf("  Ds_theory = %.6f,  Ds_sim = %.6f,  Δ=%.6f%n", DsTh, Ds, Math.abs(DsTh - Ds));
        System.out.printf("  Dq_theory = %.6f,  Dq_sim = %.6f,  Δ=%.6f%n", DqTh, Dq, Math.abs(DqTh - Dq));

    }

    // Driver
    public static void main(String[] args) {
        double lambda = 0.18;     // arrival rate
        double mu     = 0.25;     // service rate
        double rho    = lambda / mu;
        if (rho >= 1.0) {
            System.err.println("Choose lambda < mu for stability (rho<1).");
            System.exit(1);
        }

        long maxCompleted = 100000;

        new MM1_Infinite(lambda, mu, maxCompleted).run();
    }
}
