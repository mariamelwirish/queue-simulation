package MM1_w_Deadline;

import java.util.*;
import java.io.*;
import java.nio.file.*;

public class MM1wDeadline {
    // Parameters
    private final double lambda;         // arrival rate
    private final double mu;             // service rate (server capacity)
    private final double deadlineRate;   // [DEADLINE] deadline ~ Exp(deadlineRate)
    private final long   MAX_COMPLETED;  // stop after this many departures

    // Simulation state
    private double MC = 0.0;             // master clock (time)
    private double nextArrival;          // next external arrival time
    private double nextDeparture = Double.POSITIVE_INFINITY; // next service completion time
    private int n = 0;                   // # in system (in service + waiting)
    private boolean busy = false;        // server busy?

    // Track the customer currently in service
    private static final class Job {     // [DEADLINE] keep arrival + deadline together
        final double arrival;
        final double deadlineAbs;        // absolute deadline time = arrival + Exp(delta)
        Job(double arrival, double deadlineAbs) {
            this.arrival = arrival;
            this.deadlineAbs = deadlineAbs;
        }
    }
    private Job inService = null;        // [DEADLINE] job currently in service
    private Double inServiceStart = null;
    private boolean inServiceWasExpiredAtStart = false; // [DEADLINE] flag

    // FIFO queue of waiting jobs
    private final Deque<Job> Q = new ArrayDeque<>();

    // RNG (Exponential sampling)
    private final Random rng;
    private double exp(double rate) { return -Math.log(1.0 - rng.nextDouble()) / rate; }

    // Stats (areas)
    private double areaNs = 0.0;         // area under curve for Ns
    private double areaNq = 0.0;         // area under curve for Nq
    private double lastT = 0.0;          // last time we added areas

    // Delay sums
    private double sumDs = 0.0;          // sum of system delays
    private double sumDq = 0.0;          // sum of queueing delays
    private long completed = 0;          // # departures so far

    // [DEADLINE] expiry stats
    private long expiredServed = 0;      // served jobs that had expired before service start

    // For Pn estimation
    private final ArrayList<Double> timeInState = new ArrayList<>();

    private static final boolean LOG_CONSOLE = true;
    private long iterations = 0;

    // File writers
    private BufferedWriter tsWriter = null;   // time series per master-clock event

    // Constructor
    public MM1wDeadline(double lambda, double mu, double deadlineRate, long maxCompleted) {
        this.lambda = lambda;
        this.mu = mu;
        this.deadlineRate = deadlineRate; // [DEADLINE]
        this.MAX_COMPLETED = maxCompleted;
        this.rng = new Random();
    }

    // Ensure timeInState supports index n
    private void ensureStateSize(int idx) { while (timeInState.size() <= idx) timeInState.add(0.0); }

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
        double deadlineAbs = tA + exp(deadlineRate);
        Job j = new Job(tA, deadlineAbs);

        if (!busy) {
            // start service immediately
            busy = true;
            inService = j;
            inServiceStart = MC;
            inServiceWasExpiredAtStart = (inServiceStart >= inService.deadlineAbs);
            nextDeparture = MC + exp(mu);
        } else {
            // wait in FIFO queue
            Q.addLast(j);
        }
        // schedule the next Poisson arrival
        nextArrival = MC + exp(lambda);
    }

    private void onDeparture() {
        // Customer in service departs
        double tArr = inService.arrival;
        double tSrv = inServiceStart;
        double Ds = MC - tArr;     // system delay
        double Dq = tSrv - tArr;   // queueing delay

        sumDs += Ds;
        sumDq += Dq;
        completed += 1;

        // count expiry-at-start if true
        if (inServiceWasExpiredAtStart) expiredServed++;

        n -= 1;

        // Start next, if any
        if (!Q.isEmpty()) {
            inService = Q.removeFirst();
            inServiceStart = MC;
            inServiceWasExpiredAtStart = (inServiceStart >= inService.deadlineAbs);
            nextDeparture = MC + exp(mu);
            busy = true;
        } else {
            inService = null;
            inServiceStart = null;
            inServiceWasExpiredAtStart = false;
            nextDeparture = Double.POSITIVE_INFINITY;
            busy = false;
        }
    }

    // Run loop
    public void run() {
        try {
            Path resultsDir = Paths.get("src/MM1_w_Deadline/results");
            Files.createDirectories(resultsDir);

            Path tsPath  = resultsDir.resolve("mm1_timeseries.csv");

            tsWriter  = Files.newBufferedWriter(tsPath);

            tsWriter.write(
                    "event_idx,MC,completed,n,q_len,busy," +
                            "Ns_running,Nq_running,Ds_running,Dq_running," +
                            "expiredServed,pe_running\n" // [DEADLINE]
            );
            tsWriter.write("\n");

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
                double tNext = Math.min(nextArrival, nextDeparture);
                boolean isDeparture = (nextDeparture <= nextArrival);

                double dt = tNext - lastT;
                accrue(lastT, tNext);
                MC = tNext;
                lastT = MC;

                if (isDeparture) onDeparture(); else onArrival();

                iterations++;

                double Ns_running = (MC > 0) ? (areaNs / MC) : 0.0;
                double Nq_running = (MC > 0) ? (areaNq / MC) : 0.0;
                double Ds_running = (completed > 0) ? (sumDs / completed) : 0.0;
                double Dq_running = (completed > 0) ? (sumDq / completed) : 0.0;
                int qLen = Math.max(n - (busy ? 1 : 0), 0);
                double pe_running = (completed > 0) ? ((double) expiredServed / completed) : 0.0; // [DEADLINE]

                tsWriter.write(
                        iterations + "," + MC + "," + completed + "," + n + "," + qLen + "," + (busy ? 1 : 0) + "," +
                                Ns_running + "," + Nq_running + "," + Ds_running + "," + Dq_running + "," +
                                expiredServed + "," + pe_running + "\n"
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
            try { if (tsWriter  != null) tsWriter.close(); } catch (IOException ignore) {}
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
        double DqTh = NqTh / lambda;

        double pe = (completed > 0) ? ((double) expiredServed / completed) : 0.0; // [DEADLINE]

        System.out.println("==== M/M/1 Simulation Report (Fixed Completions, with Deadlines) ====");
        System.out.printf("lambda=%.6f, mu=%.6f, rho=%.6f, deadlineRate=%.6f%n", lambda, mu, rho, deadlineRate);
        System.out.printf("completed=%d, measuredTime=%.6f, lambda_hat=%.6f%n", completed, measuredTime, lambdaHat);
        System.out.println();

        System.out.println("State probabilities (Pn) from time-in-state:");
        for (int i = 0; i < timeInState.size() && i <= 20; i++) {
            double pn = (measuredTime > 0) ? (timeInState.get(i) / measuredTime) : 0.0;
            if (pn > 1e-9 || i <= 10) System.out.printf("  P[%d] = %.9f%n", i, pn);
        }
        System.out.println();

        System.out.printf("Ns (avg in system)     = %.6f%n", Ns);
        System.out.printf("Nq (avg in queue)      = %.6f%n", Nq);
        System.out.printf("Ds (avg system delay)  = %.6f%n", Ds);
        System.out.printf("Dq (avg queue delay)   = %.6f%n", Dq);
        System.out.printf("pe (served after expiry probability) = %.6f  (=%d/%d)%n",
                pe, expiredServed, completed); // [DEADLINE]
        System.out.println();

        System.out.println("Little's Theorem checks (using lambda_hat):");
        System.out.printf("  Ns = λ_hat * Ds: %.6f ≈ %.6f%n", Ns, lambdaHat * Ds);
        System.out.printf("  Nq = λ_hat * Dq: %.6f ≈ %.6f%n", Nq, lambdaHat * Dq);
        System.out.println();

        System.out.println("M/M/1 Theory vs Simulation (deadline-independent quantities):");
        System.out.printf("  Ns_theory = %.6f,  Ns_sim = %.6f,  Δ=%.6f%n", NsTh, Ns, Math.abs(NsTh - Ns));
        System.out.printf("  Nq_theory = %.6f,  Nq_sim = %.6f,  Δ=%.6f%n", NqTh, Nq, Math.abs(NqTh - Nq));
        System.out.printf("  Ds_theory = %.6f,  Ds_sim = %.6f,  Δ=%.6f%n", DsTh, Ds, Math.abs(DsTh - Ds));
        System.out.printf("  Dq_theory = %.6f,  Dq_sim = %.6f,  Δ=%.6f%n", DqTh, Dq, Math.abs(DqTh - Dq));

        System.out.println("\nCSV written to:");
        System.out.println("  results/mm1_timeseries.csv  (running averages incl. pe_running)");
        System.out.println("  results/mm1_departures.csv  (per-departure with deadlines/flags)");

        System.out.println("\nMATLAB examples:");
        System.out.println("  T = readtable('src/MM1_Infinite/results/mm1_timeseries.csv');");
        System.out.println("  D = readtable('src/MM1_Infinite/results/mm1_departures.csv');");
        System.out.println("  figure; plot(T.MC, T.Ns_running); xlabel('Time'); ylabel('N_s running avg'); grid on;");
        System.out.println("  figure; plot(T.MC, T.Nq_running); xlabel('Time'); ylabel('N_q running avg'); grid on;");
        System.out.println("  figure; plot(T.completed, T.Ds_running); xlabel('Departures'); ylabel('D_s running avg'); grid on;");
        System.out.println("  figure; plot(T.completed, T.Dq_running); xlabel('Departures'); ylabel('D_q running avg'); grid on;");
        System.out.println("  % [DEADLINE] new plots:");
        System.out.println("  figure; plot(T.completed, T.pe_running); xlabel('Departures'); ylabel('P(expired-served) running'); grid on;");
        System.out.println("  figure; stem(D.dep_idx, D.expired_at_start); xlabel('Departure index'); ylabel('Expired-at-start (0/1)'); grid on;");
        System.out.println("  % Scatter of (queueing delay) vs expiry flag:");
        System.out.println("  figure; gscatter(D.Dq, D.expired_at_start); xlabel('D_q'); ylabel('expired-at-start'); grid on;");
    }

    // Driver
    public static void main(String[] args) {
        double lambda = 0.18;     // arrival rate
        double mu     = 0.25;     // service rate
        double delta  = 0.12;     // deadline rate

        double rho    = lambda / mu;
        if (rho >= 1.0) {
            System.err.println("Choose lambda < mu for stability (rho<1).");
            System.exit(1);
        }

        long maxCompleted = 100000;

        new MM1wDeadline(lambda, mu, delta, maxCompleted).run();
    }
}
