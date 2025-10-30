package MMC_Infinite;

import java.util.*;
import java.io.*;
import java.nio.file.*;

public class MMC_Infinite {

    // Parameters
    private final double lambda;        // arrival rate
    private final double mu;            // service rate per server
    private final int c;                // number of servers
    private final long MAX_COMPLETED;   // stop after this many departures
    private final int LOG_EVERY_K;      // log every K completions

    // State
    private double MC = 0.0;                // master clock
    private double nextArrival;             // next external arrival time
    private final Deque<Double> Q = new ArrayDeque<>(); // waiting queue
    private long completed = 0;             // departures so far
    private int n = 0;                      // # in system (waiting + in service)
    private int busy = 0;                   // busy servers

    // Area integrals for time-averages
    private double areaNs = 0.0;
    private double areaNq = 0.0;
    private double lastT = 0.0;

    // Delay sums
    private double sumDs = 0.0;
    private double sumDq = 0.0;

    private final Random rng;
    private Theory theory;                  // theory results (computed once)

    // Job in service
    private static class Job {
        final double arrival;
        final double start;
        final double departure;
        Job(double arrival, double start, double departure) {
            this.arrival = arrival;
            this.start = start;
            this.departure = departure;
        }
    }

    // Min-heap by departure time
    private final PriorityQueue<Job> inService = new PriorityQueue<>(Comparator.comparingDouble(j -> j.departure));

    // RNG exponential
    private double exp(double rate) {
        return -Math.log(Math.max(1e-12, rng.nextDouble())) / rate;
    }

    // -------------- Theory (Erlang-C) --------------
    private static class Theory {
        final double P0, Pwait, Nq, Ns, Dq, Ds, rho, a;
        Theory(double lambda, double mu, int c) {
            this.a = lambda / mu;
            this.rho = lambda / (c * mu);
            if (rho >= 1.0) throw new IllegalArgumentException("Unstable: rho >= 1 (system diverges).");

            double sum = 0.0;
            for (int n = 0; n <= c - 1; n++) sum += Math.pow(a, n) / factorial(n);
            double tail = Math.pow(a, c) / factorial(c) * (1.0 / (1.0 - rho));
            this.P0 = 1.0 / (sum + tail);

            this.Pwait = (Math.pow(a, c) / factorial(c)) * (1.0 / (1.0 - rho)) * P0;
            this.Nq = (Pwait * rho) / (1.0 - rho);
            this.Ns = Nq + a;
            this.Dq = Nq / lambda;
            this.Ds = Dq + 1.0 / mu;
        }
        private static double factorial(int n) {
            double f = 1.0;
            for (int i = 2; i <= n; i++) f *= i;
            return f;
        }
    }

    // -------------- Constructor --------------
    public MMC_Infinite(double lambda, double mu, int c, long MAX_COMPLETED, int LOG_EVERY_K) {
        this.lambda = lambda;
        this.mu = mu;
        this.c = c;
        this.MAX_COMPLETED = MAX_COMPLETED;
        this.LOG_EVERY_K = LOG_EVERY_K;
        this.rng = new Random();
        this.nextArrival = exp(lambda);
        this.theory = new Theory(lambda, mu, c);
    }

    // -------------- Main simulation loop --------------
    public void run(String csvPath) throws IOException {
        try (BufferedWriter out = openCsvWriter(csvPath)) {
            // header
            out.write("MC,completed,Ns_running,Nq_running,Ds_running,Dq_running,Ns_th,Nq_th,Ds_th,Dq_th");
            out.newLine();

            while (completed < MAX_COMPLETED) {
                double nextDepartureTime = inService.isEmpty() ? Double.POSITIVE_INFINITY : inService.peek().departure;
                double tNext = Math.min(nextArrival, nextDepartureTime);

                double dt = tNext - lastT;
                areaNs += n * dt;
                areaNq += Math.max(0, n - busy) * dt;
                MC = tNext;
                lastT = tNext;

                if (nextArrival <= nextDepartureTime) {
                    // arrival
                    n++;
                    if (busy < c) {
                        busy++;
                        double svc = exp(mu);
                        inService.add(new Job(MC, MC, MC + svc));
                    } else {
                        Q.addLast(MC);
                    }
                    nextArrival = MC + exp(lambda);
                } else {
                    // departure
                    Job j = inService.poll();
                    completed++;
                    double Ds = MC - j.arrival;
                    double Dq = j.start - j.arrival;
                    sumDs += Ds;
                    sumDq += Dq;
                    n--;

                    if (!Q.isEmpty()) {
                        double arrTime = Q.removeFirst();
                        double svc = exp(mu);
                        inService.add(new Job(arrTime, MC, MC + svc));
                    } else {
                        busy--;
                    }
                }

                if (completed > 0 && (completed % LOG_EVERY_K == 0)) {
                    double Ns_run = areaNs / MC;
                    double Nq_run = areaNq / MC;
                    double Ds_run = sumDs / completed;
                    double Dq_run = sumDq / completed;

                    // write a CSV row
                    String row = String.format(
                            "%.9f,%d,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f",
                            MC, completed, Ns_run, Nq_run, Ds_run, Dq_run,
                            theory.Ns, theory.Nq, theory.Ds, theory.Dq);
                    out.write(row);
                    out.newLine();
                }
            }
            out.flush();
        }
    }

    private static BufferedWriter openCsvWriter(String csvPath) throws IOException {
        Path path = Paths.get(csvPath);
        if (!path.isAbsolute()) {
            path = Paths.get("").toAbsolutePath().resolve(path);
        }
        Path parent = path.getParent();
        if (parent != null) {
            Files.createDirectories(parent); // ensure folders exist
        }
        // Overwrite if exists; UTF-8
        return Files.newBufferedWriter(path,
                java.nio.charset.StandardCharsets.UTF_8,
                StandardOpenOption.CREATE,
                StandardOpenOption.TRUNCATE_EXISTING,
                StandardOpenOption.WRITE);
    }


    //  Report Function
    public void report() {
        double Ns_run = areaNs / MC;
        double Nq_run = areaNq / MC;
        double Ds_run = sumDs / completed;
        double Dq_run = sumDq / completed;

        System.out.println("\n===== M/M/" + c + " Queue Simulation Report =====");
        System.out.printf( "λ (arrival rate): %.4f%n", lambda);
        System.out.printf( "μ (service rate): %.4f%n", mu);
        System.out.printf( "Servers (c): %d%n", c);
        System.out.printf( "ρ (utilization): %.4f%n", theory.rho);
        System.out.printf( "Completed customers: %d%n", completed);
        System.out.printf( "Total simulation time (MC): %.4f%n", MC);
        System.out.println("-----------------------------------------------");
        System.out.println("Metric        | Simulation   | Theory       | Rel.Error%");
        System.out.println("-----------------------------------------------");
        printCompare("Ns", Ns_run, theory.Ns);
        printCompare("Nq", Nq_run, theory.Nq);
        printCompare("Ds", Ds_run, theory.Ds);
        printCompare("Dq", Dq_run, theory.Dq);
        System.out.println("-----------------------------------------------\n");
    }

    private void printCompare(String label, double sim, double th) {
        double err = 100.0 * Math.abs(sim - th) / th;
        System.out.printf( "%-12s | %-12.6f | %-12.6f | %6.2f%%%n", label, sim, th, err);
    }

    // Drive
    public static void main(String[] args) throws Exception {
        double lambda = 0.18;  // arrival rate
        double mu = 0.25;      // service rate per server
        int c = 3;            // number of servers
        long MAX_COMPLETED = 200_000;
        int LOG_EVERY_K = 1;

        MMC_Infinite sim = new MMC_Infinite(lambda, mu, c, MAX_COMPLETED, LOG_EVERY_K);
        sim.run("src/MMC_Infinite/results/mmc_timeseries.csv");
        sim.report();
    }
}
