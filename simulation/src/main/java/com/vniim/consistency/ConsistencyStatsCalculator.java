package com.vniim.consistency;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Random;

import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.optim.linear.LinearConstraintSet;
import org.apache.commons.math3.optim.linear.LinearObjectiveFunction;
import org.apache.commons.math3.optim.linear.NonNegativeConstraint;
import org.apache.commons.math3.optim.linear.Relationship;
import org.apache.commons.math3.optim.linear.SimplexSolver;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.rank.Median;

/**
 * Just some test app to illustrate the simple (linear optimization based)
 * case of http://www.sciencedirect.com/science/article/pii/S0263224115006909
 *
 * @author a.stepanov
 */
public class ConsistencyStatsCalculator {

    /**
     * minimum number of trials
     */
    public static final int NSIM_MIN = 100;

    /**
     * maximum number of trials
     */
    public static final int NSIM_MAX = 100000000;

    /**
     * default number of Monte Carlo simulation trials
     */
    public static final int DEFAULNTSIM = 1000000;

    /**
     * number of histogram bins
     */
    public static final int NHISTBINS = 20;

    /**
     * random generator seed
     */
    public static final long SEED = 1234512345L;

    private static final String ENDL = System.lineSeparator();

    /**
     * measurement data and corresponding measurement uncertainties
     */
    private ArrayList<Double> x, u;

    /**
     * number of participants (length of x). should be at least 3.
     */
    private int n;

    /**
     * internal storage for uncertainty extension coefficients, n x nSimulations
     */
    private double K[][];

    /**
     * simulated reference values and corresponding uncertainties
     */
    private double XREF[];
    private double UREF[];

    /**
     * internal storage for histogram data
     */
    private double HIST[][];

    private final String inputData;
    private final String outDir;
    private final int nSim;

    /**
     * Constructor.
     * @param in input file name
     * @param out output dir name
     * @param nSimulations number of simulation runs
     */
    public ConsistencyStatsCalculator(String in, String out, int nSimulations) {

        inputData = in;
        outDir = out;
        nSim = nSimulations;
    }

    /**
     * Load data from text file. Expected format: pairs "<value> <uncertainty>".
     * The file may contain comment lines starting from "//" or empty lines.
     * @param filePath the file path
     * @throws IOException in case of file reading issues
     */
    private void loadData(String filePath) throws IOException {

        x = new ArrayList<>();
        u = new ArrayList<>();

        int ln = 1;
        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (!(line.startsWith("//") || // start comments with "//"
                      line.isEmpty())) { 
                    String tmp[] = line.split("\\s+");
                    if (tmp.length != 2) {
                        throw new IllegalArgumentException(
                            filePath + ": invalid input data at line " + ln + ", " +
                            "correct line format is \"<value> <uncertainty>\"");
                    }
                    double xx = Double.parseDouble(tmp[0]);
                    double uu = Double.parseDouble(tmp[1]);
                    if (Double.compare(uu, 0.) <= 0) {
                        throw new IllegalArgumentException(
                            filePath + ": invalid input data: " +
                            "non-positive uncertainty at line " + ln);
                    }
                    x.add(xx);
                    u.add(uu);
                }
                ++ln;
            }
        }
        n = x.size();
        if (n < 3) {
            throw new IllegalArgumentException(
                "too few data, at least 3 participants are expected");
        }
    }

    /**
     * Check pair consistency of the measurement data.
     * @return true if the data are consistent; otherwise false
     */
    public boolean checkConsistency() {

        if (x.isEmpty()) { throw new IllegalArgumentException("empty data"); }

        for (int i = 0; i < n; ++i) {
            double ui2 = u.get(i) * u.get(i);
            for (int j = i + 1; j < n; ++j) {
                double E = 0.5 * Math.abs(x.get(i) - x.get(j)) / Math.sqrt(
                    ui2 + u.get(j) * u.get(j));
                if (Double.compare(E, 1.) > 0) { return false; }
            }
        }

        return true;
    }

    /**
     * Increase uncertainties for simulated data.
     * @param y the simulated data
     * @return the increased uncertainties
     */
    private double[] increaseUncertainties(double y[]) {

        if (y.length != n) {
            throw new IllegalArgumentException("invalid input data size");
        }

        final int n1 = (n * (n - 1)) / 2 + n;
        double of[] = new double[n1];
        for (int i = 0; i < n; ++i) { of[i] = u.get(i) * u.get(i); }
        LinearObjectiveFunction f = new LinearObjectiveFunction(of, 0.);

        Collection<LinearConstraint> c = new ArrayList<>();
        int l = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                double lc[] = new double[n1];
                lc[i] = u.get(i) * u.get(i);
                lc[j] = u.get(j) * u.get(j);
                lc[n + l] = -1.;
                c.add(new LinearConstraint(lc, Relationship.EQ,
                        0.25 * (y[i] - y[j]) * (y[i] - y[j])));
                ++l;
            }
        }

        for (int i = 0; i < n; ++i) {
            double lc[] = new double[n1];
            lc[i] = 1.;
            c.add(new LinearConstraint(lc, Relationship.GEQ, 1.));
        }

        PointValuePair solution = (new SimplexSolver()).optimize(
                new MaxIter(10000),
                f,
                new LinearConstraintSet(c),
                GoalType.MINIMIZE,
                new NonNegativeConstraint(true));

        double s[] = solution.getPoint();
        double res[] = new double[n];
        for (int i = 0; i < n; ++i) { res[i] = Math.sqrt(s[i]); }
        return res;
    }

    /**
     * Get a pair of reference value/uncertainty using the initial measured
     * values and increased (expanded) initial uncertainties.
     * @param coeffs the expansion coefficients for the initial uncertainties
     * @return the pair (ref. value, ref. uncertainty)
     */
    private double[] getRef(double coeffs[]) {

        if (coeffs.length != n) {
            throw new IllegalArgumentException("invalid number of coefficients");
        }

        double w[] = new double[n]; // weights
        double sw = 0.;
        double xRef = 0.;
        for (int i = 0; i < n; ++i) {
            double tmp = u.get(i) * coeffs[i]; // extended uncertainty
            w[i] = 1. / (tmp * tmp);
            xRef += x.get(i) * w[i];
            sw += w[i];
        }
        xRef /= sw;

        double uRef = Math.pow(sw, -0.5);
        return new double[]{xRef, uRef};
    }

    /**
     * Calculate stats for k_i, xRef and save them to the file.
     */
    private void calcAndSaveStats() {

        String sOut = outDir + File.separator + "out.txt";
        BufferedWriter w = null;
        try {
            //create a temporary file
            File f = new File(sOut);
            if (f.exists()) {
                System.out.println("WARNING: overwriting " + sOut);
            }
            w = new BufferedWriter(new FileWriter(f));

            StandardDeviation sd = new StandardDeviation();
            Mean mn = new Mean();
            Median md = new Median();

            w.write("// mean, median, stdev" + ENDL + ENDL);
            double mean, median, stdev;

            for (int i = 0; i < n; ++i) {

                mean   = mn.evaluate(K[i]);
                median = md.evaluate(K[i]);
                stdev  = sd.evaluate(K[i]);

                w.write("k_" + (i + 1) + ":\t" +
                    mean + "\t" + median + "\t" + stdev + ENDL);
            }

            mean   = mn.evaluate(XREF);
            median = md.evaluate(XREF);
            stdev  = sd.evaluate(XREF);

            w.write(ENDL + "xRef:\t" +
                mean + "\t" + median + "\t" + stdev + ENDL);

        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            try {
                if (w != null) { w.close(); }
            } catch (IOException e) {}
        }
    }

    /**
     * Run simulations.
     */
    public void run() {

        try { loadData(inputData); }
        catch (IOException e) {
            System.err.println("Error while loading data");
            e.printStackTrace();
            return;
        }

        if (checkConsistency()) {
            System.out.println("the initial data are consistent; exiting...");
            return;
        }

        System.out.println("running simulations, n = " + nSim);

        Random rnd = new Random(SEED);

        if ((nSim < NSIM_MIN) || (nSim > NSIM_MAX)) {
            throw new IllegalArgumentException("please set number of " +
                "simulations betbeen " + NSIM_MIN + " and " + NSIM_MAX);
        }

        K = new double[n][nSim];
        XREF = new double[nSim];
        UREF = new double[nSim];

        int step = nSim / 10;

        for (int i = 0; i < nSim; ++i) {

            double y[] = new double[n];
            for (int j = 0; j < n; ++j) {
                y[j] = x.get(j) + u.get(j) * rnd.nextGaussian();
            }

            double k[] = increaseUncertainties(y);
            for (int j = 0; j < n; ++j) {
                K[j][i] = k[j];
            }

            double ref[] = getRef(k);
            XREF[i] = ref[0];
            UREF[i] = ref[1];

            if (i % step == 0) {
                double progress = (100. * i) / nSim;
                System.out.println(Math.round(progress) + "% done");
            }
        }
        System.out.println("100% done");

        double h[][];
        HIST = new double[2 * n + 4][NHISTBINS];
        for (int i = 0; i < n; ++i) {
            h = histogram(K[i]);
            HIST[2 * i]     = h[0];
            HIST[2 * i + 1] = h[1];
        }
        h = histogram(XREF);
        HIST[2 * n] = h[0];
        HIST[2 * n + 1] = h[1];
        h = histogram(UREF);
        HIST[2 * n + 2] = h[0];
        HIST[2 * n + 3] = h[1];

        System.out.println("calculating/saving stats...");

        calcAndSaveStats();

        System.out.println("drawing/saving histograms...");
        new Thread() {
            @Override
            public void run() {
                javafx.application.Application.launch(PlotContainer.class);
            }
        }.start();
        PlotContainer.getInstance().showAndSavePlots(HIST, outDir);
    }

    private static double[][] histogram(double data[]) {

        if ((data == null) || (data.length < 1)) {
            throw new IllegalArgumentException("empty data");
        }

        double max = data[0], min = data[0];
        for (int i = 1; i < data.length; ++i) {
            if (data[i] > max) { max = data[i]; }
            if (data[i] < min) { min = data[i]; }
        }

        // if max == min then all the data go to the 1st bin
        double h = (max - min) / NHISTBINS;

        double res[][] = new double[2][NHISTBINS];

        for (int i = 0; i < NHISTBINS; ++i) {
            res[0][i] = min + i * h;
        }

        double s = 0.;
        for (double v: data) {

            int i = (int) Math.floor((v - min) / h);
            if (i == NHISTBINS) { --i; } // special case
            ++res[1][i];
            ++s;
        }

        // normalize
        for (int i = 0; i < NHISTBINS; ++i) { res[1][i] /= s; }
        return res;
    }


    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        int nArgs = args.length;
        if (nArgs < 2) {
            System.out.println("Usage: ConsistencyStatsCalculator "
                + "<input_file> <out_dir> [nsimulations]");
            return;
        }

        String sInputData = args[0];
        if (!(new File(sInputData)).exists()) {
            System.err.println("input file " + sInputData + " does not exist");
            return;
        }

        String sOutDir = args[1];
        if (!(new File(sOutDir)).exists()) {
            System.err.println("output directory " + sOutDir + " does not exist");
            return;
        }

        int nSimulations = DEFAULNTSIM;
        if (nArgs > 2) {
            nSimulations = Integer.parseInt(args[2]);
        }

        (new ConsistencyStatsCalculator(
            sInputData, sOutDir, nSimulations)).run();
    }
}
