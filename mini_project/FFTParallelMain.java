import java.io.File;
import java.awt.image.BufferedImage;
import javax.imageio.ImageIO;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.IntStream;

// ================================ HOW TO RUN ================================
// javac FFTParallelMain.java
// java FFTParallelMain [imagePath] [numThreads] [filterType] [cutoffDivisor]

// Example:
// java FFTParallelMain img/2048.png 4    means use 4 threads with DEFAULT lowpass filter and cutoff divisor of 8


public class FFTParallelMain {

    public static int    N              = 4096; // !!! REMEMBER change this for different image sizes (must be power of 2)!!!!!!!!!!!!!
    public static int    NUM_THREADS    = Runtime.getRuntime().availableProcessors();
    public static String FILTER_TYPE    = "lowpass";
    public static int    CUTOFF_DIVISOR = 8;

    private static ForkJoinPool pool;

    public static void main(String [] args) throws Exception {

        String imagePath = args.length > 0 ? args[0] : "img/2048.png";
        if (args.length > 1) NUM_THREADS    = Integer.parseInt(args[1]);
        if (args.length > 2) FILTER_TYPE    = args[2];
        if (args.length > 3) CUTOFF_DIVISOR = Integer.parseInt(args[3]);

        pool = new ForkJoinPool(NUM_THREADS);

        System.out.println("Threads: " + NUM_THREADS);
        System.out.print("Warming up... ");
        warmup(256, 10);
        System.out.println("Done.");

        double [][] X = new double [N][N];
        readPNG(X, imagePath, N);

        new DisplayDensity(X, N, "1. Original Image");

        double [][] CRe = new double [N][N], CIm = new double [N][N];
        for (int i = 0; i < N; i++) System.arraycopy(X[i], 0, CRe[i], 0, N);

        // --- Forward FFT SECTION ---
        long start = System.nanoTime();
        fft2dParallel(CRe, CIm, 1);
        long end = System.nanoTime();
        System.out.println("Forward FFT: " + (end - start)/1e6 + " ms");

        // Spectrum 
        double[][] CReCopy = deepCopy(CRe);
        double[][] CImCopy = deepCopy(CIm);
        new Display2dFT(CReCopy, CImCopy, N, "2. Discrete FT (Spectrum)");

        // --- Filtering SECTION ---
        applyFilter(CRe, CIm, FILTER_TYPE, CUTOFF_DIVISOR);

        // --- Inverse FFT SECTION ---
        start = System.nanoTime();
        fft2dParallel(CRe, CIm, -1);
        end = System.nanoTime();
        System.out.println("Inverse FFT: " + (end - start)/1e6 + " ms");

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) CRe[i][j] /= (N * N);
        }

        new DisplayDensity(CRe, N, "3. Reconstructed Image");

        pool.shutdown();
    }

    private static double[][] deepCopy(double[][] original) {
        double[][] copy = new double[N][N];
        for (int i = 0; i < N; i++) System.arraycopy(original[i], 0, copy[i], 0, N);
        return copy;
    }

    private static void applyFilter(double[][] re, double[][] im, String type, int div) {
        int cutoff = N / div;
        for (int k = 0; k < N; k++) {
            int kS = (k <= N / 2) ? k : k - N;
            for (int l = 0; l < N; l++) {
                int lS = (l <= N / 2) ? l : l - N;
                boolean isLowFreq = Math.abs(kS) <= cutoff && Math.abs(lS) <= cutoff;
                if ((type.equalsIgnoreCase("highpass") && isLowFreq) ||
                    (!type.equalsIgnoreCase("highpass") && !isLowFreq)) {
                    re[k][l] = 0; im[k][l] = 0;
                }
            }
        }
    }

    public static void fft2dParallel(double[][] re, double[][] im, int isgn) {
        int n = re.length;
        try {
            pool.submit(() -> IntStream.range(0, n).parallel().forEach(i -> fft1d(re[i], im[i], isgn))).get();
            transpose(re); transpose(im);
            pool.submit(() -> IntStream.range(0, n).parallel().forEach(i -> fft1d(re[i], im[i], isgn))).get();
            transpose(re); transpose(im);
        } catch (Exception e) { e.printStackTrace(); }
    }

    private static void warmup(int size, int iter) {
        double[][] r = new double[size][size], i = new double[size][size];
        for (int n = 0; n < iter; n++) { fft2dParallel(r, i, 1); fft2dParallel(r, i, -1); }
    } 

    public static void fft1d(double [] re, double [] im, int isgn) {
        final double pi = Math.PI ;
        final int N = re.length ;
        bitReverse(re, im) ;
        int ln2 = ilog2(N) ;
        for (int ilevel = 1 ; ilevel <= ln2 ; ilevel++) {
            int le = ipow(2, ilevel) ;
            int lev2 = le / 2 ;
            double uRe = 1.0, uIm = 0.0 ;
            double wRe = Math.cos(isgn * pi / lev2) ;
            double wIm = Math.sin(isgn * pi / lev2) ;
            for (int jj = 0 ; jj < lev2 ; jj++) {
                for (int ii = jj ; ii < N ; ii += le) {
                    int jndex = ii + lev2, index = ii ;
                    double tmpRe = uRe * re [jndex] - uIm * im [jndex] ;
                    double tmpIm = uRe * im [jndex] + uIm * re [jndex] ;
                    re [jndex] = re [index] - tmpRe ;
                    im [jndex] = im [index] - tmpIm ;
                    re [index] = re [index] + tmpRe ;
                    im [index] = im [index] + tmpIm ;
                }
                double tmpRe = uRe * wRe - uIm * wIm ;
                double tmpIm = uRe * wIm + uIm * wRe ;
                uRe = tmpRe ; uIm = tmpIm ;
            }
        }
    }

    public static void transpose(double [] [] a) {
        int N = a.length ;
        for (int i = 0 ; i < N ; i++) {
            for (int j = 0 ; j < i ; j++) {
                double temp = a [i] [j] ;
                a [i] [j] = a [j] [i] ;
                a [j] [i] = temp ;
            }
        }
    }

    static void bitReverse(double [] re, double [] im) {
        final int N = re.length ;
        int nm1 = N - 1, nv2 = N / 2 ;
        for (int index = 0, jndex = 0 ; index < nm1 ; index++) {
            if (jndex > index) {
                double tmpRe = re [jndex], tmpIm = im [jndex] ;
                re [jndex] = re [index] ; im [jndex] = im [index] ;
                re [index] = tmpRe ; im [index] = tmpIm ;
            }
            int m = nv2 ;
            while ((m >= 2) && (jndex >= m)) {
                jndex -= m ; m /= 2 ;
            }
            jndex += m ;
        }
    }

    static int ipow(int i, int j) {
        int tmp = 1 ;
        for (int k = 1 ; k <= j ; k++) tmp *= i ;
        return tmp ;
    }

    static int ilog2(int n) {
        int result = 0 ;
        while (n > 1) {
            result++ ; n /= 2 ;
        }
        return result ;
    }

    static void readPNG(double [] [] density, String fileName, int n) throws Exception {
        BufferedImage raw = ImageIO.read(new File(fileName)) ;
        if (raw == null) {
            System.out.println("Cannot read file: " + fileName) ;
            System.exit(1) ;
        }
        BufferedImage img = new BufferedImage(n, n, BufferedImage.TYPE_BYTE_GRAY) ;
        img.getGraphics().drawImage(raw, 0, 0, n, n, null) ;
        for (int j = 0 ; j < n ; j++) {
            for (int i = 0 ; i < n ; i++) {
                density [i] [n - 1 - j] = img.getRaster().getSample(i, j, 0) ;
            }
        }
    }
}