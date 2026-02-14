public class FFTImageFiltering_parallel {

    public static int N = 512 ;

    public static void main(String [] args) throws Exception {

        double [] [] X = new double [N] [N] ;
        ReadPGM.read(X, "OIP.pgm", N);

        DisplayDensity display =
                new DisplayDensity(X, N, "Original Image") ;

        // create array for in-place FFT, and copy original data to it
        double [] [] CRe = new double [N] [N], CIm = new double [N] [N] ;
        for(int k = 0 ; k < N ; k++) {
            for(int l = 0 ; l < N ; l++) {
                CRe [k] [l] = X [k] [l] ;
            }
        }

        // log for time 
        long startTime = System.nanoTime(); 

        fft2dParallel(CRe, CIm, 1) ;  // Fourier transform

        long endTime = System.nanoTime();
        double duration = (endTime - startTime) / 1000000.0;
        System.out.println("2D FFT (Forward) execution time: " + duration + " ms");

        Display2dFT display2 =
                new Display2dFT(CRe, CIm, N, "Discrete FT") ;

        // ======================== FILTERING ========================
        // int cutoff = N / 8;  // Define the frequency cutoff boundary
        
        // for (int k = 0; k < N; k++) {
        //     // Shift index to handle the periodic nature of FFT (centered at 0)
        //     int kSigned = (k <= N / 2) ? k : k - N;
            
        //     for (int l = 0; l < N; l++) {
        //         int lSigned = (l <= N / 2) ? l : l - N;

        //         /* * HIGH-PASS FILTER: 
        //          * Removes low-frequency components near the center (DC offset and smooth gradients).
        //          * Only frequencies outside the 'cutoff' box are preserved.
        //          */
        //         // if (Math.abs(kSigned) <= cutoff && Math.abs(lSigned) <= cutoff) {
        //         //     CRe[k][l] = 0;
        //         //     CIm[k][l] = 0;
        //         // }

        //         /* * LOW-PASS FILTER: 
        //          * Removes high-frequency components (noise and sharp edges).
        //          * Only frequencies inside the 'cutoff' box are preserved.
        //          *
        //          * */ 
        //         if (Math.abs(kSigned) > cutoff || Math.abs(lSigned) > cutoff) {
        //             CRe[k][l] = 0;
        //             CIm[k][l] = 0;
        //         }
    
        //     }
        // }
        // ======================== END FILTERING =====================

        // create array for in-place inverse FFT, and copy FT to it
        double [] [] reconRe = new double [N] [N],
                     reconIm = new double [N] [N] ;
        for(int k = 0 ; k < N ; k++) {
            for(int l = 0 ; l < N ; l++) {
                reconRe [k] [l] = CRe [k] [l] ;
                reconIm [k] [l] = CIm [k] [l] ;
            }
        }

        startTime = System.nanoTime(); 
        fft2dParallel(reconRe, reconIm, -1) ;  // Inverse Fourier transform
        endTime = System.nanoTime();
        duration = (endTime - startTime) / 1000000.0;
        System.out.println("2D FFT (Inverse) execution time: " + duration + " ms");

        DisplayDensity display3 =
                new DisplayDensity(reconRe, N, "Reconstructed Image") ;
    }

    public static void fft1d(double [] re, double [] im, int isgn) {

        // One-dimensional FFT, or inverse FFT (in-place algorithm).

        // When this method is called, the arrays re and im should contain
        // the real and imaginary parts of the input data.

        // When this method returns the values in these arrays are
        // are overwritten with the real and imaginary parts of the
        // transformed data.

        // isgn = +1 or -1 for forward or inverse transform.

        // Size of arrays should be a power or two.

        final double pi = Math.PI ;

        final int N = re.length ;  // im better be the same size

        bitReverse(re, im) ;

        int ln2   = ilog2(N)  ;  // Base 2 log of the leading dimension.

        // Danielson-Lanczos algorithm for FFT.

        for(int ilevel = 1 ; ilevel <= ln2 ; ilevel++) {
            int le   = ipow(2,ilevel) ;
            int lev2 = le / 2 ;

            double uRe = 1.0F ;
            double uIm = 0.0F ;

            double wRe = Math.cos(isgn * pi / lev2) ;
            double wIm = Math.sin(isgn * pi / lev2) ;

            for(int jj = 0 ; jj < lev2 ; jj++) {
                for(int ii = jj ; ii < N ; ii += le) {
                    int jndex = ii + lev2 ;
                    int index = ii ;

                    //tmp      = u * a(jndex) ;
                    double tmpRe = uRe * re [jndex] - uIm * im [jndex] ;
                    double tmpIm = uRe * im [jndex] + uIm * re [jndex] ;

                    //a(jndex) = a(index) - tmp ;
                    re [jndex] = re [index] - tmpRe ;
                    im [jndex] = im [index] - tmpIm ;

                    //a(index) = a(index) + tmp ;
                    re [index] = re [index] + tmpRe ;
                    im [index] = im [index] + tmpIm ;
                }
                //tmp = u * w ;
                double tmpRe = uRe * wRe - uIm * wIm ;
                double tmpIm = uRe * wIm + uIm * wRe ;

                //u   = tmp ;
                uRe   = tmpRe ;
                uIm   = tmpIm ;
            }
        }
    }

    public static void fft2dParallel(double[][] re, double[][] im, int isgn) {
        int N = re.length;

        // 1. Parallel FFT on all rows
        // IntStream.range creates a pool of tasks for the common ForkJoinPool
        java.util.stream.IntStream.range(0, N).parallel().forEach(i -> {
            fft1d(re[i], im[i], isgn);
        });

        // 2. Transpose (This is usually fast, but can be synchronized or blocked)
        transpose(re);
        transpose(im);

        // 3. Parallel FFT on all columns (now stored as rows after transpose)
        java.util.stream.IntStream.range(0, N).parallel().forEach(i -> {
            fft1d(re[i], im[i], isgn);
        });

        // 4. Transpose back
        transpose(re);
        transpose(im);
    }


    public static void transpose(double[][] a) {
        int N = a.length;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < i; j++) {
                // Swap a[i][j] and a[j][i]
                double temp = a[i][j];
                a[i][j] = a[j][i];
                a[j][i] = temp;
            }
        }
    }
    static void bitReverse(double [] re, double [] im) {

        // In place permutation to bit-reversed ordering.

        final int N = re.length ;  // im better be the same size

        int nm1   = N - 1 ;
        int nv2   = N / 2 ;

        for(int index = 0, jndex = 0 ; index < nm1 ; index++) {
            if(jndex > index) {
                
                // Swap entries
                
                double tmpRe = re [jndex] ;
                double tmpIm = im [jndex] ;

                re [jndex] = re [index] ;
                im [jndex] = im [index] ;

                re [index] = tmpRe ;
                im [index] = tmpIm ;
            }

            int m = nv2 ;
            while ((m >= 2) && (jndex >= m)) {
                jndex = jndex - m ;
                m = m / 2 ;
            }
            jndex = jndex + m ;
        }
    }

    static int ipow(int i, int j) {

        int k, tmp ;

        tmp = 1 ;
        for(k = 1 ; k <= j ; k++)
            tmp = tmp * i ;
        return tmp ;
    }

    static int ilog2(int n) {

        int i, n2, result ;

        n2     = n ;
        result = 0 ;
        for(i = 1 ; i <= n ; i++) {
            if(n2 > 1) {
                result = result + 1 ;
                n2 = n2 / 2 ;
            }
            else
                break ;
        }
        return result ;
    }
}