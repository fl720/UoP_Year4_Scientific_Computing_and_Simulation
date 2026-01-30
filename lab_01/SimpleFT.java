
import java.io.*; // Import this for saving files
public class SimpleFT {

    public static int N = 256 ;

    public static void main(String [] args) throws Exception {

        double [] [] X = new double [N] [N] ;
        ReadPGM.read(X, "wolf.pgm", N) ;
        // Add this check:
        System.out.println("Pixel at [100][100] is: " + X[100][100]);
        
        DisplayDensity display =
                new DisplayDensity(X, N, "Original Image") ;

        double [] [] CRe = new double [N] [N], CIm = new double [N] [N] ;

        for(int k = 0 ; k < N ; k++) {
            for(int l = 0 ; l < N ; l++) {
                double sumRe = 0, sumIm = 0 ;
                // Nested for loops performing sum over X elements
            //   for(int m = ...) {
            //       for(int n = ...) {
            //            double arg = ... ;
            //            double cos = ... ;
            //            double sin = ... ;
            //            sumRe += cos * X [m] [n] ;
            //            sumIm += ... ;
            //       }
            //   }
                for( int m = 0 ; m < N ; m++) {
                    for(int n = 0 ; n < N ; n++) {
                        double arg = -2 * Math.PI * ( (double)(k * m) / N +
                                                    (double)(l * n) / N ) ;
                        double cos = Math.cos(arg) ;
                        double sin = Math.sin(arg) ;
                        sumRe += cos * X [m] [n] ;
                        sumIm += sin * X [m] [n] ;
                    }
                }
                CRe [k] [l] = sumRe ;
                CIm [k] [l] = sumIm ;
            }
            System.out.println("Completed FT line " + k + " out of " + N) ;
        }

        Display2dFT display2 =
                new Display2dFT(CRe, CIm, N, "Discrete FT") ;

                  int cutoff = N/8 ;  // for example
          for(int k = 0 ; k < N ; k++) {
              int kSigned = k <= N/2 ? k : k - N ;
              for(int l = 0 ; l < N ; l++) {
                  int lSigned = l <= N/2 ? l : l - N ;
                  if(Math.abs(kSigned) > cutoff || Math.abs(lSigned) > cutoff) {
                      CRe [k] [l] = 0 ;
                      CIm [k] [l] = 0 ;
                  }
              }
          }

          Display2dFT display2a =
                  new Display2dFT(CRe, CIm, N, "Truncated FT") ;


    //   Now do inverse FT
        double [] [] reconstructed = new double [N] [N] ;

        for(int m = 0 ; m < N ; m++) {
            for(int n = 0 ; n < N ; n++) {
                double sum = 0  ;
                // Nested for loops performing sum over C elements
                for(int k = 0 ; k < N ; k++) {
                    for(int l = 0 ; l < N ; l++) {
                        double arg = 2 * Math.PI * ( (double)(m * k) / N +
                                                    (double)(n * l) / N ) ;
                        double cos = Math.cos(arg) ;
                        double sin = Math.sin(arg) ;
                        sum += cos * CRe [k] [l] - sin * CIm [k] [l] ;
                    }
                }
                reconstructed [m] [n] = sum / (N * N) ;
            //   reconstructed [m] [n] = sum ;
            }
            System.out.println("Completed inverse FT line " + m + " out of " + N) ;
        }

        DisplayDensity display3 =
                new DisplayDensity(reconstructed, N, "Reconstructed Image") ;

        //   save reconstructed image
        saveToPGM(reconstructed, "reconstructed_wolf.pgm");

    }

    // SAVE TO PGM 
    public static void saveToPGM(double[][] image, String filename) {
        try {
            int height = image.length;
            int width = image[0].length;
            PrintWriter writer = new PrintWriter(new FileWriter(filename));
            
            writer.println("P2");
            writer.println(width + " " + height);
            writer.println("255"); 
            
            for (int i = 0; i < height; i++) {
                for (int j = 0; j < width; j++) {
                    int pixel = (int) image[i][j];
                    if (pixel < 0) pixel = 0;
                    if (pixel > 255) pixel = 255;
                    writer.print(pixel + " ");
                }
                writer.println();
            }
            writer.close();
            System.out.println("Saved image to: " + filename);
        } catch (IOException e) {
            System.out.println("Error saving file: " + e.getMessage());
        }
    }
}





