  public class SimpleFT {
  
      public static int N = 256 ;
  
      public static void main(String [] args) throws Exception {
  
          double [] [] X = new double [N] [N] ;
          ReadPGM.read(X, "wolf.pgm", N) ;
  
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
      }
  }