import java.util.Arrays ;
import java.awt.* ;
import javax.swing.* ;

public class Sinogram {

    static final int N = 512 ;
    static final int CELL_SIZE = 1 ;
    static final double SCALE = 0.0045 ; 
    static final int CUTOFF = N/4 ; 

    static final float GREY_SCALE_LO = 0.95f, GREY_SCALE_HI = 1.05f ;

    public static void main(String [] args) {

        double [] [] density = new double [N] [N] ;

        // Generate the Shepp-Logan Phantom source model
        for(int i = 0 ; i < N ; i++) {
            double x = SCALE * (i - N/2) ;
            for(int j = 0 ; j < N ; j++) {
                double y = SCALE * (j - N/2) ;
                density [i] [j] = sheppLoganPhantom(x, y) ;
            }
        } 

        DisplayDensity display1 = new DisplayDensity(density, N, "Source Model",
                           GREY_SCALE_LO, GREY_SCALE_HI) ;

        // Simulate CT scan: Radon transform of density 
        double [] [] sinogram = new double [N] [N] ;

        for(int iTheta = 0 ; iTheta < N ; iTheta++) {
            double theta = (Math.PI * iTheta) / N ;
            double cos = Math.cos(theta) ;
            double sin = Math.sin(theta) ;
            for(int iR = 0 ; iR < N ; iR++) {
                double r = SCALE * (iR - N/2) ;
                double sum = 0 ;
                for(int iS = 0 ; iS < N ; iS++) {
                    double s = SCALE * (iS - N/2) ;
                    double x = r * cos + s * sin ;
                    double y = r * sin - s * cos ;
                    sum += sheppLoganPhantom(x, y) ;
                }
                sinogram [iTheta] [iR] = sum ;
            }
        }

        DisplayDensity display2 = new DisplayDensity(sinogram, N, "pure Sinogram") ;

        // Normalization factor extracted from the Radon transform 
        double normDensity = norm1(sinogram [0]) ;

        // --- Sinogram Filtering Code Start --- 

        double[][] sinogramFTRe = new double[N][N];
        double[][] sinogramFTIm = new double[N][N];

        // Copy sinogram data to the real part of the FT array 
        for (int iTheta = 0; iTheta < N; iTheta++) {
            for (int iR = 0; iR < N; iR++) {
                sinogramFTRe[iTheta][iR] = sinogram[iTheta][iR];
            }
        }

        // 1. Perform 1D FFT on each row (detector offset dimension) 
        for (int iTheta = 0; iTheta < N; iTheta++) {
            FFT.fft1d(sinogramFTRe[iTheta], sinogramFTIm[iTheta], 1);
        }

        // Display the Fourier Transform of the Sinogram
        DisplaySinogramFT display3 = new DisplaySinogramFT(sinogramFTRe, sinogramFTIm, N, "Low pass cosine Filtered Sinogram radial Fourier Transform");

        // 2. Apply the Ramp Filter |K| ======== BACK PROJECTION ========
        for (int iTheta = 0; iTheta < N; iTheta++) {
            // !!!! ramp fileter ----------------------------------------------
            // for (int iK = 0; iK < N; iK++) {
            //     // Determine the signed wave number 
            //     int kSigned = (iK <= N / 2) ? iK : iK - N;
                
            //     // Multiply FT by abs(kSigned) 
            //     sinogramFTRe[iTheta][iK] *= Math.abs(kSigned);
            //     sinogramFTIm[iTheta][iK] *= Math.abs(kSigned);
            // }
            // Low pass cosine filter (optional) --------------------------------
            for (int iK = 0; iK < N; iK++) {
                int kSigned = (iK <= N / 2) ? iK : iK - N;
                double absK = Math.abs(kSigned);

                if (absK > CUTOFF) {
                    // Zero out components greater than CUTOFF [cite: 346]
                    sinogramFTRe[iTheta][iK] = 0;
                    sinogramFTIm[iTheta][iK] = 0;
                } else {
                    // Apply RL Filter Filter 
                    double cosineFactor = Math.cos((Math.PI * kSigned) / (2.0 * CUTOFF));
                    double filter = absK * cosineFactor;

                    sinogramFTRe[iTheta][iK] *= filter;
                    sinogramFTIm[iTheta][iK] *= filter;
                }
            }
            //  ===========================================================
        }

        // 3. Perform Inverse 1D FFT on each row 
        for (int iTheta = 0; iTheta < N; iTheta++) {
            FFT.fft1d(sinogramFTRe[iTheta], sinogramFTIm[iTheta], -1);
        }

        // Display the filtered sinogram in real space 
        DisplayDensity display4 = new DisplayDensity(sinogramFTRe, N, "Low pass cosine Filtered sinogram");

        // --- Sinogram Filtering Code End ---

        double [] [] backProjection = new double [N] [N] ;
        
        // Use the FILTERED sinogram for back projection
        backProject(backProjection, sinogramFTRe) ;
        // backProject(backProjection, sinogram) ;

        // Normalize reconstruction 
        double factor = normDensity / norm2(backProjection) ;
        for(int i = 0 ; i < N ; i++) {
            for(int j = 0 ; j < N ; j++) {
                backProjection [i] [j] *= factor ;
            }
        }

        // Display reconstructed image with clipping to see brain structures 
        DisplayDensity display5 = new DisplayDensity(backProjection, N, "Low pass cosine Filtered Back Projection", GREY_SCALE_LO, GREY_SCALE_HI) ;
        // DisplayDensity display5 = new DisplayDensity(backProjection, N, "Ramp Filtered Back Projection") ;
    }

    static void backProject(double [] [] projection, double [] [] sinogram) {
        // Back Projection operation: adding signals for all rays passing through (x,y) 
        for(int i = 0 ; i < N ; i++) {
            double x = SCALE * (i - N/2) ;
            for(int j = 0 ; j < N ; j++) {
                double y = SCALE * (j - N/2) ;

                double sum = 0 ;
                for(int iTheta = 0 ; iTheta < N ; iTheta++) {
                    double theta = (Math.PI * iTheta) / N ;
                    double cos = Math.cos(theta) ;
                    double sin = Math.sin(theta) ;

                    // Ray geometry: r = x*cos(theta) + y*sin(theta) 
                    double r = x * cos + y * sin ;
                    double rBox = N/2 + r/SCALE ;

                    if(rBox < 0) continue ; 

                    int iR = (int) rBox ; 
                    double offR = rBox - iR ;
                    int iPlusR = iR + 1 ; 

                    if(iPlusR >= N) continue ; 

                    // Linear interpolation of the sinogram value 
                    double sinogramVal =
                            (1 - offR) * sinogram [iTheta] [iR] +
                            offR * sinogram [iTheta] [iPlusR] ;
                    sum += sinogramVal ;
                }
                projection [i] [j] = sum ;
            }
        }
    }

    // Shepp-Logan Phantom ellipses configuration 
    static final Ellipse [] sheppLoganEllipses = {
        new Ellipse(0.0, 0.0, 0.69, 0.92, 0, 2.0),
        new Ellipse(0.0, -0.0184, 0.6624, 0.874, 0, -0.98),
        new Ellipse(0.22, 0, 0.11, 0.31, -18.0, -0.02),
        new Ellipse(-0.22, 0, 0.16, 0.41, 18.0, -0.02),
        new Ellipse(0, 0.35, 0.21, 0.25, 0, 0.01),
        new Ellipse(0, 0.1, 0.046, 0.046, 0, 0.01),
        new Ellipse(0, -0.1, 0.046, 0.046, 0, 0.01),
        new Ellipse(-0.08, -0.605, 0.046, 0.023, 0, 0.01),
        new Ellipse(0, -0.605, 0.023, 0.023, 0, 0.01),
        new Ellipse(0.06, -0.605, 0.023, 0.046, 0, 0.01),
    } ;

    static double sheppLoganPhantom (double x, double y) {
        double total = 0 ;
        for(Ellipse ellipse : sheppLoganEllipses) {
            total += ellipse.localDensity(x, y) ;
        }
        return total ;
    }

    static class Ellipse {
        double centreX, centreY, major, minor, theta, density, cos, sin ;

        Ellipse(double centreX, double centreY,
                double major, double minor, double theta, double density) {
            this.centreX = centreX ;
            this.centreY = centreY ;
            this.major = major ;
            this.minor = minor ;
            this.theta = theta ;
            double rad = Math.PI * theta / 180 ;
            cos = Math.cos(rad) ; 
            sin = Math.sin(rad) ; 
            this.density = density;
        }

        double localDensity(double x, double y) {
            double xOff = x - centreX ;
            double yOff = y - centreY ;
            double xRot = cos * xOff - sin * yOff ;
            double yRot = sin * xOff + cos * yOff ;
            double xNorm = xRot / major ;
            double yNorm = yRot / minor ;
            return (xNorm * xNorm + yNorm * yNorm < 1) ? density : 0 ;
        }
    }

    static double norm1(double [] density) {
        double norm = 0 ;
        for(int i = 0 ; i < N ; i++) norm += density [i] ;
        return norm ;
    }

    static double norm2(double [] [] density) {
        double norm = 0 ;
        for(int i = 0 ; i < N ; i++) {
            for(int j = 0 ; j < N ; j++) {
                if(density [i] [j] > 0) norm += density [i] [j] ;
            }
        }
        return norm ;
    }
}