    static final float GREY_SCALE_LO = 0.95f, GREY_SCALE_HI = 1.05f ;
        // Clipping, for display only.  See for example Figure 1 in:
        //    http://bigwww.epfl.ch/thevenaz/shepplogan/

    public static void main(String [] args) {

        double [] [] density = new double [N] [N] ;

        for(int i = 0 ; i < N ; i++) {
            double x = SCALE * (i - N/2) ;
            for(int j = 0 ; j < N ; j++) {
                double y = SCALE * (j - N/2) ;

                density [i] [j] = sheppLoganPhantom(x, y) ;
            }
        }

        DisplayDensity display1 =
                new DisplayDensity(density, N, "Source Model",
                                   GREY_SCALE_LO, GREY_SCALE_HI) ;

        // Radon tranform of density (as measured by detectors):

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

        DisplayDensity display2 = new DisplayDensity(sinogram, N, "Sinogram") ;

        // inferred integral of density points (actually sum of density
        // points, here) for laternormalization of reconstruction

        double normDensity = norm1(sinogram [0]) ;


        // ... Insert sinogram filtering code here! ...


        double [] [] backProjection = new double [N] [N] ;
        backProject(backProjection, sinogram) ;

        // Normalize reconstruction, to have same sum as inferred for
        // original density

        double factor = normDensity / norm2(backProjection) ;
        for(int i = 0 ; i < N ; i++) {
            for(int j = 0 ; j < N ; j++) {
                backProjection [i] [j] *= factor ;
            }
        }

        DisplayDensity display5 =
                new DisplayDensity(backProjection, N,
                                   "Back projected sinogram") ;
    }