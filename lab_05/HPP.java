import java.awt.* ;
import javax.swing.* ;

public class HPP {


    final static int NX = 80, NY = 60 ;  // Lattice dimensions
    final static int q = 4 ;  // population
     
    final static int NITER = 10000 ;
    final static int DELAY = 20 ;

    final static double DENSITY = 0.5 ;  // initial state, between 0 and 1.0.

    // static Display display = new Display() ;
    static Display_grey display = new Display_grey() ;
    // static Display_average_neighbour display = new Display_average_neighbour() ;
    // static Display_average_neighbour_closed display = new Display_average_neighbour_closed() ;

    static boolean [] [] [] fin = new boolean [NX] [NY] [q] ;
    static boolean [] [] [] fout = new boolean [NX] [NY] [q] ;

    public static void main(String args []) throws Exception {

        // initialise - just two particles. 
        // fin [NX/2+1] [NY/2] [0] = true;
        // fin [NX/2-1] [NY/2] [1] = true;

        // initialise - populate a subblock of grid
        for(int i = 0; i < NX/4 ; i++) { 
            for(int j = 0; j < NY/4 ; j++) { 
                boolean [] fin_ij = fin [i] [j] ;
                for(int d = 0 ; d < q ; d++) {
                    if(Math.random() < DENSITY) {
                        fin_ij [d] = true ;
                    }
                }
            }
        }


        display.repaint() ;
        Thread.sleep(DELAY) ;

        for(int iter = 0 ; iter < NITER ; iter++) {

            // Collision

            for(int i = 0; i < NX ; i++) { 
                for(int j = 0; j < NY ; j++) { 
                    boolean [] fin_ij = fin [i] [j] ;
                    boolean [] fout_ij = fout [i] [j] ;

                    // default, no collisions case:

                    // fout_ij [0] = fin_ij [0] ;
                    // fout_ij [1] = fin_ij [1] ;
                    // fout_ij [2] = fin_ij [2] ;
                    // fout_ij [3] = fin_ij [3] ;

                    // please add collisions as per lecture!
                    // CYLINDRICAL BOUNDARY CONDITIONS:
                    if (fin_ij[0] && fin_ij[1] && !fin_ij[2] && !fin_ij[3]) {
                        fout_ij[0] = false;
                        fout_ij[1] = false;
                        fout_ij[2] = true;
                        fout_ij[3] = true;
                    } 
                    else if (!fin_ij[0] && !fin_ij[1] && fin_ij[2] && fin_ij[3]) {
                        fout_ij[0] = true;
                        fout_ij[1] = true;
                        fout_ij[2] = false;
                        fout_ij[3] = false;
                    } 
                    else {
                        fout_ij[0] = fin_ij[0];
                        fout_ij[1] = fin_ij[1];
                        fout_ij[2] = fin_ij[2];
                        fout_ij[3] = fin_ij[3];
                    }
                }
            }

            // Streaming

            for(int i = 0; i < NX ; i++) { 
                int iP1 = (i + 1) % NX ;
                int iM1 = (i - 1 + NX) % NX ;
                for(int j = 0; j < NY ; j++) { 
                    int jP1 = (j + 1) % NY ;
                    int jM1 = (j - 1 + NY) % NY ;

                    // ----- closed box-like boundary conditions: -----
                    if ( i != NX-1 ){
                        fin [i] [j] [0] = fout [i+1] [j] [0] ;
                    } else {
                        fin [i] [j] [0] = fout [NX-1] [j] [1] ;
                    }

                    if ( i != 0 ){
                        fin [i] [j] [1] = fout [i-1] [j] [1] ;
                    } else {
                        fin [i] [j] [1] = fout [0] [j] [0] ;
                    }

                    if ( j != NY-1 ){
                    fin [i] [j] [2] = fout [i] [j+1] [2] ;
                    } else {
                        fin [i] [j] [2] = fout [i] [NY-1] [3] ;
                    }

                    if ( j != 0 ){
                    fin [i] [j] [3] = fout [i] [j-1] [3] ;
                    } else {
                        fin [i] [j] [3] = fout [i] [0] [2] ;
                    }

                    // please add streaming as per lecture!
                    // ----- cylindrical boundary conditions: -----
                    // fin [i] [j] [0] = fout [iP1] [j] [0] ;
                    // fin [i] [j] [1] = fout [iM1] [j] [1] ;
                    // fin [i] [j] [2] = fout [i] [jP1] [2] ;
                    // fin [i] [j] [3] = fout [i] [jM1] [3] ;
                }
            }

            System.out.println("iter = " + iter) ;
            display.repaint() ;
      
            Thread.sleep(DELAY) ;
        }
    }

    
    static class Display extends JPanel {

        final static int CELL_SIZE = 14 ; 

        public static final int ARROW_START = 2 ;
        public static final int ARROW_END   = 7 ;
        public static final int ARROW_WIDE  = 3 ;

        Display() {

            setPreferredSize(new Dimension(CELL_SIZE * NX, CELL_SIZE * NY)) ;

            JFrame frame = new JFrame("HPP");
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.setContentPane(this);
            frame.pack();
            frame.setVisible(true);
        }

        public void paintComponent(Graphics g) {

            g.setColor(Color.WHITE) ;
            g.fillRect(0, 0, CELL_SIZE * NX, CELL_SIZE * NY) ;

            g.setColor(Color.PINK) ;
            //g.setColor(Color.LIGHT_GRAY) ;
            for(int i = 0 ; i < NX ; i++) {
                for(int j = 0 ; j < NY ; j++) {
                    int originX = CELL_SIZE * i + CELL_SIZE/2 ;
                    int originY = CELL_SIZE * j + CELL_SIZE/2 ;
                    g.fillOval(originX - 2, originY - 2, 4, 4) ;
                }
            } 

            g.setColor(Color.BLUE) ;
            int [] tri_x = new int [3], tri_y = new int [3] ;
            for(int i = 0 ; i < NX ; i++) {
                for(int j = 0 ; j < NY ; j++) {
                    boolean [] fin_ij = fin [i] [j] ;
                    int originX = CELL_SIZE * i + CELL_SIZE/2 ;
                    int originY = CELL_SIZE * j + CELL_SIZE/2 ;
                    if(fin_ij [0]) {
                        tri_x [0] = originX - ARROW_START ;
                        tri_x [1] = originX - ARROW_START ;
                        tri_x [2] = originX - ARROW_END ;
                        tri_y [0] = originY - ARROW_WIDE ;
                        tri_y [1] = originY + ARROW_WIDE ;
                        tri_y [2] = originY ;
                        //g.setColor(Color.BLUE) ;
                        g.fillPolygon(tri_x, tri_y, 3) ;
                    }
                    if(fin_ij [1]) {
                        tri_x [0] = originX + ARROW_START ;
                        tri_x [1] = originX + ARROW_START ;
                        tri_x [2] = originX + ARROW_END ;
                        tri_y [0] = originY - ARROW_WIDE ;
                        tri_y [1] = originY + ARROW_WIDE ;
                        tri_y [2] = originY ;
                        //g.setColor(Color.RED) ;
                        g.fillPolygon(tri_x, tri_y, 3) ;
                    }
                    if(fin_ij [2]) {
                        tri_x [0] = originX - ARROW_WIDE ;
                        tri_x [1] = originX + ARROW_WIDE ;
                        tri_x [2] = originX  ;
                        tri_y [0] = originY - ARROW_START ;
                        tri_y [1] = originY - ARROW_START ;
                        tri_y [2] = originY - ARROW_END ;
                        //g.setColor(Color.GREEN) ;
                        g.fillPolygon(tri_x, tri_y, 3) ;
                    }
                    if(fin_ij [3]) {
                        tri_x [0] = originX - ARROW_WIDE ;
                        tri_x [1] = originX + ARROW_WIDE ;
                        tri_x [2] = originX  ;
                        tri_y [0] = originY + ARROW_START ;
                        tri_y [1] = originY + ARROW_START ;
                        tri_y [2] = originY + ARROW_END ;
                        //g.setColor(Color.YELLOW) ;
                        g.fillPolygon(tri_x, tri_y, 3) ;
                    }
                }
            } 
        }
    }

    static class Display_grey extends JPanel {

        final static int CELL_SIZE = 14 ; 

        Display_grey() {
            setPreferredSize(new Dimension(CELL_SIZE * NX, CELL_SIZE * NY)) ;

            JFrame frame = new JFrame("HPP - Density Map");
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.setContentPane(this);
            frame.pack();
            frame.setVisible(true);
        }

        public void paintComponent(Graphics g) {
            
            for(int i = 0 ; i < NX ; i++) {
                for(int j = 0 ; j < NY ; j++) {
                    boolean [] fin_ij = fin [i] [j] ;
                    
                    int density = 0;
                    if (fin_ij[0]) density++;
                    if (fin_ij[1]) density++;
                    if (fin_ij[2]) density++;
                    if (fin_ij[3]) density++;
                    
                    // density = 0 -> grayValue = 255 (white)
                    // density = 4 -> grayValue = 0 (black)
                    int grayValue = 255 - (density * 255 / 4);
                    g.setColor(new Color(grayValue, grayValue, grayValue));
                    
                    int originX = CELL_SIZE * i ;
                    int originY = CELL_SIZE * j ;
                    g.fillRect(originX, originY, CELL_SIZE, CELL_SIZE) ;
                }
            } 
        }
    }

    static class Display_average_neighbour extends JPanel {

        final static int CELL_SIZE = 14 ; 

        Display_average_neighbour() {
            setPreferredSize(new Dimension(CELL_SIZE * NX, CELL_SIZE * NY)) ;

            JFrame frame = new JFrame("HPP - 3x3 Average Density");
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.setContentPane(this);
            frame.pack();
            frame.setVisible(true);
        }

        public void paintComponent(Graphics g) {
            
            g.setColor(Color.WHITE) ;
            g.fillRect(0, 0, CELL_SIZE * NX, CELL_SIZE * NY) ;

            for(int i = 0 ; i < NX ; i++) {
                for(int j = 0 ; j < NY ; j++) {
                    
                    int totalDensity = 0;
                    
                    for (int di = -1; di <= 1; di++) {
                        for (int dj = -1; dj <= 1; dj++) {
                            int ni = (i + di + NX) % NX;
                            int nj = (j + dj + NY) % NY;
                            
                            boolean[] f = fin[ni][nj];
                            if (f[0]) totalDensity++;
                            if (f[1]) totalDensity++;
                            if (f[2]) totalDensity++;
                            if (f[3]) totalDensity++;
                        }
                    }
                    
                    int grayValue = 255 - (totalDensity * 255 / 36);
                    
                    // boundry condition for grayValue
                    grayValue = Math.max(0, Math.min(255, grayValue));
                    
                    g.setColor(new Color(grayValue, grayValue, grayValue));
                    g.fillRect(CELL_SIZE * i, CELL_SIZE * j, CELL_SIZE, CELL_SIZE) ;
                }
            } 
        }
    }

    
    static class Display_average_neighbour_closed extends JPanel {

        final static int CELL_SIZE = 14 ; 

        Display_average_neighbour_closed() {
            setPreferredSize(new Dimension(CELL_SIZE * NX, CELL_SIZE * NY)) ;

            JFrame frame = new JFrame("HPP - 3x3 Average Density in Closed Box");
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.setContentPane(this);
            frame.pack();
            frame.setVisible(true);
        }

        public void paintComponent(Graphics g) {
            g.setColor(Color.WHITE) ;
            g.fillRect(0, 0, CELL_SIZE * NX, CELL_SIZE * NY) ;

            for(int i = 0 ; i < NX ; i++) {
                for(int j = 0 ; j < NY ; j++) {
                    
                    int totalDensity = 0;
                    int validCells = 0; 
                    
                    for (int di = -1; di <= 1; di++) {
                        for (int dj = -1; dj <= 1; dj++) {
                            
                            int ni = i + di;
                            int nj = j + dj;
                            
                            if (ni >= 0 && ni < NX && nj >= 0 && nj < NY) {
                                validCells++; 
                                
                                boolean[] f = fin[ni][nj];
                                if (f[0]) totalDensity++;
                                if (f[1]) totalDensity++;
                                if (f[2]) totalDensity++;
                                if (f[3]) totalDensity++;
                            }
                        }
                    }
                    int maxPossibleDensity = validCells * 4;
                    int grayValue = 255 - (totalDensity * 255 / maxPossibleDensity);
                    // boundry condition for grayValue
                    grayValue = Math.max(0, Math.min(255, grayValue));
                    
                    g.setColor(new Color(grayValue, grayValue, grayValue));
                    g.fillRect(CELL_SIZE * i, CELL_SIZE * j, CELL_SIZE, CELL_SIZE) ;
                }
            } 
        }
    }
}