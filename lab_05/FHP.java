import java.awt.* ;
import javax.swing.* ;

public class FHP {

    final static int NX = 80, NY = 60 ;  // Lattice dimensions
    final static int q = 6 ;  // population
     
    final static int NITER = 1000 ;
    final static int DELAY = 50 ; // Reduced delay for faster visualization, feel free to adjust

    final static double DENSITY = 0.5 ;  // initial state, between 0 and 1.0.

    // Toggle the displays below just like in HPP.java
    // static Display display = new Display() ;
    static Display_grey display = new Display_grey() ;
    // static Display_average_neighbour display = new Display_average_neighbour() ;
    // static Display_average_neighbour_closed display = new Display_average_neighbour_closed() ;

    static boolean [] [] [] fin = new boolean [NX] [NY] [q] ;
    static boolean [] [] [] fout = new boolean [NX] [NY] [q] ;

    public static void main(String args []) throws Exception {

        // initialize - populate a subblock of grid
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
                    int pop = 0 ;
                    for(int d = 0 ; d < q ; d++) {
                        if(fin_ij [d]) pop++ ;
                        fout_ij [d] = fin_ij [d] ;
                    }
                    if(pop == 2) {
                        // head on collisions
                        if(fin_ij [0] && fin_ij [1]) {
                            fout_ij [0] = false ;
                            fout_ij [1] = false ;
                            if(Math.random() < 0.5) {
                                fout_ij [2] = true ;
                                fout_ij [3] = true ;
                            }
                            else {
                                fout_ij [4] = true ;
                                fout_ij [5] = true ;
                            }
                        }    
                        if(fin_ij [2] && fin_ij [3]) {
                            fout_ij [2] = false ;
                            fout_ij [3] = false ;
                            if(Math.random() < 0.5) {
                                fout_ij [4] = true ;
                                fout_ij [5] = true ;
                            }
                            else {
                                fout_ij [0] = true ;
                                fout_ij [1] = true ;
                            }
                        }    
                        if(fin_ij [4] && fin_ij [5]) {
                            fout_ij [4] = false ;
                            fout_ij [5] = false ;
                            if(Math.random() < 0.5) {
                                fout_ij [0] = true ;
                                fout_ij [1] = true ;
                            }
                            else {
                                fout_ij [2] = true ;
                                fout_ij [3] = true ;
                            }
                        }    
                    }
                    else if(pop == 4) {
                        // double head on collisions
                        if(!fin_ij [0] && !fin_ij [1]) {
                            fout_ij [0] = true ;
                            fout_ij [1] = true ;
                            if(Math.random() < 0.5) {
                                fout_ij [2] = false ;
                                fout_ij [3] = false ;
                            }
                            else {
                                fout_ij [4] = false ;
                                fout_ij [5] = false ;
                            }
                        }    
                        if(!fin_ij [2] && !fin_ij [3]) {
                            fout_ij [2] = true ;
                            fout_ij [3] = true ;
                            if(Math.random() < 0.5) {
                                fout_ij [4] = false ;
                                fout_ij [5] = false ;
                            }
                            else {
                                fout_ij [0] = false ;
                                fout_ij [1] = false ;
                            }
                        }    
                        if(!fin_ij [4] && !fin_ij [5]) {
                            fout_ij [4] = true ;
                            fout_ij [5] = true ;
                            if(Math.random() < 0.5) {
                                fout_ij [0] = false ;
                                fout_ij [1] = false ;
                            }
                            else {
                                fout_ij [2] = false ;
                                fout_ij [3] = false ;
                            }
                        }    
                    }
                    else if(pop == 3) {
                        // three way collisions
                        if(fin_ij [0] && fin_ij [3] && fin_ij [4]) {
                            fout_ij [0] = false ;
                            fout_ij [1] = true ;
                            fout_ij [2] = true ;
                            fout_ij [3] = false ;
                            fout_ij [4] = false ;
                            fout_ij [5] = true ;
                        }
                        if(fin_ij [1] && fin_ij [2] && fin_ij [5]) {
                            fout_ij [0] = true ;
                            fout_ij [1] = false ;
                            fout_ij [2] = false ;
                            fout_ij [3] = true ;
                            fout_ij [4] = true ;
                            fout_ij [5] = false ;
                        }
                        // head on with spectator
                        if(fin_ij [0] && fin_ij [1]) {
                            fin_ij [0] = false ;
                            fin_ij [1] = false ;
                            if(fin_ij [2] || fin_ij [3]) {
                                fin_ij [4] = true ;
                                fin_ij [5] = true ;
                            }
                            else {
                                fin_ij [2] = true ;
                                fin_ij [3] = true ;
                            }
                        }
                        if(fin_ij [2] && fin_ij [3]) {
                            fin_ij [2] = false ;
                            fin_ij [3] = false ;
                            if(fin_ij [4] || fin_ij [5]) {
                                fin_ij [0] = true ;
                                fin_ij [1] = true ;
                            }
                            else {
                                fin_ij [4] = true ;
                                fin_ij [5] = true ;
                            }
                        }
                        if(fin_ij [4] && fin_ij [5]) {
                            fin_ij [4] = false ;
                            fin_ij [5] = false ;
                            if(fin_ij [0] || fin_ij [1]) {
                                fin_ij [2] = true ;
                                fin_ij [3] = true ;
                            }
                            else {
                                fin_ij [0] = true ;
                                fin_ij [1] = true ;
                            }
                        }
                    }
                }
            }
            // Streaming
            for(int i = 0; i < NX ; i++) { 
                for(int j = 0; j < NY ; j++) { 
                    
                    // ----- cylindrical boundary conditions: -----
                    int iP1 = (i + 1) % NX ;
                    int iM1 = (i - 1 + NX) % NX ;
                    int jP1 = (j + 1) % NY ;
                    int jM1 = (j - 1 + NY) % NY ;

                    // Push the particles to their new cells
                    fin [iM1] [j] [0] = fout [i] [j] [0] ;
                    fin [iP1] [j] [1] = fout [i] [j] [1] ;
                    fin [i] [jM1] [2] = fout [i] [j] [2] ;
                    fin [i] [jP1] [3] = fout [i] [j] [3] ;
                    fin [iP1] [jM1] [4] = fout [i] [j] [4] ;
                    fin [iM1] [jP1] [5] = fout [i] [j] [5] ;
                    

                    // ----- Streaming (Closed Box Boundary / Bounce-back) -----
                    // if (i < NX - 1) {
                    //     fin[i][j][0] = fout[i+1][j][0];
                    // } else {
                    //     fin[i][j][0] = fout[i][j][1];
                    // }
                    // if (i > 0) {
                    //     fin[i][j][1] = fout[i-1][j][1];
                    // } else {
                    //     fin[i][j][1] = fout[i][j][0];
                    // }
                    // if (j < NY - 1) {
                    //     fin[i][j][2] = fout[i][j+1][2];
                    // } else {
                    //     fin[i][j][2] = fout[i][j][3];
                    // }
                    // if (j > 0) {
                    //     fin[i][j][3] = fout[i][j-1][3];
                    // } else {
                    //     fin[i][j][3] = fout[i][j][2];
                    // }
                    // if (i > 0 && j < NY - 1) {
                    //     fin[i][j][4] = fout[i-1][j+1][4];
                    // } else {
                    //     fin[i][j][4] = fout[i][j][5];
                    // }
                    // if (i < NX - 1 && j > 0) {
                    //     fin[i][j][5] = fout[i+1][j-1][5];
                    // } else {
                    //     fin[i][j][5] = fout[i][j][4];
                    // }
                
                }
            }

            System.out.println("iter = " + iter) ;

            display.repaint() ;
            Thread.sleep(DELAY) ;
        }
    }

    
    static class Display extends JPanel {

        static final double ROW_HEIGHT = Math.sqrt(3) / 2 ;

        int displaySizeX = CELL_SIZE * NX ;
        int displaySizeY = (int) (ROW_HEIGHT * CELL_SIZE * NY + 0.5) ;

        final static int CELL_SIZE = 14 ; 

        public static final int ARROW_START = 2 ;
        public static final int ARROW_END   = 7 ;
        public static final int ARROW_WIDE  = 3 ;

        public static final int DIAG_X_0 = -1 ;
        public static final int DIAG_X_1 = 3 ;
        public static final int DIAG_X_2 = 4 ;
        public static final int DIAG_Y_0 = 4 ;
        public static final int DIAG_Y_1 = 0 ;
        public static final int DIAG_Y_2 = 6 ;

        Display() {

            setPreferredSize(new Dimension(displaySizeX, displaySizeY)) ;

            JFrame frame = new JFrame("FHP");
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.setContentPane(this);
            frame.pack();
            frame.setVisible(true);
        }

        public void paintComponent(Graphics g) {

            g.setColor(Color.WHITE) ;
            g.fillRect(0, 0, CELL_SIZE * NX, CELL_SIZE * NY) ;

            int [] tri_x = new int [3], tri_y = new int [3] ;
            for(int i = 0 ; i < NX ; i++) {
                for(int j = 0 ; j < NY ; j++) {
                    boolean [] f_ij = fin [i] [j] ;
                    int originX = 
                      (int) (CELL_SIZE * i + 0.5 * CELL_SIZE * j +
                             CELL_SIZE/2 + 0.5) % displaySizeX;
                    int originY =
                      (int) (ROW_HEIGHT * (CELL_SIZE * j + CELL_SIZE/2) + 0.5) ;

                    g.setColor(Color.PINK) ;
                    g.fillOval(originX - 2, originY - 2, 4, 4) ;

                    g.setColor(Color.BLUE) ;
                    if(f_ij [0]) {
                        tri_x [0] = originX - ARROW_START ; tri_x [1] = originX - ARROW_START ; tri_x [2] = originX - ARROW_END ;
                        tri_y [0] = originY - ARROW_WIDE ;  tri_y [1] = originY + ARROW_WIDE ;  tri_y [2] = originY ;
                        g.fillPolygon(tri_x, tri_y, 3) ;
                    }
                    if(f_ij [1]) {
                        tri_x [0] = originX + ARROW_START ; tri_x [1] = originX + ARROW_START ; tri_x [2] = originX + ARROW_END ;
                        tri_y [0] = originY - ARROW_WIDE ;  tri_y [1] = originY + ARROW_WIDE ;  tri_y [2] = originY ;
                        g.fillPolygon(tri_x, tri_y, 3) ;
                    }
                    if(f_ij [2]) {
                        tri_x [0] = originX - DIAG_X_0 ; tri_x [1] = originX - DIAG_X_1 ; tri_x [2] = originX - DIAG_X_2 ;
                        tri_y [0] = originY - DIAG_Y_0 ; tri_y [1] = originY - DIAG_Y_1 ; tri_y [2] = originY - DIAG_Y_2 ;
                        g.fillPolygon(tri_x, tri_y, 3) ;
                    }
                    if(f_ij [3]) {
                        tri_x [0] = originX + DIAG_X_0 ; tri_x [1] = originX + DIAG_X_1 ; tri_x [2] = originX + DIAG_X_2 ;
                        tri_y [0] = originY + DIAG_Y_0 ; tri_y [1] = originY + DIAG_Y_1 ; tri_y [2] = originY + DIAG_Y_2 ;
                        g.fillPolygon(tri_x, tri_y, 3) ;
                    }
                    if(f_ij [4]) {
                        tri_x [0] = originX + DIAG_X_0 ; tri_x [1] = originX + DIAG_X_1 ; tri_x [2] = originX + DIAG_X_2 ;
                        tri_y [0] = originY - DIAG_Y_0 ; tri_y [1] = originY - DIAG_Y_1 ; tri_y [2] = originY - DIAG_Y_2 ;
                        g.fillPolygon(tri_x, tri_y, 3) ;
                    }
                    if(f_ij [5]) {
                        tri_x [0] = originX - DIAG_X_0 ; tri_x [1] = originX - DIAG_X_1 ; tri_x [2] = originX - DIAG_X_2 ;
                        tri_y [0] = originY + DIAG_Y_0 ; tri_y [1] = originY + DIAG_Y_1 ; tri_y [2] = originY + DIAG_Y_2 ;
                        g.fillPolygon(tri_x, tri_y, 3) ;
                    }
                }
            } 
        }
    }

    static class Display_grey extends JPanel {

        static final double ROW_HEIGHT = Math.sqrt(3) / 2 ;
        final static int CELL_SIZE = 14 ; 
        
        int displaySizeX = CELL_SIZE * NX ;
        int displaySizeY = (int) (ROW_HEIGHT * CELL_SIZE * NY + 0.5) ;

        Display_grey() {
            setPreferredSize(new Dimension(displaySizeX, displaySizeY)) ;

            JFrame frame = new JFrame("FHP - Density Map");
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
                    for (int d = 0; d < q; d++) {
                        if (fin_ij[d]) density++;
                    }
                    
                    // density = 0 -> grayValue = 255 (white)
                    // density = 6 -> grayValue = 0 (black)
                    int grayValue = 255 - (density * 255 / q);
                    g.setColor(new Color(grayValue, grayValue, grayValue));
                    
                    int originX = (int) (CELL_SIZE * i + 0.5 * CELL_SIZE * j + 0.5) % displaySizeX;
                    int originY = (int) (ROW_HEIGHT * (CELL_SIZE * j) + 0.5) ;
                    g.fillRect(originX, originY, CELL_SIZE, CELL_SIZE) ;
                }
            } 
        }
    }

    static class Display_average_neighbour extends JPanel {

        static final double ROW_HEIGHT = Math.sqrt(3) / 2 ;
        final static int CELL_SIZE = 14 ; 
        
        int displaySizeX = CELL_SIZE * NX ;
        int displaySizeY = (int) (ROW_HEIGHT * CELL_SIZE * NY + 0.5) ;

        Display_average_neighbour() {
            setPreferredSize(new Dimension(displaySizeX, displaySizeY)) ;

            JFrame frame = new JFrame("FHP - 3x3 Average Density");
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.setContentPane(this);
            frame.pack();
            frame.setVisible(true);
        }

        public void paintComponent(Graphics g) {
            
            g.setColor(Color.WHITE) ;
            g.fillRect(0, 0, displaySizeX, displaySizeY) ;

            for(int i = 0 ; i < NX ; i++) {
                for(int j = 0 ; j < NY ; j++) {
                    
                    int totalDensity = 0;
                    
                    for (int di = -1; di <= 1; di++) {
                        for (int dj = -1; dj <= 1; dj++) {
                            int ni = (i + di + NX) % NX;
                            int nj = (j + dj + NY) % NY;
                            
                            boolean[] f = fin[ni][nj];
                            for(int d=0; d<q; d++) {
                                if (f[d]) totalDensity++;
                            }
                        }
                    }
                    
                    // Max density for a 3x3 block in FHP is 9 cells * 6 directions = 54
                    int grayValue = 255 - (totalDensity * 255 / (9 * q));
                    grayValue = Math.max(0, Math.min(255, grayValue));
                    
                    g.setColor(new Color(grayValue, grayValue, grayValue));
                    int originX = (int) (CELL_SIZE * i + 0.5 * CELL_SIZE * j + 0.5) % displaySizeX;
                    int originY = (int) (ROW_HEIGHT * (CELL_SIZE * j) + 0.5) ;
                    g.fillRect(originX, originY, CELL_SIZE, CELL_SIZE) ;
                }
            } 
        }
    }

    static class Display_average_neighbour_closed extends JPanel {

        static final double ROW_HEIGHT = Math.sqrt(3) / 2 ;
        final static int CELL_SIZE = 14 ; 
        
        int displaySizeX = CELL_SIZE * NX ;
        int displaySizeY = (int) (ROW_HEIGHT * CELL_SIZE * NY + 0.5) ;

        Display_average_neighbour_closed() {
            setPreferredSize(new Dimension(displaySizeX, displaySizeY)) ;

            JFrame frame = new JFrame("FHP - 3x3 Average Density in Closed Box");
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.setContentPane(this);
            frame.pack();
            frame.setVisible(true);
        }

        public void paintComponent(Graphics g) {
            g.setColor(Color.WHITE) ;
            g.fillRect(0, 0, displaySizeX, displaySizeY) ;

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
                                for(int d=0; d<q; d++) {
                                    if (f[d]) totalDensity++;
                                }
                            }
                        }
                    }
                    int maxPossibleDensity = validCells * q;
                    int grayValue = 255 - (totalDensity * 255 / maxPossibleDensity);
                    grayValue = Math.max(0, Math.min(255, grayValue));
                    
                    g.setColor(new Color(grayValue, grayValue, grayValue));
                    int originX = (int) (CELL_SIZE * i + 0.5 * CELL_SIZE * j + 0.5) % displaySizeX;
                    int originY = (int) (ROW_HEIGHT * (CELL_SIZE * j) + 0.5) ;
                    g.fillRect(originX, originY, CELL_SIZE, CELL_SIZE) ;
                }
            } 
        }
    }
}