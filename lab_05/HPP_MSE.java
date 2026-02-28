// MULTI-SPIN EDCODING HPP 
// this method uses int to represent the state of each cell, with each bit representing the presence of a particle in a particular direction.
// Using bitwise operations allows for efficient storage and manipulation of the cell states, enabling faster computations and reduced memory usage compared to using separate boolean arrays for each direction.
import java.awt.* ;
import javax.swing.* ;

public class HPP_MSE {

    final static int NX = 80, NY = 60; 
    final static int NITER = 10000;
    final static int DELAY = 20;
    final static double DENSITY = 0.5;

    final static int LEFT  = 1; // 2^0
    final static int RIGHT = 2; // 2^1
    final static int UP    = 4; // 2^2
    final static int DOWN  = 8; // 2^3

    static Display display = new Display();
    // static Display_grey display = new Display_grey();
    // static Display_average_neighbour display = new Display_average_neighbour();
    // static Display_average_neighbour_closed display = new Display_average_neighbour_closed();

    static int[][] fin = new int[NX][NY];
    static int[][] fout = new int[NX][NY];

    public static void main(String args[]) throws Exception {

        for (int i = 0; i < NX / 4; i++) {
            for (int j = 0; j < NY / 4; j++) {
                if (Math.random() < DENSITY) fin[i][j] |= LEFT;
                if (Math.random() < DENSITY) fin[i][j] |= RIGHT;
                if (Math.random() < DENSITY) fin[i][j] |= UP;
                if (Math.random() < DENSITY) fin[i][j] |= DOWN;
            }
        }

        display.repaint();
        Thread.sleep(DELAY);

        for (int iter = 0; iter < NITER; iter++) {

            // --- collision --- 
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    int s = fin[i][j];
                    if (s == (LEFT | RIGHT)) {
                        fout[i][j] = (UP | DOWN);
                    } else if (s == (UP | DOWN)) {
                        fout[i][j] = (LEFT | RIGHT);
                    } else {
                        fout[i][j] = s; 
                    }
                }
            }

            // --- Streaming  ---
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) fin[i][j] = 0;
            }

            for (int i = 0; i < NX; i++) {
                int iP1 = (i + 1) % NX;
                int iM1 = (i - 1 + NX) % NX;
                for (int j = 0; j < NY; j++) {
                    int jP1 = (j + 1) % NY;
                    int jM1 = (j - 1 + NY) % NY;

                    int s = fout[i][j];

                    // ------- closed box-like boundary conditions: -------
                    // if ((s & LEFT) != 0) {
                    //     if (i != 0) fin[iM1][j] |= LEFT;
                    //     else fin[i][j] |= RIGHT; 
                    // }
                    // if ((s & RIGHT) != 0) {
                    //     if (i != NX - 1) fin[iP1][j] |= RIGHT;
                    //     else fin[i][j] |= LEFT;
                    // }
                    // if ((s & UP) != 0) {
                    //     if (j != 0) fin[i][jM1] |= UP;
                    //     else fin[i][j] |= DOWN; 
                    // }
                    // if ((s & DOWN) != 0) {
                    //     if (j != NY - 1) fin[i][jP1] |= DOWN;
                    //     else fin[i][j] |= UP; 
                    // }

                    // ------- cylindrical boundary conditions: -------
                    if ((s & LEFT)  != 0) fin[iM1][j] |= LEFT;
                    if ((s & RIGHT) != 0) fin[iP1][j] |= RIGHT;
                    if ((s & UP)    != 0) fin[i][jM1] |= UP;
                    if ((s & DOWN)  != 0) fin[i][jP1] |= DOWN;
                }
            }

            if (iter % 10 == 0) System.out.println("iter = " + iter);
            display.repaint();
            Thread.sleep(DELAY);
        }
    }

    static class Display extends JPanel {
        final static int CELL_SIZE = 14;
        public static final int ARROW_START = 2, ARROW_END = 7, ARROW_WIDE = 3;

        Display() {
            setPreferredSize(new Dimension(CELL_SIZE * NX, CELL_SIZE * NY)) ;
            JFrame frame = new JFrame("HPP_MSE");
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.setContentPane(this);
            frame.pack();
            frame.setVisible(true);
        }

        public void paintComponent(Graphics g) {
            g.setColor(Color.WHITE);
            g.fillRect(0, 0, CELL_SIZE * NX, CELL_SIZE * NY) ;
            g.setColor(Color.PINK) ;
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    int s = fin[i][j];
                    int ox = CELL_SIZE * i + CELL_SIZE / 2;
                    int oy = CELL_SIZE * j + CELL_SIZE / 2;
                    if ((s & LEFT) != 0)  drawArrow(g, ox, oy, -1, 0);
                    if ((s & RIGHT) != 0) drawArrow(g, ox, oy, 1, 0);
                    if ((s & UP) != 0)    drawArrow(g, ox, oy, 0, -1);
                    if ((s & DOWN) != 0)  drawArrow(g, ox, oy, 0, 1);
                }
            }
        }

        private void drawArrow(Graphics g, int x, int y, int dx, int dy) {
            int[] tx = new int[3], ty = new int[3];
            g.setColor(Color.BLUE);
            if (dx != 0) {
                tx[0] = x + dx * ARROW_START; tx[1] = x + dx * ARROW_START; tx[2] = x + dx * ARROW_END;
                ty[0] = y - ARROW_WIDE;       ty[1] = y + ARROW_WIDE;       ty[2] = y;
            } else {
                tx[0] = x - ARROW_WIDE;       tx[1] = x + ARROW_WIDE;       tx[2] = x;
                ty[0] = y + dy * ARROW_START; ty[1] = y + dy * ARROW_START; ty[2] = y + dy * ARROW_END;
            }
            g.fillPolygon(tx, ty, 3);
        }
    }

    static class Display_grey extends JPanel {
        final static int CELL_SIZE = 14 ; 
        Display_grey() {
            setPreferredSize(new Dimension(CELL_SIZE * NX, CELL_SIZE * NY)) ;

            JFrame frame = new JFrame("HPP_MSE - Density Map");
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.setContentPane(this);
            frame.pack();
            frame.setVisible(true);
        }
        public void paintComponent(Graphics g) {
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    int density = Integer.bitCount(fin[i][j]);
                    int gray = 255 - (density * 255 / 4);
                    g.setColor(new Color(gray, gray, gray));
                    g.fillRect(i * 14, j * 14, 14, 14);
                }
            }
        }
    }

    static class Display_average_neighbour extends JPanel {
        final static int CELL_SIZE = 14 ; 
        Display_average_neighbour() {
            setPreferredSize(new Dimension(CELL_SIZE * NX, CELL_SIZE * NY)) ;

            JFrame frame = new JFrame("HPP_MSE - 3x3 Average Density");
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.setContentPane(this);
            frame.pack();
            frame.setVisible(true);
        }
        public void paintComponent(Graphics g) {
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    int total = 0;
                    for (int di = -1; di <= 1; di++) {
                        for (int dj = -1; dj <= 1; dj++) {
                            total += Integer.bitCount(fin[(i + di + NX) % NX][(j + dj + NY) % NY]);
                        }
                    }
                    int gray = Math.max(0, 255 - (total * 255 / 36));
                    g.setColor(new Color(gray, gray, gray));
                    g.fillRect(i * 14, j * 14, 14, 14);
                }
            }
        }
    }

    static class Display_average_neighbour_closed extends JPanel {
        
        final static int CELL_SIZE = 14 ; 
        Display_average_neighbour_closed() {
            setPreferredSize(new Dimension(CELL_SIZE * NX, CELL_SIZE * NY)) ;

            JFrame frame = new JFrame("HPP_MSE - 3x3 Average Density");
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.setContentPane(this);
            frame.pack();
            frame.setVisible(true);
        }
        public void paintComponent(Graphics g) {
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    int total = 0, cells = 0;
                    for (int di = -1; di <= 1; di++) {
                        for (int dj = -1; dj <= 1; dj++) {
                            int ni = i + di, nj = j + dj;
                            if (ni >= 0 && ni < NX && nj >= 0 && nj < NY) {
                                total += Integer.bitCount(fin[ni][nj]);
                                cells++;
                            }
                        }
                    }
                    int gray = Math.max(0, 255 - (total * 255 / (cells * 4)));
                    g.setColor(new Color(gray, gray, gray));
                    g.fillRect(i * 14, j * 14, 14, 14);
                }
            }
        }
    }

    private static void setupFrame(String title) {
        // 
    }
}