import java.awt.*;
import javax.swing.*;

/**
 * Four-state cellular automaton model for excitable media.
 * States: 0 - Resting, 1 - Recovering, 2 - Plateau, 3 - Excited (wave front)
 * Each cell has a timer that determines when to move to the next state.
 */
public class FourStateCA {

    final static int N = 100;                // grid size
    final static int CELL_SIZE = 2;         // pixel size of each cell
    final static int DELAY = 200;           // ms between frames

    // State and timer arrays
    static int[][] state = new int[N][N];
    static int[][] timer = new int[N][N];

    // Flag: does this cell have an excited neighbour? (state 2 or 3)
    static boolean[][] excitedNeighbour = new boolean[N][N];

    static Display display = new Display();

    public static void main(String args[]) throws Exception {

        // ---- Initialisation ----
        // Set bottom row to excited (state 3) with timer = 2, others resting (0)
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (j == N - 1) {
                    state[i][j] = 3;
                    timer[i][j] = 2;
                } else {
                    state[i][j] = 0;
                    timer[i][j] = 0; // not used for resting
                }
            }
        }

        display.repaint();
        pause();

        // Main simulation loop
        int iter = 0;
        while (true) {
            System.out.println("iter = " + iter++);

            // ---- Optional wave breaking (uncomment to create spiral) ----
            if (iter == N / 2) {
                // Clear left half of the grid (columns 0 .. N/2-1)
                for (int i = 0; i < N / 2; i++) {
                    for (int j = 0; j < N; j++) {
                        state[i][j] = 0;
                        timer[i][j] = 0;
                    }
                }
            }

            // ---- Compute which cells have at least one excited neighbour (state 2 or 3) in 8-neighbourhood ----
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    excitedNeighbour[i][j] = false;
                    // Check all 8 neighbours
                    for (int di = -1; di <= 1; di++) {
                        for (int dj = -1; dj <= 1; dj++) {
                            if (di == 0 && dj == 0) continue; // skip self
                            int ni = i + di;
                            int nj = j + dj;
                            if (ni >= 0 && ni < N && nj >= 0 && nj < N) {
                                int s = state[ni][nj];
                                if (s == 2 || s == 3) {
                                    excitedNeighbour[i][j] = true;
                                    break;
                                }
                            }
                        }
                        if (excitedNeighbour[i][j]) break;
                    }
                }
            }

            // ---- Update each cell's state and timer ----
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    switch (state[i][j]) {
                        case 0: // Resting
                            if (excitedNeighbour[i][j]) {
                                state[i][j] = 3;       // become excited
                                timer[i][j] = 2;       // will stay in state 3 for 2 steps
                            }
                            break;

                        case 3: // Excited (wave front)
                            timer[i][j]--;
                            if (timer[i][j] == 0) {
                                state[i][j] = 2;       // move to plateau
                                timer[i][j] = 3;       // stay in plateau for 3 steps
                            }
                            break;

                        case 2: // Plateau
                            timer[i][j]--;
                            if (timer[i][j] == 0) {
                                state[i][j] = 1;       // move to recovering
                                timer[i][j] = 3;       // stay in recovering for 3 steps
                            }
                            break;

                        case 1: // Recovering
                            timer[i][j]--;
                            if (timer[i][j] == 0) {
                                state[i][j] = 0;       // back to rest
                                // timer not used for rest
                            }
                            break;
                    }
                }
            }

            display.repaint();
            pause();
        }
    }

    // ---- Display panel ----
    static class Display extends JPanel {
        final static int WINDOW_SIZE = N * CELL_SIZE;

        Display() {
            setPreferredSize(new Dimension(WINDOW_SIZE, WINDOW_SIZE));
            JFrame frame = new JFrame("Four-State Excitable Media Model");
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.setContentPane(this);
            frame.pack();
            frame.setVisible(true);
        }

        public void paintComponent(Graphics g) {
            g.setColor(Color.WHITE);
            g.fillRect(0, 0, WINDOW_SIZE, WINDOW_SIZE);
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    if (state[i][j] != 0) {
                        switch (state[i][j]) {
                            case 3: g.setColor(Color.BLACK);       break; // excited
                            case 2: g.setColor(Color.DARK_GRAY);   break; // plateau
                            case 1: g.setColor(Color.GRAY);        break; // recovering
                        }
                        g.fillRect(CELL_SIZE * i, CELL_SIZE * j,
                                   CELL_SIZE, CELL_SIZE);
                    }
                }
            }
        }
    }

    static void pause() {
        try {
            Thread.sleep(DELAY);
        } catch (InterruptedException e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
}