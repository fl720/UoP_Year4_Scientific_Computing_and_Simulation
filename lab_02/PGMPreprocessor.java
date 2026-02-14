import java.io.*;
import java.util.Scanner;

public class PGMPreprocessor {
    public static void cropAndSave(String inputPath, String outputPath, int targetN) throws Exception {
        ReadPGM.read(X, "barbara.ascii.pgm", N);
    }

    public static void main(String[] args) throws Exception {
        cropAndSave("Square-1080x1080.pgm", "Square-1024.pgm", 1024);
    }
}