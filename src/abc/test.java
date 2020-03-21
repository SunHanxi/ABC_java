package abc;

import java.util.Arrays;
import java.util.Random;

public class test {
    public static void main(String[] args) {
        int[] Foods = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
        int n = 4, num_of_plans = 3;
        int[] sol_concurrency = new int[n];
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < num_of_plans; k++) {
                sol_concurrency[j] += Foods[j + n * k];
            }
        }
        System.out.println(Arrays.toString(sol_concurrency));
    }
}
