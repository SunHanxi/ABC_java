package abc;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class readfile {
    public static void main(String[] args) {
        List<String> list = new ArrayList<String>();

        int n = 3;
        double lb = 0;
        double ub = 3;
        Service[][] services = new Service[n][(int) (ub - lb + 1)];
        String encoding = "UTF-8";
        File file = new File("./file/数据集.txt");
        InputStreamReader read;
        {
            try {
                read = new InputStreamReader(
                        new FileInputStream(file), encoding);
                BufferedReader bufferedReader = new BufferedReader(read);
                String lineTxt;

                while ((lineTxt = bufferedReader.readLine()) != null) {
                    list.add(lineTxt);
                }
                bufferedReader.close();
                read.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        List<String> serviceList;
        serviceList = list;
        int group = 0;
        for (int i = 0; i < n * (ub - lb + 1); i++) {
            String[] tmpList = serviceList.get(i).split(" ");
            group = (int) (i / (ub - lb + 1) + 1);
            System.out.println("tmpList: " + Arrays.toString(tmpList));
            System.out.println("group: " + group);
            int[] gg = new int[2];
            gg[0] = group - 1;
            gg[1] = i - (group - 1) * (int) (ub - lb + 1);
            //System.out.println("gg:" + Arrays.toString(gg));
            services[gg[0]][gg[1]] = new Service(gg, Integer.parseInt(tmpList[0]), Integer.parseInt(tmpList[1]),
                    Integer.parseInt(tmpList[2]), Double.parseDouble(tmpList[3]), Double.parseDouble(tmpList[4]));
        }//完成初始化服务商
        for (Service[] s : services) {
            for (Service ss : s) {
                System.out.println(ss);
            }
        }
    }
}
