package abc;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class run_abc {

    public static void main(String[] args) {
        long this_env_start = System.currentTimeMillis();
        //设定单条解最长用时为600秒
        //改进了惩罚函数后，这个最长用时可能要被砍掉了。
        //double max_run_time = 600;

        //TODO 需要加一个读取配置文件的功能
        //TODO 貌似不需要了，因为配置文件就在数据集文件名里
        //获取所有数据集
        String folder_path = "./file";
        File file = new File(folder_path);
        File[] fileList = file.listFiles();
        assert fileList != null;
        List<String> fileNames = new ArrayList<>();
        for (File value : fileList) {
            if (value.getName().contains("txt")) {
                fileNames.add(value.getName());
            }
        }
        System.out.println("所有数据集：");
        for (String s : fileNames) {
            System.out.println(s);
        }
        //用来终止实验的语句，如果只是为了输出数据集的名字，即可从这儿打断实验。
        if (1 != 1) {
            return;
        }
        //新建文件流
        File f = new File("输出结果.txt");
        try {
            FileOutputStream result_file = new FileOutputStream(f);
            //输出文件采用utf-8编码
            OutputStreamWriter writer = new OutputStreamWriter(result_file, StandardCharsets.UTF_8);
            String dataset_path;
            String folder = "./file/";
            for (String fileName : fileNames) {
                dataset_path = fileName;
                System.out.println("\n本轮数据集为：" + dataset_path);
                // 初始化参数
                runFunction(writer, folder, dataset_path);
                //writer.append("\n");
            }

            writer.close();
            result_file.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
        long this_env_end = System.currentTimeMillis();
        System.out.println("此次实验用时：" + (this_env_end - this_env_start) + " ms");
    }

    private static double run_ABC(vm_abc serBee) {
        long start_time = System.currentTimeMillis();
        serBee.initial();
        serBee.MemorizeBestSource();
        for (int iter = 0; iter < serBee.maxCycle; iter++) {
            serBee.SendEmployedBees();
            serBee.CalculateProbabilities();
            serBee.SendOnlookerBees();
            serBee.MemorizeBestSource();
            serBee.SendScoutBees();
            // 每2000轮打印一下代价
//            if (iter % 2000 == 0) {
//                System.out.println(serBee.GlobalMin);
//            }
        }
        long end_time = System.currentTimeMillis();
        return (double) (end_time - start_time) / 1000;
    }

    private static void runFunction(OutputStreamWriter writer, String folder,
                                    String dataset_path) {
        //最多循环的次数
        int maxCycle = 500000;

        int n;  // 服务商的类别的数量
        int number_of_service;
        double time_want_spent; //约束条件,因为时间是5-20随机生成，取中位数12.5,共4个
        double lb;  // 随机数的上下界，此处为每个服务的数量
        double ub;  //切记，如果是10个备选服务数量，那么lb是0，ub是9
        double[] max_time = new double[2];  //随机出来的时间的范围
        int num_of_plans;
        int[] total_concurrency; //总的并发量
        //从文件名中提取参数
        String[] param = dataset_path.split("\\.txt")[0].split("_");
        //System.out.println(Arrays.toString(param));
        n = Integer.parseInt(param[0]);
        num_of_plans = Integer.parseInt(param[1]);
        lb = Integer.parseInt(param[2]);
        ub = Integer.parseInt(param[3]);
        time_want_spent = Double.parseDouble(param[4]);
        total_concurrency = new int[n];
        for (int i = 5, j = 0; i < param.length; i++) {
            total_concurrency[j] = Integer.parseInt(param[i]);
            j++;
        }

        try {
            writer.append("");
            System.out.println("控制参数为：");
            System.out.println("n: " + n + "\t路径条数：" + num_of_plans + "\t时间要求：" + time_want_spent
                    + "\t并发量：" + Arrays.toString(total_concurrency));


            int runtime = 1;  /*算法重复运行的次数，重复多次可以查看算法的稳健性*/
            //double[] result = new double[runtime];


            for (int run = 0; run < runtime; run++) {
                vm_abc serBee = new vm_abc(time_want_spent, n, num_of_plans, lb, ub,
                        folder + dataset_path, maxCycle, total_concurrency);
                double current_run_time;
                current_run_time = run_ABC(serBee);

                //输出最佳解
                int j = 0;

                //按路径依次输出每条路径的解
                //求解每条路径的cost和price

                //输出最优解
                for (int i = 0; i < serBee.GlobalParams.length; i++) {
                    System.out.println(serBee.GlobalParams[i]);
                }
                double[] total_price = new double[num_of_plans];
                double[] total_time = new double[num_of_plans];
                for (int i = 0; i < serBee.GlobalParams.length; i++) {
                    total_price[i / n] += serBee.GlobalParams[i].price;
                    total_time[i / n] += serBee.GlobalParams[i].time;
                }
                for (int i = 0; i < total_price.length; i++) {
                    System.out.println("本条路径价格： " + total_price[i] + "  " + "本条路径用时： " + total_time[i]);
                }

                System.out.print("三条路径概率分别为： ");
                for (int i = 0; i < serBee.GlobalProbability.length; i++) {
                    System.out.print(serBee.GlobalProbability[i] + "\t");
                }
                System.out.println();

                //计算解的代价和实际代价是否相同
                double real_price = Arrays.stream(total_price).sum();
                String price_flag;
                if (real_price == serBee.GlobalMin) {
                    price_flag = "=";
                } else {
                    price_flag = "≠";
                }
                System.out.println("本条解的实际代价为： " + serBee.GlobalMin + " " + price_flag + " 实际代价： " + real_price);
                System.out.println("用时：" + current_run_time + "秒");
                //result[run] = serBee.GlobalMin;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
