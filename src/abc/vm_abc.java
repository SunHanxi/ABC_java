package abc;

import java.io.*;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

class vm_abc {

    /* ABC的控制参数*/

    //不依赖于初始化的变量
    private boolean random_pro_flag = false;//是否进行等概率实验的标志位
    private int PenaltyRate = 50; //对不符合要求的虚拟机进行惩罚的倍率
    private int TimePenaltyRate = 300; //对不符合要求的虚拟机进行时间惩罚的倍率


    int maxCycle; /*实验的轮数*/
    private int NP = 50; // 种群的规模
    private int FoodNumber = NP / 2; //蜜源数量
    private double ObjValSol; //新解决方案的目标函数值
    private double FitnessSol; //新解决方案的适应度
    private double[][] probability; // 解决方案中每条路径的概率
    /*
    param2change对应于j，
    neighbour对应于等式v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) 中的i
    */
    private int neighbour, param2change;
    double GlobalMin; //ABC算法获得的最优解
    private double r; //在[0,1)范围内的随机数

    //需要初始化的参数
    private double time_want_spent; //约束条件
    private int n;  // 提供服务的虚拟机组的数量
    private double lb;  // 随机数的上下界，此处为每组虚拟机的数量
    private double ub;  //切记，如果是10台虚拟机，那么lb是0，ub是9
    private String dataset_path; //数据集路径
    private int num_of_plans;  //备选方案的条数
    private double[][] concurrency; //存储每条解的每个路径的并发量,并发量由probability计算得到。
    private int[] total_concurrency; //所有路径总的的并发量
    //依赖于初始化的参数
    private int total;// 所有节点的数量
    int D; /*要优化的问题的参数数量*/
    private Service[][] services; // 所有服务商
    /*Foods是蜜源。 Foods矩阵的每一行都是一个包含要优化的D参数的向量。
    Foods矩阵的行数等于FoodNumber*/
    private Service[][] Foods; //直接让Foods存储service对象
    private double[] f; //f是目标函数值
    private double[] fitness; //"fitness适应度
    private double[] trial;  //trail是每个蜜源的试验次数
    private double[] prob; //prob是一个保持蜜源(解决方案)概率的载体,也即轮盘赌的概率
    private Service[] solution;
    Service[] GlobalParams; //最优解的参数,直接存储service对象
    double[] GlobalProbability; // 最优解的每条路径的概率

    private long start_this_time = System.currentTimeMillis();

    //构造函数中对变量进行初始化
    vm_abc(double time_want_spent, int n, int num_of_plans, double lb, double ub,
           String dataset_path, int maxCycle, int[] total_concurrency) {
        this.time_want_spent = time_want_spent;
        this.num_of_plans = num_of_plans;
        this.n = n;
        this.lb = lb;
        this.ub = ub;
        this.maxCycle = maxCycle;
        this.total = n * num_of_plans;
        this.D = this.total;
        //this.concurrency = new double[FoodNumber][this.total];
        this.services = new Service[n][(int) (ub - lb + 1)];
        this.Foods = new Service[FoodNumber][D];
        this.f = new double[FoodNumber];
        this.fitness = new double[FoodNumber];
        this.trial = new double[FoodNumber];
        this.prob = new double[FoodNumber];
        this.solution = new Service[D];
        this.GlobalParams = new Service[D];
        this.GlobalProbability = new double[num_of_plans];
        this.dataset_path = dataset_path;
        this.total_concurrency = total_concurrency;
        this.probability = new double[FoodNumber][this.num_of_plans];
        this.concurrency = new double[FoodNumber][this.total];
    }

    /*存储最佳蜜源*/
    void MemorizeBestSource() {
        int i, j;
        for (i = 0; i < FoodNumber; i++) {
            if (f[i] < GlobalMin) {
                //存储最小的函数值
                GlobalMin = f[i];
                //存储最佳的一条解
                for (j = 0; j < D; j++) {
                    GlobalParams[j] = Foods[i][j];
                }
                //存储最优解的概率
                for (int k = 0; k < num_of_plans; k++) {
                    GlobalProbability[k] = probability[i][k];
                }
            }
        }
    }

    /*变量在[lb，ub]范围内初始化。 如果每个参数具有不同的范围，则使用数组lb [j]，ub [j]而不是lb和ub*/
    /* 蜜源的计数器也在此功能中初始化*/

    private void init(int index) {
        /*
         * 初始化一个蜜源，一个蜜源有D个参数
         *  Math.random() 生成一个[0, 1)的随机数
         *  一个蜜源的结构为：
         * [1,2,3, ... ,total]
         * 每个维度所属的类别在group_map中可以查到
         *虚拟机列表存在services[][]里面
         *
         * */
        //初始化probability[index],也即number_of_plans个概率
        Random random = new Random();
        //System.out.print("放缩前：");
        for (int i = 0; i < num_of_plans; i++) {
            probability[index][i] = random.nextDouble();
            //System.out.print(probability[index][i] + "  ");
        }
        // System.out.println();

        //放缩概率
        scaling_probability(probability[index]);
        //System.out.print("放缩后：");
        for (int i = 0; i < num_of_plans; i++) {
            //System.out.print(probability[index][i] + "  ");
        }
        //System.out.println();

        //根据probability[index]来计算并发量
        for (int i = 0, val = 0; i < this.num_of_plans; i++) {
            for (int j = 0; j < this.n; j++) {
                concurrency[index][val] = (int) (Math.ceil(total_concurrency[j] * probability[index][i]));
                val++;
            }
        }

        int pos;
        //开始初始化total个节点
        for (pos = 0; pos < total; pos++) {

            r = (Math.random() * 32767 / ((double) 32767 + (double) (1)));

            //本问题中需要是整数，所以强制类型转换,进行四舍五入后转换为int
            r = r * (ub - lb) + lb;
            BigDecimal b = new BigDecimal(r);
            b = b.setScale(0, RoundingMode.HALF_UP);
            int random_pos = b.intValue();

            //random_pos是实际需要选取的服务序号,
            //临时变量tmp_group用来存储当前节点归属的组别,
            int tmp_group = pos % n;

            // 赋值给蜜源中对应的位置

            Foods[index][pos] = services[tmp_group][random_pos];
            //此处solution是一条蜜源
            solution[pos] = Foods[index][pos];
        }
        //计算函数值,f[]是全局的
        f[index] = calculateFunction(solution);
        //如果发现约束条件不满足，就对函数值进行惩罚，将代价变得极高。
        int flag = calc_time_spent(index);
        if (flag != 0) {
            //flag如果不等于0，说明时间超出限制或者并发量不满足条件，对函数值进行惩罚。
            //惩罚时根据不满足的条件的个数，不满足的条件越多，惩罚系数越大。
            f[index] = f[index] * flag * PenaltyRate;
        }

        //计算适应度
        fitness[index] = CalculateFitness(f[index]);
        //将所有蜜源的trail置零
        trial[index] = 0;
    }

    /*所有蜜源都已初始化 */
    void initial() {

        //初始化数据集
        init_services();

        int i;
        // 初始化FoodNumber个蜜源
        for (i = 0; i < FoodNumber; i++) {
            init(i);
        }

        // 接下来将第一个蜜源设为最佳蜜源
        GlobalMin = f[0];
        for (i = 0; i < D; i++) {
            GlobalParams[i] = Foods[0][i];
        }
    }

    void SendEmployedBees() {
        int i, j;
        boolean break_flag = false;
        /*雇佣峰阶段*/
        for (i = 0; i < FoodNumber; i++) {
            /*要改变的参数是随机确定的*/
            r = (Math.random() * 32767 / ((double) (32767) + (double) (1)));
            param2change = (int) (r * D);

            /*随机选择的方案用于产生方案i的突变方案*/
            r = (Math.random() * 32767 / ((double) (32767) + (double) (1)));
            neighbour = (int) (r * FoodNumber);

            //随机选择的解决方案必须与解决方案i不同
            while (neighbour == i) {
                r = (Math.random() * 32767 / ((double) (32767) + (double) (1)));
                neighbour = (int) (r * FoodNumber);
            }

            for (j = 0; j < D; j++) {
                solution[j] = Foods[i][j];
            }
            //貌似没法进行运算，所以直接替换被随机的一个服务商即可
            solution[param2change] = Foods[neighbour][param2change];

            //这儿还得对得到的新蜜源进行惩罚，如果新蜜源不符合条件的话。
            //solution是新的解决方案
            //ObjValSol 是生成的新的解决方案的函数值
            //FitnessSol 是生成的新的解决方案的适应度
            /*for (Service s : solution) {
                System.out.println(s);
            }*/
            ObjValSol = calculateFunction(solution);
            int flag;
            flag = calc_time_spent(solution, i);
            if (flag != 0) {
                //flag如果不等于0，说明时间超出限制或者并发量不满足条件，对函数值进行惩罚。
                //惩罚时根据不满足的条件的个数，不满足的条件越多，惩罚系数越大。
                ObjValSol = ObjValSol * flag * PenaltyRate;
            }
            FitnessSol = CalculateFitness(ObjValSol);

            /*在当前解决方案i和其突变体之间应用贪婪选择*/
            if (FitnessSol > fitness[i]) {
                /*
                If the mutant solution is better than the current solution i,
                replace the solution with the mutant and reset the trial counter of solution i
                */
                trial[i] = 0;
                for (j = 0; j < D; j++) {
                    Foods[i][j] = solution[j];
                }
                f[i] = ObjValSol;
                fitness[i] = FitnessSol;
            } else {
                /*如果解决方案无法改进，增加计数器*/
                trial[i] = trial[i] + 1;
            }
        }
        /*雇佣蜂阶段结束*/
    }

    void SendOnlookerBees() {
        //跟随蜂根据轮盘赌概率进行选择
        int i, j, t;
        i = 0;
        t = 0;
        /*onlooker Bee Phase 跟随蜂*/
        while (t < FoodNumber) {
            r = (Math.random() * 32767 / ((double) (32767) + (double) (1)));
            if (r < prob[i]) /*根据其选择的概率选择蜜源*/ {
                t++;
                // 需要生成新的蜜源进行贪婪选择
                /*要改变的参数是随机确定的*/
                r = (Math.random() * 32767 / ((double) (32767) + (double) (1)));
                param2change = (int) (r * D);

                /*一个随机选择的 solution is used in producing a mutant solution of the solution i*/
                r = (Math.random() * 32767 / ((double) (32767) + (double) (1)));
                neighbour = (int) (r * FoodNumber);

                //随机选择的解决方案必须与解决方案i不同
                while (neighbour == i) {
                    r = (Math.random() * 32767 / ((double) (32767) + (double) (1)));
                    neighbour = (int) (r * FoodNumber);
                }//循环结束

                for (j = 0; j < D; j++) {
                    solution[j] = Foods[i][j];
                }

                solution[param2change] = Foods[neighbour][param2change];
                //这儿还得对得到的新蜜源进行惩罚，如果新蜜源不符合条件的话。
                //solution是新的解决方案
                //ObjValSol 是生成的新的解决方案的函数值
                //FitnessSol 是生成的新的解决方案的适应度
                ObjValSol = calculateFunction(solution);
                int flag;
                flag = calc_time_spent(solution, i);
                if (flag != 0) {
                    //flag如果不等于0，说明时间超出限制或者并发量不满足条件，对函数值进行惩罚。
                    //惩罚时根据不满足的条件的个数，不满足的条件越多，惩罚系数越大。
                    ObjValSol = ObjValSol * flag * PenaltyRate;
                }
                FitnessSol = CalculateFitness(ObjValSol);

                /*在当前解决方案i和其突变体之间应用贪婪选择*/
                if (FitnessSol > fitness[i]) {
                    /*如果突变体优于当前个体i，则用突变体替换当前个体并重置个体i的试验计数器*/
                    trial[i] = 0;
                    for (j = 0; j < D; j++) {
                        Foods[i][j] = solution[j];
                    }
                    f[i] = ObjValSol;
                    fitness[i] = FitnessSol;
                } else {   /*如果个体i未被增强，则试验计数器+1*/
                    trial[i] = trial[i] + 1;
                }

            } /*if */
            i++;
            if (i == FoodNumber) {
                i = 0;
            }
        }/*while*/
        /*跟随蜂阶段结束*/
    }

    /*
    确定试验计数器超过“limit”值的蜜源。
    在基本ABC中，每个循环中只允许一个侦察蜂。
    此处一般无需更改。
    */
    void SendScoutBees() {
        int max_trial_index = 0;
        for (int i = 1; i < FoodNumber; i++) {
            if (trial[i] > trial[max_trial_index]) {
                max_trial_index = i;
            }
        }
        /*通过“limit”试验无法改善的食物来源被其使用的蜜蜂所放弃*/
        int limit = 20;
        if (trial[max_trial_index] >= limit) {
            init(max_trial_index);
        }
    }

    /*
    计算适应度，此函数为通用函数，无需更改，
    因为是最小化问题，所以目标函数值越大，适应度越小，越不容易存活。
    */
    private double CalculateFitness(double fun) {
        double result = 0;
        if (fun >= 0) {
            result = 1 / (fun + 1);
        } else {

            result = 1 + Math.abs(fun);
        }
        return result;
    }


    /* 选择食物来源的概率与其质量成比例*/
    /*可以使用不同的方案来计算概率值，最经典的就是轮盘赌↓↓↓*/
    void CalculateProbabilities() {
        int i;
        //计算总的适应值
        double total_fit = 0;
        for (i = 0; i < FoodNumber; ++i) {
            total_fit += fitness[i];
        }
        //计算轮盘赌概率
        for (int j = 0; j < FoodNumber; j++) {
            prob[j] = fitness[j] / total_fit;
        }
    }


    /*
    此处自定义自己的目标函数
    目前目标函数设置为所需的价格，要求价格最低。
     */
    private double calculateFunction(Service[] sol) {
        double total_price = 0;
        for (int i = 0; i < total; i++) {
            total_price += sol[i].price;
        }
        return total_price;
    }

    /*
    计算随机生成的蜜源是否需要进行惩罚
    按照时间约束与并发量约束.
     */
    private int calc_time_spent(int i) {
        int flag = 0;
        double[] time_spent = new double[num_of_plans];
        //sol_concurrency[4]是此条解的四个服务的实际并发量
        //需要每个服务的三条路径的并发量之和符合要求的并发量才行
        //double[] sol_concurrency = new double[n];
        //计算每条解的每条路径是否超时
        for (int j = 0; j < Foods[i].length; j++) {
            time_spent[j / n] += Foods[i][j].time;
        }
        // 如果某条路径超时，则flag
        for (double v : time_spent) {
            if (v > time_want_spent) {
                flag += TimePenaltyRate;
            }
        }
        //计算并发量是否满足,如果某条路径的某一个服务不满足并发量的需求，则flag+1
//        for (int j = 0; j < this.n; j++) {
//            for (int k = 0; k < this.num_of_plans; k++) {
//                sol_concurrency[j] += Foods[i][j + n * k].concurrency;
//            }
//        }

        for (int j = 0; j < this.total; j++) {
            if (Foods[i][j].concurrency < this.concurrency[i][j]) {
                flag++;
            }
        }

        return flag;
    }

    /*
    第二种计算随机生成的蜜源是否需要进行惩罚的函数，用于蜜蜂进化时的交换
    需要传进来交换后生成的新蜜源，以及老蜜源在Foods里面的编号。
    */
    private int calc_time_spent(Service[] sol, int index) {
        int flag = 0;
        double[] time_spent = new double[num_of_plans];
        double[] sol_concurrency = new double[num_of_plans];
        for (int i = 0; i < sol.length; i++) {
            time_spent[(i / n)] += sol[i].time;
        }
        for (double v : time_spent) {
            if (v > time_want_spent) {
                flag += TimePenaltyRate;
            }
        }
        //计算并发量是否满足,如果某条路径的某一个服务不满足并发量的需求，则flag+1
        for (int j = 0; j < sol.length; j++) {
            if (sol[j].concurrency < concurrency[index][(j)]) {
                flag += 1;
            }
        }
        return flag;
    }

    // 如果概率和不是1，就得通过放缩变成1
    private void scaling_probability(double[] sol_pro) {
        /*
        从num_of_plans获取到一共有多少个概率，这些概率加起来必须是1
        real_probability是此条解缩放前的总的概率
        zoom_factor 是缩放倍率，默认为1，也即不缩放。
        缩放策略为倍数缩放。
        然后再对第一个进行到1的精确取整。（这个存疑，先放在这儿，等实验出结果后再考虑）
        */
        double real_probability = 0.0;
        //double zoom_factor = 0.0;
        for (int i = 0; i < num_of_plans; i++) {
            real_probability = real_probability + sol_pro[i];
        }
        if (real_probability != 1) {
            for (int i = 0; i < num_of_plans; i++) {
                sol_pro[i] = sol_pro[i] / real_probability;
            }
        }
        //此处为了等概率实验，进行概率修改。
        //如果用不到，将true改为false即可
        if (random_pro_flag) {
            for (int i = 0; i < num_of_plans; i++) {
                sol_pro[i] = (double) (1.0 / num_of_plans);
            }
        }
    }

    /*
    初始化虚拟机组
    虚拟机列表从文件中加载
    **此处功能完成**
     */
    private List<String> readFile() {
        List<String> list = new ArrayList<String>();
        String encoding = "UTF-8";
        File file = new File(dataset_path);
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
        return list;
    }

    private void init_services() {

        // 初始化group
        /*
         * num_of_plans条解决方案,group_of_plans[]存储了每个index对应的组别
         * 比如n=2,num_of_plans=2那么group=[0,0,1,1]
         * */
//        int[] group_of_plans = new int[n * num_of_plans];
//        for (int i = 0; i < n * num_of_plans; i++) {
//            int j = (int) i / n;
//            group_of_plans[i] = j;
//        }

        List<String> serviceList;
        serviceList = readFile();
        int group = 0;
        for (int i = 0; i < n * (ub - lb + 1); i++) {
            //把文件内的每一行按照空格分词，存入一个字符串数组
            String[] tmpList = serviceList.get(i).split(" ");
            //计算组别
            group = (int) (i / (ub - lb + 1) + 1);
            int[] gg = new int[2];
            //存储组别
            gg[0] = group - 1;
            //存储组内编号
            gg[1] = i - (group - 1) * (int) (ub - lb + 1);
            //初始化单个虚拟机
            services[group - 1][gg[1]] = new Service(gg, Integer.parseInt(tmpList[0]), Integer.parseInt(tmpList[1]),
                    Integer.parseInt(tmpList[2]), Double.parseDouble(tmpList[3]), Double.parseDouble(tmpList[4]));
        }//完成初始化服务商
        /*for (Service[] ss : services) {
            for (Service s : ss) {
                System.out.println(s);
            }
        }*/
    }
}