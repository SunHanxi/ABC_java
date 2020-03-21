package abc;

// 虚拟机的类，包含了组别，价格，时间。
// group[0] 代表虚拟机的组别，group[1]代表组里的第几号
public class Service {
    private int cpu;
    private int ram;
    int concurrency;
    double price;
    double time;
    int[] group;

    //初始化对象时的顺序为
    //组别  组内序号  CPU  内存  并发量  价格  响应时间
    Service(int[] group, int cpu, int ram, int concurrency, double price, double time) {
        this.group = group;
        this.cpu = cpu;
        this.ram = ram;
        this.concurrency = concurrency;
        this.price = price;
        this.time = time;
    }

    @Override
    public String toString() {
        return "group: [" + this.group[0] + " " + this.group[1] + "] CPU&RAM: " + this.cpu + "&" + this.ram
                + " time: " + this.time + " concurrency: " + this.concurrency + " price: " + this.price;
    }
}
