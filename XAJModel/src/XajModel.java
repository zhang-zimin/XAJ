import java.util.ArrayList;

public class XajModel {


    /**  新安江模型计算理论
     *
     *
     * 2）产流：进行产流计算，产流计算理论为蓄满产量，
     * 蓄满指的是包气带的含水量等于田间的持水量，在土壤没有达到田间蓄水量的时候，
     * 所有的降雨被降雨吸收，转为张力水，当土壤湿度达到田间持水量的时候，继续降雨都产流。
     *
     */



    /*** 土壤蒸发参数 */
    private double WUM;  // 上层土壤蓄水容量
    private double WLM;  // 下层土壤蓄水容量
    private double WDM; //  深层土壤蓄水容量
    private double WM;  //  土壤总的蓄水容量


    /*** 土壤蓄水量参数 */
    private double C; // 深层蒸散发系数
    private double K; // 蒸发皿折算系数:水面蒸发折算系数是指水体表面蒸发量与标准蒸发皿蒸发量之间的比值，用来表征水体表面蒸发的强度。


    /*** 产流参数 */
    private double B;       // 流域张力水蓄水容量曲线分布系数
    private double Imp;     // 不透水面积与全流域面积的比值



    /*水源划分参数*/
    private double SM;      // 流域平均自由水蓄水容量
    private double EX;      // 自由水蓄水容量曲线指数
    private double KSS;     // 自由水蓄水库对壤中流的出流系数
    private double KG;      // 自由水蓄水库对地下径流出流系数


    /*封装每一阶段计算结果*/
    RunoffGenerationResult runoffGenerationResult; // 产流计算结果
    SourcePartitionResult sourcePartitionResult;// 水源划分计算结果


    /**
     * 设置流域土壤含水量参数
     *
     * @param wum 上层土壤蓄水容量
     * @param wlm 下层土壤蓄水容量
     * @param wdm 深层土壤蓄水容量
     * @return 完成流域土壤含水量参数设置
     */
    public XajModel SetSoilWaterStorageParam(double wum, double wlm, double wdm) {
        WUM = wum;
        WLM = wlm;
        WDM = wdm;
        WM = WUM + WLM + WDM;
        return this;
    }

    /**
     * 设置流域蒸散发参数
     *
     * @param k 蒸发皿折算系数：多年平均水面蒸发量/蒸发皿的蒸发量
     * @param c 深层蒸散发系数：它决定于深根植物占流域面积的比数，同时也与 UM + LM 的值有关，
     *          此值越大，深层蒸散发越困难。一般经验，在江南湿润地区 C 值约为 0.15~0.20，而在华北半湿润地区则在 0.09~0.12 左右。
     * @return 完成流域蒸散发参数设置的模型实例
     */
    public XajModel SetEvapotranspirationParam(double k, double c) {
        K = k;
        C = c;
        return this;
    }

    /**
     * 设置流域产流计算参数
     *
     * @param b   蓄水容量曲线的指数，它反映流域上蓄水容量分布的不均匀性。一般经验，流域越大，
     *            各种地质地形配置越多样，B值也越大。在山丘区，很小面积（几平方公里）的B为0.1左右，
     *            中等面积（300平方公里以内）的B为0.2~0.3左右，较大面积（数千平方公里）的B值为0.3~0.4左右。
     * @param imp 不透水面积与全流域面积的比值
     * @return 完成流域产流计算参数设置的模型实例
     */
    public XajModel SetRunoffGenerationParam(double b, double imp) {
        B = b;
        Imp = imp;
        return this;
    }


    /**
     * 设置流域水源划分参数
     *
     * @param sm  流域平均自由水蓄水容量
     * @param ex  自由水蓄水容量曲线指数
     * @param kss 自由水蓄水库对壤中流的出流系数
     * @param kg  自由水蓄水库对地下径流的出流系数
     * @return 完成流域水源划分参数设置的模型实例
     */
    public XajModel SetSourcePartitionParam(double sm, double ex, double kss, double kg) {
        SM = sm;
        EX = ex;
        KSS = kss;
        KG = kg;
        return this;
    }


    /**
     * 执行产流计算
     *
     * @param P   降雨量序列
     * @param EI  蒸发皿蒸发量序列
     * @param wu0 初始时刻上层土壤含水量
     * @param wl0 初始时刻下层土壤含水量
     * @param wd0 初始时刻深层土壤含水量
     * @return 完成产流计算的模型实例
     */
    public XajModel ComputeRunoffGeneration(double[] P, double[] EI, double wu0, double wl0, double wd0) {
        int n = P.length;
        runoffGenerationResult = new RunoffGenerationResult(n);
        LayeredSoilParam w0 = new LayeredSoilParam(wu0, wl0, wd0);  // 时段初土壤分层含水量
        double w0Sum = w0.getSum();                   // 时段初土壤总的含水量
        double wmmax = WM * (1 + B) / (1 - Imp);      // 流域最大蓄水容量
        for (int i = 0; i < n; ++i) {
            double p = P[i];                          // 时段内降雨量
            double ep = EI[i] * K;                    // 计算时段内实际蒸发量
            double e = ComputeE(ep, p, w0.U, w0.L);   // 计算时段内土壤蒸发量
            double pe = p - e;                        // 计算时段内净雨量
            double r = ComputeR(pe, wmmax, w0Sum);    // 计算时段内产流量
            w0 = ComputeLayeredW(w0, pe, r);          // 计算时段末土壤分层含水量
            w0Sum = w0.getSum();                      // 计算时段末土壤总含水量

            runoffGenerationResult.Set(i, pe, r, w0Sum);
        }
        return this;
    }


    /**
     * 执行水源划分计算（需要先完成产流计算）
     *
     * @param s0 初始时刻流域平均自由含水量
     * @param dt 时段长度
     * @return 完成水源划分计算的模型实例
     */
    public XajModel ComputeSourcePartition(double s0, double dt) {

        int n = runoffGenerationResult.Length;
        sourcePartitionResult = new SourcePartitionResult(n);

        // 计算步长内流域自由水蓄水库的壤中流出流系数 (KSSD)和地下水出流系数(KGD)与其日模型的出流系数(KSS)和(KG)的关系
        double KSSD = (1 - Math.pow(1 - (KG + KSS), dt / 24)) / (1 + KG / KSS);     // 转化后的壤中流出流系数
        double KGD = KSSD * KG / KSS;     // 转化后的地下水出流系数
        double Smax = (1 + EX) * SM;      // 流域点自由蓄水量最大值
        for (int i = 0; i < n; ++i) {
            double pe = runoffGenerationResult.PE[i];
            double w = runoffGenerationResult.W[i];
            double r = runoffGenerationResult.R[i];


        }

        return this;
    }

    /**
     * 计算时段内的土壤蒸发量，计算思路：上层按蒸散发能力蒸发，上层含水量不够蒸发时，剩余蒸散发能力从下层蒸发;
     * 下层蒸发与剩余蒸散发能力及下层含水量成正比，与下层蓄水容量成反比。要求计算的下层蒸发量与剩余蒸散发能力之
     * 比不小于深层蒸散发系数 C，否则，不足部分由下层含水量补给，当下层水量不够补给时，用深层含水量补给。
     *
     * @param ep  时段内实际蒸发量
     * @param p   时段内降雨量
     * @param wu0 时段初上层土壤含水量
     * @param wl0 时段初下层土壤含水量
     * @return 时段内的土壤蒸发量
     */
    private double ComputeE(double ep, double p, double wu0, double wl0) {
        double eu; // eu上层流域蒸散发量
        double el; // el下层流域蒸散发量
        double ed; // eu深层流域蒸散发量
        if (p + wu0 >= ep) {                        // 上层土壤含水量+降雨量 > 实际蒸发量的时候
            eu = ep;                               // 上层土壤蒸发量 = 实际蒸发量
            el = ed = 0;                           // 下层和深层的土壤没有蒸发
        }
        else {                                    // 上层土壤含水量+降雨量 < 实际蒸发量的时候
            eu = p + wu0;                         // 上层土壤蒸发量 = 上层土壤含水量 + 降雨量
            if (wl0 >= C * WLM) {                 // 下层土壤含水量 >= 深层蒸散发系数 * 下层土壤蓄水容量
                el = (ep - eu) * wl0 / WLM;       // 下层土壤蒸发量 = (实际蒸发量 - 上层土壤蒸发量) * 下层土壤初始的含水量 / 下层土壤蓄水容量
                ed = 0;                           // 深层土壤蒸发量没有蒸发
            }
            else if (wl0 >= C * (ep - eu)) {     // 下层土壤初始含水量 >= 深层蒸散发系数 * (实际蒸发量 - 上层土壤蒸发量)
                el = C * (ep - eu);              // 下层土壤蒸发量 = 深层蒸散发系数 * (实际蒸发量 - 上层土壤蒸发量)
                ed = 0;                          // 深层土壤蒸发量没有蒸发
            }
            else {
                el = C * wl0;                   // 下层土壤蒸发量 = 深层蒸散发系数 * 下层土壤初始含水量
                ed = C * (ep - eu) - el;        // 深层土壤蒸发量 = 深层蒸散发系数 * (实际蒸发量 - 上层土壤蒸发量) - 下层土壤蒸发量
            }
        }
        return eu + el + ed;        // 返回土壤总的蒸发量
    }

    /**
     * 计算时段内的产流量
     *
     * @param pe    此时段内的净雨量
     * @param wmmax 流域点最大蓄水容量
     * @param w0    时段初始时刻的土壤总含水量
     * @return 时段内产流量
     */
    private double ComputeR(double pe, double wmmax, double w0) {
        double a = wmmax * (1 - Math.pow(1 - w0 / WM, 1 / (1 + B)));  // 流域初始平均蓄水量相应的纵坐标

        double r;  // 产流量
        if (pe <= 0) {
            r = 0;
        }
        else if (pe + a < wmmax) {
            r = pe - WM + w0 + WM * Math.pow(1 - (pe + a) / wmmax, 1 + B);
        }
        else {
            r = pe - (WM - w0);
        }

        return r;
    }


    /**
     * 计算时段末的土壤含水量
     *
     * @param w0 时段初土壤含水量
     * @param pe 时段内净雨量
     * @param r  时段内产流量
     * @return 时段末的土壤含水量
     */
    private LayeredSoilParam ComputeLayeredW(LayeredSoilParam w0, double pe, double r) {
        double wu0 = w0.U; // 时段初的上层土壤含水量
        double wl0 = w0.L; // 时段初的下层土壤含水量
        double wd0 = w0.D; // 时段初的深层土壤含水量
        double wu; // 时段末的上层土壤含水量
        double wl; // 时段末的下层土壤含水量
        double wd; // 时段末的深层土壤含水量
        double dw = pe - r;  // 时段内土壤总含水量变化量
        if (dw > 0) {
            if (wu0 + dw < WUM) {
                //  时段内上层土壤的初始含水量+水量变化量<上层土壤的蓄水容量，则下层及深层土壤含水量不变，上层土壤含水量增加
                wu = wu0 + dw;
                // 下层土壤含水量不变
                wl = wl0;
                // 深层土壤含水量不变
                wd = wd0;
            }
            else {
                //  时段内上层土壤的初始含水量+水量变化量>上层土壤的蓄水容量，首先会充满上层土壤的含水量
                wu = WUM;
                //  下层土壤含水量增加,增加的量=下层土壤的蓄水容量+水量新增量-上层土壤含水量的补给量
                wl = wl0 + dw - (WUM - wu0);
                if (wl < WLM) {
                    // 如果没有增加的水量全部被下层吃掉，深层的水量不会有变化
                    wd = wd0;
                } else {
                    // 如果下层土壤含水量已经达到下层土壤的蓄水容量，则下层土壤含水量达到初始值，深层土壤含水量增加
                    wl = WLM;
                    // 深层土壤的时刻末含水量= 深层土壤的初始含水量+ 含水量的增加值 - 上层土壤的补给量 - 下层土壤的补给量
                    wd = wd0 + dw - (WUM - wu0) - (WLM - wl0);
                }
            }
            // 如果dw<0 说明土壤水分流失，没有补给
        } else {
            // 如果流失的水量<上层初始水量，只有上层土壤含水量减小，下层和深层不会变化
            if (wu0 + dw > 0) {
                wu = wu0 + dw;
                wl = wl0;
                wd = wd0;
            } else {
                // 如果流失的水量>上层土壤的初始水量，上层土壤含水量减小到0
                wu = 0;
                // 下层含水量 =  下层土壤的初始含水量-（总的流失量-上层的流失量）
                wl = wu0 + dw + wl0;
                // 如何wl不为负，表明深层土壤水量没有损失
                if (wl > 0) {
                    wd = wd0;
                } else {
                    // 如何为负，则下层水量全部损失
                    wl = 0;
                    // 并且深层的水量也有损失，深层土壤时刻末的含水量-（总的流失量dw-上层的初始含水量-下层的初始含水量）
                    wd = wu0 + wl0 + dw + wd0;
                }
            }
        }
        return new LayeredSoilParam(wu, wl, wd);
    }


    public static void main(String[] args) {

        // 降雨序列输入
        double[] P = {10,24.1,20.4,10.1,18.3,5.5,0.6,3.1,
                1.9,4.6,5,4.8,36.2,29,6,3.6,0.4,0,0.5,3.8,0,1.8,0.2,0.3}; // 初始降雨量

        // 蒸发序列输入
        double[] E ={0.1,0,0.1,0.5,0.7,0.9,0.8,0.7,0.5,0.3,
                0.2,0.1,0,0,0.1,0.6,0.8,1,0.9,0.8,0.7,0.5,0.3,0.1};// 初始蒸发量


        // 初始化XajModel
        XajModel xajModel = new XajModel();

        // 利用链式调用设置模型实例的参数
        xajModel.SetSoilWaterStorageParam(0.1, 0.2, 0.3)  // 设置土壤分层参数：上层土壤蓄水容量WUM=0.1，下层土壤蓄水容量WLM=0.2，深层土壤蓄水容量的WDM=0.3
                .SetEvapotranspirationParam(0.5,0.1) // 设置模型蒸发参数：蒸发皿折算系数=0.5，深层蒸散发系数=0.1
                .SetRunoffGenerationParam(0.2,0.1)  // 设置模型产流初始参数：不透水面积占比imp=0.5 ，流域张力水蓄水容量曲线分布系数b=0.2
                .SetSourcePartitionParam(30,1.2,0.6,0.1); // 设置模型水源划分参数，流域平均自由水蓄水容量设置为30mm，EX自由水蓄水容量曲线指数取1.2，自由水蓄水库对壤中流的出流系数取取0.6，自由水蓄水库对壤中流的出流系数0.1


        // 执行链式计算
        xajModel.ComputeRunoffGeneration(P, E,20,30,40);  // 计算产流，上层土壤初始含水量为20，上层土壤初始含水量为30，上层土壤初始含水量为40
    }





}



/**
 * 产流计算结果
 *
 * Created by Wenxuan on 2016/1/11.
 */
class RunoffGenerationResult {
    /**
     * 序列长度
     */
    public final int Length;

    /**
     * 净雨量序列
     */
    public double[] PE;

    /**
     * 产流量序列
     */
    public double[] R;

    /**
     * 土壤含水量序列
     */
    public double[] W;

    /**
     * @param n 序列长度
     */
    public RunoffGenerationResult(int n) {
        PE = new double[n];
        R = new double[n];
        W = new double[n];
        Length = n;
    }

    /**
     * 设置序列指定位置处的产流计算结果
     *
     * @param idx 序列位置索引
     * @param pe 净雨量
     * @param r 产流量
     * @param w 土壤含水量
     */
    public void Set(int idx, double pe, double r, double w){
        PE[idx] = pe;
        R[idx] = r;
        W[idx] = w;
    }
}

/**
 * 土壤分层参数
 */
class LayeredSoilParam {
    /**
     * @param u 上层参数值
     * @param l 下层参数值
     * @param d 深层参数值
     */
    public LayeredSoilParam(double u, double l, double d) {
        Set(u, l, d);
    }

    /**
     * @param u 上层参数值
     * @param l 下层参数值
     * @param d 深层参数值
     */
    public void Set(double u, double l, double d) {
        U = u;
        L = l;
        D = d;
    }

    /**
     * @return 分层参数值总和
     */
    public double getSum() {
        return U + L + D;
    }

    /**
     * 上层参数值
     */
    public double U;

    /**
     * 下层参数值
     */
    public double L;

    /**
     * 深层参数值
     */
    public double D;
}

/** * 汇流计算结果 */
 class RunoffConcentrationResult {
    /**
     * 序列长度
     */
    public int Length;

    /**
     * 地表径流序列
     */
    public double[] QRS;

    /**
     * 壤中流序列
     */
    public double[] QRSS;

    /**
     * 地下径流序列
     */
    public double[] QRG;

    /**
     * 总径流序列
     */
    public double[] Q;

    /**
     * @param n 序列长度
     */
    public RunoffConcentrationResult(int n) {
        QRS = new double[n];
        QRSS = new double[n];
        QRG = new double[n];
        Q = new double[n];
        Length = n;
    }

    /**
     * @param QRS 地表径流序列
     * @param QRSS 壤中流序列
     * @param QRG 地下径流序列
     * @param Q 总径流序列
     */
    public RunoffConcentrationResult(double[] QRS, double[] QRSS,
                                     double[] QRG, double[] Q) {
//        this.Length = ParameterValidation.CheckIdenticalLength(QRS, QRSS, QRG, Q);

        this.QRS = QRS;
        this.QRSS = QRSS;
        this.QRG = QRG;
        this.Q = Q;
    }
}

/** * 水源地划分结果 */
class SourcePartitionResult{
    /**
     * 序列长度
     */
    public final int Length;

    /**
     * 不透水面积产流序列
     */
    public double[] RIMP;

    /**
     * 地表径流产流序列
     */
    public double[] RS;

    /**
     * 壤中流产流序列
     */
    public double[] RSS;

    /**
     * 地下径流产流序列
     */
    public double[] RG;

    /**
     * 流域平均自由含水量序列
     */
    public double[] S;

    /**
     * @param n 序列长度
     */
    public SourcePartitionResult(int n) {
        RIMP = new double[n];
        RS = new double[n];
        RSS = new double[n];
        RG = new double[n];
        S = new double[n];
        Length = n;
    }

    /**
     * 设置序列指定位置处的水源划分计算结果
     *
     * @param idx 序列位置索引
     * @param rimp 不透水面积产流
     * @param rs 地表径流产流
     * @param rss 壤中流产流
     * @param rg 地下径流产流
     * @param s 流域平均自由含水量
     */
    public void Set(int idx, double rimp, double rs, double rss, double rg, double s){
        RIMP[idx] = rimp;
        RS[idx] = rs;
        RSS[idx] = rss;
        RG[idx] = rg;
        S[idx] = s;
    }
}