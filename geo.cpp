#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define MAXLINE     65536       // 算例文件每行最大字节数
#define MAXTYPE     20          // 任务类型数量上限
#define MAXNAMELEN  32          // 模式名称最大长度
#define MININT_MOVE (-2000000000)   // 禁忌搜索中 delta 的“负无穷”初值

using namespace std;

//====================================================================
// 问题数据（由 read_instance 填充）
//====================================================================
char *instanceName;

int numBeam;          // B：波束数量
int numTask;          // N：任务数量
int numMode;          // M：服务模式数量
int numType;          // T：任务类型数量

char **modeName;      // modeName[m]          ：模式 m 的名称
int  *modeBasePower;  // modeBasePower[m]      ：模式 m 的固定基础功耗
int **typeCompatMode; // typeCompatMode[t][m]=1 表示任务类型 (t+1) 可使用模式 m

int *beamBWCap;       // beamBWCap[b]   ：波束 b 的带宽容量
int *beamPWCap;       // beamPWCap[b]   ：波束 b 的功率容量

int *taskType;        // taskType[j]    ：任务 j 的类型（从 1 起编）
int *taskProfit;      // taskProfit[j]  ：任务 j 的收益
int *taskBWDemand;    // taskBWDemand[j]：任务 j 的带宽需求
int *taskPWDemand;    // taskPWDemand[j]：任务 j 的功率需求

//====================================================================
// 当前解状态
//====================================================================
int *beamMode;        // beamMode[b]  ：波束 b 所选的模式编号
int *taskBeam;        // taskBeam[j]  ：服务任务 j 的波束（-1 表示未服务）
int *remBW;           // remBW[b]     ：波束 b 的剩余带宽
int *remPW;           // remPW[b]     ：波束 b 的剩余功率（扣除基础功耗后）
int  totalProfit;

double maxRunTime;    // 运行时间上限（秒）
double bestTime;      // 找到本阶段最优解的时刻（秒）
int    seed;          // 随机种子

//====================================================================
// 最优解快照
//====================================================================
int *bestBeamMode;    // 当前 ILS 阶段内找到的最优解
int *bestTaskBeam;
int  bestProfit;

int *globalBeamMode;  // 所有 ILS 阶段中找到的全局最优解
int *globalTaskBeam;
int  globalProfit;
double globalBestTime;

//====================================================================
// 禁忌搜索 / 搜索辅助量
//====================================================================
int *moveFreq;        // moveFreq[j]   ：任务 j 被移动的次数（用于扰动）
int **typeProfitSum;  // typeProfitSum[b][t]：波束 b 上已服务类型 t 任务的收益之和

//====================================================================
// 基于解的禁忌（参考 SBTS）：用 3 路加权哈希把整解映射到哈希表去重
//====================================================================
// 决策被选中时贡献权重、未选贡献 0（仿 SBTS）。任务-波束与波束-模式
// 两部分变量都参与哈希。每路权重为扁平数组，与现有代码风格一致。
// 持久哈希（学 MNSB-TS：全程只初始化一次，不清空），记忆贯穿整轮搜索。
// 用 char(0/1) 省内存：3*1e8 字节 ≈ 300MB/进程（与 SBTS 原版 L=1e8 对齐）。
#define HASH_L 100000000      // 哈希表槽位数（每路）
long long *taskHashW1, *taskHashW2, *taskHashW3;   // [numTask*numBeam]，下标 j*numBeam+b
long long *beamHashW1, *beamHashW2, *beamHashW3;   // [numBeam*numMode]，下标 b*numMode+m
char *hashTab1, *hashTab2, *hashTab3;        // [HASH_L]：访问过的解(全程保留，0/1)
long long Hx1, Hx2, Hx3;                     // 当前解的三路哈希和

int divStrength;      // 动态多样化强度：本次扰动要破坏的波束数（学 MNSB-TS：逃出回弱、停滞渐强）

//====================================================================
// 每步局部搜索重建的加速索引
//====================================================================
// 按波束分组的已服务任务 CSR 桶，供交换邻域跳过无关任务。
int *bucketTask;      // [numTask]   ：已服务任务编号，按波束连续存放
int *bucketStart;     // [numBeam+1] ：bucketStart[b]..bucketStart[b+1] 为波束 b 的区间
// 按（源波束，候选模式）分组的已服务任务。
int *compatBucketTask;
int *compatBucketStart;   // [numBeam*numMode+1]
// 当前模式下可接受某任务类型的波束列表。
int **typeCompatBeam;     // typeCompatBeam[t][p]：与类型 t 兼容的波束编号
int  *typeCompatBeamCount;
// insertMinPW[m][bw]：未服务、与模式 m 兼容且带宽需求<=bw 的任务的最小功率需求。
int **insertMinPW;
int   maxBeamBWCap;

static char g_line[MAXLINE];

//====================================================================
// 读入算例
//====================================================================
int find_mode_index(const char *name)
{
    for (int m = 0; m < numMode; m++)
        if (strcmp(modeName[m], name) == 0)
            return m;
    return -1;
}

// 对 val[] 做降序快速排序，同时携带 idx[]（贪心排序键）。
void qsort_desc(double *val, int *idx, int l, int r)
{
    if (l < r)
    {
        int    i = l, j = r;
        double xv = val[l];
        int    xi = idx[l];
        while (i < j)
        {
            while (i < j && val[j] <= xv) j--;
            if (i < j) { val[i] = val[j]; idx[i] = idx[j]; i++; }
            while (i < j && val[i] > xv)  i++;
            if (i < j) { val[j] = val[i]; idx[j] = idx[i]; j--; }
        }
        val[i] = xv;
        idx[i] = xi;
        qsort_desc(val, idx, l, i - 1);
        qsort_desc(val, idx, i + 1, r);
    }
}

void read_instance()
{
    ifstream FIC;
    FIC.open(instanceName);
    if (FIC.fail())
    {
        cout << "Cannot open file: " << instanceName << endl;
        exit(0);
    }

    // ---- 首行：B=5 N=50 M=3 Cbw=1578 Cpw=503 ----
    FIC.getline(g_line, MAXLINE);
    int cbw_total, cpw_total;
    sscanf(g_line, "B=%d N=%d M=%d Cbw=%d Cpw=%d",
           &numBeam, &numTask, &numMode, &cbw_total, &cpw_total);

    modeName = new char *[numMode];
    for (int m = 0; m < numMode; m++)
        modeName[m] = new char[MAXNAMELEN];
    modeBasePower = new int[numMode];

    // 空行 / "Modes" / "mode base_power"
    FIC.getline(g_line, MAXLINE);
    FIC.getline(g_line, MAXLINE);
    FIC.getline(g_line, MAXLINE);

    // ---- M 行模式："G 10" ----
    for (int m = 0; m < numMode; m++)
    {
        FIC.getline(g_line, MAXLINE);
        char tmp[MAXNAMELEN];
        sscanf(g_line, "%s %d", tmp, &modeBasePower[m]);
        strcpy(modeName[m], tmp);
    }

    // 空行 / "Task types" / "type name modes"
    FIC.getline(g_line, MAXLINE);
    FIC.getline(g_line, MAXLINE);
    FIC.getline(g_line, MAXLINE);

    // ---- 任务类型行，直到空行 ----
    static char typeLines[MAXTYPE][256];
    numType = 0;
    while (FIC.getline(g_line, MAXLINE))
    {
        if ((int)strlen(g_line) == 0) break;   // 空行结束本段
        strncpy(typeLines[numType], g_line, 255);
        typeLines[numType][255] = '\0';
        numType++;
    }

    typeCompatMode = new int *[numType];
    for (int t = 0; t < numType; t++)
    {
        typeCompatMode[t] = new int[numMode];
        for (int m = 0; m < numMode; m++) typeCompatMode[t][m] = 0;

        int  tid;
        char tname[64], tmodes[64];
        sscanf(typeLines[t], "%d %s %s", &tid, tname, tmodes);

        char *tok = strtok(tmodes, "|");
        while (tok != NULL)
        {
            int midx = find_mode_index(tok);
            if (midx >= 0) typeCompatMode[t][midx] = 1;
            tok = strtok(NULL, "|");
        }
    }

    beamBWCap = new int[numBeam];
    beamPWCap = new int[numBeam];

    // 波束带宽容量
    FIC.getline(g_line, MAXLINE);                         // 表头
    for (int b = 0; b < numBeam; b++) FIC >> beamBWCap[b];
    FIC.getline(g_line, MAXLINE);                         // 数据行结束
    FIC.getline(g_line, MAXLINE);                         // 空行

    // 波束功率容量
    FIC.getline(g_line, MAXLINE);                         // 表头
    for (int b = 0; b < numBeam; b++) FIC >> beamPWCap[b];
    FIC.getline(g_line, MAXLINE);
    FIC.getline(g_line, MAXLINE);

    taskType     = new int[numTask];
    taskProfit   = new int[numTask];
    taskBWDemand = new int[numTask];
    taskPWDemand = new int[numTask];

    // 任务类型
    FIC.getline(g_line, MAXLINE);                         // 表头
    for (int j = 0; j < numTask; j++) FIC >> taskType[j];
    FIC.getline(g_line, MAXLINE);
    FIC.getline(g_line, MAXLINE);

    // 任务收益
    FIC.getline(g_line, MAXLINE);                         // 表头
    for (int j = 0; j < numTask; j++) FIC >> taskProfit[j];
    FIC.getline(g_line, MAXLINE);
    FIC.getline(g_line, MAXLINE);

    // 任务带宽需求
    FIC.getline(g_line, MAXLINE);                         // 表头
    for (int j = 0; j < numTask; j++) FIC >> taskBWDemand[j];
    FIC.getline(g_line, MAXLINE);
    FIC.getline(g_line, MAXLINE);

    // 任务功率需求（最后一段）
    FIC.getline(g_line, MAXLINE);                         // 表头
    for (int j = 0; j < numTask; j++) FIC >> taskPWDemand[j];

    FIC.close();
}

//====================================================================
// 贪心排序键
//====================================================================
// task_score：任务 j 在波束 b、模式 m 下，单位归一化资源占用的收益。
// 带宽与功率分别按波束容量归一化，占用任一资源比例大的任务得分更低。
// 仅作贪心排序/优先级键，不作为目标函数。
double task_score(int j, int b, int m)
{
    double pwCap = (double)(beamPWCap[b] - modeBasePower[m]);
    if (pwCap <= 0.0) return -1.0;

    double bwPart = (beamBWCap[b] > 0)
                    ? (double)taskBWDemand[j] / beamBWCap[b] : 0.0;
    double pwPart = (double)taskPWDemand[j] / pwCap;
    double denom  = bwPart + pwPart;

    if (denom <= 0.0) return (double)taskProfit[j];
    return (double)taskProfit[j] / denom;
}

// task_priority：任务 j 在任意模式兼容波束上可达到的最佳 task_score。
// 若 j 没有任何兼容波束，返回 -1。
double task_priority(int j)
{
    int    t    = taskType[j] - 1;
    double best = -1.0;
    for (int b = 0; b < numBeam; b++)
    {
        if (beamMode[b] < 0)                 continue;
        if (!typeCompatMode[t][beamMode[b]]) continue;
        double s = task_score(j, b, beamMode[b]);
        if (s > best) best = s;
    }
    return best;
}

//====================================================================
// 贪心初始解
//====================================================================
// 顺序提交贪心：对每个波束，从当前未服务池中试填以选模式，
// 提交后再处理下一波束，使波束之间竞争任务。
void greedy_init()
{
    beamMode = new int[numBeam];
    taskBeam = new int[numTask];
    remBW    = new int[numBeam];
    remPW    = new int[numBeam];

    for (int j = 0; j < numTask; j++) taskBeam[j] = -1;

    int    *tmpIdx     = new int   [numTask];
    double *tmpVal     = new double[numTask];
    int    *trialAdded = new int   [numTask];
    int    *bestAdded  = new int   [numTask];

    totalProfit = 0;

    for (int b = 0; b < numBeam; b++)
    {
        int bestM      = -1;
        int bestGain   = -1;
        int nBestAdded = 0;

        for (int m = 0; m < numMode; m++)
        {
            if (modeBasePower[m] > beamPWCap[b]) continue;   // 功率不可行

            // 将当前未服务且与模式 m 兼容的任务按密度得分排序
            int nFree = 0;
            for (int j = 0; j < numTask; j++)
                if (taskBeam[j] == -1 && typeCompatMode[taskType[j] - 1][m])
                {
                    tmpIdx[nFree] = j;
                    tmpVal[nFree] = task_score(j, b, m);
                    nFree++;
                }
            if (nFree > 0) qsort_desc(tmpVal, tmpIdx, 0, nFree - 1);

            int rem_bw = beamBWCap[b];
            int rem_pw = beamPWCap[b] - modeBasePower[m];
            int nTrial = 0, gain = 0;
            for (int ki = 0; ki < nFree; ki++)
            {
                int j = tmpIdx[ki];
                if (taskBWDemand[j] <= rem_bw && taskPWDemand[j] <= rem_pw)
                {
                    rem_bw -= taskBWDemand[j];
                    rem_pw -= taskPWDemand[j];
                    gain   += taskProfit[j];
                    trialAdded[nTrial++] = j;
                }
            }

            if (gain > bestGain)
            {
                bestGain   = gain;
                bestM      = m;
                nBestAdded = nTrial;
                for (int a = 0; a < nTrial; a++) bestAdded[a] = trialAdded[a];
            }
        }

        if (bestM < 0)   // 无功率可行模式：选基础功耗最低的模式
        {
            for (int m = 0; m < numMode; m++)
                if (bestM < 0 || modeBasePower[m] < modeBasePower[bestM]) bestM = m;
            nBestAdded = 0;
        }

        // 提交所选模式及其填充（此处尚未分配 typeProfitSum）
        beamMode[b] = bestM;
        remBW[b]    = beamBWCap[b];
        remPW[b]    = beamPWCap[b] - modeBasePower[bestM];
        for (int a = 0; a < nBestAdded; a++)
        {
            int j = bestAdded[a];
            taskBeam[j]  = b;
            remBW[b]    -= taskBWDemand[j];
            remPW[b]    -= taskPWDemand[j];
            totalProfit += taskProfit[j];
        }
    }

    delete[] tmpIdx;
    delete[] tmpVal;
    delete[] trialAdded;
    delete[] bestAdded;
}

//====================================================================
// 底层解维护
//====================================================================
// feasible_on：在当前剩余资源下，任务 j 能否分配到波束 b？
// （模式兼容 + 带宽足够 + 功率足够）
int feasible_on(int j, int b)
{
    if (beamMode[b] < 0)                 return 0;   // 模式尚未确定
    int t = taskType[j] - 1;
    if (!typeCompatMode[t][beamMode[b]]) return 0;
    if (taskBWDemand[j] > remBW[b])      return 0;
    if (taskPWDemand[j] > remPW[b])      return 0;
    return 1;
}

// 重建 insertMinPW[m][bw]：未服务、与模式 m 兼容且带宽需求不超过 bw 的任务的最小功率需求。
// 用于快速判断某波束腾出空间后是否还能装入任意未服务任务。
void rebuild_insert_filter()
{
    const int INF = 1000000000;

    for (int m = 0; m < numMode; m++)
        for (int w = 0; w <= maxBeamBWCap; w++)
            insertMinPW[m][w] = INF;

    for (int j = 0; j < numTask; j++)
    {
        if (taskBeam[j] != -1) continue;
        int bw = taskBWDemand[j];
        if (bw > maxBeamBWCap) continue;

        int t = taskType[j] - 1;
        for (int m = 0; m < numMode; m++)
            if (typeCompatMode[t][m] && taskPWDemand[j] < insertMinPW[m][bw])
                insertMinPW[m][bw] = taskPWDemand[j];
    }

    // 对带宽做前缀最小，使 insertMinPW[m][w] 覆盖 BW<=w
    for (int m = 0; m < numMode; m++)
        for (int w = 1; w <= maxBeamBWCap; w++)
            if (insertMinPW[m][w - 1] < insertMinPW[m][w])
                insertMinPW[m][w] = insertMinPW[m][w - 1];
}

int beam_can_insert_unserved_with_rem(int b, int bwFree, int pwFree)
{
    if (beamMode[b] < 0) return 0;
    if (bwFree < 0 || pwFree < 0) return 0;
    if (bwFree > maxBeamBWCap) bwFree = maxBeamBWCap;
    return insertMinPW[beamMode[b]][bwFree] <= pwFree;
}

// add_task / remove_task：分配或取消任务，保持 remBW / remPW /
// taskBeam / totalProfit / typeProfitSum 一致。（此处不改变模式。）
void add_task(int j, int b)
{
    taskBeam[j]  = b;
    remBW[b]    -= taskBWDemand[j];
    remPW[b]    -= taskPWDemand[j];
    totalProfit += taskProfit[j];
    typeProfitSum[b][taskType[j] - 1] += taskProfit[j];
}

void remove_task(int j)
{
    int b = taskBeam[j];
    remBW[b]    += taskBWDemand[j];
    remPW[b]    += taskPWDemand[j];
    totalProfit -= taskProfit[j];
    typeProfitSum[b][taskType[j] - 1] -= taskProfit[j];
    taskBeam[j]  = -1;
}

//====================================================================
// 基于解的禁忌：哈希权重生成与哈希表读写
//====================================================================
// init_hash：为 3 路哈希各生成一套权重（仿 SBTS：pow(k+1,e) 后随机洗牌）。
// 任务-波束权重 numTask*numBeam 个，波束-模式权重 numBeam*numMode 个。
void init_hash()
{
    int nT = numTask * numBeam;
    int nB = numBeam * numMode;

    for (int k = 0; k < nT; k++)
    {
        taskHashW1[k] = (long long)pow((double)(k + 1), 1.2);
        taskHashW2[k] = (long long)pow((double)(k + 1), 1.6);
        taskHashW3[k] = (long long)pow((double)(k + 1), 2.0);
    }
    for (int k = 0; k < nB; k++)
    {
        beamHashW1[k] = (long long)pow((double)(nT + k + 1), 1.2);
        beamHashW2[k] = (long long)pow((double)(nT + k + 1), 1.6);
        beamHashW3[k] = (long long)pow((double)(nT + k + 1), 2.0);
    }

    // 每路各自随机两两交换打乱（仿 SBTS initial_hash 的洗牌）
    for (int k = 0; k < nT && nT > 1; k++)
    {
        int r;
        r = rand() % nT; { long long t = taskHashW1[k]; taskHashW1[k] = taskHashW1[r]; taskHashW1[r] = t; }
        r = rand() % nT; { long long t = taskHashW2[k]; taskHashW2[k] = taskHashW2[r]; taskHashW2[r] = t; }
        r = rand() % nT; { long long t = taskHashW3[k]; taskHashW3[k] = taskHashW3[r]; taskHashW3[r] = t; }
    }
    for (int k = 0; k < nB && nB > 1; k++)
    {
        int r;
        r = rand() % nB; { long long t = beamHashW1[k]; beamHashW1[k] = beamHashW1[r]; beamHashW1[r] = t; }
        r = rand() % nB; { long long t = beamHashW2[k]; beamHashW2[k] = beamHashW2[r]; beamHashW2[r] = t; }
        r = rand() % nB; { long long t = beamHashW3[k]; beamHashW3[k] = beamHashW3[r]; beamHashW3[r] = t; }
    }

    memset(hashTab1, 0, HASH_L);
    memset(hashTab2, 0, HASH_L);
    memset(hashTab3, 0, HASH_L);
}

// compute_hash：从当前解全量重算 Hx1/Hx2/Hx3（已服务任务 + 各波束模式）。
void compute_hash()
{
    Hx1 = 0; Hx2 = 0; Hx3 = 0;
    for (int j = 0; j < numTask; j++)
    {
        int b = taskBeam[j];
        if (b < 0) continue;
        int idx = j * numBeam + b;
        Hx1 += taskHashW1[idx];
        Hx2 += taskHashW2[idx];
        Hx3 += taskHashW3[idx];
    }
    for (int b = 0; b < numBeam; b++)
    {
        if (beamMode[b] < 0) continue;
        int idx = b * numMode + beamMode[b];
        Hx1 += beamHashW1[idx];
        Hx2 += beamHashW2[idx];
        Hx3 += beamHashW3[idx];
    }
}

// is_visited：给定候选解的三路哈希和，判断是否访问过（三表全中）。
int is_visited(long long h1, long long h2, long long h3)
{
    int i1 = (int)(h1 % HASH_L); if (i1 < 0) i1 += HASH_L;
    int i2 = (int)(h2 % HASH_L); if (i2 < 0) i2 += HASH_L;
    int i3 = (int)(h3 % HASH_L); if (i3 < 0) i3 += HASH_L;
    return hashTab1[i1] && hashTab2[i2] && hashTab3[i3];
}

// mark_current：把当前解（Hx1/Hx2/Hx3）三路槽位标记为已访问。
void mark_current()
{
    int i1 = (int)(Hx1 % HASH_L); if (i1 < 0) i1 += HASH_L;
    int i2 = (int)(Hx2 % HASH_L); if (i2 < 0) i2 += HASH_L;
    int i3 = (int)(Hx3 % HASH_L); if (i3 < 0) i3 += HASH_L;
    hashTab1[i1] = 1;
    hashTab2[i2] = 1;
    hashTab3[i3] = 1;
}

//====================================================================
// 快照辅助（阶段最优与全局最优）
//====================================================================
void save_best(double beginTime)
{
    bestTime   = ((double)clock() - beginTime) / CLOCKS_PER_SEC;
    bestProfit = totalProfit;
    for (int b = 0; b < numBeam; b++) bestBeamMode[b] = beamMode[b];
    for (int j = 0; j < numTask; j++) bestTaskBeam[j] = taskBeam[j];
}

// restore_best：将阶段最优解写回当前解，并重算 remBW / remPW / totalProfit / typeProfitSum。
void restore_best()
{
    for (int b = 0; b < numBeam; b++)
    {
        beamMode[b] = bestBeamMode[b];
        remBW[b]    = beamBWCap[b];
        remPW[b]    = beamPWCap[b] - modeBasePower[beamMode[b]];
        for (int t = 0; t < numType; t++) typeProfitSum[b][t] = 0;
    }
    totalProfit = 0;
    for (int j = 0; j < numTask; j++)
    {
        taskBeam[j] = bestTaskBeam[j];
        if (taskBeam[j] >= 0)
        {
            int b = taskBeam[j];
            remBW[b]    -= taskBWDemand[j];
            remPW[b]    -= taskPWDemand[j];
            totalProfit += taskProfit[j];
            typeProfitSum[b][taskType[j] - 1] += taskProfit[j];
        }
    }
}

void save_global_from_best()
{
    globalBestTime = bestTime;
    globalProfit   = bestProfit;
    for (int b = 0; b < numBeam; b++) globalBeamMode[b] = bestBeamMode[b];
    for (int j = 0; j < numTask; j++) globalTaskBeam[j] = bestTaskBeam[j];
}

void restore_global()
{
    for (int b = 0; b < numBeam; b++)
    {
        beamMode[b] = globalBeamMode[b];
        remBW[b]    = beamBWCap[b];
        remPW[b]    = beamPWCap[b] - modeBasePower[beamMode[b]];
        for (int t = 0; t < numType; t++) typeProfitSum[b][t] = 0;
    }
    totalProfit = 0;
    for (int j = 0; j < numTask; j++)
    {
        taskBeam[j] = globalTaskBeam[j];
        if (taskBeam[j] >= 0)
        {
            int b = taskBeam[j];
            remBW[b]    -= taskBWDemand[j];
            remPW[b]    -= taskPWDemand[j];
            totalProfit += taskProfit[j];
            typeProfitSum[b][taskType[j] - 1] += taskProfit[j];
        }
    }
}

void alloc_search()
{
    bestBeamMode   = new int[numBeam];
    bestTaskBeam   = new int[numTask];
    globalBeamMode = new int[numBeam];
    globalTaskBeam = new int[numTask];
    moveFreq       = new int[numTask];
    bucketTask     = new int[numTask];
    bucketStart    = new int[numBeam + 1];
    compatBucketTask  = new int[numTask * numMode];
    compatBucketStart = new int[numBeam * numMode + 1];
    typeCompatBeam      = new int*[numType];
    typeCompatBeamCount = new int[numType];
    for (int t = 0; t < numType; t++) typeCompatBeam[t] = new int[numBeam];

    maxBeamBWCap = 0;
    for (int b = 0; b < numBeam; b++)
        if (beamBWCap[b] > maxBeamBWCap) maxBeamBWCap = beamBWCap[b];
    insertMinPW = new int*[numMode];
    for (int m = 0; m < numMode; m++)
        insertMinPW[m] = new int[maxBeamBWCap + 1];

    typeProfitSum = new int*[numBeam];
    for (int b = 0; b < numBeam; b++)
    {
        typeProfitSum[b] = new int[numType];
        for (int t = 0; t < numType; t++) typeProfitSum[b][t] = 0;
    }

    // 基于解禁忌：权重表与哈希表
    taskHashW1 = new long long[numTask * numBeam];
    taskHashW2 = new long long[numTask * numBeam];
    taskHashW3 = new long long[numTask * numBeam];
    beamHashW1 = new long long[numBeam * numMode];
    beamHashW2 = new long long[numBeam * numMode];
    beamHashW3 = new long long[numBeam * numMode];
    hashTab1 = new char[HASH_L];
    hashTab2 = new char[HASH_L];
    hashTab3 = new char[HASH_L];
    init_hash();
}

//====================================================================
// 禁忌局部搜索
//====================================================================
// record_candidate：在并行（并集）邻域中，记录目前见到的最优 delta 候选移动。
// 平局时用蓄水池抽样均匀随机打破。
void record_candidate(int kind, int a, int b, int c, int delta,
                      int &bestDelta, int &numBest,
                      int &chosenKind, int &chosenA, int &chosenB, int &chosenC)
{
    if (delta < bestDelta) return;

    if (delta > bestDelta)
    {
        bestDelta  = delta;
        numBest    = 1;
        chosenKind = kind;
        chosenA    = a;
        chosenB    = b;
        chosenC    = c;
        return;
    }

    numBest++;
    if (rand() % numBest == 0)
    {
        chosenKind = kind;
        chosenA    = a;
        chosenB    = b;
        chosenC    = c;
    }
}

// apply_mode_flip：将波束 b 切换到模式 m2，仅驱逐与新模式类型不兼容的任务。
// 此处不做回填。
void apply_mode_flip(int b, int m2)
{
    int saveMode = beamMode[b];

    for (int idx = bucketStart[b]; idx < bucketStart[b + 1]; idx++)
    {
        int j = bucketTask[idx];
        if (taskBeam[j] == b && !typeCompatMode[taskType[j] - 1][m2])
        {
            remove_task(j);
            moveFreq[j]++;
        }
    }

    beamMode[b] = m2;
    remPW[b] += modeBasePower[saveMode] - modeBasePower[m2];
}

// update_hash_for_move：在提交某 move 之前调用，用提交前状态把该 move 的增量
// 加到 Hx1/Hx2/Hx3（SBTS 式增量哈希，省去全量重算）。增量公式与各邻域候选判禁忌
// 时所用的完全一致（已经自检验证逐位相符）。kind/a/bb/cc 同提交段语义。
void update_hash_for_move(int kind, int a, int bb, int cc)
{
    if (kind == 1)                       // 插入 a→bb
    {
        int idx = a * numBeam + bb;
        Hx1 += taskHashW1[idx]; Hx2 += taskHashW2[idx]; Hx3 += taskHashW3[idx];
    }
    else if (kind == 2)                  // 删除 a（在其当前波束）
    {
        int idx = a * numBeam + taskBeam[a];
        Hx1 -= taskHashW1[idx]; Hx2 -= taskHashW2[idx]; Hx3 -= taskHashW3[idx];
    }
    else if (kind == 3)                  // 交换：a 入 bb、cc 出 bb
    {
        int idi = a * numBeam + bb, idk = cc * numBeam + bb;
        Hx1 += taskHashW1[idi] - taskHashW1[idk];
        Hx2 += taskHashW2[idi] - taskHashW2[idk];
        Hx3 += taskHashW3[idi] - taskHashW3[idk];
    }
    else if (kind == 4)                  // 模式翻转：波束 a，m1->bb，驱逐不兼容任务
    {
        int bidx = a * numMode + bb, bidx0 = a * numMode + beamMode[a];
        Hx1 += beamHashW1[bidx] - beamHashW1[bidx0];
        Hx2 += beamHashW2[bidx] - beamHashW2[bidx0];
        Hx3 += beamHashW3[bidx] - beamHashW3[bidx0];
        for (int idx = bucketStart[a]; idx < bucketStart[a + 1]; idx++)
        {
            int j = bucketTask[idx];
            if (taskBeam[j] == a && !typeCompatMode[taskType[j] - 1][bb])
            {
                int idj = j * numBeam + a;
                Hx1 -= taskHashW1[idj]; Hx2 -= taskHashW2[idj]; Hx3 -= taskHashW3[idj];
            }
        }
    }
    else if (kind == 5)
    {
        int b1 = taskBeam[a];
        if (cc < 0)                      // real-dummy：a 从 b1 迁到 bb
        {
            int idi2 = a * numBeam + bb, idi1 = a * numBeam + b1;
            Hx1 += taskHashW1[idi2] - taskHashW1[idi1];
            Hx2 += taskHashW2[idi2] - taskHashW2[idi1];
            Hx3 += taskHashW3[idi2] - taskHashW3[idi1];
        }
        else                             // real-real：a:b1->b2, cc:b2->b1
        {
            int b2 = taskBeam[cc];
            int idi2 = a * numBeam + b2, idi1 = a * numBeam + b1;
            int idk1 = cc * numBeam + b1, idk2 = cc * numBeam + b2;
            Hx1 += taskHashW1[idi2] - taskHashW1[idi1] + taskHashW1[idk1] - taskHashW1[idk2];
            Hx2 += taskHashW2[idi2] - taskHashW2[idi1] + taskHashW2[idk1] - taskHashW2[idk2];
            Hx3 += taskHashW3[idi2] - taskHashW3[idi1] + taskHashW3[idk1] - taskHashW3[idk2];
        }
    }
}

void local_search(double beginTime)
{
    int ts_depth   = 300;   // 连续无改进迭代数达到此值则结束本阶段
    int nonImprove = 0;

    // 基于解禁忌：持久哈希全程不清空；这里只重算当前(扰动后)解的哈希并登记
    compute_hash();
    mark_current();

    while (nonImprove < ts_depth)
    {
        if (((double)clock() - beginTime) / CLOCKS_PER_SEC > maxRunTime)
            break;

        int bestDelta = MININT_MOVE;     // 允许 0 / 负收益移动（禁忌搜索）
        int numBest   = 0;               // 与 bestDelta 持平的候选数量
        int kind = 0, a = -1, bb = -1, cc = -1;

        // ---- 重建 CSR 桶：按波束分组已服务任务（计数排序）----
        for (int b = 0; b <= numBeam; b++) bucketStart[b] = 0;
        for (int j = 0; j < numTask; j++)
            if (taskBeam[j] >= 0) bucketStart[taskBeam[j] + 1]++;
        for (int b = 0; b < numBeam; b++) bucketStart[b + 1] += bucketStart[b];
        for (int j = 0; j < numTask; j++)
            if (taskBeam[j] >= 0)
            {
                int b = taskBeam[j];
                bucketTask[bucketStart[b]++] = j;
            }
        for (int b = numBeam; b > 0; b--) bucketStart[b] = bucketStart[b - 1];
        bucketStart[0] = 0;

        // 注：compatBucket（按 源波束×候选模式 分组）只被 kind5 用，构建延迟到
        // 下方 if(bestDelta<=0) 块内，避免在有严格改进的迭代上白白重建 O(N·M)。

        // ---- 按任务类型列出兼容波束（波束编号升序）----
        for (int t = 0; t < numType; t++) typeCompatBeamCount[t] = 0;
        for (int b = 0; b < numBeam; b++)
            if (beamMode[b] >= 0)
                for (int t = 0; t < numType; t++)
                    if (typeCompatMode[t][beamMode[b]])
                        typeCompatBeam[t][typeCompatBeamCount[t]++] = b;

        //----------------------------------------------------------------
        // (1) 插入：将未服务任务装入兼容波束
        //----------------------------------------------------------------
        for (int j = 0; j < numTask; j++)
        {
            if (taskBeam[j] != -1) continue;
            int delta = taskProfit[j];

            int t = taskType[j] - 1;
            for (int p = 0; p < typeCompatBeamCount[t]; p++)
            {
                int b = typeCompatBeam[t][p];
                if (!feasible_on(j, b)) continue;
                if (delta < bestDelta)  continue;

                // 基于解禁忌：候选解（装入 j→b）是否已访问
                int idx = j * numBeam + b;
                if (is_visited(Hx1 + taskHashW1[idx],
                               Hx2 + taskHashW2[idx],
                               Hx3 + taskHashW3[idx]))
                    continue;
                record_candidate(1, j, b, -1, delta,
                                 bestDelta, numBest, kind, a, bb, cc);
            }
        }

        //----------------------------------------------------------------
        // (2) 仅删除：移除一个已服务任务
        //----------------------------------------------------------------
        for (int j = 0; j < numTask; j++)
        {
            if (taskBeam[j] < 0) continue;
            int delta = -taskProfit[j];
            if (delta < bestDelta) continue;

            // 基于解禁忌：候选解（移除 j）是否已访问
            int idx = j * numBeam + taskBeam[j];
            if (is_visited(Hx1 - taskHashW1[idx],
                           Hx2 - taskHashW2[idx],
                           Hx3 - taskHashW3[idx]))
                continue;

            record_candidate(2, j, -1, -1, delta,
                             bestDelta, numBest, kind, a, bb, cc);
        }

        //----------------------------------------------------------------
        // (3) 交换：未服务任务 i 替换同波束上已服务任务 k
        //----------------------------------------------------------------
        for (int i = 0; i < numTask; i++)
        {
            if (taskBeam[i] != -1) continue;
            int ti = taskType[i] - 1;
            for (int p = 0; p < typeCompatBeamCount[ti]; p++)
            {
                int b = typeCompatBeam[ti][p];

                int bwFree0 = remBW[b];
                int pwFree0 = remPW[b];
                for (int idx = bucketStart[b]; idx < bucketStart[b + 1]; idx++)
                {
                    int k = bucketTask[idx];

                    int delta = taskProfit[i] - taskProfit[k];
                    if (delta < bestDelta) continue;

                    // 去掉 k、装入 i 后的可行性（同波束 b）
                    int bwFree = bwFree0 + taskBWDemand[k];
                    int pwFree = pwFree0 + taskPWDemand[k];
                    if (taskBWDemand[i] > bwFree || taskPWDemand[i] > pwFree) continue;

                    // 基于解禁忌：候选解（i 入 b、k 出 b）是否已访问
                    int idi = i * numBeam + b;
                    int idk = k * numBeam + b;
                    if (is_visited(Hx1 + taskHashW1[idi] - taskHashW1[idk],
                                   Hx2 + taskHashW2[idi] - taskHashW2[idk],
                                   Hx3 + taskHashW3[idi] - taskHashW3[idk]))
                        continue;

                    record_candidate(3, i, b, k, delta,
                                     bestDelta, numBest, kind, a, bb, cc);
                }
            }
        }

        //----------------------------------------------------------------
        // (4) 模式翻转：切换某波束模式，仅驱逐与新模式不兼容的任务
        //     （此处不评估回填）
        //----------------------------------------------------------------
        for (int b = 0; b < numBeam; b++)
        {
            if (beamMode[b] < 0) continue;

            int usedPW = (beamPWCap[b] - modeBasePower[beamMode[b]]) - remPW[b];
            for (int m2 = 0; m2 < numMode; m2++)
            {
                if (m2 == beamMode[b]) continue;
                if (modeBasePower[m2] > beamPWCap[b]) continue;

                int loss = 0;
                for (int t = 0; t < numType; t++)
                    if (!typeCompatMode[t][m2]) loss += typeProfitSum[b][t];

                int evictPW = 0;
                long long evW1 = 0, evW2 = 0, evW3 = 0;   // 被驱逐任务的权重和
                for (int idx = bucketStart[b]; idx < bucketStart[b + 1]; idx++)
                {
                    int j = bucketTask[idx];
                    if (taskBeam[j] == b && !typeCompatMode[taskType[j] - 1][m2])
                    {
                        evictPW += taskPWDemand[j];
                        int idj = j * numBeam + b;
                        evW1 += taskHashW1[idj];
                        evW2 += taskHashW2[idj];
                        evW3 += taskHashW3[idj];
                    }
                }
                if (usedPW - evictPW > beamPWCap[b] - modeBasePower[m2]) continue;

                int delta = -loss;
                if (delta < bestDelta) continue;

                // 基于解禁忌：候选解（波束 b 模式 m1->m2 并驱逐不兼容任务）是否已访问
                int bidx  = b * numMode + m2;
                int bidx0 = b * numMode + beamMode[b];
                if (is_visited(Hx1 + beamHashW1[bidx] - beamHashW1[bidx0] - evW1,
                               Hx2 + beamHashW2[bidx] - beamHashW2[bidx0] - evW2,
                               Hx3 + beamHashW3[bidx] - beamHashW3[bidx0] - evW3))
                    continue;

                record_candidate(4, b, m2, -1, delta,
                                 bestDelta, numBest, kind, a, bb, cc);
            }
        }

        //----------------------------------------------------------------
        // (5) 跨波束交换（含动态虚拟槽）。
        //     real-dummy：b1 上已服务任务 i 迁到 b2（cc == -1）。
        //     real-real ：b1 上 i 与 b2 上 k 互换波束。
        //     所有 kind=5 移动的 delta=0，故仅当尚无严格改进移动时才探索。
        //----------------------------------------------------------------
        if (bestDelta <= 0)
        {
            rebuild_insert_filter();

            // ---- 按（源波束，候选目标模式）分组已服务任务（仅 kind5 需要，延迟到此构建）----
            for (int x = 0; x <= numBeam * numMode; x++) compatBucketStart[x] = 0;
            for (int j = 0; j < numTask; j++)
                if (taskBeam[j] >= 0)
                {
                    int b = taskBeam[j];
                    int t = taskType[j] - 1;
                    for (int m = 0; m < numMode; m++)
                        if (typeCompatMode[t][m])
                            compatBucketStart[b * numMode + m + 1]++;
                }
            for (int x = 0; x < numBeam * numMode; x++)
                compatBucketStart[x + 1] += compatBucketStart[x];
            for (int j = 0; j < numTask; j++)
                if (taskBeam[j] >= 0)
                {
                    int b = taskBeam[j];
                    int t = taskType[j] - 1;
                    for (int m = 0; m < numMode; m++)
                        if (typeCompatMode[t][m])
                            compatBucketTask[compatBucketStart[b * numMode + m]++] = j;
                }
            for (int x = numBeam * numMode; x > 0; x--)
                compatBucketStart[x] = compatBucketStart[x - 1];
            compatBucketStart[0] = 0;

            for (int i = 0; i < numTask; i++)
            {
                int b1 = taskBeam[i];
                if (b1 < 0) continue;

                int ti = taskType[i] - 1;
                for (int p = 0; p < typeCompatBeamCount[ti]; p++)
                {
                    int b2 = typeCompatBeam[ti][p];
                    if (b2 == b1) continue;

                    int delta = 0;

                    // real-dummy：i 从 b1 迁到 b2 当前空位，
                    // 且 b1 腾出 i 后须能装入某个未服务任务
                    if (taskBWDemand[i] <= remBW[b2] &&
                        taskPWDemand[i] <= remPW[b2])
                    {
                        if (beam_can_insert_unserved_with_rem(
                                b1,
                                remBW[b1] + taskBWDemand[i],
                                remPW[b1] + taskPWDemand[i]))
                        {
                            // 基于解禁忌：候选解（i 从 b1 迁到 b2）是否已访问
                            int idi2 = i * numBeam + b2;
                            int idi1 = i * numBeam + b1;
                            if (!is_visited(Hx1 + taskHashW1[idi2] - taskHashW1[idi1],
                                            Hx2 + taskHashW2[idi2] - taskHashW2[idi1],
                                            Hx3 + taskHashW3[idi2] - taskHashW3[idi1]))
                                record_candidate(5, i, b2, -1, delta,
                                                 bestDelta, numBest, kind, a, bb, cc);
                        }
                    }

                    // real-real：每对跨波束组合只枚举一次
                    if (b1 > b2) continue;

                    int m1 = beamMode[b1];
                    if (m1 < 0) continue;
                    int dstStart = compatBucketStart[b2 * numMode + m1];
                    int dstEnd   = compatBucketStart[b2 * numMode + m1 + 1];
                    for (int q = dstStart; q < dstEnd; q++)
                    {
                        int k = compatBucketTask[q];

                        if (taskBWDemand[k] > remBW[b1] + taskBWDemand[i]) continue;
                        if (taskPWDemand[k] > remPW[b1] + taskPWDemand[i]) continue;
                        if (taskBWDemand[i] > remBW[b2] + taskBWDemand[k]) continue;
                        if (taskPWDemand[i] > remPW[b2] + taskPWDemand[k]) continue;

                        // 仅当交换后至少一侧能装入未服务任务时才保留
                        int b1BW = remBW[b1] + taskBWDemand[i] - taskBWDemand[k];
                        int b1PW = remPW[b1] + taskPWDemand[i] - taskPWDemand[k];
                        int b2BW = remBW[b2] + taskBWDemand[k] - taskBWDemand[i];
                        int b2PW = remPW[b2] + taskPWDemand[k] - taskPWDemand[i];
                        if (!beam_can_insert_unserved_with_rem(b1, b1BW, b1PW) &&
                            !beam_can_insert_unserved_with_rem(b2, b2BW, b2PW))
                            continue;

                        // 基于解禁忌：候选解（i:b1->b2，k:b2->b1 互换）是否已访问
                        int idi2 = i * numBeam + b2;
                        int idi1 = i * numBeam + b1;
                        int idk1 = k * numBeam + b1;
                        int idk2 = k * numBeam + b2;
                        if (is_visited(
                                Hx1 + taskHashW1[idi2] - taskHashW1[idi1] + taskHashW1[idk1] - taskHashW1[idk2],
                                Hx2 + taskHashW2[idi2] - taskHashW2[idi1] + taskHashW2[idk1] - taskHashW2[idk2],
                                Hx3 + taskHashW3[idi2] - taskHashW3[idi1] + taskHashW3[idk1] - taskHashW3[idk2]))
                            continue;

                        record_candidate(5, i, b2, k, delta,
                                         bestDelta, numBest, kind, a, bb, cc);
                    }
                }
            }
        }

        // ---- 本迭代无可行非禁忌候选：本段已无路可走，结束本段交给扰动 ----
        // （解不变则下次扫描结果相同，空转无意义；持久禁忌下此情形会随表变满更常见）
        if (numBest == 0)
            break;

        // 基于解禁忌：用提交前状态增量更新 Hx（替代提交后全量 compute_hash）
        update_hash_for_move(kind, a, bb, cc);

        // ---- 提交所选移动 ----
        if (kind == 1)                            // 插入
        {
            add_task(a, bb);
            moveFreq[a]++;
        }
        else if (kind == 2)                       // 仅删除
        {
            remove_task(a);
            moveFreq[a]++;
        }
        else if (kind == 3)                       // 交换：装入 a，移出 cc
        {
            remove_task(cc);
            add_task(a, bb);
            moveFreq[a]++; moveFreq[cc]++;
        }
        else if (kind == 4)                       // 模式翻转：波束 a 切换到模式 bb
        {
            apply_mode_flip(a, bb);
        }
        else if (kind == 5)                       // 跨波束交换
        {
            int b1 = taskBeam[a];
            if (cc < 0)                           // real-dummy：将 a 迁到 bb
            {
                remove_task(a);
                add_task(a, bb);
                moveFreq[a]++;
            }
            else                                  // real-real：a 与 cc 互换
            {
                int b2 = taskBeam[cc];
                remove_task(a);
                remove_task(cc);
                add_task(a, b2);
                add_task(cc, b1);
                moveFreq[a]++; moveFreq[cc]++;
            }
        }

        // 基于解禁忌：登记新解为已访问（Hx 已在提交前增量更新）
        mark_current();

        if (totalProfit > bestProfit)
        {
            save_best(beginTime);
            nonImprove = 0;
        }
        else
            nonImprove++;
    }
}

//====================================================================
// 扰动（破坏 + 修复），用于 ILS 阶段
//
//  破坏：清空 30% 波束，优先选本阶段局部搜索中任务移动最少的波束
//  （moveFreq 低），将搜索推向未探索区域。被清空的波束模式重置为 -1。
//
//  修复（模式）：按顺序重定每个被清空波束的模式。枚举功率可行模式，
//  从当前未服务池试填，通常保留实际收益增量最大的模式；以 20% 概率
//  选用次优模式以增加多样性。
//
//  修复（任务）：按优先级贪心重分配仍未服务的任务，兼容模式从低基础
//  功耗到高依次尝试，同档内选剩余带宽最紧的波束。
//====================================================================
void perturb(int *tmpIdx, double *tmpVal)
{
    restore_best();
    if (numBeam == 0) return;

    int numDestroy = divStrength;             // 动态多样化强度（由 ils 维护）
    if (numDestroy < 1)        numDestroy = 1;
    if (numDestroy > numBeam)  numDestroy = numBeam;

    int    *perturbBeam = new int[numBeam];      // 被选为清空的波束
    int    *prevMode    = new int[numBeam];      // 破坏前各波束的模式
    double *beamFreqKey  = new double[numBeam];   // 各波束上已服务任务 moveFreq 之和
    int    *order        = new int[numBeam];
    for (int b = 0; b < numBeam; b++) { beamFreqKey[b] = 0.0; order[b] = b; }
    for (int j = 0; j < numTask; j++)
        if (taskBeam[j] >= 0) beamFreqKey[taskBeam[j]] += moveFreq[j];

    // 按频率升序排序：对键取负后用 qsort_desc
    for (int b = 0; b < numBeam; b++) beamFreqKey[b] = -beamFreqKey[b];
    qsort_desc(beamFreqKey, order, 0, numBeam - 1);   // order[]：低频在前

    // 从最少移动的波束中抽取（池宽至少为 numDestroy）
    int poolSize = 0;
    for (int b = 0; b < numBeam; b++) if (beamFreqKey[b] == 0.0) poolSize++;
    if (poolSize < numDestroy) poolSize = numDestroy;

    // 从低频池中随机抽取 numDestroy 个互不相同的波束
    for (int k = 0; k < numDestroy; k++)             // 部分 Fisher-Yates 洗牌
    {
        int r   = k + rand() % (poolSize - k);
        int tmp = order[k]; order[k] = order[r]; order[r] = tmp;
        perturbBeam[k] = order[k];
    }

    for (int k = 0; k < numDestroy; k++)
    {
        int b = perturbBeam[k];
        prevMode[k] = beamMode[b];                  // 记录原模式
        for (int j = 0; j < numTask; j++)
            if (taskBeam[j] == b) remove_task(j);
        beamMode[b] = -1;
        remBW[b]    = beamBWCap[b];
        remPW[b]    = 0;                             // 重定模式后再设置
    }

    //---- 修复模式：按实际收益试填并顺序提交 ----
    int *trialAdded  = new int[numTask];
    int *bestAdded   = new int[numTask];
    int *secondAdded = new int[numTask];
    for (int k = 0; k < numDestroy; k++)
    {
        int b = perturbBeam[k];
        int baseProfit = totalProfit;
        int bestM      = -1;
        int bestDelta  = MININT_MOVE;
        int nBestAdded = 0;
        int secondM      = -1;
        int secondDelta  = MININT_MOVE;
        int nSecondAdded = 0;

        for (int m = 0; m < numMode; m++)
        {
            if (modeBasePower[m] > beamPWCap[b]) continue;       // 功率不可行

            beamMode[b] = m;
            remBW[b]    = beamBWCap[b];
            remPW[b]    = beamPWCap[b] - modeBasePower[m];

            int nFree = 0;
            for (int j = 0; j < numTask; j++)
                if (taskBeam[j] == -1)
                {
                    tmpIdx[nFree] = j;
                    tmpVal[nFree] = task_priority(j);
                    nFree++;
                }
            if (nFree > 0) qsort_desc(tmpVal, tmpIdx, 0, nFree - 1);

            int nTrialAdded = 0;
            for (int ki = 0; ki < nFree; ki++)
            {
                int j = tmpIdx[ki];
                if (!feasible_on(j, b)) continue;
                add_task(j, b);
                trialAdded[nTrialAdded++] = j;
            }

            int delta = totalProfit - baseProfit;
            if (delta > bestDelta)
            {
                secondDelta  = bestDelta;
                secondM      = bestM;
                nSecondAdded = nBestAdded;
                for (int a = 0; a < nBestAdded; a++) secondAdded[a] = bestAdded[a];

                bestDelta  = delta;
                bestM      = m;
                nBestAdded = nTrialAdded;
                for (int a = 0; a < nTrialAdded; a++) bestAdded[a] = trialAdded[a];
            }
            else if (delta > secondDelta)
            {
                secondDelta  = delta;
                secondM      = m;
                nSecondAdded = nTrialAdded;
                for (int a = 0; a < nTrialAdded; a++) secondAdded[a] = trialAdded[a];
            }

            for (int a = nTrialAdded - 1; a >= 0; a--) remove_task(trialAdded[a]);
        }

        if (bestM < 0)                       // 无功率可行模式：选基础功耗最低者
        {
            for (int m = 0; m < numMode; m++)
                if (bestM < 0 || modeBasePower[m] < modeBasePower[bestM]) bestM = m;
            if (bestM < 0) bestM = prevMode[k];   // 仅一种模式时无法避免
        }

        int chosenM      = bestM;
        int *chosenAdded = bestAdded;
        int nChosenAdded = nBestAdded;
        if (secondM >= 0 && rand() % 5 == 0)  // 20%：修复多样化
        {
            chosenM      = secondM;
            chosenAdded  = secondAdded;
            nChosenAdded = nSecondAdded;
        }

        beamMode[b] = chosenM;
        remBW[b] = beamBWCap[b];
        remPW[b] = beamPWCap[b] - modeBasePower[chosenM];

        for (int a = 0; a < nChosenAdded; a++)
            if (taskBeam[chosenAdded[a]] == -1 && feasible_on(chosenAdded[a], b))
                add_task(chosenAdded[a], b);
    }
    delete[] trialAdded;
    delete[] bestAdded;
    delete[] secondAdded;

    //---- 修复任务：按优先级贪心最佳适配 ----
    int nFree = 0;
    for (int j = 0; j < numTask; j++)
        if (taskBeam[j] == -1)
        {
            tmpIdx[nFree] = j;
            tmpVal[nFree] = task_priority(j);
            nFree++;
        }
    if (nFree > 0) qsort_desc(tmpVal, tmpIdx, 0, nFree - 1);

    int *modeOrder = new int[numMode];
    for (int k = 0; k < nFree; k++)
    {
        int j = tmpIdx[k];

        int bestBeam   = -1;
        int bestLeftBW = -1;

        // 任务 j 的兼容模式，按基础功耗升序（插入排序）
        int t = taskType[j] - 1;
        int numCompatMode = 0;
        for (int m = 0; m < numMode; m++)
        {
            if (!typeCompatMode[t][m]) continue;
            int pos = numCompatMode;
            while (pos > 0 && modeBasePower[m] < modeBasePower[modeOrder[pos - 1]])
            {
                modeOrder[pos] = modeOrder[pos - 1];
                pos--;
            }
            modeOrder[pos] = m;
            numCompatMode++;
        }

        for (int km = 0; km < numCompatMode && bestBeam < 0; km++)
        {
            int m = modeOrder[km];
            for (int b = 0; b < numBeam; b++)
            {
                if (beamMode[b] != m) continue;
                if (!feasible_on(j, b)) continue;
                int leftover = remBW[b] - taskBWDemand[j];
                if (bestBeam < 0 || leftover < bestLeftBW)
                {
                    bestBeam   = b;
                    bestLeftBW = leftover;
                }
            }
        }
        if (bestBeam >= 0) add_task(j, bestBeam);
    }
    delete[] modeOrder;
    delete[] perturbBeam;
    delete[] prevMode;
    delete[] beamFreqKey;
    delete[] order;
}

//====================================================================
// 迭代局部搜索主流程
//====================================================================
void ils()
{
    double beginTime = (double)clock();

    int    *tmpIdx = new int   [numTask];
    double *tmpVal = new double[numTask];

    greedy_init();
    alloc_search();

    // ILS 开始前全量重建 typeProfitSum，保证与当前解一致
    for (int b = 0; b < numBeam; b++)
        for (int t = 0; t < numType; t++) typeProfitSum[b][t] = 0;
    for (int j = 0; j < numTask; j++)
        if (taskBeam[j] >= 0)
            typeProfitSum[taskBeam[j]][taskType[j] - 1] += taskProfit[j];

    globalProfit = -1;

    // 动态多样化强度（学 MNSB-TS）：弱强度起步，刷新全局即回弱，连续停滞则渐强
    int weakDestroy = (int)(0.10 * numBeam + 0.5);
    if (weakDestroy < 1) weakDestroy = 1;
    divStrength = weakDestroy;

    double runTime = 0.0;
    while (runTime < maxRunTime)
    {
        for (int j = 0; j < numTask; j++) moveFreq[j] = 0;

        save_best(beginTime);          // 本阶段最优从当前解出发
        local_search(beginTime);

        if (bestProfit > globalProfit)
        {
            save_global_from_best();
            divStrength = weakDestroy;          // 逃出(刷新全局)→ 回到弱扰动
        }
        else
        {
            divStrength++;                       // 没逃出 → 破坏强度渐增
            if (divStrength > numBeam) divStrength = weakDestroy;  // 到顶后回弱，循环覆盖各尺度
        }

        perturb(tmpIdx, tmpVal);

        runTime = ((double)clock() - beginTime) / CLOCKS_PER_SEC;
    }

    // 将全局最优作为最终解，供校验与输出
    restore_global();
    totalProfit = globalProfit;

    delete[] tmpIdx;
    delete[] tmpVal;
}

//====================================================================
// 最终校验与输出
//====================================================================
void check_solution()
{
    int  checkProfit = 0;
    int  ok          = 1;

    int *usedBW = new int[numBeam];
    int *usedPW = new int[numBeam];
    for (int b = 0; b < numBeam; b++) { usedBW[b] = 0; usedPW[b] = 0; }

    for (int j = 0; j < numTask; j++)
    {
        if (taskBeam[j] < 0) continue;
        int b = taskBeam[j];
        int t = taskType[j] - 1;

        if (!typeCompatMode[t][beamMode[b]])
        {
            cout << "ERROR: task " << j << " (type " << (t + 1)
                 << ") incompatible with beam " << b
                 << " mode " << modeName[beamMode[b]] << endl;
            ok = 0;
        }

        usedBW[b]   += taskBWDemand[j];
        usedPW[b]   += taskPWDemand[j];
        checkProfit += taskProfit[j];
    }

    for (int b = 0; b < numBeam; b++)
    {
        int pw_avail = beamPWCap[b] - modeBasePower[beamMode[b]];

        if (usedBW[b] > beamBWCap[b])
        {
            cout << "ERROR: beam " << b << " BW overflow: used=" << usedBW[b]
                 << " cap=" << beamBWCap[b] << endl;
            ok = 0;
        }
        if (usedPW[b] > pw_avail)
        {
            cout << "ERROR: beam " << b << " PW overflow: used=" << usedPW[b]
                 << " avail=" << pw_avail
                 << " (cap=" << beamPWCap[b]
                 << " base=" << modeBasePower[beamMode[b]] << ")" << endl;
            ok = 0;
        }
    }

    if (checkProfit != totalProfit)
    {
        cout << "ERROR: profit mismatch: recorded=" << totalProfit
             << " recomputed=" << checkProfit << endl;
        ok = 0;
    }

    if (ok)
        cout << "Solution verified OK." << endl;

    delete[] usedBW;
    delete[] usedPW;
}

void print_solution()
{
    cout << "\n=== Solution Summary ===" << endl;
    cout << "Total profit : " << totalProfit << endl;
    cout << "Best time    : " << globalBestTime << " s" << endl;

    int numServed = 0;
    for (int j = 0; j < numTask; j++)
        if (taskBeam[j] >= 0) numServed++;
    cout << "Tasks served : " << numServed << " / " << numTask << endl;

    cout << "Beam details :" << endl;
    for (int b = 0; b < numBeam; b++)
    {
        int usedBW = beamBWCap[b] - remBW[b];
        int usedPW = beamPWCap[b] - modeBasePower[beamMode[b]] - remPW[b];
        int capPW  = beamPWCap[b] - modeBasePower[beamMode[b]];
        cout << "  Beam " << b
             << "  mode=" << modeName[beamMode[b]]
             << "  BW=" << usedBW << "/" << beamBWCap[b]
             << "  PW=" << usedPW << "/" << capPW
             << endl;
    }
}

void free_memory()
{
    for (int m = 0; m < numMode; m++)  delete[] modeName[m];
    delete[] modeName;
    delete[] modeBasePower;

    for (int t = 0; t < numType; t++)  delete[] typeCompatMode[t];
    delete[] typeCompatMode;

    delete[] beamBWCap;
    delete[] beamPWCap;
    delete[] taskType;
    delete[] taskProfit;
    delete[] taskBWDemand;
    delete[] taskPWDemand;

    delete[] beamMode;
    delete[] taskBeam;
    delete[] remBW;
    delete[] remPW;

    delete[] bestBeamMode;
    delete[] bestTaskBeam;
    delete[] globalBeamMode;
    delete[] globalTaskBeam;
    delete[] moveFreq;
    delete[] bucketTask;
    delete[] bucketStart;
    delete[] compatBucketTask;
    delete[] compatBucketStart;
    for (int t = 0; t < numType; t++) delete[] typeCompatBeam[t];
    delete[] typeCompatBeam;
    delete[] typeCompatBeamCount;
    for (int m = 0; m < numMode; m++) delete[] insertMinPW[m];
    delete[] insertMinPW;

    for (int b = 0; b < numBeam; b++) delete[] typeProfitSum[b];
    delete[] typeProfitSum;

    delete[] taskHashW1;
    delete[] taskHashW2;
    delete[] taskHashW3;
    delete[] beamHashW1;
    delete[] beamHashW2;
    delete[] beamHashW3;
    delete[] hashTab1;
    delete[] hashTab2;
    delete[] hashTab3;
}

//====================================================================
// 主程序
//====================================================================
int main(int argc, char **argv)
{
    if (argc < 3)
    {
        cout << "Usage: " << argv[0] << " <instance_file> <seed> [time_limit_seconds]" << endl;
        return 1;
    }
    instanceName = argv[1];
    seed         = atoi(argv[2]);
    srand(seed);

    maxRunTime = 600.0;          // 默认时间上限（秒）
    if (argc >= 4) maxRunTime = atof(argv[3]);

    read_instance();

    double t0 = (double)clock();
    ils();
    double elapsed = ((double)clock() - t0) / CLOCKS_PER_SEC;

    check_solution();
    print_solution();

    // 批处理脚本解析用的摘要行（best_profit / best_time / elapsed_time）
    cout << "ILS done.  bestProfit=" << globalProfit
         << "  bestTime=" << globalBestTime << " s" << endl;
    cout << "Elapsed time: " << elapsed << " s" << endl;

    free_memory();
    return 0;
}
