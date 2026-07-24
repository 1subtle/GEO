#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAXLINE    65536   // 算例文件单行最大字节数
#define MAXTYPE    20      // 任务类型最大数量
#define MAXNAMELEN 32      // 模式名称最大长度
#define MININT_MOVE (-2000000000)   // 禁忌搜索增量的负无穷初值
#define ORPHAN_EC_MAXLEN 4 // 未服务任务孤儿链最多连续驱逐次数
#define ORPHAN_EC_START_CAP 4  // 每轮最多尝试的链头数

using namespace std;

char *instanceName;

int numBeam;          // B：波束数量
int numTask;          // N：任务数量
int numMode;          // M：服务模式数量
int numType;          // T：任务类型数量

char **modeName;      // modeName[m]：模式 m 的名称
int  *modeBasePower;  // modeBasePower[m]：模式 m 的固定基础功率
int **typeCompatMode; // typeCompatMode[t][m]：类型 t+1 可使用模式 m 时为 1

int *beamBWCap;       // beamBWCap[b]：波束 b 的带宽容量
int *beamPWCap;       // beamPWCap[b]：波束 b 的功率容量

int *taskType;        // taskType[j]：任务 j 的类型，从 1 开始编号
int *taskProfit;      // taskProfit[j]：任务 j 的收益
int *taskBWDemand;    // taskBWDemand[j]：任务 j 的带宽需求
int *taskPWDemand;    // taskPWDemand[j]：任务 j 的功率需求

int *beamMode;   // beamMode[b]：波束 b 选择的模式下标
int *taskBeam;   // taskBeam[j]：服务任务 j 的波束，-1 表示未服务
int *remBW;      // remBW[b]：波束 b 的剩余带宽
int *remPW;      // remPW[b]：波束 b 扣除模式基础功率后的剩余功率
int  totalProfit;

double maxRunTime;   // 时间上限，单位为秒
double bestTime;     // 找到当前阶段最好解时的运行时间
int    seed;         // 随机种子

int *bestBeamMode;   // 当前 ILS 阶段最好解的模式快照
int *bestTaskBeam;
int  bestProfit;

int *globalBeamMode; // 所有 ILS 阶段全局最好解的模式快照
int *globalTaskBeam;
int  globalProfit;
double globalBestTime;

int *moveFreq;       // moveFreq[j]：任务 j 的移动次数，供扰动使用
int *tabuUntil;      // tabuUntil[j]：任务 j 的禁忌截止迭代号
int  tabuIter;       // 单次禁忌局部搜索内部的迭代计数

int  *tabuBeamMode;  // tabuBeamMode[b]：波束 b 模式的禁忌截止迭代号
int **typeProfitSum; // typeProfitSum[b][t]：波束 b 上类型 t 的已服务任务收益和

// 按波束分组存放已服务任务的 CSR 桶，每次 local_search 迭代重建，
// 使交换邻域能够跳过无关任务。
int *bucketTask;     // 长度 numTask：按波束连续存放的已服务任务编号
int *bucketStart;    // 长度 numBeam+1：波束 b 的区间为 [bucketStart[b], bucketStart[b+1])
int *mixedMinBWFree;
int **mixedMinPWFree;
int **insertMinPW;   // insertMinPW[m][bw]：模式 m 下、BW 不超过 bw 的未服务任务最小 PW
int  maxBeamBWCap;

// 孤儿收益链的工作区。链头是未服务任务；每一步用当前孤儿替换一个
// 已服务任务，新被驱逐任务成为下一代孤儿。
int *orphanWorkBW, *orphanWorkPW;
int *orphanTask, *orphanDest;
int *orphanBestTask, *orphanBestDest;
int *orphanStartMark;
int orphanBestCount;
int orphanBestDeltaProfit;
double orphanBestDeltaEval;
int orphanBestAspired;
double orphanBaseOverBW, orphanBaseOverPW;

// ---- 全程累计的性能统计计数器 ----
long long g_lsIters   = 0;   // local_search 内层总迭代数
long long g_lsMoves   = 0;   // 实际执行动作的迭代数
long long g_lsStall   = 0;   // 没有可接受候选的迭代数
long long g_tabuBlock = 0;   // 因禁忌且不满足特赦而跳过的候选数
long long g_aspire    = 0;   // 通过特赦接受的禁忌动作数
double    g_lsTime    = 0.0; // local_search 累计 CPU 时间
double    g_perturbTime = 0.0; // perturb 累计 CPU 时间
long long g_applyInsert = 0; // 已执行的插入动作数
long long g_applyRemove = 0; // 已执行的纯移除动作数
long long g_applySwap   = 0; // 已执行的同波束交换动作数
long long g_applyFlip   = 0; // 已执行的模式翻转动作数
long long g_applyCross  = 0; // 已执行的跨波束交换动作数
long long g_crossCand   = 0; // 可行的跨波束交换候选数
long long g_crossRelocCand = 0; // 第五邻域实任务-虚拟空位候选数
long long g_crossSwapCand  = 0; // 第五邻域实任务-实任务候选数
long long g_crossFiltered  = 0; // 被开放槽位过滤器拒绝的第五邻域候选数
long long g_crossRelocApplied = 0; // 已执行的实任务-虚拟空位动作数
long long g_crossSwapApplied  = 0; // 已执行的实任务-实任务动作数
long long g_orphanCalls = 0;       // 孤儿链搜索调用次数
long long g_orphanTrials = 0;      // 孤儿链完整叶子数
long long g_orphanSteps = 0;       // 孤儿链驱逐步数
long long g_orphanApplied = 0;     // 已执行的孤儿链数
long long g_orphanRelocLeaf = 0;   // 孤儿链末端重定位叶子数
long long g_orphanDropLeaf = 0;    // 孤儿链末端不服务叶子数
long long g_flipCand    = 0; // 已评价的可行模式翻转候选数
long long g_flipApplied = 0; // 已执行的模式翻转动作数
long long g_flipApplyToMode[16] = {0}; // 按目标模式统计已执行翻转数

// ---- 可行域/不可行域混合搜索状态 -----------------------
double mixedPhiBW = 1000.0;   // 归一化 BW 超载的惩罚权重
double mixedPhiPW = 1000.0;   // 归一化 PW 超载的惩罚权重
double mixedRhoBW = 0.15;     // 混合搜索当前允许的 BW 超载半径
double mixedRhoPW = 0.20;     // 混合搜索当前允许的 PW 超载半径
double mixedRhoMinBW = 0.05, mixedRhoMaxBW = 0.30;
double mixedRhoMinPW = 0.08, mixedRhoMaxPW = 0.30;
long long g_mixedIters = 0;   // 混合搜索总迭代数
long long g_mixedMoves = 0;   // 混合搜索已执行动作数
long long g_mixedFeas  = 0;   // 当前解可行的迭代数
long long g_mixedInfeas = 0;  // 当前解不可行的迭代数
long long g_mixedBestUpdate = 0; // 混合搜索中可行最好解的更新次数
long long g_mixedRhoExpand = 0;  // 自适应放宽半径扩张次数
long long g_mixedRhoShrink = 0;  // 自适应放宽半径收缩次数
double g_mixedMaxOver = 0.0;  // 访问过的最大归一化总超载
double g_mixedMaxRhoBW = 0.0; // 使用过的最大 BW 放宽半径
double g_mixedMaxRhoPW = 0.0; // 使用过的最大 PW 放宽半径
double g_mixedTime = 0.0;     // 混合搜索累计 CPU 时间

static char g_line[MAXLINE];

int find_mode_index(const char *name)
{
    for (int m = 0; m < numMode; m++)
        if (strcmp(modeName[m], name) == 0)
            return m;
    return -1;
}

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

    // 文件头示例：B=5 N=50 M=3 Cbw=1578 Cpw=503
    FIC.getline(g_line, MAXLINE);
    int cbw_total, cpw_total;
    sscanf(g_line, "B=%d N=%d M=%d Cbw=%d Cpw=%d",
           &numBeam, &numTask, &numMode, &cbw_total, &cpw_total);

    modeName     = new char *[numMode];
    for (int m = 0; m < numMode; m++)
        modeName[m] = new char[MAXNAMELEN];
    modeBasePower = new int[numMode];

    // 跳过空行、Modes 标题和 mode base_power 表头
    FIC.getline(g_line, MAXLINE);
    FIC.getline(g_line, MAXLINE);
    FIC.getline(g_line, MAXLINE);

    // 读取 M 行模式数据，例如 "G 10"
    for (int m = 0; m < numMode; m++)
    {
        FIC.getline(g_line, MAXLINE);
        char tmp[MAXNAMELEN];
        sscanf(g_line, "%s %d", tmp, &modeBasePower[m]);
        strcpy(modeName[m], tmp);
    }

    // 跳过空行、Task types 标题和 type name modes 表头
    FIC.getline(g_line, MAXLINE);
    FIC.getline(g_line, MAXLINE);
    FIC.getline(g_line, MAXLINE);

    // 读取任务类型行，直到遇到空行
    static char typeLines[MAXTYPE][256];
    numType = 0;
    while (FIC.getline(g_line, MAXLINE))
    {
        if ((int)strlen(g_line) == 0) break;   // 空行表示本节结束
        strncpy(typeLines[numType], g_line, 255);
        typeLines[numType][255] = '\0';
        numType++;
    }
    // 上面的 getline 已经读取消耗了空行

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
    FIC.getline(g_line, MAXLINE);                         // 数据行结尾
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

    // 任务功率需求，最后一节不要求尾随空行
    FIC.getline(g_line, MAXLINE);                         // 表头
    for (int j = 0; j < numTask; j++) FIC >> taskPWDemand[j];

    FIC.close();

    // 输出算例摘要
    cout << "Instance: " << instanceName << endl;
    cout << "B=" << numBeam << "  N=" << numTask
         << "  M=" << numMode << "  T=" << numType << endl;
    for (int m = 0; m < numMode; m++)
        cout << "  Mode[" << m << "] " << modeName[m]
             << "  base_pw=" << modeBasePower[m] << endl;
    cout << "Task type compatibility:" << endl;
    for (int t = 0; t < numType; t++)
    {
        cout << "  Type " << (t + 1) << ":";
        for (int m = 0; m < numMode; m++)
            if (typeCompatMode[t][m])
                cout << " " << modeName[m];
        cout << endl;
    }
}

//--------------------------------------------------------------------
// task_score：任务在波束 b、模式 m 下单位归一化资源占用的收益。
// 带宽和功率分别按波束容量归一化，其中任务功率容量等于波束功率容量
// 减去模式基础功率。占用任一资源比例越大的任务，得分越低。
// 该得分只用于贪心排序和优先级，不作为目标函数。
//--------------------------------------------------------------------
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

//--------------------------------------------------------------------
// task_priority：任务 j 在所有模式兼容波束上能够取得的最高得分。
// 使用每个波束的任务总功率容量，而非动态剩余功率，因此在一次贪心
// 扫描期间该排序键保持不变。若没有兼容波束则返回 -1。
//--------------------------------------------------------------------
double task_priority(int j)
{
    int    t    = taskType[j] - 1;
    double best = -1.0;
    for (int b = 0; b < numBeam; b++)
    {
        if (beamMode[b] < 0)                  continue;
        if (!typeCompatMode[t][beamMode[b]])  continue;
        double s = task_score(j, b, beamMode[b]);
        if (s > best) best = s;
    }
    return best;
}

void greedy_init()
{
    beamMode = new int[numBeam];
    taskBeam = new int[numTask];
    remBW    = new int[numBeam];
    remPW    = new int[numBeam];

    for (int j = 0; j < numTask; j++) taskBeam[j] = -1;

    int    *tmpIdx = new int   [numTask];
    double *tmpVal = new double[numTask];
    int    *trialAdded = new int[numTask];
    int    *bestAdded  = new int[numTask];

    totalProfit = 0;

    //================================================================
    // 顺序提交贪心，与扰动中的模式修复思路相同：
    // 从当前未服务任务池试填并决定一个波束的模式，提交装填结果后
    // 再处理下一个波束，使不同波束竞争任务。
    //================================================================
    for (int b = 0; b < numBeam; b++)
    {
        int bestM      = -1;
        int bestGain   = -1;
        int nBestAdded = 0;

        for (int m = 0; m < numMode; m++)
        {
            if (modeBasePower[m] > beamPWCap[b]) continue;   // 基础功率不可行

            // 按密度得分排序当前未服务且模式兼容的任务
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

        if (bestM < 0)   // 没有功率可行模式时选择基础功率最低的模式
        {
            for (int m = 0; m < numMode; m++)
                if (bestM < 0 || modeBasePower[m] < modeBasePower[bestM]) bestM = m;
            nBestAdded = 0;
        }

        // 提交所选模式及装填结果，此时 typeProfitSum 尚未分配
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

    cout << "Greedy init done.  totalProfit=" << totalProfit << endl;
}

//--------------------------------------------------------------------
// 输出波束模式分布和各类型服务数量，供运行诊断。
//--------------------------------------------------------------------
void probe_modes(const char *tag)
{
    int *mc = new int[numMode]; for (int m = 0; m < numMode; m++) mc[m] = 0;
    for (int b = 0; b < numBeam; b++) if (beamMode[b] >= 0) mc[beamMode[b]]++;
    int *st = new int[numType]; for (int t = 0; t < numType; t++) st[t] = 0;
    int served = 0;
    for (int j = 0; j < numTask; j++)
        if (taskBeam[j] >= 0) { st[taskType[j]-1]++; served++; }
    cout << "[PROBE " << tag << "] profit=" << totalProfit << " served=" << served << "  modes:";
    for (int m = 0; m < numMode; m++) cout << " " << modeName[m] << "=" << mc[m];
    cout << "  served_by_type:";
    for (int t = 0; t < numType; t++) cout << " T" << (t+1) << "=" << st[t];
    cout << endl;
    delete[] mc; delete[] st;
}

//--------------------------------------------------------------------
// feasible_on：判断任务 j 在当前剩余资源下能否分配到波束 b。
// 必须同时满足模式兼容、带宽充足和功率充足。
//--------------------------------------------------------------------
int feasible_on(int j, int b)
{
    if (beamMode[b] < 0)                 return 0;   // 模式尚未确定
    int t = taskType[j] - 1;
    if (!typeCompatMode[t][beamMode[b]]) return 0;
    if (taskBWDemand[j] > remBW[b])      return 0;
    if (taskPWDemand[j] > remPW[b])      return 0;
    return 1;
}

double beam_bw_over_with_rem(int b, int bwFree)
{
    if (bwFree >= 0) return 0.0;
    int denom = beamBWCap[b];
    if (denom < 1) denom = 1;
    return (double)(-bwFree) / denom;
}

double beam_pw_over_with_rem_mode(int b, int m, int pwFree)
{
    if (pwFree >= 0) return 0.0;
    int denom = beamPWCap[b] - modeBasePower[m];
    if (denom < 1) denom = 1;
    return (double)(-pwFree) / denom;
}

double beam_pw_over_with_rem(int b, int pwFree)
{
    return beam_pw_over_with_rem_mode(b, beamMode[b], pwFree);
}

double beam_over_with_rem_mode(int b, int m, int bwFree, int pwFree)
{
    return beam_bw_over_with_rem(b, bwFree) +
           beam_pw_over_with_rem_mode(b, m, pwFree);
}

double beam_over_with_rem(int b, int bwFree, int pwFree)
{
    return beam_over_with_rem_mode(b, beamMode[b], bwFree, pwFree);
}

double beam_over(int b)
{
    return beam_over_with_rem(b, remBW[b], remPW[b]);
}

double beam_bw_over(int b)
{
    return beam_bw_over_with_rem(b, remBW[b]);
}

double beam_pw_over(int b)
{
    return beam_pw_over_with_rem(b, remPW[b]);
}

double total_over()
{
    double over = 0.0;
    for (int b = 0; b < numBeam; b++)
        if (beamMode[b] >= 0) over += beam_over(b);
    return over;
}

void total_over_parts(double &bwOver, double &pwOver)
{
    bwOver = 0.0;
    pwOver = 0.0;
    for (int b = 0; b < numBeam; b++)
    {
        if (beamMode[b] < 0) continue;
        bwOver += beam_bw_over(b);
        pwOver += beam_pw_over(b);
    }
}

double clamp_double(double x, double lo, double hi)
{
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

double current_resource_tightness()
{
    double usedBW = 0.0, capBW = 0.0;
    double usedPW = 0.0, capPW = 0.0;

    for (int b = 0; b < numBeam; b++)
    {
        if (beamMode[b] < 0) continue;

        capBW  += beamBWCap[b];
        usedBW += beamBWCap[b] - remBW[b];

        int pwCap = beamPWCap[b] - modeBasePower[beamMode[b]];
        if (pwCap < 1) pwCap = 1;
        capPW  += pwCap;
        usedPW += pwCap - remPW[b];
    }

    double bwTight = capBW > 0.0 ? usedBW / capBW : 0.0;
    double pwTight = capPW > 0.0 ? usedPW / capPW : 0.0;
    return bwTight > pwTight ? bwTight : pwTight;
}

void set_mixed_rho_from_instance()
{
    double density = numBeam > 0 ? (double)numTask / numBeam : 0.0;
    double densityAdj = (density - 10.0) * 0.004;
    densityAdj = clamp_double(densityAdj, 0.0, 0.08);

    double tightness = current_resource_tightness();
    double tightAdj = 0.0;
    if (tightness > 0.92)      tightAdj = 0.03;
    else if (tightness > 0.85) tightAdj = 0.015;
    else if (tightness < 0.70) tightAdj = -0.02;

    double targetBW = clamp_double(0.12 + densityAdj + tightAdj,
                                   mixedRhoMinBW, mixedRhoMaxBW);
    double targetPW = clamp_double(targetBW + 0.04,
                                   mixedRhoMinPW, mixedRhoMaxPW);

    if (mixedRhoBW < targetBW) mixedRhoBW = targetBW;
    if (mixedRhoPW < targetPW) mixedRhoPW = targetPW;
    if (mixedRhoBW > g_mixedMaxRhoBW) g_mixedMaxRhoBW = mixedRhoBW;
    if (mixedRhoPW > g_mixedMaxRhoPW) g_mixedMaxRhoPW = mixedRhoPW;
}

void expand_mixed_rho()
{
    double oldBW = mixedRhoBW;
    double oldPW = mixedRhoPW;
    mixedRhoBW = clamp_double(mixedRhoBW * 1.10 + 0.01,
                              mixedRhoMinBW, mixedRhoMaxBW);
    mixedRhoPW = clamp_double(mixedRhoPW * 1.10 + 0.01,
                              mixedRhoMinPW, mixedRhoMaxPW);
    if (mixedRhoBW > g_mixedMaxRhoBW) g_mixedMaxRhoBW = mixedRhoBW;
    if (mixedRhoPW > g_mixedMaxRhoPW) g_mixedMaxRhoPW = mixedRhoPW;
    if (mixedRhoBW != oldBW || mixedRhoPW != oldPW) g_mixedRhoExpand++;
}

void shrink_mixed_rho()
{
    double oldBW = mixedRhoBW;
    double oldPW = mixedRhoPW;
    mixedRhoBW = clamp_double(mixedRhoBW * 0.85,
                              mixedRhoMinBW, mixedRhoMaxBW);
    mixedRhoPW = clamp_double(mixedRhoPW * 0.85,
                              mixedRhoMinPW, mixedRhoMaxPW);
    if (mixedRhoBW != oldBW || mixedRhoPW != oldPW) g_mixedRhoShrink++;
}

int max_overload_for_rho(int denom, double rho)
{
    if (denom < 1) denom = 1;
    int limit = (int)(rho * denom);
    while ((double)(limit + 1) / denom <= rho) limit++;
    while (limit > 0 && (double)limit / denom > rho) limit--;
    return limit;
}

void rebuild_relaxed_limits()
{
    for (int b = 0; b < numBeam; b++)
    {
        mixedMinBWFree[b] = -max_overload_for_rho(beamBWCap[b], mixedRhoBW);
        for (int m = 0; m < numMode; m++)
            mixedMinPWFree[b][m] =
                -max_overload_for_rho(beamPWCap[b] - modeBasePower[m], mixedRhoPW);
    }
}

int relaxed_rem_ok_mode(int b, int m, int bwFree, int pwFree)
{
    return bwFree >= mixedMinBWFree[b] && pwFree >= mixedMinPWFree[b][m];
}

int relaxed_rem_ok(int b, int bwFree, int pwFree)
{
    return relaxed_rem_ok_mode(b, beamMode[b], bwFree, pwFree);
}

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

//--------------------------------------------------------------------
// add_task / remove_task：分配或取消分配任务，并同步维护 remBW、
// remPW、taskBeam 和 totalProfit；这里不改变模式。
//--------------------------------------------------------------------
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

int tabu_tenure()
{
    return 5 + rand() % 6 + numBeam / 10;
}

// 按 RSOA 论文使用 MinP + rand(1..20) 计算阈值允许的收益下降幅度。
int threshold_delta()
{
    int minProfit = taskProfit[0];
    for (int j = 1; j < numTask; j++)
        if (taskProfit[j] < minProfit) minProfit = taskProfit[j];
    return minProfit + 1 + rand() % 20;
}

//--------------------------------------------------------------------
// save_best：将当前解保存到阶段最好解数组。
//--------------------------------------------------------------------
void save_best(double beginTime)
{
    bestTime    = ((double)clock() - beginTime) / CLOCKS_PER_SEC;
    bestProfit  = totalProfit;
    for (int b = 0; b < numBeam; b++) bestBeamMode[b] = beamMode[b];
    for (int j = 0; j < numTask; j++) bestTaskBeam[j] = taskBeam[j];
}

//--------------------------------------------------------------------
// restore_best：把阶段最好解恢复为当前解，并重新计算 remBW、
// remPW 和 totalProfit。
//--------------------------------------------------------------------
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
    bestBeamMode = new int[numBeam];
    bestTaskBeam = new int[numTask];
    globalBeamMode = new int[numBeam];
    globalTaskBeam = new int[numTask];
    moveFreq     = new int[numTask];
    tabuUntil    = new int[numTask];
    bucketTask   = new int[numTask];
    bucketStart  = new int[numBeam + 1];
    mixedMinBWFree = new int[numBeam];
    mixedMinPWFree = new int*[numBeam];
    for (int b = 0; b < numBeam; b++) mixedMinPWFree[b] = new int[numMode];
    orphanWorkBW = new int[numBeam];
    orphanWorkPW = new int[numBeam];
    orphanTask = new int[numTask];
    orphanDest = new int[numTask];
    orphanBestTask = new int[numTask];
    orphanBestDest = new int[numTask];
    orphanStartMark = new int[numTask];
    for (int j = 0; j < numTask; j++) orphanStartMark[j] = 0;

    maxBeamBWCap = 0;
    for (int b = 0; b < numBeam; b++)
        if (beamBWCap[b] > maxBeamBWCap) maxBeamBWCap = beamBWCap[b];
    insertMinPW = new int*[numMode];
    for (int m = 0; m < numMode; m++)
        insertMinPW[m] = new int[maxBeamBWCap + 1];

    tabuBeamMode = new int[numBeam];
    typeProfitSum = new int*[numBeam];
    for (int b = 0; b < numBeam; b++)
    {
        typeProfitSum[b] = new int[numType];
        for (int t = 0; t < numType; t++) typeProfitSum[b][t] = 0;
    }
}

void record_candidate(int kind, int a, int b, int c, int delta,
                      int &bestDelta, int &numBest,
                      int &chosenKind, int &chosenA, int &chosenB, int &chosenC)
{
    // FLS 接受第一个未被禁止且达到动态阈值的候选。
    if (numBest > 0 || delta < bestDelta) return;

    bestDelta  = delta;
    numBest    = 1;
    chosenKind = kind;
    chosenA    = a;
    chosenB    = b;
    chosenC    = c;
}

void record_mixed_candidate(int kind, int a, int b, int c,
                            int deltaProfit, double deltaEval,
                            double &bestDeltaEval, int &bestDeltaProfit,
                            int &numBest,
                            int &chosenKind, int &chosenA,
                            int &chosenB, int &chosenC)
{
    if (deltaEval < bestDeltaEval) return;

    if (deltaEval > bestDeltaEval)
    {
        bestDeltaEval  = deltaEval;
        bestDeltaProfit = deltaProfit;
        numBest = 1;
        chosenKind = kind;
        chosenA    = a;
        chosenB    = b;
        chosenC    = c;
        return;
    }

    numBest++;
    if (rand() % numBest == 0)
    {
        bestDeltaProfit = deltaProfit;
        chosenKind = kind;
        chosenA    = a;
        chosenB    = b;
        chosenC    = c;
    }
}

int orphan_chain_has_tabu(int count)
{
    for (int q = 0; q < count; q++)
        if (tabuIter < tabuUntil[orphanTask[q]]) return 1;
    return 0;
}

int orphan_chain_has_task(int count, int j)
{
    for (int q = 0; q < count; q++)
        if (orphanTask[q] == j) return 1;
    return 0;
}

// 保存当前参考链的一个完整叶子。tailDest 为 -1 时，尾部孤儿保持
// 未服务；否则它被放入 tailDest。严格阶段只保留正收益叶子，混合阶段
// 以现有的惩罚评价 deltaEval 选择叶子。
void orphan_consider_leaf(int count, int tailDest, int relaxed)
{
    const double EPS = 1e-12;
    int tail = orphanTask[count - 1];
    int deltaProfit = tailDest >= 0 ? taskProfit[orphanTask[0]]
                                    : taskProfit[orphanTask[0]] - taskProfit[tail];
    double newOverBW = 0.0, newOverPW = 0.0;

    g_orphanTrials++;
    if (tailDest >= 0) g_orphanRelocLeaf++;
    else               g_orphanDropLeaf++;

    for (int b = 0; b < numBeam; b++)
    {
        int bwFree = orphanWorkBW[b];
        int pwFree = orphanWorkPW[b];
        if (b == tailDest)
        {
            bwFree -= taskBWDemand[tail];
            pwFree -= taskPWDemand[tail];
        }
        newOverBW += beam_bw_over_with_rem(b, bwFree);
        newOverPW += beam_pw_over_with_rem(b, pwFree);
    }

    if (!relaxed && deltaProfit <= 0) return;

    double deltaEval = (double)deltaProfit
        - mixedPhiBW * (newOverBW - orphanBaseOverBW)
        - mixedPhiPW * (newOverPW - orphanBaseOverPW);
    if (relaxed && deltaEval <= EPS) return;

    int aspired = 0;
    if (orphan_chain_has_tabu(count))
    {
        if (newOverBW + newOverPW <= EPS && totalProfit + deltaProfit > bestProfit)
            aspired = 1;
        else
        {
            g_tabuBlock++;
            return;
        }
    }

    if (orphanBestCount > 0)
    {
        if (relaxed)
        {
            if (deltaEval < orphanBestDeltaEval) return;
            if (deltaEval == orphanBestDeltaEval &&
                deltaProfit <= orphanBestDeltaProfit) return;
        }
        else if (deltaProfit <= orphanBestDeltaProfit)
            return;
    }

    orphanBestCount = count;
    orphanBestDeltaProfit = deltaProfit;
    orphanBestDeltaEval = deltaEval;
    orphanBestAspired = aspired;
    for (int q = 0; q < count; q++)
    {
        orphanBestTask[q] = orphanTask[q];
        orphanBestDest[q] = orphanDest[q];
    }
    orphanBestDest[count - 1] = tailDest;
}

// 仅保留一个后继：严格阶段优先驱逐低收益任务，混合阶段优先选择
// 即时惩罚评价最好的替换。若当前链尚未包含禁忌任务，非禁忌后继总是
// 优先于禁忌后继，避免禁忌候选遮蔽后续可用候选。
int orphan_select_next(int orphan, int count, int relaxed, int &bestBeam)
{
    const double EPS = 1e-12;
    int bestTask = -1;
    int bestTabuGroup = 2;
    int bestEvictProfit = 0, bestSlack = 0;   // 被驱逐任务收益（勿与全局 bestProfit 混淆）
    double bestEval = -1.0e100;
    int chainTabu = orphan_chain_has_tabu(count);
    int t = taskType[orphan] - 1;

    bestBeam = -1;
    for (int b = 0; b < numBeam; b++)
    {
        if (beamMode[b] < 0) continue;
        if (!typeCompatMode[t][beamMode[b]]) continue;

        for (int idx = bucketStart[b]; idx < bucketStart[b + 1]; idx++)
        {
            int k = bucketTask[idx];
            if (orphan_chain_has_task(count, k)) continue;

            int bw2 = orphanWorkBW[b] + taskBWDemand[k] - taskBWDemand[orphan];
            int pw2 = orphanWorkPW[b] + taskPWDemand[k] - taskPWDemand[orphan];
            if (relaxed)
            {
                if (!relaxed_rem_ok(b, bw2, pw2)) continue;
            }
            else if (bw2 < 0 || pw2 < 0)
                continue;

            int tabuGroup = (!chainTabu && tabuIter < tabuUntil[k]) ? 1 : 0;
            if (tabuGroup > bestTabuGroup) continue;

            int take = 0;
            int slack = bw2 + pw2;
            if (relaxed)
            {
                double oldBW = beam_bw_over_with_rem(b, orphanWorkBW[b]);
                double oldPW = beam_pw_over_with_rem(b, orphanWorkPW[b]);
                double newBW = beam_bw_over_with_rem(b, bw2);
                double newPW = beam_pw_over_with_rem(b, pw2);
                double eval = (double)(taskProfit[orphan] - taskProfit[k])
                    - mixedPhiBW * (newBW - oldBW)
                    - mixedPhiPW * (newPW - oldPW);
                if (tabuGroup < bestTabuGroup || eval > bestEval + EPS ||
                    (eval >= bestEval - EPS && taskProfit[k] < bestEvictProfit) ||
                    (eval >= bestEval - EPS && taskProfit[k] == bestEvictProfit &&
                     slack > bestSlack))
                    take = 1;
                if (take) bestEval = eval;
            }
            else if (tabuGroup < bestTabuGroup || taskProfit[k] < bestEvictProfit ||
                     (taskProfit[k] == bestEvictProfit && slack > bestSlack))
                take = 1;

            if (take)
            {
                bestTask = k;
                bestBeam = b;
                bestTabuGroup = tabuGroup;
                bestEvictProfit = taskProfit[k];
                bestSlack = slack;
            }
        }
    }
    return bestTask;
}

void orphan_grow(int root, int relaxed)
{
    int count = 1;
    int orphan = root;

    orphanTask[0] = root;
    for (int b = 0; b < numBeam; b++)
    {
        orphanWorkBW[b] = remBW[b];
        orphanWorkPW[b] = remPW[b];
    }

    for (int depth = 0; depth < ORPHAN_EC_MAXLEN; depth++)
    {
        int source;
        int next = orphan_select_next(orphan, count, relaxed, source);
        if (next < 0) break;

        orphanWorkBW[source] += taskBWDemand[next] - taskBWDemand[orphan];
        orphanWorkPW[source] += taskPWDemand[next] - taskPWDemand[orphan];
        orphanDest[count - 1] = source;
        orphanTask[count] = next;
        count++;
        orphan = next;
        g_orphanSteps++;

        // 末端孤儿重定位。source 已被上一步前驱填补，跳过它可避免
        // 与邻域 1 的直接插入重复。
        int t = taskType[orphan] - 1;
        for (int b = 0; b < numBeam; b++)
        {
            if (b == source) continue;
            if (beamMode[b] < 0) continue;
            if (!typeCompatMode[t][beamMode[b]]) continue;

            int bw2 = orphanWorkBW[b] - taskBWDemand[orphan];
            int pw2 = orphanWorkPW[b] - taskPWDemand[orphan];
            if (relaxed)
            {
                if (!relaxed_rem_ok(b, bw2, pw2)) continue;
            }
            else if (bw2 < 0 || pw2 < 0)
                continue;
            orphan_consider_leaf(count, b, relaxed);
        }

        // 深度 1 的尾部不服务正是已有 N3 Replacement，不重复评价。
        if (count > 2) orphan_consider_leaf(count, -1, relaxed);
    }
}

// 返回是否找到正收益（严格）或正 deltaEval（混合）的完整孤儿链。
int orphan_search(int relaxed)
{
    orphanBestCount = 0;
    orphanBestDeltaProfit = MININT_MOVE;
    orphanBestDeltaEval = -1.0e100;
    orphanBestAspired = 0;
    total_over_parts(orphanBaseOverBW, orphanBaseOverPW);
    for (int j = 0; j < numTask; j++) orphanStartMark[j] = 0;

    for (int rank = 0; rank < ORPHAN_EC_START_CAP; rank++)
    {
        int root = -1;
        for (int j = 0; j < numTask; j++)
        {
            if (taskBeam[j] != -1 || orphanStartMark[j]) continue;
            if (root < 0 || taskProfit[j] > taskProfit[root]) root = j;
        }
        if (root < 0) break;
        orphanStartMark[root] = 1;
        orphan_grow(root, relaxed);
    }
    return orphanBestCount > 0;
}

void apply_orphan_best()
{
    int tenure = tabu_tenure();
    for (int q = 1; q < orphanBestCount; q++) remove_task(orphanBestTask[q]);
    for (int q = 0; q < orphanBestCount; q++)
        if (orphanBestDest[q] >= 0)
            add_task(orphanBestTask[q], orphanBestDest[q]);

    for (int q = 0; q < orphanBestCount; q++)
    {
        int j = orphanBestTask[q];
        moveFreq[j]++;
        tabuUntil[j] = tabuIter + tenure;
    }
    if (orphanBestAspired) g_aspire++;
    g_orphanApplied++;
}

//--------------------------------------------------------------------
// apply_mode_flip：把波束 b 切换到模式 m2，只驱逐与新模式不兼容
// 的任务；这里不执行回填。
//--------------------------------------------------------------------
void apply_mode_flip(int b, int m2, int tenure)
{
    int saveMode = beamMode[b];

    for (int idx = bucketStart[b]; idx < bucketStart[b + 1]; idx++)
    {
        int j = bucketTask[idx];
        if (taskBeam[j] == b && !typeCompatMode[taskType[j] - 1][m2])
        {
            remove_task(j);
            moveFreq[j]++;
            tabuUntil[j] = tabuIter + tenure;
        }
    }

    beamMode[b] = m2;
    remPW[b] += modeBasePower[saveMode] - modeBasePower[m2];
}

//--------------------------------------------------------------------
// rebuild_buckets：按波束分组重建已服务任务的 CSR 桶，供严格与混合两个
// 局部搜索的邻域扫描共用。对 taskBeam 做计数排序，复杂度 O(numTask+numBeam)。
//--------------------------------------------------------------------
void rebuild_buckets()
{
    for (int b = 0; b <= numBeam; b++) bucketStart[b] = 0;
    for (int j = 0; j < numTask; j++)
        if (taskBeam[j] >= 0) bucketStart[taskBeam[j] + 1]++;
    for (int b = 0; b < numBeam; b++) bucketStart[b + 1] += bucketStart[b];
    // 使用每个波束的局部游标填桶，并复用 bucketStart 作为游标
    for (int j = 0; j < numTask; j++)
        if (taskBeam[j] >= 0)
        {
            int b = taskBeam[j];
            bucketTask[bucketStart[b]++] = j;
        }
    // 撤销游标推进，使 bucketStart[b] 重新指向波束 b 的起点
    for (int b = numBeam; b > 0; b--) bucketStart[b] = bucketStart[b - 1];
    bucketStart[0] = 0;
}

//--------------------------------------------------------------------
// count_aspiration：若已选定动作突破了任务或波束模式禁忌，统计一次特赦。
// 两个局部搜索在执行动作前调用，特赦口径单一来源、保持一致。
//--------------------------------------------------------------------
void count_aspiration(int kind, int a, int cc)
{
    if ((kind == 1 && tabuIter < tabuUntil[a]) ||
        (kind == 2 && tabuIter < tabuUntil[a]) ||
        (kind == 3 && (tabuIter < tabuUntil[a] || tabuIter < tabuUntil[cc])) ||
        (kind == 4 && tabuIter < tabuBeamMode[a]) ||
        (kind == 5 && (tabuIter < tabuUntil[a] ||
                       (cc >= 0 && tabuIter < tabuUntil[cc]))))
        g_aspire++;
}

//--------------------------------------------------------------------
// apply_move：执行已选定的邻域动作，并同步移动频次、属性禁忌与统计。
// kind：1 插入  2 纯移除  3 同波束交换  4 模式翻转  5 跨波束交换
// (cc<0 表示实任务-虚拟空位重定位)。严格与混合两个局部搜索共用，
// 保证动作语义与禁忌写入只有一处实现，避免改一处漏另一处。
//--------------------------------------------------------------------
void apply_move(int kind, int a, int bb, int cc)
{
    if (kind == 1)                            // 插入
    {
        add_task(a, bb);
        moveFreq[a]++;
        tabuUntil[a] = tabuIter + tabu_tenure();
        g_applyInsert++;
    }
    else if (kind == 2)                       // 纯移除
    {
        remove_task(a);
        moveFreq[a]++;
        tabuUntil[a] = tabuIter + tabu_tenure();
        g_applyRemove++;
    }
    else if (kind == 3)                       // 交换：a 加入、cc 移出
    {
        remove_task(cc);
        add_task(a, bb);
        moveFreq[a]++; moveFreq[cc]++;
        tabuUntil[a]  = tabuIter + tabu_tenure();
        tabuUntil[cc] = tabuIter + tabu_tenure();
        g_applySwap++;
    }
    else if (kind == 4)                       // 把波束 a 的模式翻转为 bb
    {
        int tenure = tabu_tenure();
        apply_mode_flip(a, bb, tenure);
        tabuBeamMode[a] = tabuIter + tenure;
        g_applyFlip++;
        g_flipApplied++;
        if (bb >= 0 && bb < 16) g_flipApplyToMode[bb]++;
    }
    else if (kind == 5)                       // 带虚拟空位的跨波束交换
    {
        int b1 = taskBeam[a];
        if (cc < 0)                           // 实任务-虚拟空位：把 a 重定位到 bb
        {
            remove_task(a);
            add_task(a, bb);
            moveFreq[a]++;
            tabuUntil[a] = tabuIter + tabu_tenure();
            g_crossRelocApplied++;
        }
        else                                  // 实任务-实任务：交换 a 和 cc
        {
            int b2 = taskBeam[cc];
            remove_task(a);
            remove_task(cc);
            add_task(a, b2);
            add_task(cc, b1);
            moveFreq[a]++; moveFreq[cc]++;
            tabuUntil[a]  = tabuIter + tabu_tenure();
            tabuUntil[cc] = tabuIter + tabu_tenure();
            g_crossSwapApplied++;
        }
        g_applyCross++;
    }
}

//--------------------------------------------------------------------
void local_search(double beginTime, int *tmpIdx, double *tmpVal)
{
    (void)tmpIdx;
    (void)tmpVal;

    int ts_depth = 300;              // 连续无改进达到该次数后结束阶段并触发扰动
    int nonImprove = 0;
    int thresholdDelta = threshold_delta();
    double lsStart = (double)clock();

    for (int j = 0; j < numTask; j++) tabuUntil[j] = 0;
    for (int b = 0; b < numBeam; b++) tabuBeamMode[b] = 0;
    tabuIter = 0;

    while (nonImprove < ts_depth)
    {
        g_lsIters++;
        if (((double)clock() - beginTime) / CLOCKS_PER_SEC > maxRunTime)
            break;

        int threshold = bestProfit - thresholdDelta;
        int bestDelta = threshold - totalProfit;
        int numBest   = 0;               // 找到首个达到阈值的动作后变为 1
        int kind = 0, a = -1, bb = -1, cc = -1;

        // 重建按波束分组的已服务任务 CSR 桶。
        rebuild_buckets();

        // 邻域 1：插入一个未服务任务
        for (int j = 0; j < numTask && numBest == 0; j++)
        {
            if (taskBeam[j] != -1) continue;
            int delta = taskProfit[j];
            // 特赦：禁忌插入只有在候选收益超过阶段最好收益时才允许。
            int jTabu = (tabuIter < tabuUntil[j]);
            if (jTabu)
            {
                if (totalProfit + delta > bestProfit) { /* 满足特赦 */ }
                else { g_tabuBlock++; continue; }
            }
            if (delta < bestDelta) continue;   // delta 与波束无关，提前跳过整轮波束扫描
            for (int b = 0; b < numBeam && numBest == 0; b++)
            {
                if (!feasible_on(j, b)) continue;
                record_candidate(1, j, b, -1, delta,
                                 bestDelta, numBest, kind, a, bb, cc);
            }
        }

        // 邻域 3，第二个搜索：未服务任务 i 替换已服务任务 k。
        // 对每个未服务任务 i，只访问模式与其兼容的波束，并只扫描该波束
        // CSR 桶中的已服务任务，避免重新扫描所有任务；候选对保持不变。
        for (int i = 0; i < numTask && numBest == 0; i++)
        {
            if (taskBeam[i] != -1) continue;
            int iTabu = (tabuIter < tabuUntil[i]);
            int ti = taskType[i] - 1;
            for (int b = 0; b < numBeam && numBest == 0; b++)
            {
                if (beamMode[b] < 0) continue;
                if (!typeCompatMode[ti][beamMode[b]]) continue;

                int bwFree0 = remBW[b];
                int pwFree0 = remPW[b];
                for (int idx = bucketStart[b];
                     idx < bucketStart[b + 1] && numBest == 0; idx++)
                {
                    int k = bucketTask[idx];

                    int delta = taskProfit[i] - taskProfit[k];
                    if (delta < bestDelta) continue;

                    // 特赦：i 或 k 处于禁忌时，只有候选收益超过阶段最好
                    // 收益才允许执行交换。
                    if (iTabu || tabuIter < tabuUntil[k])
                    {
                        if (totalProfit + delta > bestProfit) { /* 满足特赦 */ }
                        else continue;
                    }

                    // 检查同一波束 b 上移除 k、加入 i 后的资源可行性
                    int bwFree = bwFree0 + taskBWDemand[k];
                    int pwFree = pwFree0 + taskPWDemand[k];
                    if (taskBWDemand[i] > bwFree || taskPWDemand[i] > pwFree) continue;

                    record_candidate(3, i, b, k, delta,
                                     bestDelta, numBest, kind, a, bb, cc);
                }
            }
        }

        // 邻域 2，第三个搜索：移除一个已服务任务
        for (int j = 0; j < numTask && numBest == 0; j++)
        {
            if (taskBeam[j] < 0) continue;
            int delta = -taskProfit[j];
            if (tabuIter < tabuUntil[j])
            {
                if (totalProfit + delta > bestProfit) { /* 满足特赦 */ }
                else { g_tabuBlock++; continue; }
            }
            if (delta < bestDelta) continue;

            record_candidate(2, j, -1, -1, delta,
                             bestDelta, numBest, kind, a, bb, cc);
        }

        // 邻域 4：翻转一个波束的模式，只驱逐与新模式不兼容的任务，
        // 这里不评价回填。
        for (int b = 0; b < numBeam && numBest == 0; b++)
        {
            if (beamMode[b] < 0) continue;

            int usedPW = (beamPWCap[b] - modeBasePower[beamMode[b]]) - remPW[b];
            for (int m2 = 0; m2 < numMode && numBest == 0; m2++)
            {
                if (m2 == beamMode[b]) continue;
                if (modeBasePower[m2] > beamPWCap[b]) continue;

                int loss = 0;
                for (int t = 0; t < numType; t++)
                    if (!typeCompatMode[t][m2]) loss += typeProfitSum[b][t];

                int evictPW = 0;
                int evictTabu = 0;
                for (int idx = bucketStart[b]; idx < bucketStart[b + 1]; idx++)
                {
                    int j = bucketTask[idx];
                    if (taskBeam[j] == b && !typeCompatMode[taskType[j] - 1][m2])
                    {
                        evictPW += taskPWDemand[j];
                        if (tabuIter < tabuUntil[j]) evictTabu = 1;
                    }
                }
                if (usedPW - evictPW > beamPWCap[b] - modeBasePower[m2]) continue;

                int delta = -loss;
                int tabuFlip = (tabuIter < tabuBeamMode[b]) || evictTabu;
                if (tabuFlip)
                {
                    if (totalProfit + delta > bestProfit) { /* 满足特赦 */ }
                    else { g_tabuBlock++; continue; }
                }

                g_flipCand++;
                if (delta < bestDelta) continue;
                record_candidate(4, b, m2, -1, delta,
                                 bestDelta, numBest, kind, a, bb, cc);
            }
        }

        // N1--N4 已有的动态阈值动作若为严格正收益，直接保留；否则
        // 先尝试未服务任务孤儿链。原有 Cross 只在孤儿链没有正收益
        // 完整叶子时作为后备。非正的动态阈值动作仍在两个邻域都失败
        // 时保留，维持原有阈值接受机制。
        int useOrphan = 0;
        int haveFallback = numBest > 0;
        int fallbackKind = kind, fallbackA = a, fallbackB = bb, fallbackC = cc;
        int fallbackDelta = bestDelta;
        if (!haveFallback || bestDelta <= 0)
        {
            g_orphanCalls++;
            if (orphan_search(0))
            {
                useOrphan = 1;
                numBest = 1;
            }
            else
            {
                // 重新从动态阈值开始记录 Cross，不能让保留的负动作
                // 遮蔽收益为 0 的重定位或交换。
                numBest = 0;
                bestDelta = threshold - totalProfit;
            }
        }

        // 邻域 5：带动态虚拟空位的跨波束交换。
        // 实任务-实任务：b1 上的已服务任务 i 与 b2 上的已服务任务 k 交换波束。
        // 实任务-虚拟空位：把 b1 上的 i 重定位到 b2，cc == -1 表示临时
        // 虚拟空位；虚拟空位不存入 taskBeam[]。
        // 只有邻域 1-4 都没有未被禁止且达到阈值的动作时才搜索该邻域。
        if (numBest == 0 && !useOrphan)
        {
            rebuild_insert_filter();

            for (int i = 0; i < numTask && numBest == 0; i++)
            {
                int b1 = taskBeam[i];
                if (b1 < 0) continue;

                int iTabu = (tabuIter < tabuUntil[i]);
                if (iTabu)
                {
                    if (totalProfit > bestProfit) { /* 满足特赦 */ }
                    else { g_tabuBlock++; continue; }
                }

                int ti = taskType[i] - 1;
                for (int b2 = 0; b2 < numBeam && numBest == 0; b2++)
                {
                    if (b2 == b1) continue;
                    if (beamMode[b2] < 0) continue;
                    if (!typeCompatMode[ti][beamMode[b2]]) continue;

                    int delta = 0;

                    // 实任务-虚拟空位：把 i 从 b1 移到 b2 的当前空位。
                    if (taskBWDemand[i] <= remBW[b2] &&
                        taskPWDemand[i] <= remPW[b2])
                    {
                        if (!beam_can_insert_unserved_with_rem(
                                b1,
                                remBW[b1] + taskBWDemand[i],
                                remPW[b1] + taskPWDemand[i]))
                        {
                            g_crossFiltered++;
                        }
                        else
                        {
                            g_crossCand++;
                            g_crossRelocCand++;
                            record_candidate(5, i, b2, -1, delta,
                                             bestDelta, numBest, kind, a, bb, cc);
                        }
                    }

                    // 实任务-实任务：每个跨波束任务对只枚举一次。
                    if (b1 > b2) continue;

                    for (int idx = bucketStart[b2];
                         idx < bucketStart[b2 + 1] && numBest == 0; idx++)
                    {
                        int k = bucketTask[idx];
                        int tk = taskType[k] - 1;

                        if (!typeCompatMode[tk][beamMode[b1]]) continue;
                        if (taskBWDemand[k] > remBW[b1] + taskBWDemand[i]) continue;
                        if (taskPWDemand[k] > remPW[b1] + taskPWDemand[i]) continue;
                        if (taskBWDemand[i] > remBW[b2] + taskBWDemand[k]) continue;
                        if (taskPWDemand[i] > remPW[b2] + taskPWDemand[k]) continue;

                        if (tabuIter < tabuUntil[k])
                        {
                            if (totalProfit > bestProfit) { /* 满足特赦 */ }
                            else { g_tabuBlock++; continue; }
                        }

                        int b1BW = remBW[b1] + taskBWDemand[i] - taskBWDemand[k];
                        int b1PW = remPW[b1] + taskPWDemand[i] - taskPWDemand[k];
                        int b2BW = remBW[b2] + taskBWDemand[k] - taskBWDemand[i];
                        int b2PW = remPW[b2] + taskPWDemand[k] - taskPWDemand[i];
                        if (!beam_can_insert_unserved_with_rem(b1, b1BW, b1PW) &&
                            !beam_can_insert_unserved_with_rem(b2, b2BW, b2PW))
                        {
                            g_crossFiltered++;
                            continue;
                        }

                        g_crossCand++;
                        g_crossSwapCand++;
                        record_candidate(5, i, b2, k, delta,
                                         bestDelta, numBest, kind, a, bb, cc);
                    }
                }
            }
        }

        if (!useOrphan && numBest == 0 && haveFallback)
        {
            bestDelta = fallbackDelta;
            numBest = 1;
            kind = fallbackKind;
            a = fallbackA;
            bb = fallbackB;
            cc = fallbackC;
        }

        // 本轮没有可接受候选时，只推进禁忌时间
        if (numBest == 0)
        {
            g_lsStall++;
            nonImprove++;
            tabuIter++;
            continue;
        }
        g_lsMoves++;
        if (useOrphan)
        {
            apply_orphan_best();
            if (totalProfit > bestProfit)
            {
                save_best(beginTime);
                nonImprove = 0;
            }
            else
                nonImprove++;
            tabuIter++;
            continue;
        }
        // 若已提交动作突破了任务或波束模式禁忌，则统计一次特赦
        count_aspiration(kind, a, cc);
        apply_move(kind, a, bb, cc);

        if (totalProfit > bestProfit)
        {
            save_best(beginTime);
            nonImprove = 0;
        }
        else
            nonImprove++;

        tabuIter++;
    }
    g_lsTime += ((double)clock() - lsStart) / CLOCKS_PER_SEC;
}

//--------------------------------------------------------------------
// local_search_mixed：可行域与不可行域混合禁忌搜索。
// BW/PW 容量约束可在小范围内放宽；模式兼容和每个任务至多分配给一个
// 波束仍是硬约束。当前解可以不可行，但只有可行解能够更新 bestProfit。
//--------------------------------------------------------------------
void local_search_mixed(double beginTime, int *tmpIdx, double *tmpVal)
{
    (void)tmpIdx;
    (void)tmpVal;

    const double EPS = 1e-12;
    const int mixedDepth = 300;
    const int phiWindow = 5;
    const double phiTau = 1.5;
    const double phiMin = 100.0;
    const double phiMax = 100000.0;

    int nonImprove = 0;
    int bwFeasibleStreak = 0, bwInfeasibleStreak = 0;
    int pwFeasibleStreak = 0, pwInfeasibleStreak = 0;
    long long localTabuBlock = 0;
    long long localCrossCand = 0;
    long long localCrossRelocCand = 0;
    long long localCrossSwapCand = 0;
    double mixedStart = (double)clock();
    double curOverBW, curOverPW;
    total_over_parts(curOverBW, curOverPW);
    double curOver = curOverBW + curOverPW;
    set_mixed_rho_from_instance();

    for (int j = 0; j < numTask; j++) tabuUntil[j] = 0;
    for (int b = 0; b < numBeam; b++) tabuBeamMode[b] = 0;
    tabuIter = 0;

    while (nonImprove < mixedDepth)
    {
        if (((double)clock() - beginTime) / CLOCKS_PER_SEC > maxRunTime)
            break;
        g_mixedIters++;

        if (curOver <= EPS)
        {
            g_mixedFeas++;
        }
        else
        {
            g_mixedInfeas++;
        }

        if (curOverBW <= EPS)
        {
            bwFeasibleStreak++;
            bwInfeasibleStreak = 0;
            if (bwFeasibleStreak >= phiWindow)
            {
                mixedPhiBW /= phiTau;
                if (mixedPhiBW < phiMin) mixedPhiBW = phiMin;
                bwFeasibleStreak = 0;
            }
        }
        else
        {
            bwInfeasibleStreak++;
            bwFeasibleStreak = 0;
            if (bwInfeasibleStreak >= phiWindow)
            {
                mixedPhiBW *= phiTau;
                if (mixedPhiBW > phiMax) mixedPhiBW = phiMax;
                bwInfeasibleStreak = 0;
            }
        }

        if (curOverPW <= EPS)
        {
            pwFeasibleStreak++;
            pwInfeasibleStreak = 0;
            if (pwFeasibleStreak >= phiWindow)
            {
                mixedPhiPW /= phiTau;
                if (mixedPhiPW < phiMin) mixedPhiPW = phiMin;
                pwFeasibleStreak = 0;
            }
        }
        else
        {
            pwInfeasibleStreak++;
            pwFeasibleStreak = 0;
            if (pwInfeasibleStreak >= phiWindow)
            {
                mixedPhiPW *= phiTau;
                if (mixedPhiPW > phiMax) mixedPhiPW = phiMax;
                pwInfeasibleStreak = 0;
            }
        }

        rebuild_relaxed_limits();

        double bestDeltaEval = -1.0e100;
        int bestDeltaProfit = MININT_MOVE;
        int numBest = 0;
        int kind = 0, a = -1, bb = -1, cc = -1;

        rebuild_buckets();

        // 邻域 1：插入一个未服务任务
        for (int j = 0; j < numTask; j++)
        {
            if (taskBeam[j] != -1) continue;
            int deltaProfit = taskProfit[j];
            int t = taskType[j] - 1;

            for (int b = 0; b < numBeam; b++)
            {
                if (beamMode[b] < 0) continue;
                if (!typeCompatMode[t][beamMode[b]]) continue;

                int bw2 = remBW[b] - taskBWDemand[j];
                int pw2 = remPW[b] - taskPWDemand[j];
                if (!relaxed_rem_ok(b, bw2, pw2)) continue;

                double newOverBW = curOverBW - beam_bw_over(b) + beam_bw_over_with_rem(b, bw2);
                double newOverPW = curOverPW - beam_pw_over(b) + beam_pw_over_with_rem(b, pw2);
                double newOver = newOverBW + newOverPW;
                if (tabuIter < tabuUntil[j])
                {
                    if (newOver <= EPS && totalProfit + deltaProfit > bestProfit) { /* 特赦 */ }
                    else { localTabuBlock++; continue; }
                }

                double deltaEval = (double)deltaProfit
                    - mixedPhiBW * (newOverBW - curOverBW)
                    - mixedPhiPW * (newOverPW - curOverPW);
                record_mixed_candidate(1, j, b, -1, deltaProfit, deltaEval,
                                       bestDeltaEval, bestDeltaProfit, numBest,
                                       kind, a, bb, cc);
            }
        }

        // 邻域 2：纯移除
        for (int j = 0; j < numTask; j++)
        {
            int b = taskBeam[j];
            if (b < 0) continue;

            int deltaProfit = -taskProfit[j];
            int bw2 = remBW[b] + taskBWDemand[j];
            int pw2 = remPW[b] + taskPWDemand[j];
            double newOverBW = curOverBW - beam_bw_over(b) + beam_bw_over_with_rem(b, bw2);
            double newOverPW = curOverPW - beam_pw_over(b) + beam_pw_over_with_rem(b, pw2);
            double newOver = newOverBW + newOverPW;

            if (tabuIter < tabuUntil[j])
            {
                if (newOver <= EPS && totalProfit + deltaProfit > bestProfit) { /* 特赦 */ }
                else { localTabuBlock++; continue; }
            }

            double deltaEval = (double)deltaProfit
                - mixedPhiBW * (newOverBW - curOverBW)
                - mixedPhiPW * (newOverPW - curOverPW);
            record_mixed_candidate(2, j, -1, -1, deltaProfit, deltaEval,
                                   bestDeltaEval, bestDeltaProfit, numBest,
                                   kind, a, bb, cc);
        }

        // 邻域 3：同波束交换
        for (int i = 0; i < numTask; i++)
        {
            if (taskBeam[i] != -1) continue;
            int ti = taskType[i] - 1;

            for (int b = 0; b < numBeam; b++)
            {
                if (beamMode[b] < 0) continue;
                if (!typeCompatMode[ti][beamMode[b]]) continue;

                for (int idx = bucketStart[b]; idx < bucketStart[b + 1]; idx++)
                {
                    int k = bucketTask[idx];
                    int deltaProfit = taskProfit[i] - taskProfit[k];
                    int bw2 = remBW[b] + taskBWDemand[k] - taskBWDemand[i];
                    int pw2 = remPW[b] + taskPWDemand[k] - taskPWDemand[i];
                    if (!relaxed_rem_ok(b, bw2, pw2)) continue;

                    double newOverBW = curOverBW - beam_bw_over(b) + beam_bw_over_with_rem(b, bw2);
                    double newOverPW = curOverPW - beam_pw_over(b) + beam_pw_over_with_rem(b, pw2);
                    double newOver = newOverBW + newOverPW;
                    if (tabuIter < tabuUntil[i] || tabuIter < tabuUntil[k])
                    {
                        if (newOver <= EPS && totalProfit + deltaProfit > bestProfit) { /* 特赦 */ }
                        else { localTabuBlock++; continue; }
                    }

                    double deltaEval = (double)deltaProfit
                        - mixedPhiBW * (newOverBW - curOverBW)
                        - mixedPhiPW * (newOverPW - curOverPW);
                    record_mixed_candidate(3, i, b, k, deltaProfit, deltaEval,
                                           bestDeltaEval, bestDeltaProfit, numBest,
                                           kind, a, bb, cc);
                }
            }
        }

        // 邻域 4：模式翻转
        for (int b = 0; b < numBeam; b++)
        {
            if (beamMode[b] < 0) continue;

            for (int m2 = 0; m2 < numMode; m2++)
            {
                if (m2 == beamMode[b]) continue;
                if (modeBasePower[m2] > beamPWCap[b]) continue;

                int loss = 0, evictBW = 0, evictPW = 0, evictTabu = 0;
                for (int idx = bucketStart[b]; idx < bucketStart[b + 1]; idx++)
                {
                    int j = bucketTask[idx];
                    if (taskBeam[j] == b && !typeCompatMode[taskType[j] - 1][m2])
                    {
                        loss += taskProfit[j];
                        evictBW += taskBWDemand[j];
                        evictPW += taskPWDemand[j];
                        if (tabuIter < tabuUntil[j]) evictTabu = 1;
                    }
                }

                int bw2 = remBW[b] + evictBW;
                int pw2 = remPW[b] + evictPW + modeBasePower[beamMode[b]] - modeBasePower[m2];
                if (!relaxed_rem_ok_mode(b, m2, bw2, pw2)) continue;

                int deltaProfit = -loss;
                double newOverBW = curOverBW - beam_bw_over(b) + beam_bw_over_with_rem(b, bw2);
                double newOverPW = curOverPW - beam_pw_over(b) + beam_pw_over_with_rem_mode(b, m2, pw2);
                double newOver = newOverBW + newOverPW;
                if (tabuIter < tabuBeamMode[b] || evictTabu)
                {
                    if (newOver <= EPS && totalProfit + deltaProfit > bestProfit) { /* 特赦 */ }
                    else { localTabuBlock++; continue; }
                }

                g_flipCand++;
                double deltaEval = (double)deltaProfit
                    - mixedPhiBW * (newOverBW - curOverBW)
                    - mixedPhiPW * (newOverPW - curOverPW);
                record_mixed_candidate(4, b, m2, -1, deltaProfit, deltaEval,
                                       bestDeltaEval, bestDeltaProfit, numBest,
                                       kind, a, bb, cc);
            }
        }

        // 混合阶段同样先比较 N1--N4；只有它们没有正 deltaEval 时，
        // 才先试孤儿链，失败后再回退到 Cross。
        int useOrphan = 0;
        if (bestDeltaEval <= EPS)
        {
            g_orphanCalls++;
            if (orphan_search(1)) useOrphan = 1;
        }

        // 邻域 5：带动态虚拟空位的跨波束交换
        if (!useOrphan && bestDeltaEval <= EPS)
        {
        for (int i = 0; i < numTask; i++)
        {
            int b1 = taskBeam[i];
            if (b1 < 0) continue;

            int ti = taskType[i] - 1;
            for (int b2 = 0; b2 < numBeam; b2++)
            {
                if (b2 == b1) continue;
                if (beamMode[b2] < 0) continue;
                if (!typeCompatMode[ti][beamMode[b2]]) continue;

                int b1BW = remBW[b1] + taskBWDemand[i];
                int b1PW = remPW[b1] + taskPWDemand[i];
                int b2BW = remBW[b2] - taskBWDemand[i];
                int b2PW = remPW[b2] - taskPWDemand[i];
                if (relaxed_rem_ok(b1, b1BW, b1PW) &&
                    relaxed_rem_ok(b2, b2BW, b2PW))
                {
                    double newOverBW = curOverBW
                        - beam_bw_over(b1) - beam_bw_over(b2)
                        + beam_bw_over_with_rem(b1, b1BW)
                        + beam_bw_over_with_rem(b2, b2BW);
                    double newOverPW = curOverPW
                        - beam_pw_over(b1) - beam_pw_over(b2)
                        + beam_pw_over_with_rem(b1, b1PW)
                        + beam_pw_over_with_rem(b2, b2PW);
                    double newOver = newOverBW + newOverPW;

                    if (tabuIter < tabuUntil[i])
                    {
                        if (newOver <= EPS && totalProfit > bestProfit) { /* 特赦 */ }
                        else { localTabuBlock++; continue; }
                    }

                    double deltaEval = -mixedPhiBW * (newOverBW - curOverBW)
                                       -mixedPhiPW * (newOverPW - curOverPW);
                    localCrossCand++;
                    localCrossRelocCand++;
                    record_mixed_candidate(5, i, b2, -1, 0, deltaEval,
                                           bestDeltaEval, bestDeltaProfit, numBest,
                                           kind, a, bb, cc);
                }

                if (b1 > b2) continue;
                for (int idx = bucketStart[b2]; idx < bucketStart[b2 + 1]; idx++)
                {
                    int k = bucketTask[idx];
                    int tk = taskType[k] - 1;
                    if (!typeCompatMode[tk][beamMode[b1]]) continue;

                    b1BW = remBW[b1] + taskBWDemand[i] - taskBWDemand[k];
                    b1PW = remPW[b1] + taskPWDemand[i] - taskPWDemand[k];
                    b2BW = remBW[b2] + taskBWDemand[k] - taskBWDemand[i];
                    b2PW = remPW[b2] + taskPWDemand[k] - taskPWDemand[i];
                    if (!relaxed_rem_ok(b1, b1BW, b1PW)) continue;
                    if (!relaxed_rem_ok(b2, b2BW, b2PW)) continue;

                    double newOverBW = curOverBW
                        - beam_bw_over(b1) - beam_bw_over(b2)
                        + beam_bw_over_with_rem(b1, b1BW)
                        + beam_bw_over_with_rem(b2, b2BW);
                    double newOverPW = curOverPW
                        - beam_pw_over(b1) - beam_pw_over(b2)
                        + beam_pw_over_with_rem(b1, b1PW)
                        + beam_pw_over_with_rem(b2, b2PW);
                    double newOver = newOverBW + newOverPW;

                    if (tabuIter < tabuUntil[i] || tabuIter < tabuUntil[k])
                    {
                        if (newOver <= EPS && totalProfit > bestProfit) { /* 特赦 */ }
                        else { localTabuBlock++; continue; }
                    }

                    double deltaEval = -mixedPhiBW * (newOverBW - curOverBW)
                                       -mixedPhiPW * (newOverPW - curOverPW);
                    localCrossCand++;
                    localCrossSwapCand++;
                    record_mixed_candidate(5, i, b2, k, 0, deltaEval,
                                           bestDeltaEval, bestDeltaProfit, numBest,
                                           kind, a, bb, cc);
                }
            }
        }
        }

        if (numBest == 0 && !useOrphan)
        {
            nonImprove++;
            if (nonImprove > 0 && nonImprove % 50 == 0)
                expand_mixed_rho();
            tabuIter++;
            continue;
        }

        g_mixedMoves++;
        if (useOrphan)
        {
            apply_orphan_best();
            total_over_parts(curOverBW, curOverPW);
            curOver = curOverBW + curOverPW;
            if (curOver > g_mixedMaxOver) g_mixedMaxOver = curOver;

            if (curOver <= EPS && totalProfit > bestProfit)
            {
                save_best(beginTime);
                shrink_mixed_rho();
                g_mixedBestUpdate++;
                nonImprove = 0;
            }
            else
            {
                nonImprove++;
                if (nonImprove > 0 && nonImprove % 50 == 0)
                    expand_mixed_rho();
            }
            tabuIter++;
            continue;
        }
        count_aspiration(kind, a, cc);
        apply_move(kind, a, bb, cc);

        total_over_parts(curOverBW, curOverPW);
        curOver = curOverBW + curOverPW;
        if (curOver > g_mixedMaxOver) g_mixedMaxOver = curOver;

        if (curOver <= EPS && totalProfit > bestProfit)
        {
            save_best(beginTime);
            shrink_mixed_rho();
            g_mixedBestUpdate++;
            nonImprove = 0;
        }
        else
        {
            nonImprove++;
            if (nonImprove > 0 && nonImprove % 50 == 0)
                expand_mixed_rho();
        }

        tabuIter++;
    }

    restore_best();       // 混合阶段始终向后续流程交付可行解
    g_tabuBlock += localTabuBlock;
    g_crossCand += localCrossCand;
    g_crossRelocCand += localCrossRelocCand;
    g_crossSwapCand += localCrossSwapCand;
    g_mixedTime += ((double)clock() - mixedStart) / CLOCKS_PER_SEC;
}

//--------------------------------------------------------------------
// perturb：破坏与修复
//
// 从当前阶段最好解出发：
//
//  破坏：选择 30% 的波束，移除其全部任务并放回未服务任务池，同时把
//  模式重置为 -1，形成部分解；未选中的波束保留模式和任务分配。
//
//  模式修复：依次为被清空波束重新选择模式。对一个波束枚举功率可行
//  模式，从当前未服务任务池临时装填，通常保留实际收益增量最大的模式；
//  以较小概率选择次优模式，避免确定性修复完全抵消破坏。修复下一个
//  波束前先提交当前模式和任务，防止模式选择时重复计算任务。
//
//  任务修复：按现有最高可达优先级排序，对仍未服务任务进行贪心回填。
//  每个任务先尝试基础功率较低的兼容波束模式；在同一模式层内使用原有
//  的最紧带宽适配规则。
//--------------------------------------------------------------------
void perturb(int *tmpIdx, double *tmpVal)
{
    restore_best();
    if (numBeam == 0) return;

    // 破坏：清空 30% 的波束，优先选择上次局部搜索中已服务任务移动次数
    // 较少的波束，推动搜索进入尚未充分探索的区域
    int numDestroy = (int)(0.60 * numBeam + 0.5);
    if (numDestroy < 1)        numDestroy = 1;
    if (numDestroy > numBeam)  numDestroy = numBeam;

    // 波束频率键等于该波束已服务任务 moveFreq 之和
    int    *perturbBeam = new int[numBeam];      // 临时保存待清空波束
    int    *prevMode    = new int[numBeam];      // 保存各清空波束破坏前的模式
    double *beamFreqKey = new double[numBeam];
    int    *order       = new int[numBeam];
    for (int b = 0; b < numBeam; b++) { beamFreqKey[b] = 0.0; order[b] = b; }
    for (int j = 0; j < numTask; j++)
        if (taskBeam[j] >= 0) beamFreqKey[taskBeam[j]] += moveFreq[j];

    // 对频率键取负后使用 qsort_desc，实现按频率升序排列
    for (int b = 0; b < numBeam; b++) beamFreqKey[b] = -beamFreqKey[b];
    qsort_desc(beamFreqKey, order, 0, numBeam - 1);   // order[] 中低频波束在前

    // 从低移动频率波束池抽取，候选池至少包含 numDestroy 个波束
    int poolSize = 0;
    for (int b = 0; b < numBeam; b++) if (beamFreqKey[b] == 0.0) poolSize++;
    if (poolSize < numDestroy) poolSize = numDestroy;

    // 从低频候选池无放回随机抽取 numDestroy 个波束
    for (int k = 0; k < numDestroy; k++)             // 部分 Fisher-Yates 洗牌
    {
        int r   = k + rand() % (poolSize - k);
        int tmp = order[k]; order[k] = order[r]; order[r] = tmp;
        perturbBeam[k] = order[k];
    }

    for (int k = 0; k < numDestroy; k++)
    {
        int b = perturbBeam[k];
        prevMode[k] = beamMode[b];                  // 保存旧模式，修复时允许重新选择
        for (int j = 0; j < numTask; j++)
            if (taskBeam[j] == b) remove_task(j);
        beamMode[b] = -1;
        remBW[b]    = beamBWCap[b];
        remPW[b]    = 0;                             // 重新确定模式后再设置
    }

    // 模式修复：按实际收益依次试探
    int *trialAdded = new int[numTask];
    int *bestAdded  = new int[numTask];
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

        if (bestM < 0)                       // 没有功率可行模式时选择
        {                                    // 基础功率最低的可用模式
            bestM = -1;
            for (int m = 0; m < numMode; m++)
            {
                if (bestM < 0 || modeBasePower[m] < modeBasePower[bestM]) bestM = m;
            }
            if (bestM < 0) bestM = prevMode[k];   // 只有一个模式时无法规避
        }

        int chosenM      = bestM;
        int *chosenAdded = bestAdded;
        int nChosenAdded = nBestAdded;
        if (secondM >= 0 && rand() % 5 == 0)  // 以 20% 概率选择次优模式以增强多样性
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

    // 任务修复：按得分进行贪心最优适配
    int nFree = 0;
    for (int j = 0; j < numTask; j++)
        if (taskBeam[j] == -1)
        {
            tmpIdx[nFree] = j;
            tmpVal[nFree] = task_priority(j);
            nFree++;
        }
    if (nFree > 0) qsort_desc(tmpVal, tmpIdx, 0, nFree - 1);

    const double FIT_EPS = 1e-12;
    for (int k = 0; k < nFree; k++)
    {
        int j = tmpIdx[k];

        int bestBeam = -1;
        double bestFit = 1.0e100;
        double bestLeftPW = 1.0e100;
        double bestLeftBW = 1.0e100;

        for (int b = 0; b < numBeam; b++)
        {
            if (!feasible_on(j, b)) continue;

            int bwCap = beamBWCap[b];
            int pwCap = beamPWCap[b] - modeBasePower[beamMode[b]];
            if (bwCap < 1) bwCap = 1;
            if (pwCap < 1) pwCap = 1;

            double leftBW = (double)(remBW[b] - taskBWDemand[j]) / bwCap;
            double leftPW = (double)(remPW[b] - taskPWDemand[j]) / pwCap;
            double fit = leftBW + leftPW;

            int better = (bestBeam < 0 || fit < bestFit - FIT_EPS);
            if (!better && fit <= bestFit + FIT_EPS)
            {
                better = leftPW < bestLeftPW - FIT_EPS ||
                    (leftPW <= bestLeftPW + FIT_EPS &&
                     leftBW < bestLeftBW - FIT_EPS);
            }

            if (better)
            {
                bestBeam   = b;
                bestFit    = fit;
                bestLeftPW = leftPW;
                bestLeftBW = leftBW;
            }
        }
        if (bestBeam >= 0) add_task(j, bestBeam);
    }
    delete[] perturbBeam;
    delete[] prevMode;
    delete[] beamFreqKey;
    delete[] order;
}

void ils()
{
    double beginTime = (double)clock();

    int    *tmpIdx = new int   [numTask];
    double *tmpVal = new double[numTask];

    greedy_init();
    alloc_search();

    // ILS 开始前完整重建 typeProfitSum，保证派生状态一致
    for (int b = 0; b < numBeam; b++)
        for (int t = 0; t < numType; t++) typeProfitSum[b][t] = 0;
    for (int j = 0; j < numTask; j++)
        if (taskBeam[j] >= 0)
            typeProfitSum[taskBeam[j]][taskType[j] - 1] += taskProfit[j];

    globalProfit = -1;

    probe_modes("after_greedy");

    double runTime = 0.0;
    long   numPhase = 0;
    while (runTime < maxRunTime)
    {
        for (int j = 0; j < numTask; j++) moveFreq[j] = 0;

        save_best(beginTime);          // 当前解作为本阶段最好解的起点
        local_search(beginTime, tmpIdx, tmpVal);
        restore_best();                // 混合搜索从最好可行点开始
        local_search_mixed(beginTime, tmpIdx, tmpVal);
        // local_search_mixed 末尾已 restore_best，此处不再重复
        numPhase++;

        if (numPhase == 1) probe_modes("after_LS1");

        if (bestProfit > globalProfit)
        {
            save_global_from_best();
            cout << "  globalProfit=" << globalProfit
                 << "  time=" << globalBestTime << " s"
                 << "  phase=" << numPhase << endl;
        }

        double pStart = (double)clock();
        perturb(tmpIdx, tmpVal);
        g_perturbTime += ((double)clock() - pStart) / CLOCKS_PER_SEC;

        runTime = ((double)clock() - beginTime) / CLOCKS_PER_SEC;
    }

    // 安装全局最好解，供最终校验和输出
    restore_global();
    totalProfit = globalProfit;

    probe_modes("final");

    delete[] tmpIdx;
    delete[] tmpVal;

    cout << "ILS done.  bestProfit=" << globalProfit
         << "  bestTime=" << globalBestTime << " s"
         << "  phases=" << numPhase << endl;
    cout << "Perturbations over the whole run: " << numPhase << endl;

    // 输出性能统计
    double totalLS = g_lsTime + g_mixedTime + g_perturbTime;
    cout << "[PROFILE] phases=" << numPhase
         << "  ls_iters=" << g_lsIters
         << "  moves=" << g_lsMoves
         << "  stalls=" << g_lsStall << endl;
    cout << "[PROFILE] iters_per_phase=" << (numPhase ? (double)g_lsIters / numPhase : 0)
         << "  move_rate=" << (g_lsIters ? 100.0 * g_lsMoves / g_lsIters : 0) << "%"
         << "  stall_rate=" << (g_lsIters ? 100.0 * g_lsStall / g_lsIters : 0) << "%" << endl;
    cout << "[PROFILE] tabu_blocks=" << g_tabuBlock
         << "  aspirations=" << g_aspire << endl;
    cout << "[PROFILE] ls_time=" << g_lsTime << "s (" << (totalLS ? 100.0 * g_lsTime / totalLS : 0) << "%)"
         << "  mixed_time=" << g_mixedTime << "s (" << (totalLS ? 100.0 * g_mixedTime / totalLS : 0) << "%)"
         << "  perturb_time=" << g_perturbTime << "s (" << (totalLS ? 100.0 * g_perturbTime / totalLS : 0) << "%)" << endl;
    cout << "[PROFILE] us_per_ls_iter=" << (g_lsIters ? 1e6 * g_lsTime / g_lsIters : 0) << endl;
    cout << "[PROFILE] mixed_iters=" << g_mixedIters
         << "  mixed_moves=" << g_mixedMoves
         << "  feasible=" << g_mixedFeas
         << " (" << (g_mixedIters ? 100.0 * g_mixedFeas / g_mixedIters : 0) << "%)"
         << "  infeasible=" << g_mixedInfeas
         << " (" << (g_mixedIters ? 100.0 * g_mixedInfeas / g_mixedIters : 0) << "%)"
         << "  best_updates=" << g_mixedBestUpdate
         << "  max_over=" << g_mixedMaxOver
         << "  final_phi_bw=" << mixedPhiBW
         << "  final_phi_pw=" << mixedPhiPW << endl;
    cout << "[PROFILE] mixed_rho final_bw=" << mixedRhoBW
         << " final_pw=" << mixedRhoPW
         << " max_bw=" << g_mixedMaxRhoBW
         << " max_pw=" << g_mixedMaxRhoPW
         << " expands=" << g_mixedRhoExpand
         << " shrinks=" << g_mixedRhoShrink << endl;
    cout << "[PROFILE] move_applied:"
         << " insert=" << g_applyInsert
         << " remove=" << g_applyRemove
         << " swap=" << g_applySwap
         << " flip=" << g_applyFlip
         << " cross=" << g_applyCross << endl;
    cout << "[PROFILE] cross_cand=" << g_crossCand
         << " reloc_cand=" << g_crossRelocCand
         << " swap_cand=" << g_crossSwapCand
         << " filtered=" << g_crossFiltered
         << " reloc_applied=" << g_crossRelocApplied
         << " swap_applied=" << g_crossSwapApplied << endl;
    cout << "[PROFILE] orphan_calls=" << g_orphanCalls
         << " trials=" << g_orphanTrials
         << " steps=" << g_orphanSteps
         << " applied=" << g_orphanApplied
         << " reloc_leaf=" << g_orphanRelocLeaf
         << " drop_leaf=" << g_orphanDropLeaf << endl;
    cout << "[PROFILE] flip_cand=" << g_flipCand
         << "  applied=" << g_flipApplied
         << " (" << (g_lsMoves ? 100.0 * g_flipApplied / g_lsMoves : 0) << "% of moves)" << endl;
    cout << "[PROFILE] flip_applied_by_mode:";
    for (int m = 0; m < numMode && m < 16; m++)
        cout << " " << modeName[m] << "=" << g_flipApplyToMode[m];
    cout << endl;
}

//--------------------------------------------------------------------
// check_solution：校验全部约束及收益一致性
//--------------------------------------------------------------------
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

        // 模式兼容性
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

//--------------------------------------------------------------------
// print_solution：输出便于阅读的解摘要
//--------------------------------------------------------------------
void print_solution()
{
    cout << "\n=== Solution Summary ===" << endl;
    cout << "Total profit : " << totalProfit << endl;

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

//--------------------------------------------------------------------
// free_memory：释放动态内存
//--------------------------------------------------------------------
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
    delete[] tabuUntil;
    delete[] bucketTask;
    delete[] bucketStart;
    delete[] mixedMinBWFree;
    for (int b = 0; b < numBeam; b++) delete[] mixedMinPWFree[b];
    delete[] mixedMinPWFree;
    delete[] orphanWorkBW;
    delete[] orphanWorkPW;
    delete[] orphanTask;
    delete[] orphanDest;
    delete[] orphanBestTask;
    delete[] orphanBestDest;
    delete[] orphanStartMark;
    for (int m = 0; m < numMode; m++) delete[] insertMinPW[m];
    delete[] insertMinPW;

    delete[] tabuBeamMode;
    for (int b = 0; b < numBeam; b++) delete[] typeProfitSum[b];
    delete[] typeProfitSum;
}

//--------------------------------------------------------------------
// main：程序入口
//--------------------------------------------------------------------
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

    maxRunTime = 600.0;          // 时间上限，单位为秒
    if (argc >= 4) maxRunTime = atof(argv[3]);

    double t0 = (double)clock();

    read_instance();
    ils();
    check_solution();
    print_solution();

    double elapsed = ((double)clock() - t0) / CLOCKS_PER_SEC;
    cout << "\nElapsed time: " << elapsed << " s" << endl;

    free_memory();
    return 0;
}
