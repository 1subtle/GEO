#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAXLINE    65536
#define MAXTYPE    20
#define MAXNAMELEN 32
#define MININT_MOVE (-2000000000)

using namespace std;

char *instanceName;

// 算例规模：波束数、任务数、模式数和任务类型数。
int numBeam;
int numTask;
int numMode;
int numType;

// 模式属性以及“任务类型-模式”兼容矩阵。
char **modeName;
int  *modeBasePower;
int **typeCompatMode;

// 每个波束的带宽容量和总功率容量。
int *beamBWCap;
int *beamPWCap;

// 每个任务的类型、收益、带宽需求和功率需求。
int *taskType;
int *taskProfit;
int *taskBWDemand;
int *taskPWDemand;

// 当前解。taskBeam[j] 为 -1 表示任务 j 未服务；
// remPW 已扣除当前模式的基础功率。
int *beamMode;       // 波束当前模式
int *taskBeam;       // 任务当前所属波束
int *remBW;          // 波束剩余带宽，可在不可行域中为负
int *remPW;          // 波束剩余任务功率，可在不可行域中为负
int  totalProfit;    // 当前解总收益

// 运行控制参数。
double maxRunTime;
double bestTime;
int    seed;

// 当前 ILS 阶段内找到的最好可行解。
int *bestBeamMode;
int *bestTaskBeam;
int  bestProfit;

// 整个运行过程中找到的全局最好可行解。
int *globalBeamMode;
int *globalTaskBeam;
int  globalProfit;
double globalBestTime;

// 禁忌搜索状态。moveFreq 同时用于扰动时识别较少探索的波束。
int *moveFreq;       // 当前 ILS 阶段中任务被移动的次数
int *tabuUntil;      // 任务禁忌截止迭代
int  tabuIter;       // 当前禁忌搜索内部迭代号

int  *tabuBeamMode;  // 波束模式翻转的禁忌截止迭代
int **typeProfitSum; // 各波束内不同任务类型的收益总和

// 邻域枚举辅助结构：按波束连续存储已服务任务，并维护快速插入筛选表。
int *bucketTask;
int *bucketStart;
int **insertMinPW;
int  maxBeamBWCap;

// 不可行域搜索参数：phi 是归一化超载惩罚，rho 是单波束最大允许超载比例。
double mixedPhiBW = 1000.0;
double mixedPhiPW = 1000.0;
double mixedRhoBW = 0.15;
double mixedRhoPW = 0.20;
double mixedRhoMinBW = 0.05, mixedRhoMaxBW = 0.30;
double mixedRhoMinPW = 0.08, mixedRhoMaxPW = 0.30;

static char g_line[MAXLINE];

// 根据输入中的模式名称返回模式下标。
int find_mode_index(const char *name)
{
    for (int m = 0; m < numMode; m++)
        if (strcmp(modeName[m], name) == 0)
            return m;
    return -1;
}

// 同步按 val 降序排列 val 和 idx，供贪心排序及扰动排序使用。
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

// 读取算例文件并分配基础数据数组。
void read_instance()
{
    ifstream FIC;
    FIC.open(instanceName);
    if (FIC.fail())
    {
        cout << "Cannot open file: " << instanceName << endl;
        exit(0);
    }

    FIC.getline(g_line, MAXLINE);
    sscanf(g_line, "B=%d N=%d M=%d Cbw=%*d Cpw=%*d",
           &numBeam, &numTask, &numMode);

    modeName     = new char *[numMode];
    for (int m = 0; m < numMode; m++)
        modeName[m] = new char[MAXNAMELEN];
    modeBasePower = new int[numMode];

    FIC.getline(g_line, MAXLINE);
    FIC.getline(g_line, MAXLINE);
    FIC.getline(g_line, MAXLINE);

    for (int m = 0; m < numMode; m++)
    {
        FIC.getline(g_line, MAXLINE);
        char tmp[MAXNAMELEN];
        sscanf(g_line, "%s %d", tmp, &modeBasePower[m]);
        strcpy(modeName[m], tmp);
    }

    FIC.getline(g_line, MAXLINE);
    FIC.getline(g_line, MAXLINE);
    FIC.getline(g_line, MAXLINE);

    static char typeLines[MAXTYPE][256];
    numType = 0;
    while (FIC.getline(g_line, MAXLINE))
    {
        if ((int)strlen(g_line) == 0) break;
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

    FIC.getline(g_line, MAXLINE);
    for (int b = 0; b < numBeam; b++) FIC >> beamBWCap[b];
    FIC.getline(g_line, MAXLINE);
    FIC.getline(g_line, MAXLINE);

    FIC.getline(g_line, MAXLINE);
    for (int b = 0; b < numBeam; b++) FIC >> beamPWCap[b];
    FIC.getline(g_line, MAXLINE);
    FIC.getline(g_line, MAXLINE);

    taskType     = new int[numTask];
    taskProfit   = new int[numTask];
    taskBWDemand = new int[numTask];
    taskPWDemand = new int[numTask];

    FIC.getline(g_line, MAXLINE);
    for (int j = 0; j < numTask; j++) FIC >> taskType[j];
    FIC.getline(g_line, MAXLINE);
    FIC.getline(g_line, MAXLINE);

    FIC.getline(g_line, MAXLINE);
    for (int j = 0; j < numTask; j++) FIC >> taskProfit[j];
    FIC.getline(g_line, MAXLINE);
    FIC.getline(g_line, MAXLINE);

    FIC.getline(g_line, MAXLINE);
    for (int j = 0; j < numTask; j++) FIC >> taskBWDemand[j];
    FIC.getline(g_line, MAXLINE);
    FIC.getline(g_line, MAXLINE);

    FIC.getline(g_line, MAXLINE);
    for (int j = 0; j < numTask; j++) FIC >> taskPWDemand[j];

    FIC.close();

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

// 贪心评分：单位归一化 BW/PW 占用能够获得的任务收益。
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

// 返回任务在当前所有兼容波束上的最高贪心评分。
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

// 逐波束枚举模式，并按任务评分贪心装填，构造初始可行解。
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

    for (int b = 0; b < numBeam; b++)
    {
        int bestM      = -1;
        int bestGain   = -1;
        int nBestAdded = 0;

        for (int m = 0; m < numMode; m++)
        {
            if (modeBasePower[m] > beamPWCap[b]) continue;

            // 在该模式下对尚未服务且类型兼容的任务排序。
            int nFree = 0;
            for (int j = 0; j < numTask; j++)
                if (taskBeam[j] == -1 && typeCompatMode[taskType[j] - 1][m])
                {
                    tmpIdx[nFree] = j;
                    tmpVal[nFree] = task_score(j, b, m);
                    nFree++;
                }
            if (nFree > 0) qsort_desc(tmpVal, tmpIdx, 0, nFree - 1);

            // 按评分从高到低装填，得到该模式对应的试探收益。
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

        // 若没有功率可行模式，则选择基础功率最低的模式并保持空波束。
        if (bestM < 0)
        {
            for (int m = 0; m < numMode; m++)
                if (bestM < 0 || modeBasePower[m] < modeBasePower[bestM]) bestM = m;
            nBestAdded = 0;
        }

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

// 判断任务 j 能否在当前解中直接插入波束 b。
int feasible_on(int j, int b)
{
    if (beamMode[b] < 0)                 return 0;
    int t = taskType[j] - 1;
    if (!typeCompatMode[t][beamMode[b]]) return 0;
    if (taskBWDemand[j] > remBW[b])      return 0;
    if (taskPWDemand[j] > remPW[b])      return 0;
    return 1;
}

// 将带宽超载量按波束带宽容量归一化；未超载时为 0。
double beam_bw_over_with_rem(int b, int bwFree)
{
    if (bwFree >= 0) return 0.0;
    int denom = beamBWCap[b];
    if (denom < 1) denom = 1;
    return (double)(-bwFree) / denom;
}

// 将功率超载量按扣除模式基础功率后的任务功率容量归一化。
double beam_pw_over_with_rem_mode(int b, int m, int pwFree)
{
    if (pwFree >= 0) return 0.0;
    int denom = beamPWCap[b] - modeBasePower[m];
    if (denom < 1) denom = 1;
    return (double)(-pwFree) / denom;
}

// 汇总当前解所有波束的归一化 BW/PW violation。
void total_over_parts(double &bwOver, double &pwOver)
{
    bwOver = 0.0;
    pwOver = 0.0;
    for (int b = 0; b < numBeam; b++)
    {
        if (beamMode[b] < 0) continue;
        bwOver += beam_bw_over_with_rem(b, remBW[b]);
        pwOver += beam_pw_over_with_rem_mode(b, beamMode[b], remPW[b]);
    }
}

// 将数值限制在闭区间 [lo, hi]。
double clamp_double(double x, double lo, double hi)
{
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

// 根据算例任务密度和当前解资源紧张度设置不可行域初始边界。
// 该函数只会把 rho 提升到目标值，不会主动缩小已有边界。
void set_mixed_rho_from_instance()
{
    double density = numBeam > 0 ? (double)numTask / numBeam : 0.0;
    double densityAdj = (density - 10.0) * 0.004;
    densityAdj = clamp_double(densityAdj, 0.0, 0.08);

    // 紧张度取当前解 BW 使用率和 PW 使用率中的较大值。
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
    double tightness = bwTight > pwTight ? bwTight : pwTight;
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
}

// 长时间没有改进时扩张不可行域边界，增强发散能力。
void expand_mixed_rho()
{
    mixedRhoBW = clamp_double(mixedRhoBW * 1.10 + 0.01,
                              mixedRhoMinBW, mixedRhoMaxBW);
    mixedRhoPW = clamp_double(mixedRhoPW * 1.10 + 0.01,
                              mixedRhoMinPW, mixedRhoMaxPW);
}

// 找到新的可行最优解后收缩不可行域边界，加强可行域附近搜索。
void shrink_mixed_rho()
{
    mixedRhoBW = clamp_double(mixedRhoBW * 0.85,
                              mixedRhoMinBW, mixedRhoMaxBW);
    mixedRhoPW = clamp_double(mixedRhoPW * 0.85,
                              mixedRhoMinPW, mixedRhoMaxPW);
}

// 检查候选动作执行后的单波束超载是否仍在 rho 硬边界内。
int relaxed_rem_ok_mode(int b, int m, int bwFree, int pwFree)
{
    if (bwFree < 0)
    {
        int denom = beamBWCap[b];
        if (denom < 1) denom = 1;
        if ((double)(-bwFree) / denom > mixedRhoBW) return 0;
    }

    if (pwFree < 0)
    {
        int denom = beamPWCap[b] - modeBasePower[m];
        if (denom < 1) denom = 1;
        if ((double)(-pwFree) / denom > mixedRhoPW) return 0;
    }

    return 1;
}

// 使用波束当前模式检查不可行域边界。
int relaxed_rem_ok(int b, int bwFree, int pwFree)
{
    return relaxed_rem_ok_mode(b, beamMode[b], bwFree, pwFree);
}

// 为“交换后是否留出可插入槽位”构建快速查询表。
// insertMinPW[m][w] 表示模式 m 下、带宽需求不超过 w 的未服务任务最小功率需求。
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

// 判断给定剩余资源是否至少能容纳一个当前未服务任务。
int beam_can_insert_unserved_with_rem(int b, int bwFree, int pwFree)
{
    if (beamMode[b] < 0) return 0;
    if (bwFree < 0 || pwFree < 0) return 0;
    if (bwFree > maxBeamBWCap) bwFree = maxBeamBWCap;

    return insertMinPW[beamMode[b]][bwFree] <= pwFree;
}

// 将任务加入波束，并同步更新资源、收益和类型收益统计。
void add_task(int j, int b)
{
    taskBeam[j]  = b;
    remBW[b]    -= taskBWDemand[j];
    remPW[b]    -= taskPWDemand[j];
    totalProfit += taskProfit[j];
    typeProfitSum[b][taskType[j] - 1] += taskProfit[j];
}

// 从当前波束移除任务，并同步恢复所有派生状态。
void remove_task(int j)
{
    int b = taskBeam[j];
    remBW[b]    += taskBWDemand[j];
    remPW[b]    += taskPWDemand[j];
    totalProfit -= taskProfit[j];
    typeProfitSum[b][taskType[j] - 1] -= taskProfit[j];
    taskBeam[j]  = -1;
}

// 禁忌期限包含随机项和随波束规模增长的项。
int tabu_tenure()
{
    return 5 + rand() % 6 + numBeam / 10;
}

// 保存当前 ILS 阶段内的最好可行解。
void save_best(double beginTime)
{
    bestTime    = ((double)clock() - beginTime) / CLOCKS_PER_SEC;
    bestProfit  = totalProfit;
    for (int b = 0; b < numBeam; b++) bestBeamMode[b] = beamMode[b];
    for (int j = 0; j < numTask; j++) bestTaskBeam[j] = taskBeam[j];
}

// 恢复阶段最好解，并从任务分配重新计算剩余资源和派生统计。
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

// 将阶段最好解提升为全局最好解。
void save_global_from_best()
{
    globalBestTime = bestTime;
    globalProfit   = bestProfit;
    for (int b = 0; b < numBeam; b++) globalBeamMode[b] = bestBeamMode[b];
    for (int j = 0; j < numTask; j++) globalTaskBeam[j] = bestTaskBeam[j];
}

// 安装全局最好解，供最终校验和输出使用。
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

// 分配禁忌搜索、快照和邻域加速结构所需内存。
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

// 记录可行域中收益增量最大的候选；并列候选使用蓄水池抽样随机打破平局。
void record_candidate(int kind, int a, int b, int c, int delta,
                      int &bestDelta, int &numBest,
                      int &chosenKind, int &chosenA, int &chosenB, int &chosenC)
{
    if (delta < bestDelta) return;

    if (delta > bestDelta)
    {
        bestDelta = delta;
        numBest   = 1;
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

// 记录不可行域中评价函数增量最大的候选，并随机处理并列最优候选。
void record_mixed_candidate(int kind, int a, int b, int c,
                            double deltaEval,
                            double &bestDeltaEval,
                            int &numBest,
                            int &chosenKind, int &chosenA,
                            int &chosenB, int &chosenC)
{
    if (deltaEval < bestDeltaEval) return;

    if (deltaEval > bestDeltaEval)
    {
        bestDeltaEval  = deltaEval;
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
        chosenKind = kind;
        chosenA    = a;
        chosenB    = b;
        chosenC    = c;
    }
}

// 翻转波束模式，只驱逐与新模式不兼容的任务。
// 被驱逐任务记为发生移动并进入禁忌；仍留在波束上的任务不记移动。
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

// 第一阶段：始终保持资源可行的五邻域禁忌搜索。
void local_search(double beginTime)
{
    int ts_depth = 300;
    int nonImprove = 0;

    // 每次搜索调用拥有独立禁忌周期。
    for (int j = 0; j < numTask; j++) tabuUntil[j] = 0;
    for (int b = 0; b < numBeam; b++) tabuBeamMode[b] = 0;
    tabuIter = 0;

    while (nonImprove < ts_depth)
    {
        if (((double)clock() - beginTime) / CLOCKS_PER_SEC > maxRunTime)
            break;

        int bestDelta = MININT_MOVE;
        int numBest   = 0;
        int kind = 0, a = -1, bb = -1, cc = -1;

        // 重建“波束 -> 已服务任务”的 CSR 索引，减少交换邻域扫描量。
        for (int b = 0; b <= numBeam; b++) bucketStart[b] = 0;
        for (int j = 0; j < numTask; j++)
            if (taskBeam[j] >= 0) bucketStart[taskBeam[j] + 1]++;
        for (int b = 0; b < numBeam; b++) bucketStart[b + 1] += bucketStart[b];
        // 邻域 1：将一个未服务任务插入某个可行波束。
        for (int j = 0; j < numTask; j++)
            if (taskBeam[j] >= 0)
            {
                int b = taskBeam[j];
                bucketTask[bucketStart[b]++] = j;
            }

        for (int b = numBeam; b > 0; b--) bucketStart[b] = bucketStart[b - 1];
        bucketStart[0] = 0;

        for (int j = 0; j < numTask; j++)
        {
            if (taskBeam[j] != -1) continue;
            int delta = taskProfit[j];

            // 特赦：即使任务处于禁忌，只要候选超过阶段最好收益仍允许执行。
            int jTabu = (tabuIter < tabuUntil[j]);
            if (jTabu)
            {
                if (totalProfit + delta > bestProfit) {  }
                else { continue; }
            }
            for (int b = 0; b < numBeam; b++)
            {
                if (!feasible_on(j, b)) continue;
                if (delta < bestDelta) continue;
                record_candidate(1, j, b, -1, delta,
                                 bestDelta, numBest, kind, a, bb, cc);
            }
        }

        // 邻域 2：移除一个已服务任务。
        for (int j = 0; j < numTask; j++)
        {
            if (taskBeam[j] < 0) continue;
            int delta = -taskProfit[j];
            if (tabuIter < tabuUntil[j])
            {
                if (totalProfit + delta > bestProfit) {  }
                else { continue; }
            }
            if (delta < bestDelta) continue;

            record_candidate(2, j, -1, -1, delta,
                             bestDelta, numBest, kind, a, bb, cc);
        }

        // 邻域 3：同一波束内用未服务任务 i 替换已服务任务 k。
        for (int i = 0; i < numTask; i++)
        {
            if (taskBeam[i] != -1) continue;
            int iTabu = (tabuIter < tabuUntil[i]);
            int ti = taskType[i] - 1;
            for (int b = 0; b < numBeam; b++)
            {
                if (beamMode[b] < 0) continue;
                if (!typeCompatMode[ti][beamMode[b]]) continue;

                int bwFree0 = remBW[b];
                int pwFree0 = remPW[b];
                for (int idx = bucketStart[b]; idx < bucketStart[b + 1]; idx++)
                {
                    int k = bucketTask[idx];

                    int delta = taskProfit[i] - taskProfit[k];
                    if (delta < bestDelta) continue;

                    if (iTabu || tabuIter < tabuUntil[k])
                    {
                        if (totalProfit + delta > bestProfit) {  }
                        else continue;
                    }

                    int bwFree = bwFree0 + taskBWDemand[k];
                    int pwFree = pwFree0 + taskPWDemand[k];
                    if (taskBWDemand[i] > bwFree || taskPWDemand[i] > pwFree) continue;

                    record_candidate(3, i, b, k, delta,
                                     bestDelta, numBest, kind, a, bb, cc);
                }
            }
        }

        // 邻域 4：翻转波束模式，并计算驱逐不兼容任务造成的收益损失。
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
                    if (totalProfit + delta > bestProfit) {  }
                    else { continue; }
                }

                if (delta < bestDelta) continue;
                record_candidate(4, b, m2, -1, delta,
                                 bestDelta, numBest, kind, a, bb, cc);
            }
        }

        // 邻域 5：跨波束重定位或交换。仅在没有正收益动作时启用，
        // 并要求动作后至少一个相关波束能插入某个未服务任务。
        if (bestDelta <= 0)
        {
            rebuild_insert_filter();

            for (int i = 0; i < numTask; i++)
            {
                int b1 = taskBeam[i];
                if (b1 < 0) continue;

                int iTabu = (tabuIter < tabuUntil[i]);
                if (iTabu)
                {
                    if (totalProfit > bestProfit) {  }
                    else { continue; }
                }

                int ti = taskType[i] - 1;
                for (int b2 = 0; b2 < numBeam; b2++)
                {
                    if (b2 == b1) continue;
                    if (beamMode[b2] < 0) continue;
                    if (!typeCompatMode[ti][beamMode[b2]]) continue;

                    int delta = 0;

                    // 实任务-虚拟空位：把任务 i 从 b1 重定位到 b2。
                    if (taskBWDemand[i] <= remBW[b2] &&
                        taskPWDemand[i] <= remPW[b2])
                    {
                        if (beam_can_insert_unserved_with_rem(
                                b1,
                                remBW[b1] + taskBWDemand[i],
                                remPW[b1] + taskPWDemand[i]))
                        {
                            record_candidate(5, i, b2, -1, delta,
                                             bestDelta, numBest, kind, a, bb, cc);
                        }
                    }

                    // 实任务-实任务：交换两个不同波束上的任务，每对只枚举一次。
                    if (b1 > b2) continue;

                    for (int idx = bucketStart[b2]; idx < bucketStart[b2 + 1]; idx++)
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
                            if (totalProfit > bestProfit) {  }
                            else { continue; }
                        }

                        int b1BW = remBW[b1] + taskBWDemand[i] - taskBWDemand[k];
                        int b1PW = remPW[b1] + taskPWDemand[i] - taskPWDemand[k];
                        int b2BW = remBW[b2] + taskBWDemand[k] - taskBWDemand[i];
                        int b2PW = remPW[b2] + taskPWDemand[k] - taskPWDemand[i];
                        if (!beam_can_insert_unserved_with_rem(b1, b1BW, b1PW) &&
                            !beam_can_insert_unserved_with_rem(b2, b2BW, b2PW))
                        {
                            continue;
                        }

                        record_candidate(5, i, b2, k, delta,
                                         bestDelta, numBest, kind, a, bb, cc);
                    }
                }
            }
        }

        // 没有可接受动作时只推进禁忌时间和停滞计数。
        if (numBest == 0)
        {
            nonImprove++;
            tabuIter++;
            continue;
        }
        // 提交本轮最优候选，并对所有实际移动的任务设置禁忌期限。
        if (kind == 1)
        {
            add_task(a, bb);
            moveFreq[a]++;
            tabuUntil[a] = tabuIter + tabu_tenure();
        }
        else if (kind == 2)
        {
            remove_task(a);
            moveFreq[a]++;
            tabuUntil[a] = tabuIter + tabu_tenure();
        }
        else if (kind == 3)
        {
            remove_task(cc);
            add_task(a, bb);
            moveFreq[a]++; moveFreq[cc]++;
            tabuUntil[a]  = tabuIter + tabu_tenure();
            tabuUntil[cc] = tabuIter + tabu_tenure();
        }
        else if (kind == 4)
        {
            int tenure = tabu_tenure();
            apply_mode_flip(a, bb, tenure);
            tabuBeamMode[a] = tabuIter + tenure;
        }
        else if (kind == 5)
        {
            int b1 = taskBeam[a];
            if (cc < 0)
            {
                remove_task(a);
                add_task(a, bb);
                moveFreq[a]++;
                tabuUntil[a] = tabuIter + tabu_tenure();
            }
            else
            {
                int b2 = taskBeam[cc];
                remove_task(a);
                remove_task(cc);
                add_task(a, b2);
                add_task(cc, b1);
                moveFreq[a]++; moveFreq[cc]++;
                tabuUntil[a]  = tabuIter + tabu_tenure();
                tabuUntil[cc] = tabuIter + tabu_tenure();
            }
        }
        // 只有严格改善阶段最好收益时才重置停滞计数。
        if (totalProfit > bestProfit)
        {
            save_best(beginTime);
            nonImprove = 0;
        }
        else
            nonImprove++;

        tabuIter++;
    }
}

// 第二阶段：允许 BW/PW 在 rho 边界内超载的混合禁忌搜索。
// 候选评价增量为：
// deltaProfit - phiBW * deltaOverBW - phiPW * deltaOverPW。
// 只有回到可行域且收益改善时，才更新阶段最好解。
void local_search_mixed(double beginTime)
{
    const double EPS = 1e-12;
    const int mixedDepth = 300;
    const int phiWindow = 5;
    const double phiTau = 1.5;
    const double phiMin = 100.0;
    const double phiMax = 100000.0;

    int nonImprove = 0;
    // BW 和 PW 分别记录连续可行/不可行状态，用于独立调整两个 phi。
    int bwFeasibleStreak = 0, bwInfeasibleStreak = 0;
    int pwFeasibleStreak = 0, pwInfeasibleStreak = 0;
    double curOverBW, curOverPW;
    total_over_parts(curOverBW, curOverPW);
    double curOver = curOverBW + curOverPW;
    set_mixed_rho_from_instance();

    // 混合搜索重新开启一个独立禁忌周期。
    for (int j = 0; j < numTask; j++) tabuUntil[j] = 0;
    for (int b = 0; b < numBeam; b++) tabuBeamMode[b] = 0;
    tabuIter = 0;

    while (nonImprove < mixedDepth)
    {
        if (((double)clock() - beginTime) / CLOCKS_PER_SEC > maxRunTime)
            break;
        // 连续 5 代 BW 可行则降低 BW 惩罚，连续 5 代不可行则提高惩罚。
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

        // PW 惩罚按同样规则独立振荡。
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

        double bestDeltaEval = -1.0e100;
        int numBest = 0;
        int kind = 0, a = -1, bb = -1, cc = -1;

        // 重建当前已服务任务的波束索引。
        for (int b = 0; b <= numBeam; b++) bucketStart[b] = 0;
        // 邻域 1：插入未服务任务，允许动作后资源在 rho 范围内超载。
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
                // rho 是候选动作执行后的硬边界，越界候选不进入评价。
                if (!relaxed_rem_ok(b, bw2, pw2)) continue;

                double newOverBW = curOverBW - beam_bw_over_with_rem(b, remBW[b])
                                   + beam_bw_over_with_rem(b, bw2);
                double newOverPW = curOverPW
                                   - beam_pw_over_with_rem_mode(b, beamMode[b], remPW[b])
                                   + beam_pw_over_with_rem_mode(b, beamMode[b], pw2);
                double newOver = newOverBW + newOverPW;
                // 混合阶段的特赦必须同时满足：候选可行且超过阶段最好收益。
                if (tabuIter < tabuUntil[j])
                {
                    if (newOver <= EPS && totalProfit + deltaProfit > bestProfit) {  }
                    else { continue; }
                }

                double deltaEval = (double)deltaProfit
                    - mixedPhiBW * (newOverBW - curOverBW)
                    - mixedPhiPW * (newOverPW - curOverPW);
                record_mixed_candidate(1, j, b, -1, deltaEval,
                                       bestDeltaEval, numBest,
                                       kind, a, bb, cc);
            }
        }

        // 邻域 2：移除已服务任务，通常用于降低收益换取减少 violation。
        for (int j = 0; j < numTask; j++)
        {
            int b = taskBeam[j];
            if (b < 0) continue;

            int deltaProfit = -taskProfit[j];
            int bw2 = remBW[b] + taskBWDemand[j];
            int pw2 = remPW[b] + taskPWDemand[j];
            double newOverBW = curOverBW - beam_bw_over_with_rem(b, remBW[b])
                               + beam_bw_over_with_rem(b, bw2);
            double newOverPW = curOverPW
                               - beam_pw_over_with_rem_mode(b, beamMode[b], remPW[b])
                               + beam_pw_over_with_rem_mode(b, beamMode[b], pw2);
            double newOver = newOverBW + newOverPW;

            if (tabuIter < tabuUntil[j])
            {
                if (newOver <= EPS && totalProfit + deltaProfit > bestProfit) {  }
                else { continue; }
            }

            double deltaEval = (double)deltaProfit
                - mixedPhiBW * (newOverBW - curOverBW)
                - mixedPhiPW * (newOverPW - curOverPW);
            record_mixed_candidate(2, j, -1, -1, deltaEval,
                                   bestDeltaEval, numBest,
                                   kind, a, bb, cc);
        }

        // 邻域 3：同波束替换。
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

                    double newOverBW = curOverBW - beam_bw_over_with_rem(b, remBW[b])
                                       + beam_bw_over_with_rem(b, bw2);
                    double newOverPW = curOverPW
                                       - beam_pw_over_with_rem_mode(b, beamMode[b], remPW[b])
                                       + beam_pw_over_with_rem_mode(b, beamMode[b], pw2);
                    double newOver = newOverBW + newOverPW;
                    if (tabuIter < tabuUntil[i] || tabuIter < tabuUntil[k])
                    {
                        if (newOver <= EPS && totalProfit + deltaProfit > bestProfit) {  }
                        else { continue; }
                    }

                    double deltaEval = (double)deltaProfit
                        - mixedPhiBW * (newOverBW - curOverBW)
                        - mixedPhiPW * (newOverPW - curOverPW);
                    record_mixed_candidate(3, i, b, k, deltaEval,
                                           bestDeltaEval, numBest,
                                           kind, a, bb, cc);
                }
            }
        }

        // 邻域 4：模式翻转，候选资源状态按新模式的功率容量计算。
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
                double newOverBW = curOverBW - beam_bw_over_with_rem(b, remBW[b])
                                   + beam_bw_over_with_rem(b, bw2);
                double newOverPW = curOverPW
                                   - beam_pw_over_with_rem_mode(b, beamMode[b], remPW[b])
                                   + beam_pw_over_with_rem_mode(b, m2, pw2);
                double newOver = newOverBW + newOverPW;
                if (tabuIter < tabuBeamMode[b] || evictTabu)
                {
                    if (newOver <= EPS && totalProfit + deltaProfit > bestProfit) {  }
                    else { continue; }
                }

                double deltaEval = (double)deltaProfit
                    - mixedPhiBW * (newOverBW - curOverBW)
                    - mixedPhiPW * (newOverPW - curOverPW);
                record_mixed_candidate(4, b, m2, -1, deltaEval,
                                       bestDeltaEval, numBest,
                                       kind, a, bb, cc);
            }
        }

        // 邻域 5：跨波束重定位或交换。该阶段不要求动作后留出插入槽位。
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
                // 实任务-虚拟空位重定位。
                if (relaxed_rem_ok(b1, b1BW, b1PW) &&
                    relaxed_rem_ok(b2, b2BW, b2PW))
                {
                    double newOverBW = curOverBW
                        - beam_bw_over_with_rem(b1, remBW[b1])
                        - beam_bw_over_with_rem(b2, remBW[b2])
                        + beam_bw_over_with_rem(b1, b1BW)
                        + beam_bw_over_with_rem(b2, b2BW);
                    double newOverPW = curOverPW
                        - beam_pw_over_with_rem_mode(b1, beamMode[b1], remPW[b1])
                        - beam_pw_over_with_rem_mode(b2, beamMode[b2], remPW[b2])
                        + beam_pw_over_with_rem_mode(b1, beamMode[b1], b1PW)
                        + beam_pw_over_with_rem_mode(b2, beamMode[b2], b2PW);
                    double newOver = newOverBW + newOverPW;

                    if (tabuIter < tabuUntil[i])
                    {
                        if (newOver <= EPS && totalProfit > bestProfit) {  }
                        else { continue; }
                    }

                    double deltaEval = -mixedPhiBW * (newOverBW - curOverBW)
                                       -mixedPhiPW * (newOverPW - curOverPW);
                    record_mixed_candidate(5, i, b2, -1, deltaEval,
                                           bestDeltaEval, numBest,
                                           kind, a, bb, cc);
                }

                // 实任务-实任务交换。
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
                        - beam_bw_over_with_rem(b1, remBW[b1])
                        - beam_bw_over_with_rem(b2, remBW[b2])
                        + beam_bw_over_with_rem(b1, b1BW)
                        + beam_bw_over_with_rem(b2, b2BW);
                    double newOverPW = curOverPW
                        - beam_pw_over_with_rem_mode(b1, beamMode[b1], remPW[b1])
                        - beam_pw_over_with_rem_mode(b2, beamMode[b2], remPW[b2])
                        + beam_pw_over_with_rem_mode(b1, beamMode[b1], b1PW)
                        + beam_pw_over_with_rem_mode(b2, beamMode[b2], b2PW);
                    double newOver = newOverBW + newOverPW;

                    if (tabuIter < tabuUntil[i] || tabuIter < tabuUntil[k])
                    {
                        if (newOver <= EPS && totalProfit > bestProfit) {  }
                        else { continue; }
                    }

                    double deltaEval = -mixedPhiBW * (newOverBW - curOverBW)
                                       -mixedPhiPW * (newOverPW - curOverPW);
                    record_mixed_candidate(5, i, b2, k, deltaEval,
                                           bestDeltaEval, numBest,
                                           kind, a, bb, cc);
                }
            }
        }

        // 无候选或持续无改进时，每 50 代扩张一次 rho。
        if (numBest == 0)
        {
            nonImprove++;
            if (nonImprove > 0 && nonImprove % 50 == 0)
                expand_mixed_rho();
            tabuIter++;
            continue;
        }

        // 提交评价函数增量最大的候选，并设置任务/模式禁忌。
        if (kind == 1)
        {
            add_task(a, bb);
            moveFreq[a]++;
            tabuUntil[a] = tabuIter + tabu_tenure();
        }
        else if (kind == 2)
        {
            remove_task(a);
            moveFreq[a]++;
            tabuUntil[a] = tabuIter + tabu_tenure();
        }
        else if (kind == 3)
        {
            remove_task(cc);
            add_task(a, bb);
            moveFreq[a]++; moveFreq[cc]++;
            tabuUntil[a]  = tabuIter + tabu_tenure();
            tabuUntil[cc] = tabuIter + tabu_tenure();
        }
        else if (kind == 4)
        {
            int tenure = tabu_tenure();
            apply_mode_flip(a, bb, tenure);
            tabuBeamMode[a] = tabuIter + tenure;
        }
        else if (kind == 5)
        {
            int b1 = taskBeam[a];
            if (cc < 0)
            {
                remove_task(a);
                add_task(a, bb);
                moveFreq[a]++;
                tabuUntil[a] = tabuIter + tabu_tenure();
            }
            else
            {
                int b2 = taskBeam[cc];
                remove_task(a);
                remove_task(cc);
                add_task(a, b2);
                add_task(cc, b1);
                moveFreq[a]++; moveFreq[cc]++;
                tabuUntil[a]  = tabuIter + tabu_tenure();
                tabuUntil[cc] = tabuIter + tabu_tenure();
            }
        }

        // 动作执行后重新计算总 violation。
        total_over_parts(curOverBW, curOverPW);
        curOver = curOverBW + curOverPW;
        // 只保存可行改进；成功回到更好的可行解后收缩 rho。
        if (curOver <= EPS && totalProfit > bestProfit)
        {
            save_best(beginTime);
            shrink_mixed_rho();
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

    // 不可行的末状态不向后传递，阶段结束统一恢复最好可行解。
    restore_best();
}

// ILS 扰动：从阶段最好解出发，执行“破坏波束 + 重选模式 + 贪心修复”。
void perturb(int *tmpIdx, double *tmpVal)
{
    restore_best();
    if (numBeam == 0) return;

    // 每次破坏约 30% 的波束，至少破坏一个。
    int numDestroy = (int)(0.30 * numBeam + 0.5);
    if (numDestroy < 1)        numDestroy = 1;
    if (numDestroy > numBeam)  numDestroy = numBeam;

    // 按波束内任务移动频率之和排序，优先破坏探索较少的波束。
    int    *perturbBeam = new int[numBeam];
    int    *prevMode    = new int[numBeam];
    double *beamFreqKey = new double[numBeam];
    int    *order       = new int[numBeam];
    for (int b = 0; b < numBeam; b++) { beamFreqKey[b] = 0.0; order[b] = b; }
    for (int j = 0; j < numTask; j++)
        if (taskBeam[j] >= 0) beamFreqKey[taskBeam[j]] += moveFreq[j];

    for (int b = 0; b < numBeam; b++) beamFreqKey[b] = -beamFreqKey[b];
    qsort_desc(beamFreqKey, order, 0, numBeam - 1);

    int poolSize = 0;
    for (int b = 0; b < numBeam; b++) if (beamFreqKey[b] == 0.0) poolSize++;
    if (poolSize < numDestroy) poolSize = numDestroy;

    // 从低频候选池中无放回随机选择待破坏波束。
    for (int k = 0; k < numDestroy; k++)
    {
        int r   = k + rand() % (poolSize - k);
        int tmp = order[k]; order[k] = order[r]; order[r] = tmp;
        perturbBeam[k] = order[k];
    }

    // 破坏：移除选中波束上的全部任务，并暂时取消其模式。
    for (int k = 0; k < numDestroy; k++)
    {
        int b = perturbBeam[k];
        prevMode[k] = beamMode[b];
        for (int j = 0; j < numTask; j++)
            if (taskBeam[j] == b) remove_task(j);
        beamMode[b] = -1;
        remBW[b]    = beamBWCap[b];
        remPW[b]    = 0;
    }

    // 模式修复：逐个波束枚举模式，并试探性贪心装填当前未服务任务。
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
            if (modeBasePower[m] > beamPWCap[b]) continue;

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

        if (bestM < 0)
        {
            bestM = -1;
            for (int m = 0; m < numMode; m++)
            {
                if (bestM < 0 || modeBasePower[m] < modeBasePower[bestM]) bestM = m;
            }
            if (bestM < 0) bestM = prevMode[k];
        }

        int chosenM      = bestM;
        int *chosenAdded = bestAdded;
        int nChosenAdded = nBestAdded;
        // 以 20% 概率选用次优模式，避免修复过程完全确定化。
        if (secondM >= 0 && rand() % 5 == 0)
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

    // 任务修复：对仍未服务的任务做一次全局贪心回填。
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

        // 先尝试基础功率较低的兼容模式，同一模式内选择剩余 BW 最紧的波束。
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

// 迭代局部搜索主流程：
// 贪心初解 -> 可行域禁忌搜索 -> 混合可行/不可行禁忌搜索 -> 扰动。
void ils()
{
    double beginTime = (double)clock();

    int    *tmpIdx = new int   [numTask];
    double *tmpVal = new double[numTask];

    greedy_init();
    alloc_search();

    // 根据贪心初解重建类型收益统计，保证派生状态一致。
    for (int b = 0; b < numBeam; b++)
        for (int t = 0; t < numType; t++) typeProfitSum[b][t] = 0;
    for (int j = 0; j < numTask; j++)
        if (taskBeam[j] >= 0)
            typeProfitSum[taskBeam[j]][taskType[j] - 1] += taskProfit[j];

    globalProfit = -1;

    double runTime = 0.0;
    while (runTime < maxRunTime)
    {
        // moveFreq 在一个完整 ILS 阶段内累计，供阶段末扰动使用。
        for (int j = 0; j < numTask; j++) moveFreq[j] = 0;

        // 两段搜索共享阶段最好可行解，但各自重新初始化禁忌周期。
        save_best(beginTime);
        local_search(beginTime);
        restore_best();
        local_search_mixed(beginTime);
        restore_best();

        if (bestProfit > globalProfit)
        {
            save_global_from_best();
            cout << "  globalProfit=" << globalProfit
                 << "  time=" << globalBestTime << " s" << endl;
        }

        // 从阶段最好可行解出发产生下一阶段起点。
        perturb(tmpIdx, tmpVal);

        runTime = ((double)clock() - beginTime) / CLOCKS_PER_SEC;
    }

    // 时间结束后恢复全局最好可行解。
    restore_global();
    totalProfit = globalProfit;

    delete[] tmpIdx;
    delete[] tmpVal;

    cout << "ILS done.  bestProfit=" << globalProfit
         << "  bestTime=" << globalBestTime << " s" << endl;

}

// 重新计算最终解的兼容性、BW/PW 容量和总收益，防止增量状态失配。
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

// 输出最终收益、服务任务数以及每个波束的资源使用情况。
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

// 释放运行期间分配的全部动态内存。
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
    for (int m = 0; m < numMode; m++) delete[] insertMinPW[m];
    delete[] insertMinPW;

    delete[] tabuBeamMode;
    for (int b = 0; b < numBeam; b++) delete[] typeProfitSum[b];
    delete[] typeProfitSum;
}

// 命令行参数：算例路径、随机种子、可选时间上限（秒）。
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

    maxRunTime = 600.0;
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
