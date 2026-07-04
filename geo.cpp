#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAXLINE    65536   // max bytes per line in instance file
#define MAXTYPE    20      // max number of task types
#define MAXNAMELEN 32      // max length of mode name
#define MININT_MOVE (-2000000000)   // "minus infinity" delta seed for tabu search

using namespace std;

char *instanceName;

int numBeam;          // B: number of beams
int numTask;          // N: number of tasks
int numMode;          // M: number of service modes
int numType;          // T: number of task types

char **modeName;      // modeName[m]             : name string of mode m
int  *modeBasePower;  // modeBasePower[m]         : fixed power consumed by mode m
int **typeCompatMode; // typeCompatMode[t][m] = 1 if task-type (t+1) can use mode m

int *beamBWCap;       // beamBWCap[b]   : bandwidth capacity of beam b
int *beamPWCap;       // beamPWCap[b]   : power capacity of beam b

int *taskType;        // taskType[j]    : type of task j (1-indexed)
int *taskProfit;      // taskProfit[j]  : profit of task j
int *taskBWDemand;    // taskBWDemand[j]: bandwidth demand of task j
int *taskPWDemand;    // taskPWDemand[j]: power demand of task j

int *beamMode;   // beamMode[b]   : mode index chosen for beam b
int *taskBeam;   // taskBeam[j]   : beam serving task j  (-1 = not served)
int *remBW;      // remBW[b]      : remaining bandwidth of beam b
int *remPW;      // remPW[b]      : remaining power of beam b (after base power)
int  totalProfit;

double maxRunTime;   // time limit in seconds
double bestTime;     // time (s) at which the phase-best solution was found
int    seed;         // random seed

int *bestBeamMode;   // snapshot of the best solution found in the current ILS phase
int *bestTaskBeam;
int  bestProfit;

int *globalBeamMode; // snapshot of the best solution found across all ILS phases
int *globalTaskBeam;
int  globalProfit;
double globalBestTime;

int *moveFreq;       // moveFreq[j]  : how often task j has been moved (for perturbation)
int *tabuUntil;      // tabuUntil[j] : iteration index until which task j is tabu
int  tabuIter;       // iteration counter within one tabu local search call

int  *tabuBeamMode;  // tabuBeamMode[b]: iteration until which beam b's mode is tabu
int **typeProfitSum; // typeProfitSum[b][t]: sum of profit of served tasks of type t on beam b

// CSR buckets of served tasks grouped by beam, rebuilt each local_search
// iteration so the swap neighbourhood can skip unrelated tasks.
int *bucketTask;     // [numTask]   : served task ids, contiguous per beam
int *bucketStart;    // [numBeam+1] : bucketStart[b]..bucketStart[b+1] = beam b's slice
int **insertMinPW;   // insertMinPW[m][bw]: min PW of an unserved task compatible with mode m
int  maxBeamBWCap;

// ---- profiling counters (accumulated over the whole run) ----
long long g_lsIters   = 0;   // total local_search inner iterations
long long g_lsMoves   = 0;   // iterations that actually applied a move
long long g_lsStall   = 0;   // iterations with no admissible candidate
long long g_tabuBlock = 0;   // candidate (task) skips due to tabu (no aspiration)
long long g_aspire    = 0;   // tabu moves admitted via aspiration
double    g_lsTime    = 0.0; // cumulative CPU seconds spent inside local_search
double    g_perturbTime = 0.0; // cumulative CPU seconds spent inside perturb
long long g_applyInsert = 0; // committed insert moves
long long g_applyRemove = 0; // committed remove-only moves
long long g_applySwap   = 0; // committed swap moves
long long g_applyFlip   = 0; // committed mode-flip moves
long long g_applyCross  = 0; // committed cross-beam exchange moves
long long g_crossCand   = 0; // feasible cross-beam exchange candidates
long long g_crossRelocCand = 0; // kind=5 real-dummy candidates
long long g_crossSwapCand  = 0; // kind=5 real-real candidates
long long g_crossFiltered  = 0; // feasible kind=5 moves rejected by open-slot filter
long long g_crossRelocApplied = 0; // committed real-dummy moves
long long g_crossSwapApplied  = 0; // committed real-real moves
long long g_flipCand    = 0; // feasible mode-flip candidates evaluated
long long g_flipApplied = 0; // committed mode-flip moves
long long g_flipApplyToMode[16] = {0}; // committed flips counted by target mode index

// ---- mixed feasible/infeasible search state -----------------------
double mixedPhi = 1000.0;     // penalty weight in eval = profit - phi * over
long long g_mixedIters = 0;   // total mixed-search iterations
long long g_mixedMoves = 0;   // committed mixed-search moves
long long g_mixedFeas  = 0;   // iterations whose current solution is feasible
long long g_mixedInfeas = 0;  // iterations whose current solution is infeasible
long long g_mixedBestUpdate = 0; // feasible best updates found in mixed search
double g_mixedMaxOver = 0.0;  // largest normalized total overflow visited
double g_mixedTime = 0.0;     // cumulative CPU seconds spent inside mixed search

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

    // ---- Header line: B=5 N=50 M=3 Cbw=1578 Cpw=503 ----
    FIC.getline(g_line, MAXLINE);
    int cbw_total, cpw_total;
    sscanf(g_line, "B=%d N=%d M=%d Cbw=%d Cpw=%d",
           &numBeam, &numTask, &numMode, &cbw_total, &cpw_total);

    modeName     = new char *[numMode];
    for (int m = 0; m < numMode; m++)
        modeName[m] = new char[MAXNAMELEN];
    modeBasePower = new int[numMode];

    // blank / "Modes" / "mode base_power"
    FIC.getline(g_line, MAXLINE);
    FIC.getline(g_line, MAXLINE);
    FIC.getline(g_line, MAXLINE);

    // ---- M mode lines: "G 10" ----
    for (int m = 0; m < numMode; m++)
    {
        FIC.getline(g_line, MAXLINE);
        char tmp[MAXNAMELEN];
        sscanf(g_line, "%s %d", tmp, &modeBasePower[m]);
        strcpy(modeName[m], tmp);
    }

    // blank / "Task types" / "type name modes"
    FIC.getline(g_line, MAXLINE);
    FIC.getline(g_line, MAXLINE);
    FIC.getline(g_line, MAXLINE);

    // ---- Task type lines until blank line ----
    static char typeLines[MAXTYPE][256];
    numType = 0;
    while (FIC.getline(g_line, MAXLINE))
    {
        if ((int)strlen(g_line) == 0) break;   // blank line ends this section
        strncpy(typeLines[numType], g_line, 255);
        typeLines[numType][255] = '\0';
        numType++;
    }
    // blank line was already consumed by the breaking getline above

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

    // Beam bandwidth capacity
    FIC.getline(g_line, MAXLINE);                         // header
    for (int b = 0; b < numBeam; b++) FIC >> beamBWCap[b];
    FIC.getline(g_line, MAXLINE);                         // end of data line
    FIC.getline(g_line, MAXLINE);                         // blank line

    // Beam power capacity
    FIC.getline(g_line, MAXLINE);                         // header
    for (int b = 0; b < numBeam; b++) FIC >> beamPWCap[b];
    FIC.getline(g_line, MAXLINE);
    FIC.getline(g_line, MAXLINE);

    taskType     = new int[numTask];
    taskProfit   = new int[numTask];
    taskBWDemand = new int[numTask];
    taskPWDemand = new int[numTask];

    // Task types
    FIC.getline(g_line, MAXLINE);                         // header
    for (int j = 0; j < numTask; j++) FIC >> taskType[j];
    FIC.getline(g_line, MAXLINE);
    FIC.getline(g_line, MAXLINE);

    // Task profits
    FIC.getline(g_line, MAXLINE);                         // header
    for (int j = 0; j < numTask; j++) FIC >> taskProfit[j];
    FIC.getline(g_line, MAXLINE);
    FIC.getline(g_line, MAXLINE);

    // Task bandwidth demands
    FIC.getline(g_line, MAXLINE);                         // header
    for (int j = 0; j < numTask; j++) FIC >> taskBWDemand[j];
    FIC.getline(g_line, MAXLINE);
    FIC.getline(g_line, MAXLINE);

    // Task power demands (last section, no trailing blank needed)
    FIC.getline(g_line, MAXLINE);                         // header
    for (int j = 0; j < numTask; j++) FIC >> taskPWDemand[j];

    FIC.close();

    // ---- Print instance summary ----
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
// task_score: profit per unit of normalized resource usage on beam b
// under mode m.  Bandwidth and power are each normalized by the beam's
// capacity (power capacity = beamPWCap minus the mode base power), so a
// task that eats a large fraction of either resource scores lower.
// Used only as a greedy ordering / priority key, never as the objective.
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
// task_priority: the best score task j can attain on any mode-compatible
// beam (using each beam's total task power capacity, not the dynamic
// remaining power, so the key is static during one greedy pass).
// Returns -1 if j has no compatible beam at all.
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
    // Sequential-commit greedy (same idea as perturb's repair-modes):
    // decide each beam's mode by trial-filling from the current unserved
    // pool, then commit that fill before moving on so beams compete for tasks.
    //================================================================
    for (int b = 0; b < numBeam; b++)
    {
        int bestM      = -1;
        int bestGain   = -1;
        int nBestAdded = 0;

        for (int m = 0; m < numMode; m++)
        {
            if (modeBasePower[m] > beamPWCap[b]) continue;   // power infeasible

            // rank current unserved, mode-compatible tasks by density score
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

        if (bestM < 0)   // no power-feasible mode: take the cheapest base-power
        {
            for (int m = 0; m < numMode; m++)
                if (bestM < 0 || modeBasePower[m] < modeBasePower[bestM]) bestM = m;
            nBestAdded = 0;
        }

        // commit chosen mode + its fill (inline: typeProfitSum not yet allocated)
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
// Print beam-mode distribution + served-by-type for run diagnostics.
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
// feasible_on: can task j be assigned to beam b given the current
// remaining resources?  (mode-compatible + BW + PW all sufficient)
//--------------------------------------------------------------------
int feasible_on(int j, int b)
{
    if (beamMode[b] < 0)                 return 0;   // mode not yet decided
    int t = taskType[j] - 1;
    if (!typeCompatMode[t][beamMode[b]]) return 0;
    if (taskBWDemand[j] > remBW[b])      return 0;
    if (taskPWDemand[j] > remPW[b])      return 0;
    return 1;
}

double beam_over_with_rem_mode(int b, int m, int bwFree, int pwFree)
{
    double over = 0.0;

    if (bwFree < 0)
    {
        int denom = beamBWCap[b];
        if (denom < 1) denom = 1;
        over += (double)(-bwFree) / denom;
    }

    if (pwFree < 0)
    {
        int denom = beamPWCap[b] - modeBasePower[m];
        if (denom < 1) denom = 1;
        over += (double)(-pwFree) / denom;
    }

    return over;
}

double beam_over_with_rem(int b, int bwFree, int pwFree)
{
    return beam_over_with_rem_mode(b, beamMode[b], bwFree, pwFree);
}

double beam_over(int b)
{
    return beam_over_with_rem(b, remBW[b], remPW[b]);
}

double total_over()
{
    double over = 0.0;
    for (int b = 0; b < numBeam; b++)
        if (beamMode[b] >= 0) over += beam_over(b);
    return over;
}

int relaxed_rem_ok_mode(int b, int m, int bwFree, int pwFree)
{
    const double rhoBW = 0.08;
    const double rhoPW = 0.12;

    if (bwFree < 0)
    {
        int denom = beamBWCap[b];
        if (denom < 1) denom = 1;
        if ((double)(-bwFree) / denom > rhoBW) return 0;
    }

    if (pwFree < 0)
    {
        int denom = beamPWCap[b] - modeBasePower[m];
        if (denom < 1) denom = 1;
        if ((double)(-pwFree) / denom > rhoPW) return 0;
    }

    return 1;
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
// add_task / remove_task: assign or unassign a task, keeping remBW /
// remPW / taskBeam / totalProfit consistent.  (No mode change here.)
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

//--------------------------------------------------------------------
// save_best: snapshot current solution into the phase-best arrays.
//--------------------------------------------------------------------
void save_best(double beginTime)
{
    bestTime    = ((double)clock() - beginTime) / CLOCKS_PER_SEC;
    bestProfit  = totalProfit;
    for (int b = 0; b < numBeam; b++) bestBeamMode[b] = beamMode[b];
    for (int j = 0; j < numTask; j++) bestTaskBeam[j] = taskBeam[j];
}

//--------------------------------------------------------------------
// restore_best: copy the phase-best solution back into the current
// solution and recompute remBW / remPW / totalProfit from scratch.
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

//--------------------------------------------------------------------
// apply_mode_flip: switch beam b to mode m2 and evict only tasks whose
// types are incompatible with the new mode.  No refill is performed here.
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
void local_search(double beginTime, int *tmpIdx, double *tmpVal)
{
    (void)tmpIdx;
    (void)tmpVal;

    int ts_depth = 300;              // non-improving iterations before a phase ends (triggers perturb)
    int nonImprove = 0;
    double lsStart = (double)clock();

    for (int j = 0; j < numTask; j++) tabuUntil[j] = 0;
    for (int b = 0; b < numBeam; b++) tabuBeamMode[b] = 0;
    tabuIter = 0;

    while (nonImprove < ts_depth)
    {
        g_lsIters++;
        if (((double)clock() - beginTime) / CLOCKS_PER_SEC > maxRunTime)
            break;

        int bestDelta = MININT_MOVE;     // allow 0 / negative moves (tabu search)
        int numBest   = 0;               // number of candidates tied at bestDelta
        int kind = 0, a = -1, bb = -1, cc = -1;

        // ---- rebuild CSR buckets: served tasks grouped by beam ----------
        // counting sort over taskBeam, O(numTask + numBeam) per iteration.
        for (int b = 0; b <= numBeam; b++) bucketStart[b] = 0;
        for (int j = 0; j < numTask; j++)
            if (taskBeam[j] >= 0) bucketStart[taskBeam[j] + 1]++;
        for (int b = 0; b < numBeam; b++) bucketStart[b + 1] += bucketStart[b];
        // fill: walk a local cursor per beam (reuse the bucket as cursor base)
        for (int j = 0; j < numTask; j++)
            if (taskBeam[j] >= 0)
            {
                int b = taskBeam[j];
                bucketTask[bucketStart[b]++] = j;
            }
        // undo the cursor advance so bucketStart[b] points at beam b's start again
        for (int b = numBeam; b > 0; b--) bucketStart[b] = bucketStart[b - 1];
        bucketStart[0] = 0;

        // (1) insert unserved task ----------------------------------
        for (int j = 0; j < numTask; j++)
        {
            if (taskBeam[j] != -1) continue;
            int delta = taskProfit[j];
            // Aspiration: a tabu insert is admitted only if it would push the
            // current solution past the phase-best profit.
            int jTabu = (tabuIter < tabuUntil[j]);
            if (jTabu)
            {
                if (totalProfit + delta > bestProfit) { /* aspiration: admit */ }
                else { g_tabuBlock++; continue; }
            }
            for (int b = 0; b < numBeam; b++)
            {
                if (!feasible_on(j, b)) continue;
                if (delta < bestDelta) continue;
                record_candidate(1, j, b, -1, delta,
                                 bestDelta, numBest, kind, a, bb, cc);
            }
        }

        // (2) remove-only: drop a served task ------------------------
        for (int j = 0; j < numTask; j++)
        {
            if (taskBeam[j] < 0) continue;
            int delta = -taskProfit[j];
            if (tabuIter < tabuUntil[j])
            {
                if (totalProfit + delta > bestProfit) { /* aspiration: admit */ }
                else { g_tabuBlock++; continue; }
            }
            if (delta < bestDelta) continue;

            record_candidate(2, j, -1, -1, delta,
                             bestDelta, numBest, kind, a, bb, cc);
        }

        // (3) swap: unserved i replaces served k on the same beam -------
        //     For each unserved i, only visit beams whose mode is compatible
        //     with i, and within such a beam only its served tasks (bucket),
        //     instead of rescanning every task.  Same (i,k) pairs as before.
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

                    // Aspiration: a swap whose i or k is tabu is admitted only
                    // if it would beat the phase-best profit.
                    if (iTabu || tabuIter < tabuUntil[k])
                    {
                        if (totalProfit + delta > bestProfit) { /* aspiration: admit */ }
                        else continue;
                    }

                    // feasibility after removing k, adding i (same beam b)
                    int bwFree = bwFree0 + taskBWDemand[k];
                    int pwFree = pwFree0 + taskPWDemand[k];
                    if (taskBWDemand[i] > bwFree || taskPWDemand[i] > pwFree) continue;

                    record_candidate(3, i, b, k, delta,
                                     bestDelta, numBest, kind, a, bb, cc);
                }
            }
        }

        // (4) mode-flip: switch one beam's mode, evicting only tasks
        //     incompatible with the new mode.  No refill is evaluated here.
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
                    if (totalProfit + delta > bestProfit) { /* aspiration: admit */ }
                    else { g_tabuBlock++; continue; }
                }

                g_flipCand++;
                if (delta < bestDelta) continue;
                record_candidate(4, b, m2, -1, delta,
                                 bestDelta, numBest, kind, a, bb, cc);
            }
        }

        // (5) cross-beam exchange with a dynamic dummy slot --------------
        //     real-real:  served i on b1 swaps beam with served k on b2.
        //     real-dummy: served i on b1 relocates to b2; cc == -1 records
        //     the temporary dummy(b2).  Dummy is not stored in taskBeam[].
        //     Since every kind=5 move has delta=0, skip this neighbourhood
        //     when a strictly improving move is already available.
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
                    if (totalProfit > bestProfit) { /* aspiration: admit */ }
                    else { g_tabuBlock++; continue; }
                }

                int ti = taskType[i] - 1;
                for (int b2 = 0; b2 < numBeam; b2++)
                {
                    if (b2 == b1) continue;
                    if (beamMode[b2] < 0) continue;
                    if (!typeCompatMode[ti][beamMode[b2]]) continue;

                    int delta = 0;

                    // real-dummy: move i from b1 to the current free slot on b2.
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

                    // real-real: enumerate each cross-beam pair once.
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
                            if (totalProfit > bestProfit) { /* aspiration: admit */ }
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

        // ---- no admissible candidate in this iteration: advance tabu time ----
        if (numBest == 0)
        {
            g_lsStall++;
            nonImprove++;
            tabuIter++;
            continue;
        }
        g_lsMoves++;
        // count an aspiration whenever the committed move broke a task/beam tabu
        if ((kind == 1 && tabuIter < tabuUntil[a]) ||
            (kind == 2 && tabuIter < tabuUntil[a]) ||
            (kind == 3 && (tabuIter < tabuUntil[a] || tabuIter < tabuUntil[cc])) ||
            (kind == 4 && tabuIter < tabuBeamMode[a]) ||
            (kind == 5 && (tabuIter < tabuUntil[a] ||
                           (cc >= 0 && tabuIter < tabuUntil[cc]))))
            g_aspire++;

        if (kind == 1)                            // insert
        {
            add_task(a, bb);
            moveFreq[a]++;
            tabuUntil[a] = tabuIter + tabu_tenure();
            g_applyInsert++;
        }
        else if (kind == 2)                       // remove-only
        {
            remove_task(a);
            moveFreq[a]++;
            tabuUntil[a] = tabuIter + tabu_tenure();
            g_applyRemove++;
        }
        else if (kind == 3)                       // swap a in, cc out
        {
            remove_task(cc);
            add_task(a, bb);
            moveFreq[a]++; moveFreq[cc]++;
            tabuUntil[a]  = tabuIter + tabu_tenure();
            tabuUntil[cc] = tabuIter + tabu_tenure();
            g_applySwap++;
        }
        else if (kind == 4)                       // mode-flip beam a to mode bb
        {
            int tenure = tabu_tenure();
            apply_mode_flip(a, bb, tenure);
            tabuBeamMode[a] = tabuIter + tenure;
            g_applyFlip++;
            g_flipApplied++;
            if (bb >= 0 && bb < 16) g_flipApplyToMode[bb]++;
        }
        else if (kind == 5)                       // cross-beam exchange with dummy
        {
            int b1 = taskBeam[a];
            if (cc < 0)                           // real-dummy: relocate a to bb
            {
                remove_task(a);
                add_task(a, bb);
                moveFreq[a]++;
                tabuUntil[a] = tabuIter + tabu_tenure();
                g_crossRelocApplied++;
            }
            else                                  // real-real: swap a and cc
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
// local_search_mixed: feasible + infeasible tabu search.
// BW/PW capacity constraints are relaxed within a small boundary.  Mode
// compatibility and one-beam-per-task remain hard constraints.  The current
// solution may be infeasible, but only feasible solutions update bestProfit.
//--------------------------------------------------------------------
void local_search_mixed(double beginTime, int *tmpIdx, double *tmpVal)
{
    (void)tmpIdx;
    (void)tmpVal;

    const double EPS = 1e-12;
    const int mixedDepth = 300;
    const int phiWindow = 5;
    const double phiTau = 2.0;
    const double phiMin = 100.0;
    const double phiMax = 100000.0;

    int nonImprove = 0;
    int feasibleStreak = 0;
    int infeasibleStreak = 0;
    double mixedStart = (double)clock();
    double curOver = total_over();

    for (int j = 0; j < numTask; j++) tabuUntil[j] = 0;
    for (int b = 0; b < numBeam; b++) tabuBeamMode[b] = 0;
    tabuIter = 0;

    while (nonImprove < mixedDepth)
    {
        g_mixedIters++;
        if (((double)clock() - beginTime) / CLOCKS_PER_SEC > maxRunTime)
            break;

        if (curOver <= EPS)
        {
            g_mixedFeas++;
            feasibleStreak++;
            infeasibleStreak = 0;
            if (feasibleStreak >= phiWindow)
            {
                mixedPhi /= phiTau;
                if (mixedPhi < phiMin) mixedPhi = phiMin;
                feasibleStreak = 0;
            }
        }
        else
        {
            g_mixedInfeas++;
            infeasibleStreak++;
            feasibleStreak = 0;
            if (infeasibleStreak >= phiWindow)
            {
                mixedPhi *= phiTau;
                if (mixedPhi > phiMax) mixedPhi = phiMax;
                infeasibleStreak = 0;
            }
        }

        double bestDeltaEval = -1.0e100;
        int bestDeltaProfit = MININT_MOVE;
        int numBest = 0;
        int kind = 0, a = -1, bb = -1, cc = -1;

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

        // (1) insert unserved task ----------------------------------
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

                double newOver = curOver - beam_over(b) + beam_over_with_rem(b, bw2, pw2);
                if (tabuIter < tabuUntil[j])
                {
                    if (newOver <= EPS && totalProfit + deltaProfit > bestProfit) { /* aspiration */ }
                    else { g_tabuBlock++; continue; }
                }

                double deltaEval = (double)deltaProfit - mixedPhi * (newOver - curOver);
                record_mixed_candidate(1, j, b, -1, deltaProfit, deltaEval,
                                       bestDeltaEval, bestDeltaProfit, numBest,
                                       kind, a, bb, cc);
            }
        }

        // (2) remove-only -------------------------------------------
        for (int j = 0; j < numTask; j++)
        {
            int b = taskBeam[j];
            if (b < 0) continue;

            int deltaProfit = -taskProfit[j];
            int bw2 = remBW[b] + taskBWDemand[j];
            int pw2 = remPW[b] + taskPWDemand[j];
            double newOver = curOver - beam_over(b) + beam_over_with_rem(b, bw2, pw2);

            if (tabuIter < tabuUntil[j])
            {
                if (newOver <= EPS && totalProfit + deltaProfit > bestProfit) { /* aspiration */ }
                else { g_tabuBlock++; continue; }
            }

            double deltaEval = (double)deltaProfit - mixedPhi * (newOver - curOver);
            record_mixed_candidate(2, j, -1, -1, deltaProfit, deltaEval,
                                   bestDeltaEval, bestDeltaProfit, numBest,
                                   kind, a, bb, cc);
        }

        // (3) same-beam swap ----------------------------------------
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

                    double newOver = curOver - beam_over(b) + beam_over_with_rem(b, bw2, pw2);
                    if (tabuIter < tabuUntil[i] || tabuIter < tabuUntil[k])
                    {
                        if (newOver <= EPS && totalProfit + deltaProfit > bestProfit) { /* aspiration */ }
                        else { g_tabuBlock++; continue; }
                    }

                    double deltaEval = (double)deltaProfit - mixedPhi * (newOver - curOver);
                    record_mixed_candidate(3, i, b, k, deltaProfit, deltaEval,
                                           bestDeltaEval, bestDeltaProfit, numBest,
                                           kind, a, bb, cc);
                }
            }
        }

        // (4) mode-flip ---------------------------------------------
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
                double newOver = curOver - beam_over(b) + beam_over_with_rem_mode(b, m2, bw2, pw2);
                if (tabuIter < tabuBeamMode[b] || evictTabu)
                {
                    if (newOver <= EPS && totalProfit + deltaProfit > bestProfit) { /* aspiration */ }
                    else { g_tabuBlock++; continue; }
                }

                g_flipCand++;
                double deltaEval = (double)deltaProfit - mixedPhi * (newOver - curOver);
                record_mixed_candidate(4, b, m2, -1, deltaProfit, deltaEval,
                                       bestDeltaEval, bestDeltaProfit, numBest,
                                       kind, a, bb, cc);
            }
        }

        // (5) cross-beam exchange with a dynamic dummy slot ----------
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
                    double newOver = curOver
                        - beam_over(b1) - beam_over(b2)
                        + beam_over_with_rem(b1, b1BW, b1PW)
                        + beam_over_with_rem(b2, b2BW, b2PW);

                    if (tabuIter < tabuUntil[i])
                    {
                        if (newOver <= EPS && totalProfit > bestProfit) { /* aspiration */ }
                        else { g_tabuBlock++; continue; }
                    }

                    double deltaEval = -mixedPhi * (newOver - curOver);
                    g_crossCand++;
                    g_crossRelocCand++;
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

                    double newOver = curOver
                        - beam_over(b1) - beam_over(b2)
                        + beam_over_with_rem(b1, b1BW, b1PW)
                        + beam_over_with_rem(b2, b2BW, b2PW);

                    if (tabuIter < tabuUntil[i] || tabuIter < tabuUntil[k])
                    {
                        if (newOver <= EPS && totalProfit > bestProfit) { /* aspiration */ }
                        else { g_tabuBlock++; continue; }
                    }

                    double deltaEval = -mixedPhi * (newOver - curOver);
                    g_crossCand++;
                    g_crossSwapCand++;
                    record_mixed_candidate(5, i, b2, k, 0, deltaEval,
                                           bestDeltaEval, bestDeltaProfit, numBest,
                                           kind, a, bb, cc);
                }
            }
        }

        if (numBest == 0)
        {
            nonImprove++;
            tabuIter++;
            continue;
        }

        g_mixedMoves++;
        if ((kind == 1 && tabuIter < tabuUntil[a]) ||
            (kind == 2 && tabuIter < tabuUntil[a]) ||
            (kind == 3 && (tabuIter < tabuUntil[a] || tabuIter < tabuUntil[cc])) ||
            (kind == 4 && tabuIter < tabuBeamMode[a]) ||
            (kind == 5 && (tabuIter < tabuUntil[a] ||
                           (cc >= 0 && tabuIter < tabuUntil[cc]))))
            g_aspire++;

        if (kind == 1)
        {
            add_task(a, bb);
            moveFreq[a]++;
            tabuUntil[a] = tabuIter + tabu_tenure();
            g_applyInsert++;
        }
        else if (kind == 2)
        {
            remove_task(a);
            moveFreq[a]++;
            tabuUntil[a] = tabuIter + tabu_tenure();
            g_applyRemove++;
        }
        else if (kind == 3)
        {
            remove_task(cc);
            add_task(a, bb);
            moveFreq[a]++; moveFreq[cc]++;
            tabuUntil[a]  = tabuIter + tabu_tenure();
            tabuUntil[cc] = tabuIter + tabu_tenure();
            g_applySwap++;
        }
        else if (kind == 4)
        {
            int tenure = tabu_tenure();
            apply_mode_flip(a, bb, tenure);
            tabuBeamMode[a] = tabuIter + tenure;
            g_applyFlip++;
            g_flipApplied++;
            if (bb >= 0 && bb < 16) g_flipApplyToMode[bb]++;
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
                g_crossRelocApplied++;
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
                g_crossSwapApplied++;
            }
            g_applyCross++;
        }

        curOver = total_over();
        if (curOver > g_mixedMaxOver) g_mixedMaxOver = curOver;

        if (curOver <= EPS && totalProfit > bestProfit)
        {
            save_best(beginTime);
            g_mixedBestUpdate++;
            nonImprove = 0;
        }
        else
            nonImprove++;

        tabuIter++;
    }

    restore_best();       // mixed phase always hands a feasible solution forward
    g_mixedTime += ((double)clock() - mixedStart) / CLOCKS_PER_SEC;
}

//--------------------------------------------------------------------
// perturb  (destroy + repair)
//
// Starting from the current phase-best solution:
//
//  Destroy: pick 35% of the beams, strip every
//  task off them (back to the unserved pool) and reset their mode to -1
//  ("empty"), leaving a partial solution.  Untouched beams keep both
//  their mode and their assigned tasks.
//
//  Repair (modes): redecide emptied beams sequentially.  For one beam,
//  enumerate its power-feasible modes, temporarily fill that beam from
//  the current unserved pool, and usually keep the mode whose committed
//  repair gives the largest actual profit increase.  With small
//  probability the second-best mode is used to avoid deterministic
//  repair undoing the destroy step.  The chosen mode/tasks are committed
//  before repairing the next emptied beam, so tasks cannot be counted by
//  multiple beams during mode choice.
//
//  Repair (tasks): greedily reassign the still-unserved tasks, ordered
//  by the existing best-attainable priority key.  For each task, try
//  compatible beam modes from low base-power to high base-power first;
//  within the chosen mode tier, keep the original tightest-BW fit rule.
//--------------------------------------------------------------------
void perturb(int *tmpIdx, double *tmpVal)
{
    restore_best();
    if (numBeam == 0) return;

    //---- destroy: empty 30% of the beams, preferring those whose served
    //     tasks were moved least during the last local search (low
    //     moveFreq) so the search is pushed toward unexplored regions ----
    int numDestroy = (int)(0.30 * numBeam + 0.5);
    if (numDestroy < 1)        numDestroy = 1;
    if (numDestroy > numBeam)  numDestroy = numBeam;

    // per-beam frequency key = sum of moveFreq over its served tasks
    int    *perturbBeam = new int[numBeam];      // scratch: beams chosen to empty
    int    *prevMode    = new int[numBeam];      // mode each emptied beam had before destroy
    double *beamFreqKey = new double[numBeam];
    int    *order       = new int[numBeam];
    for (int b = 0; b < numBeam; b++) { beamFreqKey[b] = 0.0; order[b] = b; }
    for (int j = 0; j < numTask; j++)
        if (taskBeam[j] >= 0) beamFreqKey[taskBeam[j]] += moveFreq[j];

    // sort beams ascending by frequency: qsort_desc on the negated key
    for (int b = 0; b < numBeam; b++) beamFreqKey[b] = -beamFreqKey[b];
    qsort_desc(beamFreqKey, order, 0, numBeam - 1);   // order[]: low-freq first

    // pool of least-moved beams to draw from (at least numDestroy wide)
    int poolSize = 0;
    for (int b = 0; b < numBeam; b++) if (beamFreqKey[b] == 0.0) poolSize++;
    if (poolSize < numDestroy) poolSize = numDestroy;

    // draw numDestroy distinct beams at random from the low-frequency pool
    for (int k = 0; k < numDestroy; k++)             // partial Fisher-Yates
    {
        int r   = k + rand() % (poolSize - k);
        int tmp = order[k]; order[k] = order[r]; order[r] = tmp;
        perturbBeam[k] = order[k];
    }

    for (int k = 0; k < numDestroy; k++)
    {
        int b = perturbBeam[k];
        prevMode[k] = beamMode[b];                  // remember old mode; repair may choose it again
        for (int j = 0; j < numTask; j++)
            if (taskBeam[j] == b) remove_task(j);
        beamMode[b] = -1;
        remBW[b]    = beamBWCap[b];
        remPW[b]    = 0;                             // set once mode is redecided
    }

    //---- repair modes: sequential actual-profit trial ----------------
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
            if (modeBasePower[m] > beamPWCap[b]) continue;       // power infeasible

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

        if (bestM < 0)                       // no power-feasible mode: take the cheapest
        {                                    // base-power mode available
            bestM = -1;
            for (int m = 0; m < numMode; m++)
            {
                if (bestM < 0 || modeBasePower[m] < modeBasePower[bestM]) bestM = m;
            }
            if (bestM < 0) bestM = prevMode[k];   // only one mode exists: unavoidable
        }

        int chosenM      = bestM;
        int *chosenAdded = bestAdded;
        int nChosenAdded = nBestAdded;
        if (secondM >= 0 && rand() % 5 == 0)  // 20%: diversify repair
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

    //---- repair tasks: greedy best-fit by score ----------------------
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

void ils()
{
    double beginTime = (double)clock();

    int    *tmpIdx = new int   [numTask];
    double *tmpVal = new double[numTask];

    greedy_init();
    alloc_search();

    // full rebuild of typeProfitSum to guarantee consistency before ILS
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

        save_best(beginTime);          // phase-best starts from the current solution
        local_search(beginTime, tmpIdx, tmpVal);
        restore_best();                // mixed search starts from the best feasible point
        local_search_mixed(beginTime, tmpIdx, tmpVal);
        restore_best();                // perturb and global update receive a feasible solution
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

    // install the global best as the final solution for checking / printing
    restore_global();
    totalProfit = globalProfit;

    probe_modes("final");

    delete[] tmpIdx;
    delete[] tmpVal;

    cout << "ILS done.  bestProfit=" << globalProfit
         << "  bestTime=" << globalBestTime << " s"
         << "  phases=" << numPhase << endl;
    cout << "Perturbations over the whole run: " << numPhase << endl;

    // ---- profiling report ----
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
         << "  final_phi=" << mixedPhi << endl;
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
    cout << "[PROFILE] flip_cand=" << g_flipCand
         << "  applied=" << g_flipApplied
         << " (" << (g_lsMoves ? 100.0 * g_flipApplied / g_lsMoves : 0) << "% of moves)" << endl;
    cout << "[PROFILE] flip_applied_by_mode:";
    for (int m = 0; m < numMode && m < 16; m++)
        cout << " " << modeName[m] << "=" << g_flipApplyToMode[m];
    cout << endl;

}

//--------------------------------------------------------------------
// check_solution: verify all constraints and profit consistency
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

        // mode-compatibility
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
// print_solution: human-readable summary
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
// free_memory
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
    for (int m = 0; m < numMode; m++) delete[] insertMinPW[m];
    delete[] insertMinPW;

    delete[] tabuBeamMode;
    for (int b = 0; b < numBeam; b++) delete[] typeProfitSum[b];
    delete[] typeProfitSum;
}

//--------------------------------------------------------------------
// main
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

    maxRunTime = 600.0;          // time limit in seconds
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
