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

// CSR buckets of served tasks grouped by beam, rebuilt each local_search
// iteration so the swap neighbourhood can skip unrelated tasks.
int *bucketTask;     // [numTask]   : served task ids, contiguous per beam
int *bucketStart;    // [numBeam+1] : bucketStart[b]..bucketStart[b+1] = beam b's slice


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

    //================================================================
    // Phase 1: determine mode for each beam
    //================================================================
    for (int b = 0; b < numBeam; b++)
    {
        int bestMode      = 0;
        int bestEstProfit = -1;

        for (int m = 0; m < numMode; m++)
        {
            int avail_pw = beamPWCap[b] - modeBasePower[m];
            if (avail_pw < 0) continue;   // cannot even power on this mode

            // collect tasks compatible with mode m, ranked by resource-density score
            int numComp = 0;
            for (int j = 0; j < numTask; j++)
            {
                int t = taskType[j] - 1;
                if (typeCompatMode[t][m])
                {
                    tmpIdx[numComp] = j;
                    tmpVal[numComp] = task_score(j, b, m);
                    numComp++;
                }
            }
            if (numComp > 0)
                qsort_desc(tmpVal, tmpIdx, 0, numComp - 1);

            // greedy fill to estimate profit
            int rem_bw = beamBWCap[b];
            int rem_pw = avail_pw;
            int est    = 0;
            for (int k = 0; k < numComp; k++)
            {
                int j = tmpIdx[k];
                if (taskBWDemand[j] <= rem_bw && taskPWDemand[j] <= rem_pw)
                {
                    rem_bw -= taskBWDemand[j];
                    rem_pw -= taskPWDemand[j];
                    est    += taskProfit[j];
                }
            }

            if (est > bestEstProfit)
            {
                bestEstProfit = est;
                bestMode      = m;
            }
        }

        beamMode[b] = bestMode;
        remBW[b]    = beamBWCap[b];
        remPW[b]    = beamPWCap[b] - modeBasePower[bestMode];
    }

    //================================================================
    // Phase 2: task assignment
    //================================================================

    // Sort tasks by best attainable score (descending)
    for (int j = 0; j < numTask; j++)
    {
        tmpIdx[j] = j;
        tmpVal[j] = task_priority(j);
    }
    qsort_desc(tmpVal, tmpIdx, 0, numTask - 1);

    totalProfit = 0;

    for (int ki = 0; ki < numTask; ki++)
    {
        int j = tmpIdx[ki];
        int t = taskType[j] - 1;

        int bestBeam    = -1;
        int bestLeftBW  = -1;   // tightest BW fit: smallest leftover BW after add

        for (int b = 0; b < numBeam; b++)
        {
            if (!typeCompatMode[t][beamMode[b]])       continue;
            if (taskBWDemand[j] > remBW[b])            continue;
            if (taskPWDemand[j] > remPW[b])            continue;

            int leftover = remBW[b] - taskBWDemand[j];
            if (bestBeam < 0 || leftover < bestLeftBW)
            {
                bestBeam   = b;
                bestLeftBW = leftover;
            }
        }

        if (bestBeam >= 0)
        {
            taskBeam[j]      = bestBeam;
            remBW[bestBeam] -= taskBWDemand[j];
            remPW[bestBeam] -= taskPWDemand[j];
            totalProfit     += taskProfit[j];
        }
    }

    delete[] tmpIdx;
    delete[] tmpVal;

    cout << "Greedy init done.  totalProfit=" << totalProfit << endl;
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
}

void remove_task(int j)
{
    int b = taskBeam[j];
    remBW[b]    += taskBWDemand[j];
    remPW[b]    += taskPWDemand[j];
    totalProfit -= taskProfit[j];
    taskBeam[j]  = -1;
}

int tabu_tenure()
{
    return 5 + rand() % 6 + (int)(0.1 * numTask);
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
        }
    }
}

//--------------------------------------------------------------------
// alloc_search: allocate ILS bookkeeping arrays.
//--------------------------------------------------------------------
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

void local_search(double beginTime, int *tmpIdx, double *tmpVal)
{
    (void)tmpIdx;
    (void)tmpVal;

    while (1)
    {
        if (((double)clock() - beginTime) / CLOCKS_PER_SEC > maxRunTime)
            break;

        int bestDelta = 0;               // pure descent accepts improving moves only
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
            for (int b = 0; b < numBeam; b++)
            {
                if (!feasible_on(j, b)) continue;
                if (delta <= 0 || delta < bestDelta) continue;
                record_candidate(1, j, b, -1, delta,
                                 bestDelta, numBest, kind, a, bb, cc);
            }
        }

        // (2) swap: unserved i replaces served k on the same beam -------
        //     Same (i,k) pairs as the task-by-task scan, but grouped by beam.
        for (int i = 0; i < numTask; i++)
        {
            if (taskBeam[i] != -1) continue;
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
                    if (delta <= 0 || delta < bestDelta) continue;

                    // feasibility after removing k, adding i (same beam b)
                    int bwFree = bwFree0 + taskBWDemand[k];
                    int pwFree = pwFree0 + taskPWDemand[k];
                    if (taskBWDemand[i] > bwFree || taskPWDemand[i] > pwFree) continue;

                    record_candidate(2, i, b, k, delta,
                                     bestDelta, numBest, kind, a, bb, cc);
                }
            }
        }

        // ---- no improving candidate -> local optimum ----------------
        if (numBest == 0)
            break;

        if (kind == 1)                            // insert
        {
            add_task(a, bb);
            moveFreq[a]++;
        }
        else if (kind == 2)                       // swap a in, cc out
        {
            remove_task(cc);
            add_task(a, bb);
            moveFreq[a]++; moveFreq[cc]++;
        }

        if (totalProfit > bestProfit)
            save_best(beginTime);
    }
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

    //---- destroy: empty 35% of the beams, preferring those whose served
    //     tasks were moved least during the last local search (low
    //     moveFreq) so the search is pushed toward unexplored regions ----
    int numDestroy = (int)(0.35 * numBeam + 0.5);
    if (numDestroy < 1)        numDestroy = 1;
    if (numDestroy > numBeam)  numDestroy = numBeam;

    // per-beam frequency key = sum of moveFreq over its served tasks
    int    perturbBeam[numBeam];                 // scratch: beams chosen to empty
    double beamFreqKey[numBeam];
    int    order[numBeam];
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
        {
            bestM = 0;
            for (int m = 1; m < numMode; m++)
                if (modeBasePower[m] < modeBasePower[bestM]) bestM = m;
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
}

//--------------------------------------------------------------------
// iterated_descent_search
//
// Iterated Local Search driver:
//   greedy_init  -> pure descent local_search -> perturb -> local_search -> ...
// keeping the global best, until the time limit is reached.
//--------------------------------------------------------------------
void iterated_descent_search()
{
    double beginTime = (double)clock();

    int    *tmpIdx = new int   [numTask];
    double *tmpVal = new double[numTask];

    greedy_init();
    alloc_search();

    globalProfit = -1;

    double runTime = 0.0;
    long   numPhase = 0;
    while (runTime < maxRunTime)
    {
        for (int j = 0; j < numTask; j++) moveFreq[j] = 0;

        save_best(beginTime);          // phase-best starts from the current solution
        local_search(beginTime, tmpIdx, tmpVal);
        numPhase++;

        if (bestProfit > globalProfit)
        {
            save_global_from_best();
            cout << "  globalProfit=" << globalProfit
                 << "  time=" << globalBestTime << " s"
                 << "  phase=" << numPhase << endl;
        }

        perturb(tmpIdx, tmpVal);

        runTime = ((double)clock() - beginTime) / CLOCKS_PER_SEC;
    }

    // install the global best as the final solution for checking / printing
    restore_global();
    totalProfit = globalProfit;

    delete[] tmpIdx;
    delete[] tmpVal;

    cout << "Iterated descent done.  bestProfit=" << globalProfit
         << "  bestTime=" << globalBestTime << " s"
         << "  phases=" << numPhase << endl;
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
    iterated_descent_search();
    check_solution();
    print_solution();

    double elapsed = ((double)clock() - t0) / CLOCKS_PER_SEC;
    cout << "\nElapsed time: " << elapsed << " s" << endl;

    free_memory();
    return 0;
}
