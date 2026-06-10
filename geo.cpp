#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unordered_set>

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
int *tabuUntil;      // tabuUntil[j] : iteration index until which task j is frozen
int *tabuBeamMode;   // tabuBeamMode[b]: iteration until which beam b's mode is frozen
int  tabuIter;       // iteration counter within one local search call

// ---- solution-based tabu (SBTS-style) state ----------------------------
// Instead of forbidding move *attributes* (which over-constrains and stalls on
// plateaus), we fingerprint the whole solution and forbid re-visiting it.
// Each (variable, value) pair gets a random 64-bit weight; the fingerprint H is
// the sum of the weights of every variable's current value, wrapping naturally
// in unsigned 64-bit arithmetic (no modulo needed).  A single weight set + an
// unordered_set membership test replaces MSBTS's three hash bitmaps: there are
// no bitmap bit-collision false positives; 64-bit fingerprint collision remains
// possible in theory but is negligible in these runs.
unsigned long long **taskW; // taskW[j][b]: weight of taskBeam[j]==b; b==numBeam means "unserved"
unsigned long long **beamV; // beamV[b][m]: weight of beamMode[b]==m; m==numMode means "empty(-1)"
unsigned long long  curHash;                 // running fingerprint of the current solution
unordered_set<unsigned long long> visited;   // fingerprints seen in the current local_search phase

// CSR buckets of served tasks grouped by beam, rebuilt each local_search
// iteration so the swap neighbourhood can skip unrelated tasks.
int *bucketTask;     // [numTask]   : served task ids, contiguous per beam
int *bucketStart;    // [numBeam+1] : bucketStart[b]..bucketStart[b+1] = beam b's slice

static char g_line[MAXLINE];

// ---- independent RNG for the hash-weight tables only -------------------
// Kept separate from srand()/rand() so that filling the weight tables never
// consumes from the rand() stream that drives move selection (tie-breaking,
// tenure, perturbation).  This keeps an A/B comparison fair: the baseline and
// this SBTS variant walk the *same* rand() sequence under the same seed.
static unsigned long long g_wstate;          // weight-RNG state

void wrng_seed(unsigned long long s)
{
    g_wstate = s ? s : 0x9E3779B97F4A7C15ULL; // avoid the all-zero fixed point
}

unsigned long long wrng_next()
{
    unsigned long long x = g_wstate;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    g_wstate = x;
    return x;
}

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


void alloc_solution()
{
    beamMode = new int[numBeam];
    taskBeam = new int[numTask];
    remBW    = new int[numBeam];
    remPW    = new int[numBeam];
}

void reset_solution()
{
    totalProfit = 0;
    for (int b = 0; b < numBeam; b++)
    {
        beamMode[b] = -1;
        remBW[b]    = beamBWCap[b];
        remPW[b]    = 0;
    }
    for (int j = 0; j < numTask; j++) taskBeam[j] = -1;
}

void greedy_init()
{
    reset_solution();

    int    *tmpIdx = new int   [numTask];
    double *tmpVal = new double[numTask];
    int    *trialAdded = new int[numTask];
    int    *bestAdded  = new int[numTask];

    totalProfit = 0;

    //================================================================
    // Sequential-commit greedy (same idea as perturb's repair-modes):
    // decide each beam's mode by trial-filling from the *current* unserved
    // pool, then commit that fill before moving on.  Because the pool shrinks
    // as beams commit, beams compete for scarce high-density tasks instead of
    // each independently assuming it can grab them all — which is what made
    // the old per-beam-independent estimate collapse every beam onto the
    // low-base-power mode.
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

        // commit chosen mode + its fill
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

void greedy_randomized_init()
{
    reset_solution();

    int    *tmpIdx = new int   [numTask];
    double *tmpVal = new double[numTask];
    int    *trialAdded = new int[numTask];
    int    *bestAdded  = new int[numTask];
    int    *secondAdded = new int[numTask];
    int    *beamOrder = new int[numBeam];

    for (int b = 0; b < numBeam; b++) beamOrder[b] = b;
    for (int b = numBeam - 1; b > 0; b--)
    {
        int r = rand() % (b + 1);
        int tmp = beamOrder[b]; beamOrder[b] = beamOrder[r]; beamOrder[r] = tmp;
    }

    for (int ob = 0; ob < numBeam; ob++)
    {
        int b = beamOrder[ob];
        int bestM        = -1;
        int bestGain     = -1;
        int nBestAdded   = 0;
        int secondM      = -1;
        int secondGain   = -1;
        int nSecondAdded = 0;

        for (int m = 0; m < numMode; m++)
        {
            if (modeBasePower[m] > beamPWCap[b]) continue;

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
                secondGain   = bestGain;
                secondM      = bestM;
                nSecondAdded = nBestAdded;
                for (int a = 0; a < nBestAdded; a++) secondAdded[a] = bestAdded[a];

                bestGain   = gain;
                bestM      = m;
                nBestAdded = nTrial;
                for (int a = 0; a < nTrial; a++) bestAdded[a] = trialAdded[a];
            }
            else if (gain > secondGain)
            {
                secondGain   = gain;
                secondM      = m;
                nSecondAdded = nTrial;
                for (int a = 0; a < nTrial; a++) secondAdded[a] = trialAdded[a];
            }
        }

        if (bestM < 0)
        {
            for (int m = 0; m < numMode; m++)
                if (bestM < 0 || modeBasePower[m] < modeBasePower[bestM]) bestM = m;
            nBestAdded = 0;
        }

        int chosenM      = bestM;
        int *chosenAdded = bestAdded;
        int nChosenAdded = nBestAdded;
        if (secondM >= 0 && rand() % 5 == 0)
        {
            chosenM      = secondM;
            chosenAdded  = secondAdded;
            nChosenAdded = nSecondAdded;
        }

        beamMode[b] = chosenM;
        remBW[b]    = beamBWCap[b];
        remPW[b]    = beamPWCap[b] - modeBasePower[chosenM];
        for (int a = 0; a < nChosenAdded; a++)
        {
            int j = chosenAdded[a];
            if (taskBeam[j] == -1 && taskBWDemand[j] <= remBW[b] && taskPWDemand[j] <= remPW[b])
            {
                taskBeam[j]  = b;
                remBW[b]    -= taskBWDemand[j];
                remPW[b]    -= taskPWDemand[j];
                totalProfit += taskProfit[j];
            }
        }
    }

    delete[] tmpIdx;
    delete[] tmpVal;
    delete[] trialAdded;
    delete[] bestAdded;
    delete[] secondAdded;
    delete[] beamOrder;
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


void alloc_search()
{
    bestBeamMode = new int[numBeam];
    bestTaskBeam = new int[numTask];
    globalBeamMode = new int[numBeam];
    globalTaskBeam = new int[numTask];
    moveFreq     = new int[numTask];
    tabuUntil    = new int[numTask];
    tabuBeamMode = new int[numBeam];
    bucketTask   = new int[numTask];
    bucketStart  = new int[numBeam + 1];

    // ---- hash-weight tables for solution-based tabu --------------------
    // own seed derived from `seed` (so different seeds -> different weights,
    // same seed -> reproducible) but drawn from the independent weight-RNG.
    wrng_seed((unsigned long long)seed * 2654435761ULL + 12345ULL);

    taskW = new unsigned long long *[numTask];
    for (int j = 0; j < numTask; j++)
    {
        taskW[j] = new unsigned long long[numBeam + 1];   // [numBeam] = "unserved"
        for (int b = 0; b <= numBeam; b++) taskW[j][b] = wrng_next();
    }

    beamV = new unsigned long long *[numBeam];
    for (int b = 0; b < numBeam; b++)
    {
        beamV[b] = new unsigned long long[numMode + 1];   // [numMode] = "empty(-1)"
        for (int m = 0; m <= numMode; m++) beamV[b][m] = wrng_next();
    }
}

//--------------------------------------------------------------------
// compute_hash: full fingerprint of the current solution from scratch.
// O(numTask + numBeam); called once on entering local_search and after
// perturb (which does not maintain the running hash incrementally).
//--------------------------------------------------------------------
unsigned long long compute_hash()
{
    unsigned long long h = 0;
    for (int j = 0; j < numTask; j++)
        h += taskW[j][ taskBeam[j] == -1 ? numBeam : taskBeam[j] ];
    for (int b = 0; b < numBeam; b++)
        h += beamV[b][ beamMode[b] == -1 ? numMode : beamMode[b] ];
    return h;
}

void record_candidate(int kind, int a, int b, int c, int delta,
                      unsigned long long cand,
                      int &bestDelta, int &numBest,
                      int &chosenKind, int &chosenA, int &chosenB, int &chosenC,
                      unsigned long long &chosenHash)
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
        chosenHash = cand;
        return;
    }

    numBest++;
    if (rand() % numBest == 0)
    {
        chosenKind = kind;
        chosenA    = a;
        chosenB    = b;
        chosenC    = c;
        chosenHash = cand;
    }
}

//--------------------------------------------------------------------
// apply_mode_flip: commit a pure mode switch of beam b to mode m2 — evict
// only the served tasks whose type is incompatible with m2, switch the base
// power, and stop.  No refill: the freed capacity is left for subsequent
// insert/swap moves.  Matches the candidate block's delta = -eviction_loss.
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
    tabuBeamMode[b] = tabuIter + tenure;
    remPW[b] += modeBasePower[saveMode] - modeBasePower[m2];
    // (remPW >= 0 guaranteed by the power check that admitted this move)
}

//--------------------------------------------------------------------
void local_search(double beginTime)
{
    int ts_depth = 300;              // non-improving iterations before this restart ends
    int nonImprove = 0;

    tabuIter = 0;
    for (int j = 0; j < numTask; j++) tabuUntil[j] = 0;
    for (int b = 0; b < numBeam; b++) tabuBeamMode[b] = 0;

    // ---- solution-based tabu: start a fresh visited set for this phase ----
    // curHash tracks the current solution's fingerprint incrementally; visited
    // holds every fingerprint reached in this phase so a move that would land on
    // an already-seen solution fingerprint is forbidden.
    curHash = compute_hash();
    visited.clear();
    visited.insert(curHash);

    while (nonImprove < ts_depth)
    {
        if (((double)clock() - beginTime) / CLOCKS_PER_SEC > maxRunTime)
            break;

        int bestDelta = MININT_MOVE;     // allow 0 / negative moves (tabu search)
        int numBest   = 0;               // number of candidates tied at bestDelta
        int kind = 0, a = -1, bb = -1, cc = -1;
        int timedTabuBlocked = 0;         // true if waiting may reopen candidates
        unsigned long long chosenHash = 0;   // fingerprint of the chosen move's result

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
            if (tabuIter < tabuUntil[j]) { timedTabuBlocked = 1; continue; }
            int delta = taskProfit[j];
            for (int b = 0; b < numBeam; b++)
            {
                if (!feasible_on(j, b)) continue;
                if (delta < bestDelta) continue;
                // dH: taskBeam[j] moves from "unserved"(numBeam) to b
                unsigned long long cand =
                    curHash + taskW[j][b] - taskW[j][numBeam];
                if (visited.count(cand)) continue;       // solution-level tabu
                record_candidate(1, j, b, -1, delta, cand,
                                 bestDelta, numBest, kind, a, bb, cc, chosenHash);
            }
        }

        // (2) remove-only: drop a served task ------------------------
        for (int j = 0; j < numTask; j++)
        {
            if (taskBeam[j] < 0) continue;
            if (tabuIter < tabuUntil[j]) { timedTabuBlocked = 1; continue; }
            int delta = -taskProfit[j];
            if (delta < bestDelta) continue;

            // dH: taskBeam[j] moves from its beam to "unserved"(numBeam)
            unsigned long long cand =
                curHash + taskW[j][numBeam] - taskW[j][taskBeam[j]];
            if (visited.count(cand)) continue;           // solution-level tabu
            record_candidate(2, j, -1, -1, delta, cand,
                             bestDelta, numBest, kind, a, bb, cc, chosenHash);
        }

        // (3) swap: unserved i replaces served k on the same beam -------
        //     For each unserved i, only visit beams whose mode is compatible
        //     with i, and within such a beam only its served tasks (bucket),
        //     instead of rescanning every task.  Same (i,k) pairs as before.
        for (int i = 0; i < numTask; i++)
        {
            if (taskBeam[i] != -1) continue;
            if (tabuIter < tabuUntil[i]) { timedTabuBlocked = 1; continue; }
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
                    if (tabuIter < tabuUntil[k]) { timedTabuBlocked = 1; continue; }

                    int delta = taskProfit[i] - taskProfit[k];
                    if (delta < bestDelta) continue;

                    // feasibility after removing k, adding i (same beam b)
                    int bwFree = bwFree0 + taskBWDemand[k];
                    int pwFree = pwFree0 + taskPWDemand[k];
                    if (taskBWDemand[i] > bwFree || taskPWDemand[i] > pwFree) continue;

                    // dH: i moves unserved->b, k moves b->unserved
                    unsigned long long cand = curHash
                        + taskW[i][b] - taskW[i][numBeam]
                        + taskW[k][numBeam] - taskW[k][b];
                    if (visited.count(cand)) continue;   // solution-level tabu
                    record_candidate(3, i, b, k, delta, cand,
                                     bestDelta, numBest, kind, a, bb, cc, chosenHash);
                }
            }
        }

        // (4) mode-flip: switch beam b's mode, evicting ONLY the served tasks
        //     whose type is incompatible with the new mode.  No refill here:
        //     the freed room is left for later insert/swap moves to fill, so
        //     the move's delta is exactly the eviction loss (<= 0).  This keeps
        //     a flip honest about its real local cost (it no longer looks like
        //     a gain by greedily refilling) and lets insert/swap compete for
        //     the freed capacity on equal footing.
        for (int b = 0; b < numBeam; b++)
        {
            if (beamMode[b] < 0) continue;
            if (tabuIter < tabuBeamMode[b]) { timedTabuBlocked = 1; continue; }
            int usedPW_b = (beamPWCap[b] - modeBasePower[beamMode[b]]) - remPW[b];
            for (int m2 = 0; m2 < numMode; m2++)
            {
                if (m2 == beamMode[b]) continue;
                if (modeBasePower[m2] > beamPWCap[b]) continue;  // mode power-infeasible

                // one pass over beam b's served tasks: those whose type is
                // incompatible with m2 are the ones the switch would evict.  Sum
                // their profit (= the move's loss) and power, and accumulate the
                // hash delta of each evicted task (b -> unserved) on the same pass
                // so the flip's fingerprint stays O(bucket) with no extra scan.
                int loss = 0, evictPW = 0;
                int blockedByTaskTabu = 0;
                unsigned long long dHevict = 0;
                for (int idx = bucketStart[b]; idx < bucketStart[b + 1]; idx++)
                {
                    int j = bucketTask[idx];
                    if (!typeCompatMode[taskType[j] - 1][m2])
                    {
                        if (tabuIter < tabuUntil[j])
                        {
                            timedTabuBlocked = 1;
                            blockedByTaskTabu = 1;
                            break;
                        }
                        loss    += taskProfit[j];
                        evictPW += taskPWDemand[j];
                        dHevict += taskW[j][numBeam] - taskW[j][b];
                    }
                }
                if (blockedByTaskTabu) continue;
                if (usedPW_b - evictPW > beamPWCap[b] - modeBasePower[m2]) continue;

                int delta = -loss;                 // pure switch: exact, always <= 0
                if (delta < bestDelta) continue;

                // dH: beamMode[b] changes m->m2, plus every evicted task's move
                unsigned long long cand = curHash
                    + beamV[b][m2] - beamV[b][beamMode[b]] + dHevict;
                if (visited.count(cand)) continue;   // solution-level tabu
                record_candidate(4, b, m2, -1, delta, cand,
                                 bestDelta, numBest, kind, a, bb, cc, chosenHash);
            }
        }

        // ---- no admissible candidate this iteration -------------------
        // Attribute tabu has a time axis, so advance it instead of ending the
        // restart immediately; candidates frozen this iteration may reopen.
        if (numBest == 0)
        {
            if (!timedTabuBlocked) break;
            nonImprove++;
            tabuIter++;
            continue;
        }

        int tenure = tabu_tenure();
        if (kind == 1)                            // insert
        {
            add_task(a, bb);
            moveFreq[a]++;
            tabuUntil[a] = tabuIter + tenure;
        }
        else if (kind == 2)                       // remove-only
        {
            remove_task(a);
            moveFreq[a]++;
            tabuUntil[a] = tabuIter + tenure;
        }
        else if (kind == 3)                       // swap a in, cc out
        {
            remove_task(cc);
            add_task(a, bb);
            moveFreq[a]++; moveFreq[cc]++;
            tabuUntil[a]  = tabuIter + tenure;
            tabuUntil[cc] = tabuIter + tenure;
        }
        else if (kind == 4)                       // mode-flip beam a to mode bb
        {
            apply_mode_flip(a, bb, tenure);
        }

        // commit the new fingerprint and mark it visited for this phase
        curHash = chosenHash;
        visited.insert(curHash);

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

//--------------------------------------------------------------------
// perturb  (destroy + repair)
//
// Legacy ILS perturbation kept for comparison experiments.  The current
// MSBTS-GEO driver does not call it; diversity comes from randomized restarts.
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
    int    *prevMode = new int[numBeam];         // mode each emptied beam had before destroy
    double *beamFreqKey = new double[numBeam];
    int    *order = new int[numBeam];
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
        prevMode[k] = beamMode[b];                  // remember the old mode: must not be reselected
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

void msbts_geo()
{
    double beginTime = (double)clock();

    alloc_solution();
    alloc_search();

    globalProfit = -1;

    double runTime = 0.0;
    long   numRestart = 0;
    while (runTime < maxRunTime)
    {
        greedy_randomized_init();
        for (int j = 0; j < numTask; j++) moveFreq[j] = 0;

        save_best(beginTime);          // restart-best starts from the constructed solution
        local_search(beginTime);
        numRestart++;

        if (bestProfit > globalProfit)
        {
            save_global_from_best();
            cout << "  globalProfit=" << globalProfit
                 << "  time=" << globalBestTime << " s"
                 << "  restart=" << numRestart << endl;
        }

        runTime = ((double)clock() - beginTime) / CLOCKS_PER_SEC;
    }

    if (globalProfit < 0)
    {
        greedy_randomized_init();
        save_best(beginTime);
        save_global_from_best();
    }

    // install the global best as the final solution for checking / printing
    restore_global();
    totalProfit = globalProfit;

    cout << "MSBTS-GEO done.  bestProfit=" << globalProfit
         << "  bestTime=" << globalBestTime << " s"
         << "  restarts=" << numRestart << endl;

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
    delete[] tabuBeamMode;
    delete[] bucketTask;
    delete[] bucketStart;

    for (int j = 0; j < numTask; j++) delete[] taskW[j];
    delete[] taskW;
    for (int b = 0; b < numBeam; b++) delete[] beamV[b];
    delete[] beamV;
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
    msbts_geo();
    check_solution();
    print_solution();

    double elapsed = ((double)clock() - t0) / CLOCKS_PER_SEC;
    cout << "\nElapsed time: " << elapsed << " s" << endl;

    free_memory();
    return 0;
}
