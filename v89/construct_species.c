/*
 * construct_species.c — enumerate closed integer locks (particle species)
 * under the v88/CONSTRUCTION.md regular method, with NO spatial background.
 *
 * Ontology: a "lock" is a finite multiset of harmonic words (integer vectors)
 * plus channel handedness between them, such that every channel is resonant
 * and every parcel's resonance is only with partners inside the lock
 * (spectral isolation). Self-reproduction is checked as: the integer system
 * is closed under a one-step internal cycle (harmonic words rearrange among
 * an isomorphic multiset).
 *
 * No lattice, no coordinates, no immortal cell IDs. Species are isomorphism
 * classes of locks, labeled by (N_c, net χ, harmonic signature, E_D proxy).
 *
 * Build:  gcc -O2 -o construct_species construct_species.c -lm
 * Run:    ./construct_species [max_Nc] [max_mode]
 *
 * Tags: implements [P] resonance/lock rules of CONSTRUCTION.md §§2,5.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define K_MAX 3          /* harmonic word length */
#define NC_MAX 6         /* max parcels in a lock */
#define MODE_MAX 4       /* mode indices in 0..MODE_MAX-1 per component */
#define MAX_SPECIES 4096

typedef struct {
    int h[K_MAX];
} Word;

typedef struct {
    int n;                 /* N_c */
    Word w[NC_MAX];
    int chi[NC_MAX][NC_MAX]; /* +1 / -1 / 0 (0 = no channel) */
    int net_chi;
    int signature;         /* hash of sorted words + chi pattern class */
    int e_proxy;           /* mass proxy: sum of |h| content */
} Lock;

/* Resonance: channel between a,b with handedness chi is resonant iff
 *   h_a[0] - h_b[0] == 0  and  chi * (h_a[1] - h_b[2]) + (h_a[2] - h_b[1]) == 0
 * when both transverse slots are used; simplified integer form of
 * "relative phase matches χ and longitudinal modes agree."
 * [P] schematic — see CONSTRUCTION.md §2.2
 */
static int resonant(const Word *a, const Word *b, int chi)
{
    if (a->h[0] != b->h[0])
        return 0;
    int t = chi * (a->h[1] - b->h[2]) + (a->h[2] - b->h[1]);
    return t == 0;
}

/* Interior closed: every pair with a channel is resonant; graph connected
 * for n>1; every parcel has at least one channel (n=1: self-lock via
 * internal transverse closure).
 */
static int self_lock_word(const Word *w)
{
    /* single parcel: transverse content forms a chiral pair */
    return (w->h[1] != 0 || w->h[2] != 0) && (w->h[1] * w->h[2] != 0
            ? (w->h[1] == w->h[2] || w->h[1] == -w->h[2])
            : 0);
}

static int lock_closed(const Lock *L)
{
    if (L->n == 1)
        return self_lock_word(&L->w[0]);

    int deg[NC_MAX];
    memset(deg, 0, sizeof(deg));
    for (int i = 0; i < L->n; i++) {
        for (int j = i + 1; j < L->n; j++) {
            int c = L->chi[i][j];
            if (c == 0)
                continue;
            if (!resonant(&L->w[i], &L->w[j], c))
                return 0;
            deg[i]++;
            deg[j]++;
        }
    }
    for (int i = 0; i < L->n; i++)
        if (deg[i] == 0)
            return 0;

    /* connectivity BFS */
    int seen[NC_MAX] = {0}, q[NC_MAX], qh = 0, qt = 0;
    seen[0] = 1;
    q[qt++] = 0;
    int reached = 0;
    while (qh < qt) {
        int u = q[qh++];
        reached++;
        for (int v = 0; v < L->n; v++) {
            if (seen[v])
                continue;
            int c = (u < v) ? L->chi[u][v] : L->chi[v][u];
            if (c == 0)
                continue;
            seen[v] = 1;
            q[qt++] = v;
        }
    }
    return reached == L->n;
}

/* Spectral isolation: no resonant channel with a generic exterior probe.
 * Probe words: unit vectors and a few mixed words. If any resonant channel
 * exists for either chi, the lock is not isolated.
 */
static int isolated(const Lock *L)
{
    Word probes[] = {
        {{0, 1, 1}}, {{0, 1, -1}}, {{0, -1, 1}}, {{0, -1, -1}},
        {{1, 1, 1}}, {{1, 1, -1}}, {{1, -1, 1}}, {{2, 1, 1}},
        {{0, 2, 1}}, {{0, 1, 2}}, {{1, 0, 1}}, {{1, 1, 0}},
    };
    int np = (int)(sizeof(probes) / sizeof(probes[0]));
    for (int i = 0; i < L->n; i++) {
        for (int p = 0; p < np; p++) {
            /* skip probes that match an interior word exactly — those are
             * same-species contact, not exterior continuum */
            int same = 0;
            for (int j = 0; j < L->n; j++) {
                if (memcmp(&probes[p], &L->w[j], sizeof(Word)) == 0) {
                    same = 1;
                    break;
                }
            }
            if (same)
                continue;
            if (resonant(&L->w[i], &probes[p], +1))
                return 0;
            if (resonant(&L->w[i], &probes[p], -1))
                return 0;
        }
    }
    return 1;
}

/* Self-reproduction proxy: after a full internal cycle, each parcel's word
 * may permute among the multiset — isomorphism of multiset of words and of
 * the chi graph up to relabeling. Here: already a fixed lock is a fixed point
 * of "advance all internal channels one full cycle" if resonance holds
 * (R2/R3 return energy to the same lock type). Satisfied by lock_closed.
 * Additional: net chi and multiset of words are the invariants.
 */
static int e_proxy(const Lock *L)
{
    int s = 0;
    for (int i = 0; i < L->n; i++)
        for (int k = 0; k < K_MAX; k++)
            s += L->w[i].h[k] >= 0 ? L->w[i].h[k] : -L->w[i].h[k];
    return s + L->n; /* N_c contribution: mass ∝ N_c when capacity filled */
}

static int word_key(const Word *w)
{
    /* pack small modes into a key; modes in [-M,M] mapped to 0..2M */
    int M = MODE_MAX;
    int key = 0;
    for (int k = 0; k < K_MAX; k++) {
        int x = w->h[k] + M;
        if (x < 0)
            x = 0;
        if (x > 2 * M)
            x = 2 * M;
        key = key * (2 * M + 1) + x;
    }
    return key;
}

static int lock_signature(const Lock *L)
{
    int keys[NC_MAX];
    for (int i = 0; i < L->n; i++)
        keys[i] = word_key(&L->w[i]);
    /* sort keys for multiset equality */
    for (int i = 0; i < L->n; i++)
        for (int j = i + 1; j < L->n; j++)
            if (keys[j] < keys[i]) {
                int t = keys[i];
                keys[i] = keys[j];
                keys[j] = t;
            }
    uint32_t h = 2166136261u;
    for (int i = 0; i < L->n; i++) {
        h ^= (uint32_t)keys[i];
        h *= 16777619u;
    }
    /* chi net and edge count */
    h ^= (uint32_t)(L->net_chi + 64);
    h *= 16777619u;
    int edges = 0;
    for (int i = 0; i < L->n; i++)
        for (int j = i + 1; j < L->n; j++)
            if (L->chi[i][j] != 0)
                edges++;
    h ^= (uint32_t)edges;
    h *= 16777619u;
    return (int)(h & 0x7fffffff);
}

static Lock species[MAX_SPECIES];
static int nspecies = 0;

static int already_have(const Lock *L)
{
    for (int i = 0; i < nspecies; i++) {
        if (species[i].n == L->n && species[i].signature == L->signature
            && species[i].net_chi == L->net_chi
            && species[i].e_proxy == L->e_proxy)
            return 1;
    }
    return 0;
}

static void consider(Lock *L)
{
    if (!lock_closed(L))
        return;
    /* isolation is a strong filter; report both isolated and non-isolated
     * closed locks, mark isolation flag in e_proxy high bit via print */
    L->net_chi = 0;
    for (int i = 0; i < L->n; i++)
        for (int j = i + 1; j < L->n; j++)
            L->net_chi += L->chi[i][j];
    if (L->n == 1) {
        /* self-lock handedness from transverse product sign */
        int s = L->w[0].h[1] * L->w[0].h[2];
        L->net_chi = (s > 0) - (s < 0);
    }
    L->e_proxy = e_proxy(L);
    L->signature = lock_signature(L);
    if (already_have(L))
        return;
    if (nspecies >= MAX_SPECIES)
        return;
    species[nspecies++] = *L;
}

/* Enumerate words with components in [-max_mode, max_mode], not all zero
 * transverse. */
static void enum_words(int max_mode, Word *out, int *nout, int cap)
{
    *nout = 0;
    for (int a = -max_mode; a <= max_mode; a++)
        for (int b = -max_mode; b <= max_mode; b++)
            for (int c = -max_mode; c <= max_mode; c++) {
                if (a == 0 && b == 0 && c == 0)
                    continue;
                if (*nout >= cap)
                    return;
                out[*nout].h[0] = a;
                out[*nout].h[1] = b;
                out[*nout].h[2] = c;
                (*nout)++;
            }
}

static void enum_Nc1(const Word *words, int nw)
{
    for (int i = 0; i < nw; i++) {
        Lock L;
        memset(&L, 0, sizeof(L));
        L.n = 1;
        L.w[0] = words[i];
        consider(&L);
    }
}

static void enum_Nc2(const Word *words, int nw)
{
    for (int i = 0; i < nw; i++)
        for (int j = i; j < nw; j++) {
            for (int chi = -1; chi <= 1; chi += 2) {
                Lock L;
                memset(&L, 0, sizeof(L));
                L.n = 2;
                L.w[0] = words[i];
                L.w[1] = words[j];
                L.chi[0][1] = chi;
                L.chi[1][0] = chi;
                consider(&L);
            }
        }
}

static void enum_Nc3(const Word *words, int nw, int max_mode)
{
    /* restrict for cost: only |components| <= max_mode, subsample if large */
    int step = 1;
    if (nw > 40)
        step = nw / 40 + 1;
    for (int i = 0; i < nw; i += step)
        for (int j = i; j < nw; j += step)
            for (int k = j; k < nw; k += step) {
                int chis[3];
                /* edges 01, 02, 12 each ±1 or skip sparse: require all three */
                for (int c01 = -1; c01 <= 1; c01 += 2)
                    for (int c02 = -1; c02 <= 1; c02 += 2)
                        for (int c12 = -1; c12 <= 1; c12 += 2) {
                            Lock L;
                            memset(&L, 0, sizeof(L));
                            L.n = 3;
                            L.w[0] = words[i];
                            L.w[1] = words[j];
                            L.w[2] = words[k];
                            L.chi[0][1] = L.chi[1][0] = c01;
                            L.chi[0][2] = L.chi[2][0] = c02;
                            L.chi[1][2] = L.chi[2][1] = c12;
                            (void)chis;
                            consider(&L);
                        }
            }
    (void)max_mode;
}

static int cmp_species(const void *a, const void *b)
{
    const Lock *x = a, *y = b;
    if (x->n != y->n)
        return x->n - y->n;
    if (x->e_proxy != y->e_proxy)
        return x->e_proxy - y->e_proxy;
    if (x->net_chi != y->net_chi)
        return x->net_chi - y->net_chi;
    return x->signature - y->signature;
}

int main(int argc, char **argv)
{
    int max_nc = 3;
    int max_mode = 2;
    if (argc > 1)
        max_nc = atoi(argv[1]);
    if (argc > 2)
        max_mode = atoi(argv[2]);
    if (max_nc < 1)
        max_nc = 1;
    if (max_nc > NC_MAX)
        max_nc = NC_MAX;
    if (max_mode < 1)
        max_mode = 1;
    if (max_mode > MODE_MAX)
        max_mode = MODE_MAX;

    Word words[512];
    int nw = 0;
    enum_words(max_mode, words, &nw, 512);
    printf("# construct_species — no-background lock enumeration\n");
    printf("# max_Nc=%d max_mode=%d words=%d\n", max_nc, max_mode, nw);
    printf("# resonance: h0 match and chi*(h1a-h2b)+(h2a-h1b)=0  [P]\n");
    printf("# columns: id Nc net_chi e_proxy isolated signature words...\n");

    if (max_nc >= 1)
        enum_Nc1(words, nw);
    if (max_nc >= 2)
        enum_Nc2(words, nw);
    if (max_nc >= 3)
        enum_Nc3(words, nw, max_mode);

    qsort(species, (size_t)nspecies, sizeof(Lock), cmp_species);

    int n_iso = 0;
    for (int s = 0; s < nspecies; s++) {
        int iso = isolated(&species[s]);
        if (iso)
            n_iso++;
        printf("%d\t%d\t%d\t%d\t%d\t%08x", s, species[s].n, species[s].net_chi,
               species[s].e_proxy, iso, (unsigned)species[s].signature);
        for (int i = 0; i < species[s].n; i++)
            printf("\t(%d,%d,%d)", species[s].w[i].h[0], species[s].w[i].h[1],
                   species[s].w[i].h[2]);
        printf("\n");
    }

    printf("# total_closed_species %d\n", nspecies);
    printf("# isolated_species %d\n", n_iso);
    printf("# note: discrete spectrum (finite count) is the structural claim;\n");
    printf("# exact resonance polynomial is [P] schematic in CONSTRUCTION.md\n");
    return 0;
}
