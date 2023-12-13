/* ----------------------------------------------------------------------- */
/* Tunable parameters                                                      */
/* ----------------------------------------------------------------------- */

/* Frequency of comprehensive reports-- lower values will provide more
 * info while slowing down the simulation. Higher values will give less
 * frequent updates. */
/* This is also the frequency of screen refreshes if SDL is enabled. */
#define REPORT_FREQUENCY 2000000

/* Mutation rate -- range is from 0 (none) to 0xffffffff (all mutations!) */
/* To get it from a float probability from 0.0 to 1.0, multiply it by
 * 4294967295 (0xffffffff) and round. */
#define MUTATION_RATE 5000

/* How frequently should random cells / energy be introduced?
 * Making this too high makes things very chaotic. Making it too low
 * might not introduce enough energy. */
#define INFLOW_FREQUENCY 100

/* Base amount of energy to introduce per INFLOW_FREQUENCY ticks */
#define INFLOW_RATE_BASE 600

/* A random amount of energy between 0 and this is added to
 * INFLOW_RATE_BASE when energy is introduced. Comment this out for
 * no variation in inflow rate. */
#define INFLOW_RATE_VARIATION 1000

/* Size of pond in X and Y dimensions. */
#define POND_SIZE_X 800
#define POND_SIZE_Y 600

/* Depth of pond in four-bit codons -- this is the maximum
 * genome size. This *must* be a multiple of 16! */
#define POND_DEPTH 1024

/* This is the divisor that determines how much energy is taken
 * from cells when they try to KILL a viable cell neighbor and
 * fail. Higher numbers mean lower penalties. */
#define FAILED_KILL_PENALTY 3

/* Define this to use SDL. To use SDL, you must have SDL headers
 * available and you must link with the SDL library when you compile. */
/* Comment this out to compile without SDL visualization support. */
//#define USE_SDL 1

/* Define this to use threads, and how many threads to create */
// #define USE_PTHREADS_COUNT 4

/* ----------------------------------------------------------------------- */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


#define BUFFER_SIZE 1000  // Size of the circular buffer

/* Pond depth in machine-size words.  This is calculated from
 * POND_DEPTH and the size of the machine word. (The multiplication
 * by two is due to the fact that there are two four-bit values in
 * each eight-bit byte.) */
#define POND_DEPTH_SYSWORDS (POND_DEPTH / (sizeof(uintptr_t) * 2))

/* Number of bits in a machine-size word */
#define SYSWORD_BITS (sizeof(uintptr_t) * 8)

/* Constants representing neighbors in the 2D grid. */
#define N_LEFT 0
#define N_RIGHT 1
#define N_UP 2
#define N_DOWN 3

/* Word and bit at which to start execution */
/* This is after the "logo" */
#define EXEC_START_WORD 0
#define EXEC_START_BIT 4

/* Number of bits set in binary numbers 0000 through 1111 */
static const uintptr_t BITS_IN_FOURBIT_WORD[16] = { 0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4 };

static uintptr_t buffer[BUFFER_SIZE];
static int in = 0;
static uintptr_t last_random_number;

volatile uint64_t prngState[2];

static inline uintptr_t getRandomPre(int rollback)
{
	// https://en.wikipedia.org/wiki/Xorshift#xorshift.2B
	uint64_t x = prngState[0];
	const uint64_t y = prngState[1];
	prngState[0] = prngState[0] * !rollback + rollback * y;
	x ^= x << 23;
	const uint64_t z = x ^ y ^ (x >> 17) ^ (y >> 26);
	prngState[1] = prngState[1] * !rollback + rollback * z;
	return (uintptr_t)(z + y);
}
void precalculate_random_numbers() {
    for (int i = 0; i < BUFFER_SIZE; i++) {
        buffer[i] = getRandomPre(1);
    }
}
static inline uintptr_t getRandomRollback(uintptr_t rollback) {
    uintptr_t num = buffer[in];
    last_random_number = num;  // Store the last random number
    uintptr_t new_num = getRandomPre(rollback);
    buffer[in] = (new_num & -rollback) | (num & ~-rollback);
    in = (((in + 1) & -rollback) | (in & ~-rollback)) % BUFFER_SIZE;
    return num;
}

/**
 * Structure for a cell in the pond
 */
struct Cell
        {
	/* Globally unique cell ID */
	uint64_t ID;
	
	/* ID of the cell's parent */
	uint64_t parentID;
	
	/* Counter for original lineages -- equal to the cell ID of
	 * the first cell in the line. */
	uint64_t lineage;
	
	/* Generations start at 0 and are incremented from there. */
	uintptr_t generation;
	
	/* Energy level of this cell */
	uintptr_t energy;

	/* Memory space for cell genome (genome is stored as four
	 * bit instructions packed into machine size words) */
	uintptr_t genome[POND_DEPTH_SYSWORDS];

#ifdef USE_PTHREADS_COUNT
	pthread_mutex_t lock;
#endif
};

/* The pond is a 2D array of cells */
static struct Cell pond[POND_SIZE_X][POND_SIZE_Y];

/* This is used to generate unique cell IDs */
static volatile uint64_t cellIdCounter = 0;

volatile struct {
	/* Counts for the number of times each instruction was
	 * executed since the last report. */
	double instructionExecutions[16];
	
	/* Number of cells executed since last report */
	double cellExecutions;
	
	/* Number of viable cells replaced by other cells' offspring */
	uintptr_t viableCellsReplaced;
	
	/* Number of viable cells KILLed */
	uintptr_t viableCellsKilled;
	
	/* Number of successful SHARE operations */
	uintptr_t viableCellShares;
} statCounters;

static void doReport(const uint64_t clock)
{
	static uint64_t lastTotalViableReplicators = 0;
	
	uintptr_t x,y;
	
	uint64_t totalActiveCells = 0;
	uint64_t totalEnergy = 0;
	uint64_t totalViableReplicators = 0;
	uintptr_t maxGeneration = 0;
	
	for(x=0;x<POND_SIZE_X;++x) {
		for(y=0;y<POND_SIZE_Y;++y) {
			struct Cell *const c = &pond[x][y];
			if (c->energy) {
				++totalActiveCells;
				totalEnergy += (uint64_t)c->energy;
				if (c->generation > 2)
					++totalViableReplicators;
				if (c->generation > maxGeneration)
					maxGeneration = c->generation;
			}
		}
	}
	
	/* Look here to get the columns in the CSV output */
	
	/* The first five are here and are self-explanatory */
	printf("%lu,%lu,%lu,%lu,%lu,%lu,%lu,%lu",
		(uint64_t)clock,
		(uint64_t)totalEnergy,
		(uint64_t)totalActiveCells,
		(uint64_t)totalViableReplicators,
		(uint64_t)maxGeneration,
		(uint64_t)statCounters.viableCellsReplaced,
		(uint64_t)statCounters.viableCellsKilled,
		(uint64_t)statCounters.viableCellShares
		);
	
	/* The next 16 are the average frequencies of execution for each
	 * instruction per cell execution. */
	double totalMetabolism = 0.0;
	for(x=0;x<16;++x) {
		totalMetabolism += statCounters.instructionExecutions[x];
		printf(",%.4f",(statCounters.cellExecutions > 0.0) ? (statCounters.instructionExecutions[x] / statCounters.cellExecutions) : 0.0);
	}
	
	/* The last column is the average metabolism per cell execution */
	printf(",%.4f\n",(statCounters.cellExecutions > 0.0) ? (totalMetabolism / statCounters.cellExecutions) : 0.0);
	fflush(stdout);
	
	if ((lastTotalViableReplicators > 0)&&(totalViableReplicators == 0))
		fprintf(stderr,"[EVENT] Viable replicators have gone extinct. Please reserve a moment of silence.\n");
	else if ((lastTotalViableReplicators == 0)&&(totalViableReplicators > 0))
		fprintf(stderr,"[EVENT] Viable replicators have appeared!\n");
	
	lastTotalViableReplicators = totalViableReplicators;
	
	/* Reset per-report stat counters */
	for(x=0;x<sizeof(statCounters);++x)
		((uint8_t *)&statCounters)[x] = (uint8_t)0;
}

static inline struct Cell *getNeighbor(const uintptr_t x, const uintptr_t y, const uintptr_t dir)
{
    /* Define the changes in the x and y coordinates for each direction */
    int dx[] = {-1, 1, 0, 0}; // Changes in x for N_LEFT, N_RIGHT, N_UP, N_DOWN
    int dy[] = {0, 0, -1, 1}; // Changes in y for N_LEFT, N_RIGHT, N_UP, N_DOWN

    /* Calculate the new coordinates */
    uintptr_t newX = (x + dx[dir] + POND_SIZE_X) % POND_SIZE_X;
    uintptr_t newY = (y + dy[dir] + POND_SIZE_Y) % POND_SIZE_Y;

    return &pond[newX][newY];
}

static inline int accessAllowed(struct Cell *const c2, const uintptr_t c1guess, int sense, int rollback)
{
	uintptr_t random = (uintptr_t)(getRandomRollback(rollback) & 0xf);
	return ((((random >= BITS_IN_FOURBIT_WORD[(c2->genome[0] & 0xf) ^ (c1guess & 0xf)]) || !c2->parentID) & sense) | (((random <= BITS_IN_FOURBIT_WORD[(c2->genome[0] & 0xf) ^ (c1guess & 0xf)]) || !c2->parentID) & ~sense));
}

static void *run(uintptr_t *ptr_wordPtr, uintptr_t *ptr_shiftPtr, uintptr_t *reg, 
uintptr_t *loopStackPtr, uintptr_t *wordPtr, uintptr_t *shiftPtr, uintptr_t *facing, 
uintptr_t *falseLoopDepth, int *stop, int *skip, int *access_neg_used, 
int *access_pos_used, uintptr_t *outputBuf, uintptr_t *clock, struct Cell *pptr, struct Cell *tmpptr)
{
    uintptr_t x,y,i;

    if (!(clock % INFLOW_FREQUENCY)) {
        x = getRandomRollback(1) % POND_SIZE_X;
        y = getRandomRollback(1) % POND_SIZE_Y;
        pptr = &pond[x][y];
        pptr->ID = cellIdCounter;
        pptr->parentID = 0;
        pptr->lineage = cellIdCounter;
        pptr->generation = 0;
#ifdef INFLOW_RATE_VARIATION
        pptr->energy += INFLOW_RATE_BASE + (getRandomRollback(1) % INFLOW_RATE_VARIATION);
#else
        pptr->energy += INFLOW_RATE_BASE;
#endif /* INFLOW_RATE_VARIATION */
        for(i=0;i<POND_DEPTH_SYSWORDS;++i) 
            pptr->genome[i] = getRandomRollback(1);
        ++cellIdCounter;
		}

		/* Pick a random cell to execute */
		i = getRandomRollback(1);
		x = i % POND_SIZE_X;
		y = ((i / POND_SIZE_X) >> 1) % POND_SIZE_Y;
		pptr = &pond[x][y];

		/* Reset the state of the VM prior to execution */
		for(i=0;i<POND_DEPTH_SYSWORDS;++i)
			outputBuf[i] = ~((uintptr_t)0); /* ~0 == 0xfffff... */
	
		ptr_wordPtr = 0;
		ptr_shiftPtr = 0;
		reg = 0;
		loopStackPtr = 0;
		wordPtr = EXEC_START_WORD;
		shiftPtr = EXEC_START_BIT;
		facing = 0;
		falseLoopDepth = 0;
		stop = 0;
        skip=0;
		access_neg_used = 0;
		access_pos_used = 0;
		access_neg = 0;
		access_pos = 0;

		/* We use a currentWord buffer to hold the word we're
		 * currently working on.  This speeds things up a bit
		 * since it eliminates a pointer dereference in the
		 * inner loop. We have to be careful to refresh this
		 * whenever it might have changed... take a look at
		 * the code. :) */
		currentWord = pptr->genome[0];

		/* Keep track of how many cells have been executed */
		statCounters.cellExecutions += 1.0;

		/* Core execution loop */
		while ((pptr->energy)&&(!stop)) {
			/* Get the next instruction */
			inst = (currentWord >> shiftPtr) & 0xf;
            skip=0;

			/* Randomly frob either the instruction or the register with a
			 * probability defined by MUTATION_RATE. This introduces variation,
			 * and since the variation is introduced into the state of the VM
			 * it can have all manner of different effects on the end result of
			 * replication: insertions, deletions, duplications of entire
			 * ranges of the genome, etc. */
			if ((getRandomRollback(1) & 0xffffffff) < MUTATION_RATE) {
				tmp = getRandomRollback(1); // Call getRandom() only once for speed 
				if (tmp & 0x80) // Check for the 8th bit to get random boolean 
					inst = tmp & 0xf; // Only the first four bits are used here 
				else reg = tmp & 0xf;
			}
			
            /*
			uintptr_t mutation_occurred = (getRandomRollback(1) & 0xffffffff) < MUTATION_RATE;
			uintptr_t tmp = getRandomRollback(mutation_occurred) * mutation_occurred;
			uintptr_t is_inst = (tmp & 0x80) >> 7; // Shift right by 7 to get a 1 or 0
			uintptr_t is_reg = ~is_inst & 0x1; // Invert is_inst and mask with 0x1 to get a 1 or 0
			inst = (tmp & 0xf) * is_inst + inst * (!is_inst); // Update inst only if is_inst is 1
			reg = (tmp & 0xf) * is_reg + reg * (!is_reg); // Update reg only if is_reg is 1
			*/

			/* Each instruction processed costs one unit of energy */
			--pptr->energy;

			/* Execute the instruction */
			if (falseLoopDepth) {
				/* Skip forward to matching REP if we're in a false loop. */
				if (inst == 0x9) /* Increment false LOOP depth */
					++falseLoopDepth;
				else if (inst == 0xa) /* Decrement on REP */
					--falseLoopDepth;
			} else {

				ptr_shiftPtr = (inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xd || inst == 0xe || inst == 0xf) * (ptr_shiftPtr) +((inst == 0x0)*0)+((inst == 0x1)*((ptr_shiftPtr+4)*((ptr_shiftPtr+4)<SYSWORD_BITS)))+((inst == 0x2)*(((ptr_shiftPtr==0)*SYSWORD_BITS)+ptr_shiftPtr-4)); 
                ptr_wordPtr = (inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xd || inst == 0xe || inst == 0xf) * (ptr_wordPtr) +((inst == 0x0)*0)+((inst == 0x1)*(((ptr_wordPtr*(ptr_shiftPtr!=0||((ptr_wordPtr+1)<POND_DEPTH_SYSWORDS))+(ptr_shiftPtr==0)*((ptr_wordPtr+1)<POND_DEPTH_SYSWORDS)))))+((inst == 0x2)*(((ptr_wordPtr==0&&ptr_shiftPtr==(SYSWORD_BITS-4))*(POND_DEPTH_SYSWORDS))+ptr_wordPtr-(ptr_shiftPtr==(SYSWORD_BITS-4))));
				wordPtr=(inst==0x0||inst==0x1||inst==0x2||inst==0x3||inst==0x4||inst==0x5|| inst == 0x6 || inst==0x7||inst==0x8||inst==0x9||inst==0xb||inst==0xd||inst==0xe||inst==0xf)*(wordPtr)+((inst==0xa)*(wordPtr*!(reg&&loopStackPtr)+(loopStack_wordPtr[loopStackPtr-1]*(reg&&loopStackPtr))))+((inst==0xc)*(wordPtr*((shiftPtr+4<SYSWORD_BITS)||(wordPtr+1<POND_DEPTH_SYSWORDS))+((shiftPtr+4>=SYSWORD_BITS)&&(wordPtr+1<POND_DEPTH_SYSWORDS))+EXEC_START_WORD*((wordPtr+1>=POND_DEPTH_SYSWORDS)&&(shiftPtr+4>=SYSWORD_BITS)))); 
                shiftPtr=(inst==0x0||inst==0x1||inst==0x2||inst==0x3||inst==0x4||inst==0x5|| inst == 0x6 || inst==0x7||inst==0x8||inst==0x9||inst==0xb||inst==0xd||inst==0xe||inst==0xf)*(shiftPtr)+((inst==0xa)*(shiftPtr*!(reg&&loopStackPtr)+(loopStack_shiftPtr[loopStackPtr-1]*(reg&&loopStackPtr))))+((inst==0xc)*((shiftPtr+4)+(shiftPtr+4>=SYSWORD_BITS)*(-shiftPtr-4)));
                skip=(reg&&loopStackPtr)*(inst==0xa);
				facing=(inst==0x1||inst==0x2||inst==0x3||inst==0x4||inst==0x5||inst==0x6||inst==0x7||inst==0x8||inst==0x9||inst==0xa||inst==0xc||inst==0xd||inst==0xe||inst==0xf)*(facing) + ((inst==0x0)*0)+((inst==0xb)*(reg & 3));
				pptr->genome[ptr_wordPtr]=(inst==0x0||inst==0x1||inst==0x2||inst==0x3||inst==0x4||inst==0x5||inst==0x7||inst==0x8||inst==0x9||inst==0xa||inst==0xb||inst==0xc||inst==0xd||inst==0xe||inst==0xf)*(pptr->genome[ptr_wordPtr])+((inst==0x6)*((pptr->genome[ptr_wordPtr]&~(((uintptr_t)0xf)<<ptr_shiftPtr))|reg<<ptr_shiftPtr)); 
				outputBuf[ptr_wordPtr]=(inst==0x0||inst==0x1||inst==0x2||inst==0x3||inst==0x4||inst==0x5|| inst == 0x6 || inst==0x7||inst==0x9||inst==0xa||inst==0xb||inst==0xc||inst==0xd||inst==0xe||inst==0xf)*(outputBuf[ptr_wordPtr])+((inst==0x8)*((outputBuf[ptr_wordPtr]&~(((uintptr_t)0xf) << ptr_shiftPtr))|reg << ptr_shiftPtr));
				currentWord=(inst==0x0||inst==0x1||inst==0x2||inst==0x3||inst==0x4||inst==0x5||inst==0x7||inst==0x8||inst==0x9|| inst==0xb || inst==0xd||inst==0xe||inst==0xf)*(currentWord)+((inst==0x6)*(pptr->genome[wordPtr]))+((inst==0xa)*(currentWord*!(reg&&loopStackPtr)+(pptr->genome[wordPtr])*(reg&&loopStackPtr)))+((inst == 0xc)*(pptr->genome[wordPtr]));
				stop=(inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xd || inst == 0xe)*(stop)+((inst == 0x9)*(stop*!(reg&&(loopStackPtr>=POND_DEPTH))+(reg&&(loopStackPtr>=POND_DEPTH))))+((inst == 0xf)*(1));
				loopStack_wordPtr[loopStackPtr]=(inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xd || inst == 0xe || inst == 0xf)*(loopStack_wordPtr[loopStackPtr])+((inst == 0x9)*(loopStack_wordPtr[loopStackPtr]*(!reg||(loopStackPtr>=POND_DEPTH))+(wordPtr*(reg&&(loopStackPtr<POND_DEPTH)))));
				loopStack_shiftPtr[loopStackPtr]=(inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xd || inst == 0xe || inst == 0xf)*(loopStack_shiftPtr[loopStackPtr])+((inst == 0x9) * (loopStack_shiftPtr[loopStackPtr]*(!reg||(loopStackPtr>=POND_DEPTH))+(shiftPtr*(reg&&(loopStackPtr<POND_DEPTH)))));
				loopStackPtr=(inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0xb || inst == 0xc || inst == 0xd || inst == 0xe || inst == 0xf)*(loopStackPtr)+((inst == 0x9)*(loopStackPtr + (reg&&(loopStackPtr<POND_DEPTH))))+((inst == 0xa)*(loopStackPtr-!!loopStackPtr));
				falseLoopDepth=(inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xd || inst == 0xe || inst == 0xf)*(falseLoopDepth)+((inst == 0x9)*(falseLoopDepth + (!reg)));
				tmpptr = getNeighbor(x,y,facing);
				access_neg_used = 0;
				access_pos_used = 0;
				access_pos_used = (inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xd || inst == 0xf)*(access_pos_used)+((inst == 0xe)*(1));
				access_neg_used = (inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xe|| inst == 0xf)*(access_neg_used)+((inst == 0xd)*(1));
				access_neg = accessAllowed(tmpptr,reg,0, access_neg_used);
				access_pos = accessAllowed(tmpptr,reg,1, access_pos_used);
				statCounters.viableCellsKilled=(inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xf)*(statCounters.viableCellsKilled)+((inst == 0xd)*(statCounters.viableCellsKilled+(access_neg)*(tmpptr->generation>2)))+((inst == 0xe)*(statCounters.viableCellsKilled+(access_pos)*(tmpptr->generation>2)));
				tmpptr->genome[0]=(inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xe || inst == 0xf)*(tmpptr->genome[0])+((inst == 0xd)*(tmpptr->genome[0]*!(access_neg)+(access_neg)*~((uintptr_t)0)));
				tmpptr->genome[1]=(inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xe || inst == 0xf)*(tmpptr->genome[1])+((inst == 0xd)*(tmpptr->genome[0]*!(access_neg)+(access_neg)*~((uintptr_t)0)));
				tmpptr->ID=(inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xe || inst == 0xf)*(tmpptr->ID)+((inst == 0xd)*(tmpptr->ID * !(access_neg)+ (access_neg)*cellIdCounter));
				tmpptr->parentID=(inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xe || inst == 0xf)*(tmpptr->parentID)+((inst == 0xd)*(tmpptr->parentID * !(access_neg)));
				tmpptr->lineage=(inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xe || inst == 0xf)*(tmpptr->lineage)+((inst == 0xd)*(tmpptr->lineage * !(access_neg) + (access_neg)*cellIdCounter));
				cellIdCounter=(inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xe || inst == 0xf)*(cellIdCounter)+((inst == 0xd)*(cellIdCounter * !(access_neg) + (access_neg)* cellIdCounter));
				tmp = (inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xf)*(tmp)+((inst == 0xc)*(reg))+((inst == 0xd)*((access_neg) + (tmpptr->generation>2)*!(access_neg)*(pptr->energy / FAILED_KILL_PENALTY)))+((inst == 0xe)* (pptr->energy + tmpptr->energy));
				reg = (inst == 0x1 || inst == 0x2 || inst == 0x6 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xd ||inst == 0xe || inst == 0xf) * (reg) + ((inst==0x0)*0) + ((inst==0x3)*((reg + 1) & 0xf)) +((inst==0x4)*((reg - 1) & 0xf)) +((inst==0x5)*((pptr->genome[ptr_wordPtr] >> ptr_shiftPtr) & 0xf)) +((inst==0x7)*((outputBuf[ptr_wordPtr] >> ptr_shiftPtr) & 0xf)) +((inst==0xc)*((pptr->genome[wordPtr] >> shiftPtr) & 0xf));
				pptr->genome[wordPtr]= (inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xd || inst == 0xe || inst == 0xf)*(pptr->genome[wordPtr])+((inst == 0xc)* (((pptr->genome[wordPtr]&~(((uintptr_t)0xf) << shiftPtr))|tmp << shiftPtr)));
				pptr->energy = (inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xf)*(pptr->energy)+((inst == 0xd)*(pptr->energy+!(access_neg)*(tmpptr->generation>2)*(-pptr->energy) + !(access_neg)*(tmpptr->generation>2)*(pptr->energy-tmp)))+((inst == 0xe)*((access_pos * (tmp - (access_pos * (tmp / 2) + (1 - access_pos) * tmpptr->energy)) + (1 - access_pos) * pptr->energy)));
				tmpptr->generation = (inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xe || inst == 0xf)*(tmpptr->generation)+((inst == 0xd)*(tmpptr->generation * (access_neg)));
				tmpptr->energy=(inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xd || inst == 0xf)*(tmpptr->energy)+((inst == 0xe)*((access_pos * (tmp / 2) + (1 - access_pos) * tmpptr->energy)));
				
				/* Keep track of execution frequencies for each instruction */
				statCounters.instructionExecutions[inst] += 1.0;
			}
            wordPtr=(wordPtr*((shiftPtr+4<SYSWORD_BITS)||(wordPtr+1<POND_DEPTH_SYSWORDS))+((shiftPtr+4>=SYSWORD_BITS)&&(wordPtr+1<POND_DEPTH_SYSWORDS))+EXEC_START_WORD*((wordPtr+1>=POND_DEPTH_SYSWORDS)&&(shiftPtr+4>=SYSWORD_BITS)))*!skip + wordPtr*skip;
            currentWord=(currentWord*(shiftPtr+4<SYSWORD_BITS)+(pptr->genome[wordPtr])*(shiftPtr+4>=SYSWORD_BITS))*!skip+ currentWord*skip;
            shiftPtr=((shiftPtr+4)+(shiftPtr+4>=SYSWORD_BITS)*(-shiftPtr-4))*!skip+shiftPtr*skip;
        }   
		if ((outputBuf[0] & 0xff) != 0xff) {
			tmpptr = getNeighbor(x,y,facing);

			//printf("%lu\n", tmpptr->energy);
			if ((tmpptr->energy)&&accessAllowed(tmpptr,reg,0,1)) {
				/* Log it if we're replacing a viable cell */
				if (tmpptr->generation > 2)
					++statCounters.viableCellsReplaced;
				
				tmpptr->ID = ++cellIdCounter;
				tmpptr->parentID = pptr->ID;
				tmpptr->lineage = pptr->lineage; /* Lineage is copied in offspring */
				tmpptr->generation = pptr->generation + 1;

				for(i=0;i<POND_DEPTH_SYSWORDS;++i)
					tmpptr->genome[i] = outputBuf[i];
			}
		}
	}

	return (void *)0;
}

/**
 * Main method
 *
 * @param argc Number of args
 * @param argv Argument array
 */
int main()
{
	uintptr_t outputBuf[POND_DEPTH_SYSWORDS];
	uintptr_t currentWord,wordPtr,shiftPtr,inst,tmp;
	struct Cell *pptr,*tmpptr;
	uintptr_t ptr_wordPtr;
	uintptr_t ptr_shiftPtr;
	uintptr_t reg;
	uintptr_t facing;
	uintptr_t loopStack_wordPtr[POND_DEPTH];
	uintptr_t loopStack_shiftPtr[POND_DEPTH];
	uintptr_t loopStackPtr;
	uintptr_t falseLoopDepth;
    uintptr_t clock;
    int skip;
    int access_neg_used;
    int access_pos_used;
    int access_neg;
    int access_pos; 
	int stop;

	/* Seed and init the random number generator */
	prngState[0] = 0; //(uint64_t)time(NULL);
	srand(13);
	prngState[1] = (uint64_t)rand();
	
	precalculate_random_numbers();

	/* Reset per-report stat counters */
	for(x=0;x<sizeof(statCounters);++x)
		((uint8_t *)&statCounters)[x] = (uint8_t)0;
 
	/* Clear the pond and initialize all genomes
	 * to 0xffff... */
	for(x=0;x<POND_SIZE_X;++x) {
		for(y=0;y<POND_SIZE_Y;++y) {
			pond[x][y].ID = 0;
			pond[x][y].parentID = 0;
			pond[x][y].lineage = 0;
			pond[x][y].generation = 0;
			pond[x][y].energy = 0;
			for(i=0;i<POND_DEPTH_SYSWORDS;++i)
				pond[x][y].genome[i] = ~((uintptr_t)0);
		}	
	}

    volatile int exitNow = 0;
    /* Main loop */
    clock = 0;
	while (!exitNow) {
		/* Increment clock and run reports periodically */
		/* Clock is incremented at the start, so it starts at 1 */
		++clock;
		if (clock == 1000000)
        {
            exitNow = 1;
        }
        if (!(clock % REPORT_FREQUENCY)) {
			doReport(clock);
		}
        /*
        static void *run(uintptr_t *ptr_wordPtr, uintptr_t *ptr_shiftPtr, uintptr_t *reg, 
        uintptr_t *loopStackPtr, uintptr_t *wordPtr, uintptr_t *shiftPtr, uintptr_t *facing, 
        uintptr_t *falseLoopDepth, int *stop, int *skip, int *access_neg_used, 
        int *access_pos_used, uintptr_t *outputBuf, uintptr_t *clock, struct Cell *pptr, struct Cell *tmpptr)
        */
	    run(ptr_wordPtr, ptr_shiftPtr, reg, loopStackPtr, wordPtr, shiftPtr, facing, falseLoopDepth, stop, skip, access_neg_used, access_pos_used, outputBuf, clock, pptr, tmpptr);
    }
	return 0;
    
}
