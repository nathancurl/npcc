/* *********************************************************************** */
/*                                                                         */
/* Nanopond version 2.0 -- A teeny tiny artificial life virtual machine    */
/* Copyright (C) Adam Ierymenko                                            */
/* MIT license -- see LICENSE.txt                                          */
/*                                                                         */
/* *********************************************************************** */

/*
 * Changelog:
 *
 * 1.0 - Initial release
 * 1.1 - Made empty cells get initialized with 0xffff... instead of zeros
 *       when the simulation starts. This makes things more consistent with
 *       the way the output buf is treated for self-replication, though
 *       the initial state rapidly becomes irrelevant as the simulation
 *       gets going.  Also made one or two very minor performance fixes.
 * 1.2 - Added statistics for execution frequency and metabolism, and made
 *       the visualization use 16bpp color.
 * 1.3 - Added a few other statistics.
 * 1.4 - Replaced SET with KILL and changed EAT to SHARE. The SHARE idea
 *       was contributed by Christoph Groth (http://www.falma.de/). KILL
 *       is a variation on the original EAT that is easier for cells to
 *       make use of.
 * 1.5 - Made some other instruction changes such as XCHG and added a
 *       penalty for failed KILL attempts. Also made access permissions
 *       stochastic.
 * 1.6 - Made cells all start facing in direction 0. This removes a bit
 *       of artificiality and requires cells to evolve the ability to
 *       turn in various directions in order to reproduce in anything but
 *       a straight line. It also makes pretty graphics.
 * 1.7 - Added more statistics, such as original lineage, and made the
 *       genome dump files CSV files as well.
 * 1.8 - Fixed LOOP/REP bug reported by user Sotek.  Thanks!  Also
 *       reduced the default mutation rate a bit.
 * 1.9 - Added a bunch of changes suggested by Christoph Groth: a better
 *       coloring algorithm, right click to switch coloring schemes (two
 *       are currently supported), and a few speed optimizations. Also
 *       changed visualization so that cells with generations less than 2
 *       are no longer shown.
 * 2.0 - Ported to SDL2 by Charles Huber, and added simple pthread based
 *       threading to make it take advantage of modern machines.
 */

/*
 * Nanopond is just what it says: a very very small and simple artificial
 * life virtual machine.
 *
 * It is a "small evolving program" based artificial life system of the same
 * general class as Tierra, Avida, and Archis.  It is written in very tight
 * and efficient C code to make it as fast as possible, and is so small that
 * it consists of only one .c file.
 *
 * How Nanopond works:
 *
 * The Nanopond world is called a "pond."  It is an NxN two dimensional
 * array of Cell structures, and it wraps at the edges (it's toroidal).
 * Each Cell structure consists of a few attributes that are there for
 * statistics purposes, an energy level, and an array of POND_DEPTH
 * four-bit values.  (The four-bit values are actually stored in an array
 * of machine-size words.)  The array in each cell contains the genome
 * associated with that cell, and POND_DEPTH is therefore the maximum
 * allowable size for a cell genome.
 *
 * The first four bit value in the genome is called the "logo." What that is
 * for will be explained later. The remaining four bit values each code for
 * one of 16 instructions. Instruction zero (0x0) is NOP (no operation) and
 * instruction 15 (0xf) is STOP (stop cell execution). Read the code to see
 * what the others are. The instructions are exceptionless and lack fragile
 * operands. This means that *any* arbitrary sequence of instructions will
 * always run and will always do *something*. This is called an evolvable
 * instruction set, because programs coded in an instruction set with these
 * basic characteristics can mutate. The instruction set is also
 * Turing-complete, which means that it can theoretically do anything any
 * computer can do. If you're curious, the instruciton set is based on this:
 * http://www.muppetlabs.com/~breadbox/bf/
 *
 * At the center of Nanopond is a core loop. Each time this loop executes,
 * a clock counter is incremented and one or more things happen:
 *
 * - Every REPORT_FREQUENCY clock ticks a line of comma seperated output
 *   is printed to STDOUT with some statistics about what's going on.
 * - Every INFLOW_FREQUENCY clock ticks a random x,y location is picked,
 *   energy is added (see INFLOW_RATE_MEAN and INFLOW_RATE_DEVIATION)
 *   and it's genome is filled with completely random bits.  Statistics
 *   are also reset to generation==0 and parentID==0 and a new cell ID
 *   is assigned.
 * - Every tick a random x,y location is picked and the genome inside is
 *   executed until a STOP instruction is encountered or the cell's
 *   energy counter reaches zero. (Each instruction costs one unit energy.)
 *
 * The cell virtual machine is an extremely simple register machine with
 * a single four bit register, one memory pointer, one spare memory pointer
 * that can be exchanged with the main one, and an output buffer. When
 * cell execution starts, this output buffer is filled with all binary 1's
 * (0xffff....). When cell execution is finished, if the first byte of
 * this buffer is *not* 0xff, then the VM says "hey, it must have some
 * data!". This data is a candidate offspring; to reproduce cells must
 * copy their genome data into the output buffer.
 *
 * When the VM sees data in the output buffer, it looks at the cell
 * adjacent to the cell that just executed and checks whether or not
 * the cell has permission (see below) to modify it. If so, then the
 * contents of the output buffer replace the genome data in the
 * adjacent cell. Statistics are also updated: parentID is set to the
 * ID of the cell that generated the output and generation is set to
 * one plus the generation of the parent.
 *
 * A cell is permitted to access a neighboring cell if:
 *    - That cell's energy is zero
 *    - That cell's parentID is zero
 *    - That cell's logo (remember?) matches the trying cell's "guess"
 *
 * Since randomly introduced cells have a parentID of zero, this allows
 * real living cells to always replace them or eat them.
 *
 * The "guess" is merely the value of the register at the time that the
 * access attempt occurs.
 *
 * Permissions determine whether or not an offspring can take the place
 * of the contents of a cell and also whether or not the cell is allowed
 * to EAT (an instruction) the energy in it's neighbor.
 *
 * If you haven't realized it yet, this is why the final permission
 * criteria is comparison against what is called a "guess." In conjunction
 * with the ability to "eat" neighbors' energy, guess what this permits?
 *
 * Since this is an evolving system, there have to be mutations. The
 * MUTATION_RATE sets their probability. Mutations are random variations
 * with a frequency defined by the mutation rate to the state of the
 * virtual machine while cell genomes are executing. Since cells have
 * to actually make copies of themselves to replicate, this means that
 * these copies can vary if mutations have occurred to the state of the
 * VM while copying was in progress.
 *
 * What results from this simple set of rules is an evolutionary game of
 * "corewar." In the beginning, the process of randomly generating cells
 * will cause self-replicating viable cells to spontaneously emerge. This
 * is something I call "random genesis," and happens when some of the
 * random gak turns out to be a program able to copy itself. After this,
 * evolution by natural selection takes over. Since natural selection is
 * most certainly *not* random, things will start to get more and more
 * ordered and complex (in the functional sense). There are two commodities
 * that are scarce in the pond: space in the NxN grid and energy. Evolving
 * cells compete for access to both.
 *
 * If you want more implementation details such as the actual instruction
 * set, read the source. It's well commented and is not that hard to
 * read. Most of it's complexity comes from the fact that four-bit values
 * are packed into machine size words by bit shifting. Once you get that,
 * the rest is pretty simple.
 *
 * Nanopond, for it's simplicity, manifests some really interesting
 * evolutionary dynamics. While I haven't run the kind of multiple-
 * month-long experiment necessary to really see this (I might!), it
 * would appear that evolution in the pond doesn't get "stuck" on just
 * one or a few forms the way some other simulators are apt to do.
 * I think simplicity is partly reponsible for this along with what
 * biologists call embeddedness, which means that the cells are a part
 * of their own world.
 *
 * Run it for a while... the results can be... interesting!
 *
 * Running Nanopond:
 *
 * Nanopond can use SDL (Simple Directmedia Layer) for screen output. If
 * you don't have SDL, comment out USE_SDL below and you'll just see text
 * statistics and get genome data dumps. (Turning off SDL will also speed
 * things up slightly.)
 *
 * After looking over the tunable parameters below, compile Nanopond and
 * run it. Here are some example compilation commands from Linux:
 *
 * For Pentiums:
 *  gcc -O6 -march=pentium -funroll-loops -fomit-frame-pointer -s
 *   -o nanopond nanopond.c -lSDL
 *
 * For Athlons with gcc 4.0+:
 *  gcc -O6 -msse -mmmx -march=athlon -mtune=athlon -ftree-vectorize
 *   -funroll-loops -fomit-frame-pointer -o nanopond nanopond.c -lSDL
 *
 * The second line is for gcc 4.0 or newer and makes use of GCC's new
 * tree vectorizing feature. This will speed things up a bit by
 * compiling a few of the loops into MMX/SSE instructions.
 *
 * This should also work on other Posix-compliant OSes with relatively
 * new C compilers. (Really old C compilers will probably not work.)
 * On other platforms, you're on your own! On Windows, you will probably
 * need to find and download SDL if you want pretty graphics and you
 * will need a compiler. MinGW and Borland's BCC32 are both free. I
 * would actually expect those to work better than Microsoft's compilers,
 * since MS tends to ignore C/C++ standards. If stdint.h isn't around,
 * you can fudge it like this:
 *
 * #define uintptr_t unsigned long (or whatever your machine size word is)
 * #define uint8_t unsigned char
 * #define uint16_t unsigned short
 * #define uint64_t unsigned long long (or whatever is your 64-bit int)
 *
 * When Nanopond runs, comma-seperated stats (see doReport() for
 * the columns) are output to stdout and various messages are output
 * to stderr. For example, you might do:
 *
 * ./nanopond >>stats.csv 2>messages.txt &
 *
 * To get both in seperate files.
 *
 * Have fun!
 */

/* ----------------------------------------------------------------------- */
/* Tunable parameters                                                      */
/* ----------------------------------------------------------------------- */

/* Frequency of comprehensive reports-- lower values will provide more
 * info while slowing down the simulation. Higher values will give less
 * frequent updates. */
/* This is also the frequency of screen refreshes if SDL is enabled. */
#define REPORT_FREQUENCY 200000

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
#include <cuda.h>

#define BUFFER_SIZE 1000  // Size of the circular buffer
__managed__ uintptr_t buffer[BUFFER_SIZE];
static int in = 0;
static uintptr_t last_random_number;

__managed__ volatile uint64_t prngState[2];

__device__ void getRandomPre(int rollback, uintptr_t &ret)
{
	// https://en.wikipedia.org/wiki/Xorshift#xorshift.2B
	uint64_t x = prngState[0];
	const uint64_t y = prngState[1];
	prngState[0] = prngState[0] * !rollback + rollback * y;
	x ^= x << 23;
	const uint64_t z = x ^ y ^ (x >> 17) ^ (y >> 26);
	prngState[1] = prngState[1] * !rollback + rollback * z;
	ret = (uintptr_t)(z + y);
}
__global__ void precalculate_random_numbers() {
    for (int i = 0; i < BUFFER_SIZE; i++) {
        buffer[i] = getRandomPre(1);
    }
}
__device__ uintptr_t getRandomRollback(uintptr_t rollback, int &ret) {
    uintptr_t num = buffer[in];
    last_random_number = num;  // Store the last random number
    uintptr_t new_num = getRandomPre(rollback);
    buffer[in] = (new_num & -rollback) | (num & ~-rollback);
    in = (((in + 1) & -rollback) | (in & ~-rollback)) % BUFFER_SIZE;
    ret = num;
}
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
    /* Access permission is more probable if they are more similar in sense 0,
	 * and more probable if they are different in sense 1. Sense 0 is used for
	 * "negative" interactions and sense 1 for "positive" ones. */
	//return sense ? (((getRandomRollback(1) & 0xf) >= BITS_IN_FOURBIT_WORD[(c2->genome[0] & 0xf) ^ (c1guess & 0xf)])||(!c2->parentID)) : (((getRandomRollback(1) & 0xf) <= BITS_IN_FOURBIT_WORD[(c2->genome[0] & 0xf) ^ (c1guess & 0xf)])||(!c2->parentID));
	return ((((random >= BITS_IN_FOURBIT_WORD[(c2->genome[0] & 0xf) ^ (c1guess & 0xf)]) || !c2->parentID) & sense) | (((random <= BITS_IN_FOURBIT_WORD[(c2->genome[0] & 0xf) ^ (c1guess & 0xf)]) || !c2->parentID) & ~sense));
}

volatile int exitNow = 0;

static void *run(void *targ)
{
	const uintptr_t threadNo = (uintptr_t)targ;
	uintptr_t x,y,i;
	uintptr_t clock = 0;

	/* Buffer used for execution output of candidate offspring */
	uintptr_t outputBuf[POND_DEPTH_SYSWORDS];

	/* Miscellaneous variables used in the loop */
	uintptr_t currentWord,wordPtr,shiftPtr,inst,tmp;
	struct Cell *pptr,*tmpptr;
	
	/* Virtual machine memory pointer register (which
	 * exists in two parts... read the code below...) */
	uintptr_t ptr_wordPtr;
	uintptr_t ptr_shiftPtr;
	
	/* The main "register" */
	uintptr_t reg;
	
	/* Which way is the cell facing? */
	uintptr_t facing;
	
	/* Virtual machine loop/rep stack */
	uintptr_t loopStack_wordPtr[POND_DEPTH];
	uintptr_t loopStack_shiftPtr[POND_DEPTH];
	uintptr_t loopStackPtr;
	
	/* If this is nonzero, we're skipping to matching REP */
	/* It is incremented to track the depth of a nested set
	 * of LOOP/REP pairs in false state. */
	uintptr_t falseLoopDepth;


	/* If this is nonzero, cell execution stops. This allows us
	 * to avoid the ugly use of a goto to exit the loop. :) */
	int stop;

	/* Main loop */
	while (!exitNow) {
		/* Increment clock and run reports periodically */
		/* Clock is incremented at the start, so it starts at 1 */
		++clock;
		if (clock == 1000000)
        {
            exitNow = 1;
        }
        if ((threadNo == 0)&&(!(clock % REPORT_FREQUENCY))) {
			doReport(clock);
		}

		/* Introduce a random cell somewhere with a given energy level */
		/* This is called seeding, and introduces both energy and
		 * entropy into the substrate. This happens every INFLOW_FREQUENCY
		 * clock ticks. */
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
        int skip=0;
		int access_neg_used = 0;
		int access_pos_used = 0;

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
				/*
				* ptr_shiftPtr
				*/
				ptr_shiftPtr = 
				(inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xd || inst == 0xe || inst == 0xf) * (ptr_shiftPtr) +
				((inst == 0x0)*0)+
				((inst == 0x1)*((ptr_shiftPtr+4)*((ptr_shiftPtr+4)<SYSWORD_BITS)))+
				((inst == 0x2)*(((ptr_shiftPtr==0)*SYSWORD_BITS)+ptr_shiftPtr-4));
				/*
				* ptr_wordPtr
				* set in 0x0, 0x1, 0x2
				*/				
				ptr_wordPtr =
				(inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xd || inst == 0xe || inst == 0xf) * (ptr_wordPtr) +
				((inst == 0x0)*0)+
				((inst == 0x1)*(((ptr_wordPtr*(ptr_shiftPtr!=0||((ptr_wordPtr+1)<POND_DEPTH_SYSWORDS))+(ptr_shiftPtr==0)*((ptr_wordPtr+1)<POND_DEPTH_SYSWORDS)))))+
				((inst == 0x2)*(((ptr_wordPtr==0&&ptr_shiftPtr==(SYSWORD_BITS-4))*(POND_DEPTH_SYSWORDS))+ptr_wordPtr-(ptr_shiftPtr==(SYSWORD_BITS-4))));
				
               /*
                * wordPtr
                * set in 0xc
                */ 
                wordPtr=
				(inst==0x0||inst==0x1||inst==0x2||inst==0x3||inst==0x4||inst==0x5|| inst == 0x6 || inst==0x7||inst==0x8||inst==0x9||inst==0xb||inst==0xd||inst==0xe||inst==0xf)*(wordPtr)+
                ((inst==0xa)*(wordPtr*!(reg&&loopStackPtr)+(loopStack_wordPtr[loopStackPtr-1]*(reg&&loopStackPtr))))+
                ((inst==0xc)*(wordPtr*((shiftPtr+4<SYSWORD_BITS)||(wordPtr+1<POND_DEPTH_SYSWORDS))+((shiftPtr+4>=SYSWORD_BITS)&&(wordPtr+1<POND_DEPTH_SYSWORDS))+EXEC_START_WORD*((wordPtr+1>=POND_DEPTH_SYSWORDS)&&(shiftPtr+4>=SYSWORD_BITS))));

               /*
                * shiftPtr
                * set in 0xc
                */ 
                shiftPtr=
				(inst==0x0||inst==0x1||inst==0x2||inst==0x3||inst==0x4||inst==0x5|| inst == 0x6 || inst==0x7||inst==0x8||inst==0x9||inst==0xb||inst==0xd||inst==0xe||inst==0xf)*(shiftPtr)+
                ((inst==0xa)*(shiftPtr*!(reg&&loopStackPtr)+(loopStack_shiftPtr[loopStackPtr-1]*(reg&&loopStackPtr))))+
                ((inst==0xc)*((shiftPtr+4)+(shiftPtr+4>=SYSWORD_BITS)*(-shiftPtr-4)));


                /*
                 * skip from 0xc 
                 */
                skip=(reg&&loopStackPtr)*(inst==0xa);
                
                
				/*facing is called in 0x0 and 0xb
				* facing is used to determine which direction the cell is facing
				*/
				facing=
				(inst==0x1||inst==0x2||inst==0x3||inst==0x4||inst==0x5||inst==0x6||inst==0x7||inst==0x8||inst==0x9||inst==0xa||inst==0xc||inst==0xd||inst==0xe||inst==0xf)*(facing) + 
				((inst==0x0)*0)+
				((inst==0xb)*(reg & 3));

				/* pptr->genome[ptr_wordPtr]*/
				pptr->genome[ptr_wordPtr]=
				(inst==0x0||inst==0x1||inst==0x2||inst==0x3||inst==0x4||inst==0x5||inst==0x7||inst==0x8||inst==0x9||inst==0xa||inst==0xb||inst==0xc||inst==0xd||inst==0xe||inst==0xf)
				*(pptr->genome[ptr_wordPtr])+((inst==0x6)*((pptr->genome[ptr_wordPtr]&~(((uintptr_t)0xf)<<ptr_shiftPtr))|reg<<ptr_shiftPtr)); 
				/*
				* wordPtr
				* set in 0xa
				*/
				/*
				* outputBuf[ptr_wordPtr]
				* set in 0x8
				*/
				outputBuf[ptr_wordPtr]=
				(inst==0x0||inst==0x1||inst==0x2||inst==0x3||inst==0x4||inst==0x5|| inst == 0x6 || inst==0x7||inst==0x9||inst==0xa||inst==0xb||inst==0xc||inst==0xd||inst==0xe||inst==0xf)*(outputBuf[ptr_wordPtr])+
				((inst==0x8)*((outputBuf[ptr_wordPtr]&~(((uintptr_t)0xf) << ptr_shiftPtr))|reg << ptr_shiftPtr));
				
                
                
                /*
				* currentWord
				* set in 0x6, 0xa, 0xc
				* TODO: 0xa, 0xc
				*/
				currentWord=
				(inst==0x0||inst==0x1||inst==0x2||inst==0x3||inst==0x4||inst==0x5||inst==0x7||inst==0x8||inst==0x9|| inst==0xb || inst==0xd||inst==0xe||inst==0xf)*(currentWord)+
				((inst==0x6)*(pptr->genome[wordPtr]))+
				((inst==0xa)*(currentWord*!(reg&&loopStackPtr)+(pptr->genome[wordPtr])*(reg&&loopStackPtr)))+
				((inst == 0xc)*(pptr->genome[wordPtr]));
				
                /*
				* stop
				* set in 0x9 & 0xf
				*/
				stop=
				(inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xd || inst == 0xe)*(stop)+
				((inst == 0x9)*(stop*!(reg&&(loopStackPtr>=POND_DEPTH))+(reg&&(loopStackPtr>=POND_DEPTH))))+
				((inst == 0xf)*(1));
				/*
				* loopStack_wordPtr[loopStackPtr]
				* loopStack_wordPtr[loopStackPtr] set in 0x9
				*/
				loopStack_wordPtr[loopStackPtr]=
				(inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xd || inst == 0xe || inst == 0xf)*(loopStack_wordPtr[loopStackPtr])+
				((inst == 0x9)*(loopStack_wordPtr[loopStackPtr]*(!reg||(loopStackPtr>=POND_DEPTH))+(wordPtr*(reg&&(loopStackPtr<POND_DEPTH)))));
				/*
				* loopStack_shiftPtr[loopStackPtr]
				* set in 0x9
				*/
				loopStack_shiftPtr[loopStackPtr]=
				(inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xd || inst == 0xe || inst == 0xf)*(loopStack_shiftPtr[loopStackPtr])+
				((inst == 0x9) * (loopStack_shiftPtr[loopStackPtr]*(!reg||(loopStackPtr>=POND_DEPTH))+(shiftPtr*(reg&&(loopStackPtr<POND_DEPTH)))));
				/*
				* loopStackPtr
				* set in 0x9
				*/
				loopStackPtr=
				(inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0xb || inst == 0xc || inst == 0xd || inst == 0xe || inst == 0xf)*(loopStackPtr)+
				((inst == 0x9)*(loopStackPtr + (reg&&(loopStackPtr<POND_DEPTH))))+
                ((inst == 0xa)*(loopStackPtr-!!loopStackPtr));
				/*
				* falseLoopDepth
				* set in 0x9,
				*/
				falseLoopDepth=
				(inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xd || inst == 0xe || inst == 0xf)*(falseLoopDepth)+
				((inst == 0x9)*(falseLoopDepth + (!reg)));
				/*
				* tmpptr
				* set is 0xd, 0xe
				*/
				tmpptr = getNeighbor(x,y,facing);
				access_neg_used = 0;
				access_pos_used = 0;
				
				access_pos_used =
				(inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xd || inst == 0xf)*(access_pos_used)+
				((inst == 0xe)*(1));

				access_neg_used =
				(inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xe|| inst == 0xf)*(access_neg_used)+
				((inst == 0xd)*(1));
				
				int access_neg = accessAllowed(tmpptr,reg,0, access_neg_used);
				int access_pos = accessAllowed(tmpptr,reg,1, access_pos_used);

				statCounters.viableCellsKilled=
				(inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xf)*(statCounters.viableCellsKilled)+
				((inst == 0xd)*(statCounters.viableCellsKilled+(access_neg)*(tmpptr->generation>2)))+
				((inst == 0xe)*(statCounters.viableCellsKilled+(access_pos)*(tmpptr->generation>2)));

				tmpptr->genome[0]=
				(inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xe || inst == 0xf)*(tmpptr->genome[0])+
				((inst == 0xd)*(tmpptr->genome[0]*!(access_neg)+(access_neg)*~((uintptr_t)0)));

				tmpptr->genome[1]=
				(inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xe || inst == 0xf)*(tmpptr->genome[1])+
				((inst == 0xd)*(tmpptr->genome[0]*!(access_neg)+(access_neg)*~((uintptr_t)0)));

				tmpptr->ID=
				(inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xe || inst == 0xf)*(tmpptr->ID)+
				((inst == 0xd)*(tmpptr->ID * !(access_neg)+ (access_neg)*cellIdCounter));

				tmpptr->parentID=
				(inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xe || inst == 0xf)*(tmpptr->parentID)+
				((inst == 0xd)*(tmpptr->parentID * !(access_neg)));

				tmpptr->lineage=
				(inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xe || inst == 0xf)*(tmpptr->lineage)+
				((inst == 0xd)*(tmpptr->lineage * !(access_neg) + (access_neg)*cellIdCounter));

				cellIdCounter=
				(inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xe || inst == 0xf)*(cellIdCounter)+
				((inst == 0xd)*(cellIdCounter * !(access_neg) + (access_neg)* cellIdCounter));

				tmp =
				(inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xf)*(tmp)+
				((inst == 0xc)*(reg))+
				((inst == 0xd)*((access_neg) + (tmpptr->generation>2)*!(access_neg)*(pptr->energy / FAILED_KILL_PENALTY)))+
				((inst == 0xe)* (pptr->energy + tmpptr->energy));

				reg=
				(inst == 0x1 || inst == 0x2 || inst == 0x6 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xd ||inst == 0xe || inst == 0xf) * (reg) + 
				((inst==0x0)*0) + 
				((inst==0x3)*((reg + 1) & 0xf)) +
				((inst==0x4)*((reg - 1) & 0xf)) +
				((inst==0x5)*((pptr->genome[ptr_wordPtr] >> ptr_shiftPtr) & 0xf)) +
				((inst==0x7)*((outputBuf[ptr_wordPtr] >> ptr_shiftPtr) & 0xf))+
				((inst==0xc)*((pptr->genome[wordPtr] >> shiftPtr) & 0xf));

				pptr->genome[wordPtr]=
				(inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xd || inst == 0xe || inst == 0xf)*(pptr->genome[wordPtr])+
				((inst == 0xc)* (((pptr->genome[wordPtr]&~(((uintptr_t)0xf) << shiftPtr))|tmp << shiftPtr)));

				pptr->energy=
				(inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xf)*(pptr->energy)+
				((inst == 0xd)*(pptr->energy+!(access_neg)*(tmpptr->generation>2)*(-pptr->energy) + !(access_neg)*(tmpptr->generation>2)*(pptr->energy-tmp)))+
				((inst == 0xe)*((access_pos * (tmp - (access_pos * (tmp / 2) + (1 - access_pos) * tmpptr->energy)) + (1 - access_pos) * pptr->energy)));

				tmpptr->generation=
				(inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xe || inst == 0xf)*(tmpptr->generation)+
				((inst == 0xd)*(tmpptr->generation * (access_neg)));

				tmpptr->energy=
				(inst == 0x0 || inst == 0x1 || inst == 0x2 || inst == 0x3 || inst == 0x4 || inst == 0x5 || inst == 0x6 || inst == 0x7 || inst == 0x8 || inst == 0x9 || inst == 0xa || inst == 0xb || inst == 0xc || inst == 0xd || inst == 0xf)*(tmpptr->energy)+
				((inst == 0xe)*((access_pos * (tmp / 2) + (1 - access_pos) * tmpptr->energy)));
				
				/* Keep track of execution frequencies for each instruction */
				statCounters.instructionExecutions[inst] += 1.0;

			}
			
			/* Advance the shift and word pointers, and loop around
			 * to the beginning at the end of the genome. */

            // increment wordptr by 1 if the shift Ptr is going to go
            // beyond the current word it is reading.
            // Set the wordptr to EXEC_START_WORD if the end of the
            // POND_DEPTH_SYSWORDS has been reached.
            wordPtr=(wordPtr*((shiftPtr+4<SYSWORD_BITS)||(wordPtr+1<POND_DEPTH_SYSWORDS))
                +
                ((shiftPtr+4>=SYSWORD_BITS)&&(wordPtr+1<POND_DEPTH_SYSWORDS))
                +
                EXEC_START_WORD*((wordPtr+1>=POND_DEPTH_SYSWORDS)&&(shiftPtr+4>=SYSWORD_BITS)))*!skip + wordPtr*skip;

            //currentWord gets incremented when the shiftptr is greater than
            //SYSWORD_BITS, and it's time to move to the next word
            currentWord=(currentWord*(shiftPtr+4<SYSWORD_BITS)
                +
                (pptr->genome[wordPtr])*(shiftPtr+4>=SYSWORD_BITS))*!skip
                + currentWord*skip;
            
            //shiftPtr shifts the current nibble being read by the machine
            //It incrememnts four bits until it gets past SYSWORD_BITS, the 
            //number of bits in a word, and then resets at either 0 or
            //EXEC_START_BIT
            shiftPtr=((shiftPtr+4)
                    +
                    (shiftPtr+4>=SYSWORD_BITS)*(-shiftPtr-4))*!skip+shiftPtr*skip;
            //+
                //(EXEC_START_BIT)*(wordPtr+1>=POND_DEPTH_SYSWORDS)
                //*(shiftPtr+4>=SYSWORD_BITS);
            /*
            if ((shiftPtr += 4) >= SYSWORD_BITS) {
				if (++wordPtr >= POND_DEPTH_SYSWORDS) {
					wordPtr = EXEC_START_WORD;
					shiftPtr = EXEC_START_BIT;
				} else shiftPtr = 0;
				currentWord = pptr->genome[wordPtr];
			}
            */
        }   

		/* Copy outputBuf into neighbor if access is permitted and there
		 * is energy there to make something happen. There is no need
		 * to copy to a cell with no energy, since anything copied there
		 * would never be executed and then would be replaced with random
		 * junk eventually. See the seeding code in the main loop above. */
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

	uintptr_t i,x,y;

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
	run((void *)0);

	return 0;
}
