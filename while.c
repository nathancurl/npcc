#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>



/* Core execution loop */
__global__ void whileLoop(){
volatile int exitNow = 0;
uintptr_t x,y,i;

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

		/* If this is nonzero, we're skipping to matching REP */
/* It is incremented to track the depth of a nested set
 * of LOOP/REP pairs in false state. */
uintptr_t falseLoopDepth;

/* If this is nonzero, cell execution stops. This allows us
 * to avoid the ugly use of a goto to exit the loop. :) */
int stop;

	while ((pptr->energy)&&(!stop)) {
		/* Get the next instruction */
		inst = (currentWord >> shiftPtr) & 0xf;

		/* Randomly frob either the instruction or the register with a
		 * probability defined by MUTATION_RATE. This introduces variation,
		 * and since the variation is introduced into the state of the VM
		 * it can have all manner of different effects on the end result of
		 * replication: insertions, deletions, duplications of entire
		 * ranges of the genome, etc. */
		if ((getRandom() & 0xffffffff) < MUTATION_RATE) {
			tmp = getRandom(); /* Call getRandom() only once for speed */
			if (tmp & 0x80) /* Check for the 8th bit to get random boolean */
				inst = tmp & 0xf; /* Only the first four bits are used here */
			else reg = tmp & 0xf;
		}

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
			/* If we're not in a false LOOP/REP, execute normally */
			
			/* Keep track of execution frequencies for each instruction */
			statCounters.instructionExecutions[inst] += 1.0;
			
			switch(inst) {
				case 0x0: /* ZERO: Zero VM state registers */
					reg = 0;
					ptr_wordPtr = 0;
					ptr_shiftPtr = 0;
					facing = 0;
					break;
				case 0x1: /* FWD: Increment the pointer (wrap at end) */
					if ((ptr_shiftPtr += 4) >= SYSWORD_BITS) {
						if (++ptr_wordPtr >= POND_DEPTH_SYSWORDS)
							ptr_wordPtr = 0;
						ptr_shiftPtr = 0;
					}
					break;
				case 0x2: /* BACK: Decrement the pointer (wrap at beginning) */
					if (ptr_shiftPtr)
						ptr_shiftPtr -= 4;
					else {
						if (ptr_wordPtr)
							--ptr_wordPtr;
						else ptr_wordPtr = POND_DEPTH_SYSWORDS - 1;
						ptr_shiftPtr = SYSWORD_BITS - 4;
					}
					break;
				case 0x3: /* INC: Increment the register */
					reg = (reg + 1) & 0xf;
					break;
				case 0x4: /* DEC: Decrement the register */
					reg = (reg - 1) & 0xf;
					break;
				case 0x5: /* READG: Read into the register from genome */
					reg = (pptr->genome[ptr_wordPtr] >> ptr_shiftPtr) & 0xf;
					break;
				case 0x6: /* WRITEG: Write out from the register to genome */
					pptr->genome[ptr_wordPtr] &= ~(((uintptr_t)0xf) << ptr_shiftPtr);
					pptr->genome[ptr_wordPtr] |= reg << ptr_shiftPtr;
					currentWord = pptr->genome[wordPtr]; /* Must refresh in case this changed! */
					break;
				case 0x7: /* READB: Read into the register from buffer */
					reg = (outputBuf[ptr_wordPtr] >> ptr_shiftPtr) & 0xf;
					break;
				case 0x8: /* WRITEB: Write out from the register to buffer */
					outputBuf[ptr_wordPtr] &= ~(((uintptr_t)0xf) << ptr_shiftPtr);
					outputBuf[ptr_wordPtr] |= reg << ptr_shiftPtr;
					break;
				case 0x9: /* LOOP: Jump forward to matching REP if register is zero */
					if (reg) {
						if (loopStackPtr >= POND_DEPTH)
							stop = 1; /* Stack overflow ends execution */
						else {
							loopStack_wordPtr[loopStackPtr] = wordPtr;
							loopStack_shiftPtr[loopStackPtr] = shiftPtr;
							++loopStackPtr;
						}
					} else falseLoopDepth = 1;
					break;
				case 0xa: /* REP: Jump back to matching LOOP if register is nonzero */
					if (loopStackPtr) {
						--loopStackPtr;
						if (reg) {
							wordPtr = loopStack_wordPtr[loopStackPtr];
							shiftPtr = loopStack_shiftPtr[loopStackPtr];
							currentWord = pptr->genome[wordPtr];
							/* This ensures that the LOOP is rerun */
							continue;
						}
					}
					break;
				case 0xb: /* TURN: Turn in the direction specified by register */
					facing = reg & 3;
					break;
				case 0xc: /* XCHG: Skip next instruction and exchange value of register with it */
					if ((shiftPtr += 4) >= SYSWORD_BITS) {
						if (++wordPtr >= POND_DEPTH_SYSWORDS) {
							wordPtr = EXEC_START_WORD;
							shiftPtr = EXEC_START_BIT;
						} else shiftPtr = 0;
					}
					tmp = reg;
					reg = (pptr->genome[wordPtr] >> shiftPtr) & 0xf;
					pptr->genome[wordPtr] &= ~(((uintptr_t)0xf) << shiftPtr);
					pptr->genome[wordPtr] |= tmp << shiftPtr;
					currentWord = pptr->genome[wordPtr];
					break;
				case 0xd: /* KILL: Blow away neighboring cell if allowed with penalty on failure */
					tmpptr = getNeighbor(x,y,facing);
					if (accessAllowed(tmpptr,reg,0)) {
						if (tmpptr->generation > 2)
							++statCounters.viableCellsKilled;

						/* Filling first two words with 0xfffff... is enough */
						tmpptr->genome[0] = ~((uintptr_t)0);
						tmpptr->genome[1] = ~((uintptr_t)0);
						tmpptr->ID = cellIdCounter;
						tmpptr->parentID = 0;
						tmpptr->lineage = cellIdCounter;
						tmpptr->generation = 0;
						++cellIdCounter;
					} else if (tmpptr->generation > 2) {
						tmp = pptr->energy / FAILED_KILL_PENALTY;
						if (pptr->energy > tmp)
							pptr->energy -= tmp;
						else pptr->energy = 0;
					}
					break;
				case 0xe: /* SHARE: Equalize energy between self and neighbor if allowed */
					tmpptr = getNeighbor(x,y,facing);
					if (accessAllowed(tmpptr,reg,1)) {
#ifdef USE_PTHREADS_COUNT
						pthread_mutex_lock(&(tmpptr->lock));
#endif
						if (tmpptr->generation > 2)
							++statCounters.viableCellShares;
						tmp = pptr->energy + tmpptr->energy;
						tmpptr->energy = tmp / 2;
						pptr->energy = tmp - tmpptr->energy;
#ifdef USE_PTHREADS_COUNT
						pthread_mutex_unlock(&(tmpptr->lock));
#endif
					}
					break;
				case 0xf: /* STOP: End execution */
					stop = 1;
					break;
			}
		}
		
		/* Advance the shift and word pointers, and loop around
		 * to the beginning at the end of the genome. */
		if ((shiftPtr += 4) >= SYSWORD_BITS) {
			if (++wordPtr >= POND_DEPTH_SYSWORDS) {
				wordPtr = EXEC_START_WORD;
				shiftPtr = EXEC_START_BIT;
			} else shiftPtr = 0;
			currentWord = pptr->genome[wordPtr];
		}
	}
}

int run{
	whileLoop<<<1,1>>>();
	return 0;
}
