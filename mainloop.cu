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
	/* other variables */
	int skip;
	int access_neg_used;
	int access_pos_used;
	int access_neg ;
	int access_pos;
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
			* uintptr_t mutation_occurred = (getRandomRollback(1) & 0xffffffff) < MUTATION_RATE;
			* uintptr_t tmp = getRandomRollback(mutation_occurred) * mutation_occurred;
			* uintptr_t is_inst = (tmp & 0x80) >> 7; // Shift right by 7 to get a 1 or 0
			* uintptr_t is_reg = ~is_inst & 0x1; // Invert is_inst and mask with 0x1 to get a 1 or 0
			* inst = (tmp & 0xf) * is_inst + inst * (!is_inst); // Update inst only if is_inst is 1
			* reg = (tmp & 0xf) * is_reg + reg * (!is_reg); // Update reg only if is_reg is 1
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
			
            wordPtr=(wordPtr*((shiftPtr+4<SYSWORD_BITS)||(wordPtr+1<POND_DEPTH_SYSWORDS)) + ((shiftPtr+4>=SYSWORD_BITS)&&(wordPtr+1<POND_DEPTH_SYSWORDS)) + EXEC_START_WORD*((wordPtr+1>=POND_DEPTH_SYSWORDS)&&(shiftPtr+4>=SYSWORD_BITS)))*!skip + wordPtr*skip;

            currentWord=(currentWord*(shiftPtr+4<SYSWORD_BITS)+(pptr->genome[wordPtr])*(shiftPtr+4>=SYSWORD_BITS))*!skip+ currentWord*skip;

            shiftPtr=((shiftPtr+4)+(shiftPtr+4>=SYSWORD_BITS)*(-shiftPtr-4))*!skip+shiftPtr*skip;

        }   
3
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
__global__ void executionLoop(
    uintptr_t *d_threadNo, uintptr_t *d_x, uintptr_t *d_y, uintptr_t *d_i, uintptr_t *d_clock, 
    uintptr_t *d_outputBuf, uintptr_t *d_currentWord, uintptr_t *d_wordPtr, uintptr_t *d_shiftPtr, 
    uintptr_t *d_inst, uintptr_t *d_tmp, uintptr_t *d_ptr_wordPtr, uintptr_t *d_ptr_shiftPtr, 
    uintptr_t *d_reg, uintptr_t *d_facing, uintptr_t *d_loopStack_wordPtr, uintptr_t *d_loopStack_shiftPtr, 
    uintptr_t *d_loopStackPtr, uintptr_t *d_falseLoopDepth, int *d_stop, int *d_skip, 
    int *d_access_neg_used, int *d_access_pos_used, int *d_access_neg, int *d_access_pos, 
    struct Cell *d_pptr, struct Cell *d_tmpptr) 
    {
        
        uintptr_t threadNo = *d_threadNo;
        uintptr_t x = *d_x;
        uintptr_t y = *d_y;
        uintptr_t i = *d_i;
        uintptr_t clock = *d_clock;
        uintptr_t outputBuf[POND_DEPTH_SYSWORDS];
        memcpy(outputBuf, d_outputBuf, POND_DEPTH_SYSWORDS * sizeof(uintptr_t));
        uintptr_t currentWord = *d_currentWord;
        uintptr_t wordPtr = *d_wordPtr;
        uintptr_t shiftPtr = *d_shiftPtr;
        uintptr_t inst = *d_inst;
        uintptr_t tmp = *d_tmp;
        uintptr_t ptr_wordPtr = *d_ptr_wordPtr;
        uintptr_t ptr_shiftPtr = *d_ptr_shiftPtr;
        uintptr_t reg = *d_reg;
        uintptr_t facing = *d_facing;
        uintptr_t loopStack_wordPtr[POND_DEPTH];
        memcpy(loopStack_wordPtr, d_loopStack_wordPtr, POND_DEPTH * sizeof(uintptr_t));
        uintptr_t loopStack_shiftPtr[POND_DEPTH];
        memcpy(loopStack_shiftPtr, d_loopStack_shiftPtr, POND_DEPTH * sizeof(uintptr_t));
        uintptr_t loopStackPtr = *d_loopStackPtr;
        uintptr_t falseLoopDepth = *d_falseLoopDepth;
        int stop = *d_stop;
        int skip = *d_skip;
        int access_neg_used = *d_access_neg_used;
        int access_pos_used = *d_access_pos_used;
        int access_neg = *d_access_neg;
        int access_pos = *d_access_pos;
        struct Cell pptr = *d_pptr;
        struct Cell tmpptr = *d_tmpptr;


}
    


int cudaMain(){

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
	/* other variables */
	int skip;
	int access_neg_used;
	int access_pos_used;
	int access_neg ;
	int access_pos;

	// Declare pointers for the variables
	uintptr_t *d_threadNo, *d_x, *d_y, *d_i, *d_clock, *d_outputBuf, *d_currentWord, *d_wordPtr, *d_shiftPtr, *d_inst, *d_tmp;
	uintptr_t *d_ptr_wordPtr, *d_ptr_shiftPtr, *d_reg, *d_facing, *d_loopStack_wordPtr, *d_loopStack_shiftPtr, *d_loopStackPtr, *d_falseLoopDepth;
	int *d_stop, *d_skip, *d_access_neg_used, *d_access_pos_used, *d_access_neg, *d_access_pos;
	struct Cell *d_pptr, *d_tmpptr;

	// Allocate memory on the GPU
	cudaMalloc(&d_threadNo, sizeof(uintptr_t));
	cudaMalloc(&d_x, sizeof(uintptr_t));
	cudaMalloc(&d_y, sizeof(uintptr_t));
	cudaMalloc(&d_i, sizeof(uintptr_t));
	cudaMalloc(&d_clock, sizeof(uintptr_t));
	cudaMalloc(&d_outputBuf, POND_DEPTH_SYSWORDS * sizeof(uintptr_t));
	cudaMalloc(&d_currentWord, sizeof(uintptr_t));
	cudaMalloc(&d_wordPtr, sizeof(uintptr_t));
	cudaMalloc(&d_shiftPtr, sizeof(uintptr_t));
	cudaMalloc(&d_inst, sizeof(uintptr_t));
	cudaMalloc(&d_tmp, sizeof(uintptr_t));
	cudaMalloc(&d_ptr_wordPtr, sizeof(uintptr_t));
	cudaMalloc(&d_ptr_shiftPtr, sizeof(uintptr_t));
	cudaMalloc(&d_reg, sizeof(uintptr_t));
	cudaMalloc(&d_facing, sizeof(uintptr_t));
	cudaMalloc(&d_loopStack_wordPtr, POND_DEPTH * sizeof(uintptr_t));
	cudaMalloc(&d_loopStack_shiftPtr, POND_DEPTH * sizeof(uintptr_t));
	cudaMalloc(&d_loopStackPtr, sizeof(uintptr_t));
	cudaMalloc(&d_falseLoopDepth, sizeof(uintptr_t));
	cudaMalloc(&d_stop, sizeof(int));
	cudaMalloc(&d_skip, sizeof(int));
	cudaMalloc(&d_access_neg_used, sizeof(int));
	cudaMalloc(&d_access_pos_used, sizeof(int));
	cudaMalloc(&d_access_neg, sizeof(int));
	cudaMalloc(&d_access_pos, sizeof(int));
	cudaMalloc(&d_pptr, sizeof(struct Cell));
	cudaMalloc(&d_tmpptr, sizeof(struct Cell));

	// Call the kernel function
	run<<<1, 1>>>(d_outputBuf, d_loopStack_wordPtr, d_loopStack_shiftPtr, d_stop, d_skip, d_access_neg_used, d_access_pos_used);

	// Copy data from host to device
	cudaMemcpy(d_threadNo, &threadNo, sizeof(uintptr_t), cudaMemcpyHostToDevice);
	cudaMemcpy(d_x, &x, sizeof(uintptr_t), cudaMemcpyHostToDevice);
	cudaMemcpy(d_y, &y, sizeof(uintptr_t), cudaMemcpyHostToDevice);
	cudaMemcpy(d_i, &i, sizeof(uintptr_t), cudaMemcpyHostToDevice);
	cudaMemcpy(d_clock, &clock, sizeof(uintptr_t), cudaMemcpyHostToDevice);
	cudaMemcpy(d_outputBuf, outputBuf, POND_DEPTH_SYSWORDS * sizeof(uintptr_t), cudaMemcpyHostToDevice);
	cudaMemcpy(d_currentWord, &currentWord, sizeof(uintptr_t), cudaMemcpyHostToDevice);
	cudaMemcpy(d_wordPtr, &wordPtr, sizeof(uintptr_t), cudaMemcpyHostToDevice);
	cudaMemcpy(d_shiftPtr, &shiftPtr, sizeof(uintptr_t), cudaMemcpyHostToDevice);
	cudaMemcpy(d_inst, &inst, sizeof(uintptr_t), cudaMemcpyHostToDevice);
	cudaMemcpy(d_tmp, &tmp, sizeof(uintptr_t), cudaMemcpyHostToDevice);
	cudaMemcpy(d_ptr_wordPtr, &ptr_wordPtr, sizeof(uintptr_t), cudaMemcpyHostToDevice);
	cudaMemcpy(d_ptr_shiftPtr, &ptr_shiftPtr, sizeof(uintptr_t), cudaMemcpyHostToDevice);
	cudaMemcpy(d_reg, &reg, sizeof(uintptr_t), cudaMemcpyHostToDevice);
	cudaMemcpy(d_facing, &facing, sizeof(uintptr_t), cudaMemcpyHostToDevice);
	cudaMemcpy(d_loopStack_wordPtr, loopStack_wordPtr, POND_DEPTH * sizeof(uintptr_t), cudaMemcpyHostToDevice);
	cudaMemcpy(d_loopStack_shiftPtr, loopStack_shiftPtr, POND_DEPTH * sizeof(uintptr_t), cudaMemcpyHostToDevice);
	cudaMemcpy(d_loopStackPtr, &loopStackPtr, sizeof(uintptr_t), cudaMemcpyHostToDevice);
	cudaMemcpy(d_falseLoopDepth, &falseLoopDepth, sizeof(uintptr_t), cudaMemcpyHostToDevice);
	cudaMemcpy(d_stop, &stop, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_skip, &skip, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_access_neg_used, &access_neg_used, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_access_pos_used, &access_pos_used, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_access_neg, &access_neg, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_access_pos, &access_pos, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_pptr, pptr, sizeof(struct Cell), cudaMemcpyHostToDevice);
	cudaMemcpy(d_tmpptr, tmpptr, sizeof(struct Cell), cudaMemcpyHostToDevice);

	// Free memory on the GPU
	cudaFree(d_outputBuf);
	cudaFree(d_loopStack_wordPtr);
	cudaFree(d_loopStack_shiftPtr);
	cudaFree(d_stop);
	cudaFree(d_skip);
	cudaFree(d_access_neg_used);
	cudaFree(d_access_pos_used);
}