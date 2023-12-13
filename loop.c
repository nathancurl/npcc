


switch(inst) {
					case 0x0: /* ZERO: Zero VM state registers */
						reg = 0;
						ptr_wordPtr = 0;
						ptr_shiftPtr = 0;
						facing = 0;
						break;
					case 0x1: /* FWD: Increment the pointer (wrap at end) */
                        ptr_shiftPtr=(ptr_shiftPtr+4)*((ptr_shiftPtr+4)<SYSWORD_BITS);
                        ptr_wordPtr=(ptr_wordPtr*(ptr_shiftPtr!=0||((ptr_wordPtr+1)<POND_DEPTH_SYSWORDS))+(ptr_shiftPtr==0)*((ptr_wordPtr+1)<POND_DEPTH_SYSWORDS));
                        break;
					case 0x2: /* BACK: Decrement the pointer (wrap at beginning) */ 
                        ptr_shiftPtr=((ptr_shiftPtr==0)*SYSWORD_BITS)+ptr_shiftPtr-4;
                        ptr_wordPtr=((ptr_wordPtr==0&&ptr_shiftPtr==(SYSWORD_BITS-4))*(POND_DEPTH_SYSWORDS))+ptr_wordPtr-(ptr_shiftPtr==(SYSWORD_BITS-4));
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
                        outputBuf[ptr_wordPtr] = (outputBuf[ptr_wordPtr] & ~(((uintptr_t)0xf) << ptr_shiftPtr)) | (reg << ptr_shiftPtr);
                        break;
					case 0x9: /* LOOP: Jump forward to matching REP if register is zero */
                        stop=stop*!(reg&&(loopStackPtr>=POND_DEPTH))+(reg&&(loopStackPtr>=POND_DEPTH));
                        loopStack_wordPtr[loopStackPtr]=loopStack_wordPtr[loopStackPtr]*(!reg||(loopStackPtr>=POND_DEPTH))+(wordPtr*(reg&&(loopStackPtr<POND_DEPTH)));
                        loopStack_shiftPtr[loopStackPtr]=loopStack_shiftPtr[loopStackPtr]*(!reg||(loopStackPtr>=POND_DEPTH))+(shiftPtr*(reg&&(loopStackPtr<POND_DEPTH)));
                        loopStackPtr = loopStackPtr + (reg&&(loopStackPtr<POND_DEPTH));
                        falseLoopDepth = !reg; 
                        break;
					case 0xa: /* REP: Jump back to matching LOOP if register is nonzero */
                                 if (loopStackPtr) {
							--loopStackPtr;
							if (reg) {
								wordPtr = loopStack_wordPtr[loopStackPtr];
								shiftPtr = loopStack_shiftPtr[loopStackPtr];
								currentWord = pptr->genome[wordPtr];
								//This ensures that the LOOP is rerun 
								continue;
							}
						}
						break;
					case 0xb: /* TURN: Turn in the direction specified by register */
						facing = reg & 3;
						break;
					case 0xc: /* XCHG: Skip next instruction and exchange value of register with it */
                        wordPtr=wordPtr*((shiftPtr+4<SYSWORD_BITS)||(wordPtr+1<POND_DEPTH_SYSWORDS))+((shiftPtr+4>=SYSWORD_BITS)&&(wordPtr+1<POND_DEPTH_SYSWORDS))+EXEC_START_WORD*((wordPtr+1>=POND_DEPTH_SYSWORDS)&&(shiftPtr+4>=SYSWORD_BITS));
                        shiftPtr=(shiftPtr+4)+(shiftPtr+4>=SYSWORD_BITS)*(-shiftPtr-4);
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
						int access = accessAllowed(tmpptr,reg,1);
						int generationCondition = tmpptr->generation > 2;
						tmp = pptr->energy + tmpptr->energy;
						statCounters.viableCellShares += access * generationCondition;
						uintptr_t newEnergyNeighbor = access * (tmp / 2) + (1 - access) * tmpptr->energy;
						uintptr_t newEnergySelf = access * (tmp - newEnergyNeighbor) + (1 - access) * pptr->energy;
						tmpptr->energy = newEnergyNeighbor;
						pptr->energy = newEnergySelf;
						break; 
					case 0xf: /* STOP: End execution */
						stop = 1;
						break;
				}