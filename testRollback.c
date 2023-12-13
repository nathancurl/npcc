#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


#define BUFFER_SIZE 1000  // Size of the circular buffer
static uintptr_t buffer[BUFFER_SIZE];
static int in = 0;
static int last_random_number = 0;


volatile uint64_t prngState[2];
volatile uint64_t prngStateOG[2];

static inline uintptr_t getRandomPreOG()
{
	// https://en.wikipedia.org/wiki/Xorshift#xorshift.2B
	uint64_t x = prngStateOG[0];
	const uint64_t y = prngStateOG[1];
	prngStateOG[0] = y;
	x ^= x << 23;
	const uint64_t z = x ^ y ^ (x >> 17) ^ (y >> 26);
	prngStateOG[1] = z;
	return (uintptr_t)(z + y);
}
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
/*
static inline uintptr_t getRandom() {
    uintptr_t num = buffer[in];
    last_random_number = num;  // Store the last random number
	buffer[in] = getRandomPre();  // Generate a new random number and add it to the buffer
    in = (in + 1) % BUFFER_SIZE;  // Wrap around to 0 when index reaches BUFFER_SIZE
    return num;
}
*/
/*
static inline uintptr_t getRandomRollback(uintptr_t rollback) {
    uintptr_t num = buffer[in];
    last_random_number = num;  // Store the last random number
    
    buffer[in] = getRandomPre(rollback);  // Generate a new random number and add it to the buffer

    in = ((in + 1) % BUFFER_SIZE) * rollback + in * (!rollback);  // Roll back if rollback is zero
    
    return num;
}
*/

static inline uintptr_t getRandomRollback(uintptr_t rollback) {
    uintptr_t num = buffer[in];
    last_random_number = num;  // Store the last random number
    uintptr_t new_num = getRandomPre(rollback);
    buffer[in] = (new_num & -rollback) | (num & ~-rollback);
    in = ((in + 1) & -rollback) | (in & ~-rollback) % BUFFER_SIZE;
    return num;
}

/*
static inline uintptr_t getRandomRollback(uintptr_t rollback) {
    uintptr_t num = buffer[in];
    buffer[in] = getRandomPre();  // Generate a new random number and add it to the buffer
    in = next_in;  // Use the next index
    next_in = ((next_in + 1) % BUFFER_SIZE) * rollback + next_in * (!rollback);  // Update the next index only if rollback is not zero
    return num;
}
*/
int main() {
    // Initialize the PRNG state
    prngState[0] = 13;
    prngState[1] = 0xDEADBEEF;

    prngStateOG[0] = 13;
    prngStateOG[1] = 0xDEADBEEF;

    // Precalculate random numbers
    precalculate_random_numbers();

    uintptr_t numrollback = getRandomRollback(1);   // match numroll
    uintptr_t numrollback2 = getRandomRollback(0); // nothing
    uintptr_t numrollback3 = getRandomRollback(0);  // match numroll2
    uintptr_t numrollback4 = getRandomRollback(1);  // match numroll3 
    
    uintptr_t numroll = getRandomPreOG();   // one number
    uintptr_t numroll2 = getRandomPreOG();   //two number
    uintptr_t numroll3 = getRandomPreOG();   //three number
    uintptr_t numroll4 = getRandomPreOG();   //four number


    printf("numrollback: %lu, numrollback2: %lu, numrollback3 = %lu, numrollback4 = %lu\n", numrollback, numrollback2, numrollback3, numrollback4);
    printf("numroll: %lu, numroll2: %lu, numroll3 = %lu, numroll4 = %lu\n", numroll, numroll2, numroll3, numroll4);

    

    
    return 0;
}