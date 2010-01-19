/*! 	\file pfda_memory.h
	\brief Definitions for memory scheme.
*/
#ifndef PFDA_MEMORY_H_
#define PFDA_MEMORY_H_

/*! \defgroup memory Memory Management Functions

these functions are to manage the memory pools, the dp and ip parameters in nearly every function, that are passed into the code.

the reson for this is that the dynamically allocated memory was causeing a bottlenecck for performance.  
Also passing auxilliary parameters is too cumbersome.  The solution was to have pools of memory that were
divided up.

The functions use pointers to be able to incriment the poitner for the allocated memory, 
but not affect the memory pool pointer beyond the calling function.

*/
/*\@{*/
extern int * pfdaAlloc_i(int size, int ** pool);
extern double * pfdaAlloc_d(int size, double ** pool);
extern void pfdaFree(void*p);
/*\@}*/

#endif /* PFDA_MEMORY_H_ */
