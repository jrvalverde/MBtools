/*
 *  P_SORT.H
 *
 *	Include file with the definitions needed for module
 *  P_Sort.C
 *
 *	Designed by:
 *	    J. R. Valverde	    8  - apr - 1990
 */

extern void DirSort(info_t **, counter, int (*)());

extern void BinSort(info_t **, counter, int (*)());

extern void SelSort(info_t **, counter, int (*)());

extern void BubSort(info_t **, counter, int (*)());

extern void ShkSort(info_t **, counter, int (*)());

extern void ShellSort(info_t **, counter, int (*)());

extern void HeapSort(info_t **, counter, int (*)());

extern void RQkSort(info_t **, counter, counter, int (*)());

extern void RQuickSort(info_t **, counter, int (*)());

extern void IQuickSort(info_t **, counter, int (*)());

extern void DirMerge(info_t **, counter, int (*)());

/*
 * J. R.
 */

