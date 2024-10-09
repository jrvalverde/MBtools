
/* Part 1 of torture-test of the memory controller */

# include <stdio.h>
/*
# include <console.h>
*/
void torture(void);
main()
{
    int i;
    for (i = 0; i < 100; i++)
	torture();
}



