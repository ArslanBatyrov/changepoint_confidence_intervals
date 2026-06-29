#include <R.h>

/* below is just a fake function to test that R can call into C and get a
   result back. It does not do anything useful. This is just a test */
void PELTC_checklist_fake(int *positions, double *likes, int *n){
    int i;
    int fill = 3;
    if (fill > *n) {fill = *n; } /*something to study further for me, was done for safety */

    for(i = 0; i < fill; i++){
        positions[i] = i + 1;  /* fake candidates */
        likes[i] = 10 * (i + 1); /*fake score paired per candidate*/
    }
}
