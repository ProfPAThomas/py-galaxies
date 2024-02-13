/* dynamics.c
Routines related to the dynamics of galaxies within halos and subhalos within halos.
*/

#include "all_headers.h"

void F_set_central_galaxy(struct struct_sub *sub) {
    /*
    Determine the central galaxy in a subhalo.

    Will eventually have fancy code to determine which, if any, galaxy is the central one.
    For now, just make that the first (and only) galaxy.

    Arguments
    ---------
    sub : obj : struct struct_sub
       Properties of the subhalo
    */
    
    (*sub).gal_central_sid = (*sub).gal_start_sid;
    return;
}
