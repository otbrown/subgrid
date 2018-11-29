#include <stdio.h>
#include "d3q15.h"

void report(Lattice *);

int main (void) {
  int nx = 10;
  int ny = 10;
  int nz = 50;
  double tau = 0.5;
  double amplitude = 0.01;
  const int numel = (nx+2) * (ny+2) * (nz+2);
  int i, j, k;

  // Initialisation
  printf("Creating lattice with nx = %d, ny = %d, nz = %d, and tau = %0.1f\n",
          nx, ny, nz, tau);
  Lattice * lat = d3q15_init(nx, ny, nz, tau, tau);
  bc_pbc_init(lat);
  lat->bc_func = bc_pbc_update;
  force_none_init(lat);
  lat->force_func = force_none_calc;

  // Initialise density and velocity
  for (i = 0; i < numel; i++) lat->rho_ptr[i] = 1.0;
  for (i = 0; i < numel * DQ_d; i++) lat->u_ptr[i] = 0.0;

  // Set enhanced density
  for (i = 1; i <= nx; i++) {
    for (j = 1; j <= ny; j++) {
      DQ_rho_get(lat, i, j, nz/2) += amplitude;
    }
  }
  double z_mean = 0.0;
  for (k = 1; k <= nz; k++) {
    z_mean += DQ_rho_get(lat, nx/2, ny/2, k);
  }
  z_mean /= nz;

  // renormalise
  for (i = 1; i <= nx; i++) {
    for (j = 1; j <= ny; j++) {
      for (k = 1; k <= nz; k++) {
        DQ_rho_get(lat, i, j, k) /= z_mean;
      }
    }
  }

  // initFromHydroVars()
  for (i = 1; i <= nx; i++) {
    for (j = 1; j <= ny; j++) {
      for (k = 1; k <= nz; k++) {
        calc_equil(DQ_rho_get(lat,i,j,k), &DQ_u_get(lat,i,j,k,0), &DQ_f_get(lat,i,j,k,0));
      }
    }
  }

  // calculation
  report(lat);
  d3q15_iterate(lat, 100);
  calc_hydro(lat);
  report(lat);
  d3q15_iterate(lat, 100);
  calc_hydro(lat);
  report(lat);

  // Finalisation
  printf("Cleaning up.\n");
  d3q15_destroy(lat);

  return 0;
}

void report(Lattice * lat) {
  printf("Time step %d\n", lat->time_step);
  printf("Density:\n[ ");
  for (int i = 1; i <= lat->nz; i++) {
    printf("%g ", DQ_rho_get(lat, lat->nx/2, lat->ny/2, i) - 1.0);
  }
  printf("]\n");
  return;
}
