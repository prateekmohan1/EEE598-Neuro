/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * neuron_new.h
 *
 * Code generation for function 'neuron_new'
 *
 */

#ifndef NEURON_NEW_H
#define NEURON_NEW_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "neuron_new_types.h"

/* Function Declarations */
extern double neuron_new(const boolean_T spike_exists_data[], const int
  spike_exists_size[2], const boolean_T synapses_data[], const int
  synapses_size[2], const double syn_weights_data[], const int syn_weights_size
  [2], double Vin, double ref_in, double lambda);

#endif

/* End of code generation (neuron_new.h) */
