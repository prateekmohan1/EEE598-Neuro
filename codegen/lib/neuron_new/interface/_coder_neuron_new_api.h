/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_neuron_new_api.h
 *
 * Code generation for function '_coder_neuron_new_api'
 *
 */

#ifndef _CODER_NEURON_NEW_API_H
#define _CODER_NEURON_NEW_API_H

/* Include files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_neuron_new_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern real_T neuron_new(boolean_T spike_exists_data[], int32_T
  spike_exists_size[2], boolean_T synapses_data[], int32_T synapses_size[2],
  real_T syn_weights_data[], int32_T syn_weights_size[2], real_T Vin, real_T
  ref_in, real_T lambda);
extern void neuron_new_api(const mxArray * const prhs[6], const mxArray *plhs[1]);
extern void neuron_new_atexit(void);
extern void neuron_new_initialize(void);
extern void neuron_new_terminate(void);
extern void neuron_new_xil_terminate(void);

#endif

/* End of code generation (_coder_neuron_new_api.h) */
