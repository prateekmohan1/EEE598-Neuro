/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * neuron_new.c
 *
 * Code generation for function 'neuron_new'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "neuron_new.h"

/* Function Definitions */

/*
 * function [Vout] = neuron_new (spike_exists, synapses, syn_weights, Vin,  ...
 *                                       ref_in, lambda)
 */
double neuron_new(const boolean_T spike_exists_data[], const int
                  spike_exists_size[2], const boolean_T synapses_data[], const
                  int synapses_size[2], const double syn_weights_data[], const
                  int syn_weights_size[2], double Vin, double ref_in, double
                  lambda)
{
  double Vout;
  double temp;
  int i;
  (void)spike_exists_size;
  (void)syn_weights_size;

  /* Scaling factor for increasing membrane potentials, the larger this the */
  /* larger the membrane increase */
  /* 'neuron_new:6' sc_factor = 1.0; */
  /* 'neuron_new:8' if ~ref_in */
  if (!(ref_in != 0.0)) {
    /* Vout = Vin - (Vin/(R*C)) + (I/C); */
    /* 'neuron_new:10' temp = 0; */
    temp = 0.0;

    /* 'neuron_new:11' for i = 1:size(synapses,2) */
    for (i = 0; i < synapses_size[1]; i++) {
      /* The spike_exists, synapses, and syn_weights are only for the layer that you are interested in */
      /* For example, if you are looking at a neuron in L2, the spike_exists contains if spike exists in neurons */
      /* in L1 */
      /* 'neuron_new:15' if (spike_exists(i) && synapses(i)) */
      if (spike_exists_data[i] && synapses_data[i]) {
        /* 'neuron_new:16' temp = temp + syn_weights(i)*sc_factor; */
        temp += syn_weights_data[i];
      }
    }

    /* 'neuron_new:19' Vout = Vin + temp - lambda; */
    Vout = (Vin + temp) - lambda;
  } else {
    /* 'neuron_new:20' else */
    /* 'neuron_new:21' Vout = 0; */
    Vout = 0.0;

    /*  reset voltage */
  }

  return Vout;
}

/* End of code generation (neuron_new.c) */
