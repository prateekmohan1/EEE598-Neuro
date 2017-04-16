/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * main.c
 *
 * Code generation for function 'main'
 *
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/
/* Include files */
#include "rt_nonfinite.h"
#include "neuron_new.h"
#include "main.h"
#include "neuron_new_terminate.h"
#include "neuron_new_initialize.h"

/* Function Declarations */
static void argInit_1xd784_boolean_T(boolean_T result_data[], int result_size[2]);
static void argInit_1xd784_real_T(double result_data[], int result_size[2]);
static boolean_T argInit_boolean_T(void);
static double argInit_real_T(void);
static void main_neuron_new(void);

/* Function Definitions */
static void argInit_1xd784_boolean_T(boolean_T result_data[], int result_size[2])
{
  int idx1;

  /* Set the size of the array.
     Change this size to the value that the application requires. */
  result_size[0] = 1;
  result_size[1] = 2;

  /* Loop over the array to initialize each element. */
  for (idx1 = 0; idx1 < 2; idx1++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result_data[idx1] = argInit_boolean_T();
  }
}

static void argInit_1xd784_real_T(double result_data[], int result_size[2])
{
  int idx1;

  /* Set the size of the array.
     Change this size to the value that the application requires. */
  result_size[0] = 1;
  result_size[1] = 2;

  /* Loop over the array to initialize each element. */
  for (idx1 = 0; idx1 < 2; idx1++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result_data[idx1] = argInit_real_T();
  }
}

static boolean_T argInit_boolean_T(void)
{
  return false;
}

static double argInit_real_T(void)
{
  return 0.0;
}

static void main_neuron_new(void)
{
  boolean_T spike_exists_data[784];
  int spike_exists_size[2];
  boolean_T synapses_data[784];
  int synapses_size[2];
  double syn_weights_data[784];
  int syn_weights_size[2];
  double Vout;

  /* Initialize function 'neuron_new' input arguments. */
  /* Initialize function input argument 'spike_exists'. */
  argInit_1xd784_boolean_T(spike_exists_data, spike_exists_size);

  /* Initialize function input argument 'synapses'. */
  argInit_1xd784_boolean_T(synapses_data, synapses_size);

  /* Initialize function input argument 'syn_weights'. */
  argInit_1xd784_real_T(syn_weights_data, syn_weights_size);

  /* Call the entry-point 'neuron_new'. */
  Vout = neuron_new(spike_exists_data, spike_exists_size, synapses_data,
                    synapses_size, syn_weights_data, syn_weights_size,
                    argInit_real_T(), argInit_real_T(), argInit_real_T());
}

int main(int argc, const char * const argv[])
{
  (void)argc;
  (void)argv;

  /* Initialize the application.
     You do not need to do this more than one time. */
  neuron_new_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_neuron_new();

  /* Terminate the application.
     You do not need to do this more than one time. */
  neuron_new_terminate();
  return 0;
}

/* End of code generation (main.c) */
