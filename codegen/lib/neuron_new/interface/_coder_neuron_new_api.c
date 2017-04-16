/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_neuron_new_api.c
 *
 * Code generation for function '_coder_neuron_new_api'
 *
 */

/* Include files */
#include "tmwtypes.h"
#include "_coder_neuron_new_api.h"
#include "_coder_neuron_new_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131450U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "neuron_new",                        /* fFunctionName */
  NULL,                                /* fRTCallStack */
  false,                               /* bDebugMode */
  { 2045744189U, 2170104910U, 2743257031U, 4284093946U },/* fSigWrd */
  NULL                                 /* fSigMem */
};

/* Function Declarations */
static void b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, boolean_T **y_data, int32_T y_size[2]);
static void c_emlrt_marshallIn(const mxArray *syn_weights, const char_T
  *identifier, real_T **y_data, int32_T y_size[2]);
static void d_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, real_T **y_data, int32_T y_size[2]);
static real_T e_emlrt_marshallIn(const mxArray *Vin, const char_T *identifier);
static void emlrt_marshallIn(const mxArray *spike_exists, const char_T
  *identifier, boolean_T **y_data, int32_T y_size[2]);
static const mxArray *emlrt_marshallOut(const real_T u);
static real_T f_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId);
static void g_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, boolean_T **ret_data, int32_T ret_size[2]);
static void h_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, real_T **ret_data, int32_T ret_size[2]);
static real_T i_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId);

/* Function Definitions */
static void b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, boolean_T **y_data, int32_T y_size[2])
{
  g_emlrt_marshallIn(emlrtAlias(u), parentId, y_data, y_size);
  emlrtDestroyArray(&u);
}

static void c_emlrt_marshallIn(const mxArray *syn_weights, const char_T
  *identifier, real_T **y_data, int32_T y_size[2])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  d_emlrt_marshallIn(emlrtAlias(syn_weights), &thisId, y_data, y_size);
  emlrtDestroyArray(&syn_weights);
}

static void d_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, real_T **y_data, int32_T y_size[2])
{
  h_emlrt_marshallIn(emlrtAlias(u), parentId, y_data, y_size);
  emlrtDestroyArray(&u);
}

static real_T e_emlrt_marshallIn(const mxArray *Vin, const char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = f_emlrt_marshallIn(emlrtAlias(Vin), &thisId);
  emlrtDestroyArray(&Vin);
  return y;
}

static void emlrt_marshallIn(const mxArray *spike_exists, const char_T
  *identifier, boolean_T **y_data, int32_T y_size[2])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  b_emlrt_marshallIn(emlrtAlias(spike_exists), &thisId, y_data, y_size);
  emlrtDestroyArray(&spike_exists);
}

static const mxArray *emlrt_marshallOut(const real_T u)
{
  const mxArray *y;
  const mxArray *m0;
  y = NULL;
  m0 = emlrtCreateDoubleScalar(u);
  emlrtAssign(&y, m0);
  return y;
}

static real_T f_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId)
{
  real_T y;
  y = i_emlrt_marshallIn(emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static void g_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, boolean_T **ret_data, int32_T ret_size[2])
{
  static const int32_T dims[2] = { 1, 784 };

  const boolean_T bv0[2] = { false, true };

  int32_T iv0[2];
  emlrtCheckVsBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "logical", false, 2U,
    dims, &bv0[0], iv0);
  ret_size[0] = iv0[0];
  ret_size[1] = iv0[1];
  *ret_data = (boolean_T *)mxGetData(src);
  emlrtDestroyArray(&src);
}

static void h_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, real_T **ret_data, int32_T ret_size[2])
{
  static const int32_T dims[2] = { 1, 784 };

  const boolean_T bv1[2] = { false, true };

  int32_T iv1[2];
  emlrtCheckVsBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", false, 2U,
    dims, &bv1[0], iv1);
  ret_size[0] = iv1[0];
  ret_size[1] = iv1[1];
  *ret_data = (real_T *)mxGetData(src);
  emlrtDestroyArray(&src);
}

static real_T i_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId)
{
  real_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", false, 0U,
    &dims);
  ret = *(real_T *)mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

void neuron_new_api(const mxArray * const prhs[6], const mxArray *plhs[1])
{
  boolean_T (*spike_exists_data)[784];
  int32_T spike_exists_size[2];
  boolean_T (*synapses_data)[784];
  int32_T synapses_size[2];
  real_T (*syn_weights_data)[784];
  int32_T syn_weights_size[2];
  real_T Vin;
  real_T ref_in;
  real_T lambda;

  /* Marshall function inputs */
  emlrt_marshallIn(emlrtAlias(prhs[0]), "spike_exists", (boolean_T **)
                   &spike_exists_data, spike_exists_size);
  emlrt_marshallIn(emlrtAlias(prhs[1]), "synapses", (boolean_T **)&synapses_data,
                   synapses_size);
  c_emlrt_marshallIn(emlrtAlias(prhs[2]), "syn_weights", (real_T **)
                     &syn_weights_data, syn_weights_size);
  Vin = e_emlrt_marshallIn(emlrtAliasP(prhs[3]), "Vin");
  ref_in = e_emlrt_marshallIn(emlrtAliasP(prhs[4]), "ref_in");
  lambda = e_emlrt_marshallIn(emlrtAliasP(prhs[5]), "lambda");

  /* Invoke the target function */
  Vin = neuron_new(*spike_exists_data, spike_exists_size, *synapses_data,
                   synapses_size, *syn_weights_data, syn_weights_size, Vin,
                   ref_in, lambda);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(Vin);
}

void neuron_new_atexit(void)
{
  mexFunctionCreateRootTLS();
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  neuron_new_xil_terminate();
}

void neuron_new_initialize(void)
{
  mexFunctionCreateRootTLS();
  emlrtClearAllocCountR2012b(emlrtRootTLSGlobal, false, 0U, 0);
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

void neuron_new_terminate(void)
{
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (_coder_neuron_new_api.c) */
