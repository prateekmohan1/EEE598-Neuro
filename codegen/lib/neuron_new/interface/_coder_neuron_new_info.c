/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_neuron_new_info.c
 *
 * Code generation for function 'neuron_new'
 *
 */

/* Include files */
#include "_coder_neuron_new_info.h"

/* Function Definitions */
const mxArray *emlrtMexFcnResolvedFunctionsInfo(void)
{
  const mxArray *nameCaptureInfo;
  const char * data[7] = {
    "789ced58cb6ed340149da2b6824551251e628360811052a54c1ba5af5d93d42115495a9a5484562838ce6ded761ee9781c520989f21f7c004b962c59b265c18a"
    "9fc1713dad1dd538c8252d5546b2eca313f9dc7b737d7c3d686cad3c8610baed1ef50f087d7f80bc35e5e3e913886ea0f0eae7c7fcf37800f778b5265ca687a7",
    "02fc2b177ff2b1c19984ae3c01c46250716813840b984ee1f4362d4e2da633593b6a03126073d28196c7ec5a046a1685120f80a2e5025a0850a7a047f5aef326"
    "180755872261da67e192204081fa7c8bc87f3cfcf3c8fabc8ea8cfb48f15bfa3bdc126a780db42970007949b3ac3ab601f48dec6d9ea56a3c01d21cd461528d8",
    "1204d6346d7e79096f08be0f86c45497446f366cd96a63068ee0acc1e05d8af6e5f33622dec901f3099e83f9dc423703f8f18ad23b4e58bf8731f553bcc15b20"
    "5296db5482e924257989ef59864ebcbf7b80fce3e2e95f51f1a8a5f4be24cc7f23464ff13b6ba5badb428e2d30e16edeb89cad95b239bc999e9d5bd4b1e49c34",
    "79170325de31e3d50bcfa882e170c1527478fd727ce7fd8f9fd9e1f6e7b09f87cbd3eb46dc6fd0febb1fa1a7fa4ff115a22ff3627da1bc5d685bbb5bdb1a677b"
    "2c771687ead3289db83850009f97e7bfbeffe8393e3ffe70df3d5d49eab34a6722b25e131eb66c772c40c97d5de94d0670586fd2c32dee34099ce97d4ea8978f",
    "d43bc18aff8b7e2056d39f023010f7027b35f2278021faf92afa35f2f3abeae7f722f454df29fe30ebacd3cc513d33bb29d24e7e0e963bb482ae8f9f8f9edf3f"
    "f5db930b9bdf1fc5d449f17df3bb65e71c8bc835e67e1682b08c4bf3f9af09f53663f2577cd2f77e7fc1d487df107d3f3ff2fdabebfb83cef1eb9dc37c76b7d8",
    "59aabdcc1c191a65e9c236295e1fdf1fcdf1e7c71feebb672b17d5777763eaa5f83eff0721786f33f0c27c3f6a9f72dadfa7344c3da497743e781e93b7e293f6"
    "895728d7ec87e8f31fdb30f2f9ffdde79d8cb6554f3b2fba76c5b473e5fdea823d0fd760bfe637d63648c7",
    "" };

  nameCaptureInfo = NULL;
  emlrtNameCaptureMxArrayR2016a(data, 6368U, &nameCaptureInfo);
  return nameCaptureInfo;
}

mxArray *emlrtMexFcnProperties(void)
{
  mxArray *xResult;
  mxArray *xEntryPoints;
  const char * fldNames[4] = { "Name", "NumberOfInputs", "NumberOfOutputs",
    "ConstantInputs" };

  mxArray *xInputs;
  const char * b_fldNames[4] = { "Version", "ResolvedFunctions", "EntryPoints",
    "CoverageInfo" };

  xEntryPoints = emlrtCreateStructMatrix(1, 1, 4, fldNames);
  xInputs = emlrtCreateLogicalMatrix(1, 6);
  emlrtSetField(xEntryPoints, 0, "Name", mxCreateString("neuron_new"));
  emlrtSetField(xEntryPoints, 0, "NumberOfInputs", mxCreateDoubleScalar(6.0));
  emlrtSetField(xEntryPoints, 0, "NumberOfOutputs", mxCreateDoubleScalar(1.0));
  emlrtSetField(xEntryPoints, 0, "ConstantInputs", xInputs);
  xResult = emlrtCreateStructMatrix(1, 1, 4, b_fldNames);
  emlrtSetField(xResult, 0, "Version", mxCreateString("9.2.0.538062 (R2017a)"));
  emlrtSetField(xResult, 0, "ResolvedFunctions", (mxArray *)
                emlrtMexFcnResolvedFunctionsInfo());
  emlrtSetField(xResult, 0, "EntryPoints", xEntryPoints);
  return xResult;
}

/* End of code generation (_coder_neuron_new_info.c) */
