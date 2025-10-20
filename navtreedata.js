/*
 @licstart  The following is the entire license notice for the JavaScript code in this file.

 The MIT License (MIT)

 Copyright (C) 1997-2020 by Dimitri van Heesch

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software
 and associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
 BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 @licend  The above is the entire license notice for the JavaScript code in this file
*/
var NAVTREE =
[
  [ "SG++-Doxygen-Documentation", "index.html", [
    [ "SG++: General Sparse Grid Toolbox", "index.html", "index" ],
    [ "Copyright", "copyright.html", null ],
    [ "Developer Manual", "development.html", [
      [ "Overview", "development.html#development_overview", null ],
      [ "Eclipse", "development.html#development_eclipse", [
        [ "Setup", "development.html#development_eclipse_setup", null ],
        [ "Fix C++11 header indexing", "development.html#development_eclipse_header", null ],
        [ "Configure SCons build", "development.html#development_eclipse_scons", [
          [ "Overwrite Eclipse Build", "development.html#development_eclipse_scons_overwrite", null ],
          [ "Install SConsolidator", "development.html#development_eclipse_scons_sconsolidator", null ]
        ] ]
      ] ],
      [ "Testing", "development.html#development_testing", null ],
      [ "Documentations, Styleguide, and Doxygen", "development.html#development_doxygen", [
        [ "Usage", "development.html#development_doxygen_usage", null ]
      ] ],
      [ "Coding", "development.html#development_coding", [
        [ "Comments", "development.html#development_coding_comments", null ],
        [ "Style", "development.html#development_coding_style", null ],
        [ "Naming", "development.html#development_coding_naming", null ],
        [ "Don'ts", "development.html#development_coding_donts", null ]
      ] ],
      [ "Creating new modules", "development.html#development_newmodules", null ]
    ] ],
    [ "Usage Examples", "examples.html", [
      [ "C++ Examples", "examples_cpp.html", [
        [ "Module sgpp::base", "examples_cpp.html#examples_cpp_module_base", null ],
        [ "Module sgpp::combigrid", "examples_cpp.html#examples_cpp_module_combigrid", null ],
        [ "Module sgpp::datadriven", "examples_cpp.html#examples_cpp_module_datadriven", null ],
        [ "Module sgpp::optimization", "examples_cpp.html#examples_cpp_module_optimization", null ],
        [ "Module sgpp::solver", "examples_cpp.html#examples_cpp_module_solver", null ],
        [ "benchmark_gridInteraction.cpp", "example_benchmark_gridInteraction_cpp.html", null ],
        [ "Using the DataMatrix object", "example_dataMatrixSerializeDemo_cpp.html", null ],
        [ "Using the DataVector object", "example_dataVectorSerializeDemo_cpp.html", null ],
        [ "Detect the configuration of OpenCL platforms", "example_detectPlatformConfiguration_cpp.html", null ],
        [ "Interaction-Term aware sparse grids.", "example_gridInteractionExample_cpp.html", null ],
        [ "Generalised Sparse Grids", "example_gridTExample_cpp.html", null ],
        [ "Using JSON", "example_json_cpp.html", null ],
        [ "Spatially-Dimension-Adaptive Refinement in C++", "example_predictiveRefinement_cpp.html", null ],
        [ "Quadrature in C++", "example_quadrature_cpp.html", null ],
        [ "Refinement Example", "example_refinement_cpp.html", null ],
        [ "tutorial.cpp (Start Here)", "example_tutorial_cpp.html", null ],
        [ "Grid unserialization", "example_unserializeGrid_cpp.html", null ],
        [ "List of different Grid Types", "GridTypes.html", null ],
        [ "Combigrid Example Dimensional Adaptivity (C++)", "example_combigrid_adaptive_cpp.html", null ],
        [ "Combigrid Example (C++)", "example_combigrid_cpp.html", null ],
        [ "examplePCE.cpp", "example_examplePCE_cpp.html", null ],
        [ "Learner Classification Test", "example_learnerClassificationTest_cpp.html", null ],
        [ "Regression Learner", "example_learnerRegressionTest_cpp.html", null ],
        [ "Learner SGDE OnOff", "example_learnerSGDEOnOffTest_cpp.html", null ],
        [ "learner SGDE", "example_learnerSGDETest_cpp.html", null ],
        [ "Learner SGD", "example_learnerSGDTest_cpp.html", null ],
        [ "new_sgde.cpp", "example_new_sgde_cpp.html", null ],
        [ "optimize_kde_bandwidth.cpp", "example_optimize_kde_bandwidth_cpp.html", null ],
        [ "Constrained Optimization", "example_constrainedOptimization_cpp.html", null ],
        [ "Fuzzy Extension Principle (C++)", "example_fuzzy_cpp.html", null ],
        [ "Optimization Example (C++)", "example_optimization_cpp.html", null ],
        [ "FISTA Solver", "example_fistaExample_cpp.html", null ]
      ] ],
      [ "Python Examples", "examples_py.html", [
        [ "Module sgpp::base", "examples_py.html#examples_py_module_base", null ],
        [ "Module sgpp::combigrid", "examples_py.html#examples_py_module_combigrid", null ],
        [ "Module sgpp::datadriven", "examples_py.html#examples_py_module_datadriven", null ],
        [ "Module sgpp::optimization", "examples_py.html#examples_py_module_optimization", null ],
        [ "Module sgpp::pde", "examples_py.html#examples_py_module_pde", null ],
        [ "Using the DataMatrix object", "example_dataMatrixSerializeDemo_py.html", null ],
        [ "Using the DataVector object", "example_dataVectorSerializeDemo_py.html", null ],
        [ "Generalised Sparse Grids", "example_gridTExample_py.html", null ],
        [ "Spatially-Dimension-Adaptive Refinement of ANOVA Components in Python", "example_predictiveANOVARefinement_py.html", null ],
        [ "Spatially-Dimension-Adaptive Refinement in Python", "example_predictiveRefinement_py.html", null ],
        [ "Quadrature in Python", "example_quadrature_py.html", null ],
        [ "refinement.py", "example_refinement_py.html", null ],
        [ "Dimension-Adaptive Refinement in Python", "example_subspaceRefinement_py.html", null ],
        [ "tutorial.py (Start Here)", "example_tutorial_py.html", null ],
        [ "Combigrid Example Dimensional Adaptivity (Python)", "example_combigrid_adaptive_py.html", null ],
        [ "Combigrid Example (Python)", "example_combigrid_py.html", null ],
        [ "Generalised Sparse Grids", "example_generalisedGridsest_py.html", null ],
        [ "learnerExample.py", "example_learnerExample_py.html", null ],
        [ "learnerSGDETest.py", "example_learnerSGDETest_py.html", null ],
        [ "positive_density.py", "example_positive_density_py.html", null ],
        [ "test_Rosenblatt.py", "example_test_Rosenblatt_py.html", null ],
        [ "Fuzzy Extension Principle (Python)", "example_fuzzy_py.html", null ],
        [ "Optimization Example (Python)", "example_optimization_py.html", null ],
        [ "splineResponseSurface_example.py", "example_splineResponseSurface_example_py.html", null ],
        [ "LTwoDotTest.py", "example_LTwoDotTest_py.html", null ]
      ] ],
      [ "Java Examples", "examples_java.html", [
        [ "Module sgpp::base", "examples_java.html#examples_java_module_base", null ],
        [ "Module sgpp::datadriven", "examples_java.html#examples_java_module_datadriven", null ],
        [ "Module sgpp::optimization", "examples_java.html#examples_java_module_optimization", null ],
        [ "Refinement Example", "example_refinement_java.html", null ],
        [ "tutorial.java (Start Here)", "example_tutorial_java.html", null ],
        [ "Learner SGDE", "example_example_learnerSGDE_java.html", null ],
        [ "Optimization Example (Java)", "example_optimization_java.html", null ]
      ] ],
      [ "MATLAB Examples", "examples_m.html", [
        [ "Module sgpp::base", "examples_m.html#examples_m_module_base", null ],
        [ "Module sgpp::optimization", "examples_m.html#examples_m_module_optimization", null ],
        [ "tutorial.m (Start Here)", "example_tutorial_m.html", null ],
        [ "Optimization Example (MATLAB)", "example_optimization_m.html", null ]
      ] ]
    ] ],
    [ "Integrate Dakota", "install_dakota.html", null ],
    [ "SGDE Miner", "example_SGDEMinerFromConfigFile_py.html", null ],
    [ "Deprecated List", "deprecated.html", null ],
    [ "Todo List", "todo.html", null ],
    [ "Namespaces", "namespaces.html", [
      [ "Namespace List", "namespaces.html", "namespaces_dup" ],
      [ "Namespace Members", "namespacemembers.html", [
        [ "All", "namespacemembers.html", "namespacemembers_dup" ],
        [ "Functions", "namespacemembers_func.html", "namespacemembers_func" ],
        [ "Variables", "namespacemembers_vars.html", "namespacemembers_vars" ],
        [ "Typedefs", "namespacemembers_type.html", null ],
        [ "Enumerations", "namespacemembers_enum.html", null ],
        [ "Enumerator", "namespacemembers_eval.html", null ]
      ] ]
    ] ],
    [ "Classes", "annotated.html", [
      [ "Class List", "annotated.html", "annotated_dup" ],
      [ "Class Hierarchy", "hierarchy.html", "hierarchy" ],
      [ "Class Members", "functions.html", [
        [ "All", "functions.html", "functions_dup" ],
        [ "Functions", "functions_func.html", "functions_func" ],
        [ "Variables", "functions_vars.html", "functions_vars" ],
        [ "Typedefs", "functions_type.html", null ],
        [ "Enumerations", "functions_enum.html", null ],
        [ "Enumerator", "functions_eval.html", null ],
        [ "Related Symbols", "functions_rela.html", null ]
      ] ]
    ] ],
    [ "Files", "files.html", [
      [ "File List", "files.html", "files_dup" ],
      [ "File Members", "globals.html", [
        [ "All", "globals.html", null ],
        [ "Functions", "globals_func.html", null ],
        [ "Variables", "globals_vars.html", null ],
        [ "Typedefs", "globals_type.html", null ],
        [ "Macros", "globals_defs.html", null ]
      ] ]
    ] ]
  ] ]
];

var NAVTREEINDEX =
[
"ANOVAHashRefinement_8cpp.html",
"DBMatDMSOrthoAdapt_8cpp.html",
"DensityEstimationConfiguration_8hpp.html",
"GridFactory_8cpp.html",
"LaplaceUpPrewavelet_8cpp.html",
"NakPBsplineBasis_8hpp.html#aa9a5d095d317c6a00687331efe460055",
"OperationEvalModLinear_8cpp.html",
"OperationLaplaceModPoly_8cpp.html",
"OperationMultipleEvalSubspaceSimple_8cpp.html#a0ab892f09d0d501b4cba6e0873f85dd9",
"OperationTestPoly_8hpp.html",
"SamplerTypes_8hpp.html#a4a0d17759496b89f8d523bfded8ecceba27abb4b8e57a88a6a15db17f673ee25f",
"VizualizerDummy_8cpp.html",
"classjson_1_1IDNode.html#aaaeafb872d409c682837a61a5b14cc99",
"classpython_1_1controller_1_1CheckpointController_1_1CheckpointController.html#a0186ca412b5240f4e0a15cf624fd9307",
"classpython_1_1learner_1_1LearnerBuilder_1_1LearnerBuilder.html#a6c31d4bcb38566e53138964b89203bfb",
"classpython_1_1learner_1_1TrainingStopPolicy_1_1TrainingStopPolicy.html#af2f157694e4039e43e9eba2470862696",
"classpython_1_1uq_1_1analysis_1_1asgc_1_1ASGCAnalysisBuilder_1_1ASGCAnalysisBuilder.html#a306e0637cf85b0faa319c8dd8442b385",
"classpython_1_1uq_1_1dists_1_1Beta_1_1Beta.html#a55dd018edca05f4ec9111c2d8c6b6682",
"classpython_1_1uq_1_1dists_1_1LibAGFDist_1_1LibAGFDist.html#a3e56281a89255f8ec5d2454fdb936cb2",
"classpython_1_1uq_1_1dists_1_1TLognormal_1_1TLognormal.html#a5867b3ae3d058e135925f32e6c30cea3",
"classpython_1_1uq_1_1learner_1_1Learner_1_1Learner.html#a64b0db6d2172236f874b47e1690ef7cc",
"classpython_1_1uq_1_1learner_1_1builder_1_1RegressorSpecificationDescriptor_1_1FoldingDescriptor.html",
"classpython_1_1uq_1_1models_1_1Model_1_1Model.html",
"classpython_1_1uq_1_1operations_1_1forcePositivity_1_1localFullGridSearch_1_1LocalFullGridCandidates.html#ae42c8f4e674818884d26cd805497e1a2",
"classpython_1_1uq_1_1parameters_1_1ParameterBuilder_1_1ParameterBuilder.html",
"classpython_1_1uq_1_1quadrature_1_1HashQuadrature_1_1HashQuadratureMap.html#a07242899dece99878ffdc7941348a9cb",
"classpython_1_1uq_1_1refinement_1_1RefinementManagerDescriptor_1_1RefineCurrentNodesDescriptor.html#abef871198a80d341295c94d8d114326d",
"classpython_1_1uq_1_1refinement_1_1RefinementStrategy_1_1VarianceOptRanking.html#a05bf14caf1ee21c17a58539117d6b96b",
"classpython_1_1uq_1_1sampler_1_1asgc_1_1ASGCSamplerSpecification_1_1ASGCSamplerSpecification.html#a8b2de6f93a1365e08d7c1f752e165998",
"classpython_1_1uq_1_1uq__setting_1_1UQSettingManager_1_1UQSettingManager.html#a24e707b5f9c2f9df9821bcc365e38c09",
"classpython_1_1uq_1_1uq__setting_1_1samplingresult_1_1Samplingresult.html#ab4595159dd450f582cc1b762dd8c0e11",
"classsgpp_1_1base_1_1BoundingBox.html#a6a2b256c642761425a663fe585f2ab54",
"classsgpp_1_1base_1_1BsplineModifiedClenshawCurtisBasis.html#a4de75b50de92c24b78963e7c618e901e",
"classsgpp_1_1base_1_1DataMatrix.html#a02759810456bd088f67e21b968f77a89",
"classsgpp_1_1base_1_1DataMatrixSP.html#adf3ad0d5c79543f886c5697fb8a78744",
"classsgpp_1_1base_1_1DataVectorSP.html#afda70671d3f3cdfdd86cde55a84501f0",
"classsgpp_1_1base_1_1DehierarchisationPolyBoundary.html#a15bd069212175961b3398efe6b0164ab",
"classsgpp_1_1base_1_1DistributionsVector.html#aedf40ccf6a1cd34bc07e561001c89b04",
"classsgpp_1_1base_1_1FundamentalSplineModifiedBasis.html#a16b15f2cb3b3aa5eb1b1d1dd96073f44",
"classsgpp_1_1base_1_1Grid.html#af615aa3ebfcadbfabb175052149c9482",
"classsgpp_1_1base_1_1HashGridIterator.html#ab56178e3bfe2c9ac32c4c25cca4c2d21",
"classsgpp_1_1base_1_1HashRefinement.html#a3659fd43ed26f12da3075d713cd44fa5",
"classsgpp_1_1base_1_1HierarchisationModLinearClenshawCurtis.html#ab77711a838bf90911e019de3f5811a6c",
"classsgpp_1_1base_1_1HierarchisationSLE.html#adfbad40f73db959d2a2d1adb3c1a3aa7",
"classsgpp_1_1base_1_1L0BoundaryGridGenerator.html",
"classsgpp_1_1base_1_1LinearModifiedClenshawCurtisBasis.html#aa0f7c0592cb9ba8939601cd17c9b77c7",
"classsgpp_1_1base_1_1ModPolyClenshawCurtisGrid.html#a09799571e4600fd65be46549905eb24d",
"classsgpp_1_1base_1_1NakBsplineExtendedBasis.html#a571e4a89b10463161e76f998e7ebd065",
"classsgpp_1_1base_1_1OCLClonedBuffer.html#abbe14159a4acb142122138f3c2131955",
"classsgpp_1_1base_1_1OperationConvertPrewavelet.html#a0c538ba910cb66398ccb4a302e044af3",
"classsgpp_1_1base_1_1OperationEvalGradientModFundamentalSplineNaive.html#a0d1830ef8f44b2c23d19ede376ca86f4",
"classsgpp_1_1base_1_1OperationEvalHessian.html#aca80c9a0bd821c8b6bfe06edf2b2eae5",
"classsgpp_1_1base_1_1OperationEvalHessianWaveletNaive.html",
"classsgpp_1_1base_1_1OperationEvalModLinearNaive.html#ab98adf0563e5de77471ca96271879b43",
"classsgpp_1_1base_1_1OperationEvalPartialDerivativeBsplineNaive.html#a8633c6a9d3d428215cd1cfd5b24b6ebd",
"classsgpp_1_1base_1_1OperationEvalPolyBoundary.html#a0daee92ea6c3850cc7b5c4400189579e",
"classsgpp_1_1base_1_1OperationFirstMomentModPolyClenshawCurtis.html",
"classsgpp_1_1base_1_1OperationHierarchisationModPoly.html#aa4730ee55a9a268b7c464a328efaba1e",
"classsgpp_1_1base_1_1OperationMultipleEvalLinearBoundary.html#a9e003d4952fd97b5515fc3726d347c9e",
"classsgpp_1_1base_1_1OperationMultipleEvalNakBsplineModifiedNaive.html#a9b60af532ad49cf831485ec839a34502",
"classsgpp_1_1base_1_1OperationQuadratureMC.html#a4cea599d5e4aeb41da9e293b26fa48dd",
"classsgpp_1_1base_1_1OperationSecondMomentModBspline.html",
"classsgpp_1_1base_1_1OperationWeightedQuadratureNakBsplineModified.html#af417a6ca1450996c0cb037ea2dcf0ad0",
"classsgpp_1_1base_1_1PolyGrid.html#a17d168df4830d426d92922e31c7fc6bf",
"classsgpp_1_1base_1_1Printer.html#a9e56f510e95aa56cb7217990bc0566b5",
"classsgpp_1_1base_1_1ScaledScalarFunctionGradient.html#a6ccac4192482d0985240749fc80e9c4c",
"classsgpp_1_1base_1_1StencilHierarchisationModLinear.html#a751c5db43f582d9a31495b576bf884f7",
"classsgpp_1_1base_1_1VectorFunctionGradient.html#a812d75ded199767c298c714da1b30f66",
"classsgpp_1_1base_1_1WeaklyFundamentalNakSplineModifiedBasisDeriv2.html#a488bffddf625ca24d7762a16268980ee",
"classsgpp_1_1base_1_1not__implemented__exception.html#aa5cf5e97d6a97ba277d5b99d989929c5",
"classsgpp_1_1combigrid_1_1CombinationGrid.html#acf3c0f0dec52aad69274e53dc4e9265a",
"classsgpp_1_1combigrid_1_1OperationPole.html#af18108fb3f956bdcb602d6a6acf92471",
"classsgpp_1_1datadriven_1_1AlgorithmAdaBoostBase.html#ac156f3c66500adff03be8b94d127ca26",
"classsgpp_1_1datadriven_1_1ClassificationRefinementFunctor.html#a0d91d8b148c4935cb5b1b1046f7a0649",
"classsgpp_1_1datadriven_1_1DBMatOffline.html#a60f5101fbbece0d114ecd689dd941f1d",
"classsgpp_1_1datadriven_1_1DBMatOnlineDE.html#a60527aafece10a30d4965b1e23f507be",
"classsgpp_1_1datadriven_1_1DMWeightMatrix.html",
"classsgpp_1_1datadriven_1_1DataShufflingFunctorCrossValidation.html#a3d92df49c86337c058444c5da746cd17",
"classsgpp_1_1datadriven_1_1DensityDifferenceSystemMatrix.html",
"classsgpp_1_1datadriven_1_1DiscreteParameter.html#a7330e07b1faeae0d84d5bc4a800fbec0",
"classsgpp_1_1datadriven_1_1GridFactory.html#a7cc4e8f146c12ff4465ac8015ecf0063",
"classsgpp_1_1datadriven_1_1KernelDensityEstimator.html#af20a9ae1cce434c3435ee5e1f5b50e51",
"classsgpp_1_1datadriven_1_1LearnerSGDE.html#a0ec4cfd39e0623940a977250dc454f14",
"classsgpp_1_1datadriven_1_1LearnerSVM.html#a66e6ebc2035c1975366903b502559a60",
"classsgpp_1_1datadriven_1_1ModelFittingBase.html#a658d0e678584bb4f3f232529e15d831c",
"classsgpp_1_1datadriven_1_1ModelFittingDensityEstimationCombi.html",
"classsgpp_1_1datadriven_1_1MultiSurplusRefinementFunctor.html#afb7e9f19bf58872d10ca55086f8f8da7",
"classsgpp_1_1datadriven_1_1OperationInverseRosenblattTransformation1DLinear.html#a1b647abf5d230be56acd0256e366d67c",
"classsgpp_1_1datadriven_1_1OperationInverseRosenblattTransformationPoly.html#a090f2e3eeaf7c0821372df5b57b1efd7",
"classsgpp_1_1datadriven_1_1OperationMultiEvalCuda.html#a205678f5f8bd43e6628d5bf4ef251be7",
"classsgpp_1_1datadriven_1_1OperationMultiEvalStreamingModOCLFastMultiPlatform.html#a524beb5481b3b2b6655da0431fe1b018",
"classsgpp_1_1datadriven_1_1OperationMultipleEvalDistributed.html#a73ce51812d66f74798d348a4d6bbe6cd",
"classsgpp_1_1datadriven_1_1OperationRosenblattTransformation1DPolyClenshawCurtisBoundary.html#a19eeb1d2020a1a21a3cb3151242d969e",
"classsgpp_1_1datadriven_1_1OperationTestLinear.html",
"classsgpp_1_1datadriven_1_1PolynomialChaosExpansion.html#a4e812d60ae4aa1a7ca4295cf635829b6",
"classsgpp_1_1datadriven_1_1ScorerFactory.html#a3ea671c7d6de2abeecd9626ba926c25c",
"classsgpp_1_1datadriven_1_1StreamingBSplineOCLKernelSourceBuilder.html#a601406c44ecb554799e5b17a5367af1c",
"classsgpp_1_1datadriven_1_1SubspaceNodeCombined.html#a1fad7f511e02aaf95b7a9d3b642af42c",
"classsgpp_1_1datadriven_1_1VisualizerDensityEstimation.html#ab827b0b7fda7b052fd719e2595d534ee",
"classsgpp_1_1datadriven_1_1clusteringmpi_1_1MPIWorkerPackageBase.html#af784e2e0e3072540077e3106740072e3",
"classsgpp_1_1optimization_1_1FuzzyExtensionPrincipleViaTransformation.html#aa01a0f263038432be5cb9f695b0f1a00",
"classsgpp_1_1optimization_1_1IterativeGridGeneratorLinearSurplus.html",
"classsgpp_1_1optimization_1_1OperationMultipleHierarchisationLinearClenshawCurtisBoundary.html#a6e3e350f912bb18a69275f40605858ca",
"classsgpp_1_1optimization_1_1QuasiGaussianFuzzyNumber.html",
"classsgpp_1_1optimization_1_1optimizer_1_1AdaptiveGradientDescent.html#af7cb90348d31dfe18056233f5283a7af",
"classsgpp_1_1optimization_1_1optimizer_1_1GradientDescent.html#a87af226146e72f2d7d1cb29ed31da587",
"classsgpp_1_1optimization_1_1optimizer_1_1NelderMead.html#a3f7ac46dd28c87c45829abc3fa9798db",
"classsgpp_1_1optimization_1_1test__problems_1_1Ackley.html#a04d1d8ec9d9adaa8996ad13a529de94b",
"classsgpp_1_1optimization_1_1test__problems_1_1FloudasObjective.html#aadcfebff4721c8bc3030e02d72937113",
"classsgpp_1_1optimization_1_1test__problems_1_1G06Objective.html",
"classsgpp_1_1optimization_1_1test__problems_1_1G12Objective.html",
"classsgpp_1_1optimization_1_1test__problems_1_1Mladineo.html#aec432c735db7960d1743ea5bf67145bb",
"classsgpp_1_1optimization_1_1test__problems_1_1SolandInequalityConstraint.html#a8b4801fa9c8298584680ad019d3b4428",
"classsgpp_1_1pde_1_1HeatEquationSolverWithStretching.html#a4b842e84c97d964e4cdd65efcbf63d0f",
"classsgpp_1_1pde_1_1OperationLTwoDotProductLinearBoundary.html#aae7eaf13f2dfba6d8e77e70ad2308ef8",
"classsgpp_1_1pde_1_1OperationLaplacePolyClenshawCurtisBoundary.html#a8e280f8ee2fc54f102b1ed4c54ddf319",
"classsgpp_1_1pde_1_1OperationMatrixLTwoDotModLinear.html#af0cca8fa1f2bed5d6db8acfb6afe9216",
"classsgpp_1_1pde_1_1PhiPhiDownModLinear.html#aaac53892be12cd3816da85fc93a94ace",
"classsgpp_1_1pde_1_1UpDownFourOpDims.html#af58426e27e8136901e2f67e076b7208b",
"classsgpp_1_1quadrature_1_1OperationQuadratureMCAdvanced.html#af474430edf1d73df656750e6ca4ffa6b",
"classsgpp_1_1solver_1_1SGSolver.html#a25e172e5f3db16e9bf7b4d6f11a7fea6",
"datadriven_2python_2data_2____init_____8py.html",
"dir_cc58c086ea718db481d79e6432c2a345.html",
"friedman_8mod_8py.html#a6801ba8a9f103c083cf618824c33a5b2",
"namespaceSGDEMinerFromConfigFile.html#ac7c3ec3af7e40bcaf2d0a37dd6387ad3",
"namespacepython_1_1classifier.html#ae32c7bfcf85448d7f817255569c48471",
"namespacepython_1_1uq_1_1jsonLib.html#a3d8f4ce89e1c6922c937bf26a871715a",
"namespacepython_1_1uq_1_1sampler.html",
"namespacesgpp_1_1base.html#aee7c83d256e58a2d9de502fe390e9fdaa576fda5c28603004964d20937bcf803f",
"namespacesgpp_1_1datadriven_1_1PermutationUtil.html",
"plot1d_8py.html#a766e1c2c62d3dae1d2cb3b3d17c64da5",
"structsgpp_1_1base_1_1GeneralGridConfiguration.html#a8b0a54d8083289ed34f2e8612cd5eddd",
"structsgpp_1_1datadriven_1_1MergeGridNetworkMessage.html",
"toolsExtended_8py.html#a0bed280ba5dcab74f89911b4f05eb860"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';