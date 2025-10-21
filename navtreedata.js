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
"LaplaceUpGradientPrewavelet_8cpp.html",
"NakBsplineModifiedBasis_8hpp.html#aeaeef09d115821a9f5cba93c27c321d9",
"OperationEvalModLinearNaive_8cpp.html",
"OperationLaplaceModPolyClenshawCurtis_8cpp.html",
"OperationMultipleEvalSubspaceSimpleParameters_8hpp.html#ab848a0ca31255d248f794a89a8ef77fb",
"OperationTestModWavelet_8hpp.html",
"SamplerTypes_8hpp.html",
"Visualizer_8cpp.html",
"classjson_1_1IDNode.html#aa890554fba906d1bf87c4b41a21e8d06",
"classpython_1_1controller_1_1CheckpointController_1_1CheckpointController.html",
"classpython_1_1learner_1_1LearnerBuilder_1_1LearnerBuilder.html#a6b797690343912f1e51745a6b49bc64c",
"classpython_1_1learner_1_1TrainingStopPolicy_1_1TrainingStopPolicy.html#ad31d954b2c986dc9f4e1ac7f9ac4612c",
"classpython_1_1uq_1_1analysis_1_1asgc_1_1ASGCAnalysisBuilder_1_1ASGCAnalysisBuilder.html#a1691ae72126876cc40991298d15968a2",
"classpython_1_1uq_1_1dists_1_1Beta_1_1Beta.html#a4ca7149e7a14c81497f9a1891b3e0671",
"classpython_1_1uq_1_1dists_1_1LibAGFDist_1_1LibAGFDist.html#a3dc7d1fff377a0276e79ef6fbe37cfc0",
"classpython_1_1uq_1_1dists_1_1TLognormal_1_1TLognormal.html#a41cf21cae67ec7a05ab996e5d887042d",
"classpython_1_1uq_1_1learner_1_1Learner_1_1Learner.html#a62fb1891f8bf7c7bd733c0fbb2eb718a",
"classpython_1_1uq_1_1learner_1_1builder_1_1LearnerBuilder_1_1LearnerBuilder.html#ad07496e243fac38f20becaf98e7b233b",
"classpython_1_1uq_1_1models_1_1AnalyticModels_1_1Parabola.html#a4e17edbde8d0aa26edcffc1451804d93",
"classpython_1_1uq_1_1operations_1_1forcePositivity_1_1localFullGridSearch_1_1LocalFullGridCandidates.html#adae7db32f6840133502dc2ae5c56d1cc",
"classpython_1_1uq_1_1parameters_1_1ParameterBuilder_1_1GeneralParameterBuilder.html#af2d74f0280a64f4439efaa2f4174deea",
"classpython_1_1uq_1_1quadrature_1_1HashQuadrature_1_1HashQuadratureMap.html",
"classpython_1_1uq_1_1refinement_1_1RefinementManagerDescriptor_1_1RefineCurrentNodesDescriptor.html#abc0ed6b0d40b03874c01566cf76aede5",
"classpython_1_1uq_1_1refinement_1_1RefinementStrategy_1_1VarianceOptRanking.html#a03cb8fd7b8debc2eba258b1444890543",
"classpython_1_1uq_1_1sampler_1_1asgc_1_1ASGCSamplerSpecification_1_1ASGCSamplerSpecification.html#a76fb6ab2566f417ac560a796b007cd93",
"classpython_1_1uq_1_1uq__setting_1_1UQSettingManager_1_1UQSettingManager.html#a1c5651f21beeaa650038c802fb4c6204",
"classpython_1_1uq_1_1uq__setting_1_1samplingresult_1_1Samplingresult.html#aab1c0d8d6d8df94fc692e9879f16217a",
"classsgpp_1_1base_1_1BoundingBox.html#a66c6c477b29f219a0450f081f3a6719b",
"classsgpp_1_1base_1_1BsplineModifiedClenshawCurtisBasis.html#a460d65a0cdabfdc4916f7312c9a09d4f",
"classsgpp_1_1base_1_1DataMatrix.html#a01b4c9775403ca3c7a12597806776aec",
"classsgpp_1_1base_1_1DataMatrixSP.html#adeefc1e7cd74d31952441f4e79cebdcf",
"classsgpp_1_1base_1_1DataVectorSP.html#afbed7c4b858a891f2fe9310c3a679f83",
"classsgpp_1_1base_1_1DehierarchisationPolyBoundary.html",
"classsgpp_1_1base_1_1DistributionsVector.html#ae3034d1def40eea62744023822f28986",
"classsgpp_1_1base_1_1FundamentalSplineModifiedBasis.html#a0c88c1b11bcb9e601a3163630eaea941",
"classsgpp_1_1base_1_1Grid.html#aef696e83f78c80da9332df27369b48b6",
"classsgpp_1_1base_1_1HashGridIterator.html#ab0b2d7b65dec1f1ac383d819e019aa66",
"classsgpp_1_1base_1_1HashRefinement.html#a2f59c8994044b294bb1b5ed0980c882c",
"classsgpp_1_1base_1_1HierarchisationModLinearClenshawCurtis.html#a8bb50e2d6644d1e25202b1ab5a4a52f7",
"classsgpp_1_1base_1_1HierarchisationSLE.html#ad69cd5e92ae2c845fd400e07d134ae10",
"classsgpp_1_1base_1_1KernelSourceBuilderBase.html#ae9d3fb3a1d8ff09135b3d3fdd155ea4a",
"classsgpp_1_1base_1_1LinearModifiedClenshawCurtisBasis.html#aa0234cb133a0775aa299d1b560023178",
"classsgpp_1_1base_1_1ModPolyClenshawCurtisGrid.html",
"classsgpp_1_1base_1_1NakBsplineExtendedBasis.html#a34e08e31dae0e9b947d9d279705a12c4",
"classsgpp_1_1base_1_1OCLClonedBuffer.html#ab1e4713352f3112657e01617835ee56b",
"classsgpp_1_1base_1_1OperationConvertPrewavelet.html",
"classsgpp_1_1base_1_1OperationEvalGradientModFundamentalSplineNaive.html",
"classsgpp_1_1base_1_1OperationEvalHessian.html#a184e2fc5920cddce4f9458517c44d7f0",
"classsgpp_1_1base_1_1OperationEvalHessianWaveletBoundaryNaive.html#aed49aad3470ea83b25d04db0ff040d84",
"classsgpp_1_1base_1_1OperationEvalModLinearNaive.html#aadb03cda0bad083ca002a01d54682d11",
"classsgpp_1_1base_1_1OperationEvalPartialDerivativeBsplineNaive.html#a854e5a5288a5561041c1c25a2b6081fa",
"classsgpp_1_1base_1_1OperationEvalPolyBoundary.html#a0a28e2771d8e3a0541350f524a624846",
"classsgpp_1_1base_1_1OperationFirstMomentModPoly.html#ad14930e9c591b0b7604fa691efbbb111",
"classsgpp_1_1base_1_1OperationHierarchisationModPoly.html#a88194a7f681a434c07463a8032da3819",
"classsgpp_1_1base_1_1OperationMultipleEvalLinearBoundary.html#a895b6bf5fb4c80890f31ff26e1748693",
"classsgpp_1_1base_1_1OperationMultipleEvalNakBsplineModifiedNaive.html#a4729d912574ae29ef1d8668a260b624d",
"classsgpp_1_1base_1_1OperationQuadratureMC.html#a1ef6ab42d522f797f42ad40bbd3f5f05",
"classsgpp_1_1base_1_1OperationSecondMomentLinearBoundary.html#ad6060e40a9bf10a43661840b71cb1697",
"classsgpp_1_1base_1_1OperationWeightedQuadratureNakBsplineModified.html#ac529c49bb7274fb6bb9f473cca45fe98",
"classsgpp_1_1base_1_1PolyGrid.html",
"classsgpp_1_1base_1_1Printer.html#a95f7da9d1d317f3908d9129d080d6183",
"classsgpp_1_1base_1_1ScaledScalarFunctionGradient.html#a5fd73130b08e40ef0cd714f71af0ea20",
"classsgpp_1_1base_1_1StencilHierarchisationModLinear.html#a7450537d0631c72834c90fe87d486e2a",
"classsgpp_1_1base_1_1VectorFunctionGradient.html#a67db96f358bc4776a74912e0f9f6f091",
"classsgpp_1_1base_1_1WeaklyFundamentalNakSplineModifiedBasisDeriv2.html#a36e4e5ca9072eca56c2b3a69d4b8aa6f",
"classsgpp_1_1base_1_1not__implemented__exception.html#a922873059fe070f6278915a44d947ec7",
"classsgpp_1_1combigrid_1_1CombinationGrid.html#acaeec79668a7a8645061e7566b8309dd",
"classsgpp_1_1combigrid_1_1OperationPole.html#aaca3a8e6320786d8c8e511dfdabc6321",
"classsgpp_1_1datadriven_1_1AlgorithmAdaBoostBase.html#abe6ef4d7cabc8199a9684b40bef49a1d",
"classsgpp_1_1datadriven_1_1ClassificationRefinementFunctor.html#a0d59a0c4af398a2e7112e4e4c1597187",
"classsgpp_1_1datadriven_1_1DBMatOffline.html#a5feef6f8a489406f6ec2d295bd5d70fc",
"classsgpp_1_1datadriven_1_1DBMatOnlineDE.html#a5f51892c11a96dc0838402db47c60798",
"classsgpp_1_1datadriven_1_1DMSystemMatrixDRE.html#adfb9c8ea88fe2f74709f9870e87387f1",
"classsgpp_1_1datadriven_1_1DataShufflingFunctorCrossValidation.html",
"classsgpp_1_1datadriven_1_1DensityDifferenceEstimationMinerFactory.html#ae6caccaf2692654b08abbf11289c09db",
"classsgpp_1_1datadriven_1_1DiscreteParameter.html#a4bb0b1695f7351e1cc74506db42ddc5b",
"classsgpp_1_1datadriven_1_1GridFactory.html#a780539c2cc9fc90deb573fb80d0a15cb",
"classsgpp_1_1datadriven_1_1KernelDensityEstimator.html#aee85643bbd7e7f6a72ebb07a1dd39905",
"classsgpp_1_1datadriven_1_1LearnerSGDE.html#a041e6817818a21799b84abec8d5d40b8",
"classsgpp_1_1datadriven_1_1LearnerSVM.html#a5472dd46079341725d9dc9ac3d33a426",
"classsgpp_1_1datadriven_1_1ModelFittingBase.html#a5edebc687deb1021e147b708a9b83fcd",
"classsgpp_1_1datadriven_1_1ModelFittingDensityEstimationCG.html#af5e583b310145dd781552eeaa13bf4e4",
"classsgpp_1_1datadriven_1_1MultiSurplusRefinementFunctor.html#af008bea7f3fd2c78a71979c0c6319bbb",
"classsgpp_1_1datadriven_1_1OperationInverseRosenblattTransformation1DLinear.html#a0ab77446789d80d371bd7c4dec1532ef",
"classsgpp_1_1datadriven_1_1OperationInverseRosenblattTransformationPoly.html#a07474e2dafe30ee28600f4b788138626",
"classsgpp_1_1datadriven_1_1OperationMultiEvalCuda.html#a1a7e2dfbfdc9cf1cf0f18de1ef8bb811",
"classsgpp_1_1datadriven_1_1OperationMultiEvalStreamingModOCLFastMultiPlatform.html#a4e40265669a417ae91fc4b29e2d2164a",
"classsgpp_1_1datadriven_1_1OperationMultipleEvalDistributed.html#a1d43f3bdb82ad137012c803c5bf95158",
"classsgpp_1_1datadriven_1_1OperationRosenblattTransformation1DPolyClenshawCurtisBoundary.html",
"classsgpp_1_1datadriven_1_1OperationTest.html#ae15296bdaf7d20e6422440175f40d6b2",
"classsgpp_1_1datadriven_1_1PolynomialChaosExpansion.html#a46a9c80a828ca70ec382b07e5d71d3b1",
"classsgpp_1_1datadriven_1_1ScorerFactory.html",
"classsgpp_1_1datadriven_1_1StreamingBSplineOCLKernelSourceBuilder.html#a5280f7c3491be6c9e2f532e2add6e29f",
"classsgpp_1_1datadriven_1_1SubspaceNodeCombined.html#a14c8303c1fe020aa88888771fe9746fd",
"classsgpp_1_1datadriven_1_1VisualizerDensityEstimation.html#ab077c01c34e929b77c1684ed34a03a18",
"classsgpp_1_1datadriven_1_1clusteringmpi_1_1MPIWorkerPackageBase.html#ae9febee2ec6823409ef74c292e3b7a79",
"classsgpp_1_1optimization_1_1FuzzyExtensionPrincipleViaTransformation.html#a8fc52d21eb997f1101814ee726cf083f",
"classsgpp_1_1optimization_1_1IterativeGridGeneratorFullAdaptiveRitterNovak.html#aaad8c0fef60d9222e3a308b14165af93",
"classsgpp_1_1optimization_1_1OperationMultipleHierarchisationLinearBoundary.html#a22129edd29d12bc968a7a100aadb58b2",
"classsgpp_1_1optimization_1_1OperationMultipleHierarchisationWaveletBoundary.html#ab38e9f0400059a19c0fc63f8bfa193cb",
"classsgpp_1_1optimization_1_1TriangularFuzzyInterval.html#af8659629a3bcb6019d0b45afafedac12",
"classsgpp_1_1optimization_1_1optimizer_1_1DifferentialEvolution.html#a8e43959430d1c50beaf6ce29d128d52a",
"classsgpp_1_1optimization_1_1optimizer_1_1NLCG.html#a49bc84278659a522fe0b13c5ac5acec9",
"classsgpp_1_1optimization_1_1optimizer_1_1UnconstrainedOptimizer.html#ad948db088fc279fbc0e51dbcb2495b0d",
"classsgpp_1_1optimization_1_1test__problems_1_1Floudas.html#a6f6ca6ad79a5b4ab03bc03ea17dbeb80",
"classsgpp_1_1optimization_1_1test__problems_1_1G06.html#a067282d862ec47914e2b05ffda6ec950",
"classsgpp_1_1optimization_1_1test__problems_1_1G12.html#a043afd9b3356a37ef6164b1dab2f5647",
"classsgpp_1_1optimization_1_1test__problems_1_1IncreasingPowerObjective.html#a61ae35b9d7fc4d480bec85acbdf4e8e6",
"classsgpp_1_1optimization_1_1test__problems_1_1SimionescuObjective.html#aad9bdce8642927b5194d8782700b1b93",
"classsgpp_1_1pde_1_1HeatEquationParabolicPDESolverSystemParallelOMP.html#aefcd2635ba0b6c162e4070c637634cfb",
"classsgpp_1_1pde_1_1OperationEllipticPDESolverSystemDirichlet.html#a9916f636e46c1f83efcfe55b29513d4a",
"classsgpp_1_1pde_1_1OperationLaplaceModPoly.html#aef19ae301b9c5cf6219bbf32a0656fc1",
"classsgpp_1_1pde_1_1OperationMatrixLTwoDotExplicitPolyClenshawCurtisBoundary.html",
"classsgpp_1_1pde_1_1PhiPhiDownBBLinearBoundary.html#a02048b9864fec3936ea72567a15d23a7",
"classsgpp_1_1pde_1_1UpDownFourOpDims.html#a98fd2ccf4fdae3c4e50c46694edb40f5",
"classsgpp_1_1quadrature_1_1NaiveSampleGenerator.html#ae5d948dd1d670afc52e0007c802fdf86",
"classsgpp_1_1solver_1_1OperationParabolicPDESolverSystem.html#acf2801531d702b8f2d989e4493f4d90f",
"create__scripts_8py.html#a7d20e2108a25a24aecbad65826774d7d",
"dir_a45ba08a15e1109ebad7c1b4903995b5.html",
"examples_java.html#examples_java_module_datadriven",
"learner_8cpp.html#a8d4de72b320fa113a4cd4775d5ce2d58a0dc58ae2bdfe414fe21af0475304b851",
"namespacepython_1_1classifier.html#a2d735a4fc0c280699c64210a8369b163",
"namespacepython_1_1uq_1_1dists_1_1Dist.html",
"namespacepython_1_1uq_1_1quadrature_1_1bilinearform_1_1bilinear__form.html#a0e6b6e193baea2974787d735f3c341d5",
"namespacesgpp_1_1base.html#a64916dd79050cf0b9094d35bc94007bdaa4ef7dc8a562f95768dd2bfc4d97ec8f",
"namespacesgpp_1_1datadriven.html#ad622a0a106238178faf9e1477bbf2026aeeceef5d9dbc1e62a270af389f4e15b0",
"parabolasimple_8py.html#a477e6feb211373118632334b148e6308",
"structsgpp_1_1base_1_1AdaptivityConfiguration.html#adcdaf4e7774e0c4edca304d6b8d2f1f0",
"structsgpp_1_1datadriven_1_1LearnerConfiguration.html#afefbbc63415f8ea8f74a1dbcf579ebc1",
"structsgpp_1_1solver_1_1SLESolverConfiguration.html"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';