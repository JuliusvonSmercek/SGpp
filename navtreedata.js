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
        [ "RitterNovakExample.cpp", "example_RitterNovakExample_cpp.html", null ],
        [ "RitterNovakHyperparameterOptimisation.cpp", "example_RitterNovakHyperparameterOptimisation_cpp.html", null ],
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
"LaplaceEnhancedUpBBLinear_8hpp.html",
"NakBsplineModifiedBasis_8hpp.html",
"OperationEvalModLinearClenshawCurtisNaive_8hpp.html",
"OperationLaplaceModLinear_8hpp.html",
"OperationMultipleEvalSubspaceSimpleParameters_8hpp.html#a9d70d34494742a192a71f116b7cdf51c",
"OperationTestModWavelet_8cpp.html",
"Sample_8py.html",
"VisualizerDummy_8hpp.html",
"classjson_1_1IDNode.html#aa3fc0ea699281cf7a0e1b6cdb5213dfd",
"classjson_1_1json__exception.html#aedeb3d4a00245aa6d95dfa0adace911a",
"classpython_1_1learner_1_1LearnerBuilder_1_1LearnerBuilder.html#a45c9f5a5ca3770025a6a5f8d1fe34e6c",
"classpython_1_1learner_1_1TrainingStopPolicy_1_1TrainingStopPolicy.html#ace02169fa34a9f718f72cfc6371b7146",
"classpython_1_1uq_1_1analysis_1_1asgc_1_1ASGCAnalysisBuilder_1_1ASGCAnalysisBuilder.html#a028e0f54e2a64a5886bad7d1bed47ee5",
"classpython_1_1uq_1_1dists_1_1Beta_1_1Beta.html#a324ff0cafacab387d1cd49a9669adafa",
"classpython_1_1uq_1_1dists_1_1LibAGFDist_1_1LibAGFDist.html#a2db2e44bd5eb6ca6a98e411a6f5cf97b",
"classpython_1_1uq_1_1dists_1_1TLognormal_1_1TLognormal.html#a41bf47913190245a9600b15d8f2f9cc4",
"classpython_1_1uq_1_1learner_1_1Learner_1_1Learner.html#a5e08946a7cb14caf18064cacce88979e",
"classpython_1_1uq_1_1learner_1_1builder_1_1LearnerBuilder_1_1LearnerBuilder.html#ac9901a9414f83d807e100497b258b96e",
"classpython_1_1uq_1_1models_1_1AnalyticModels_1_1Parabola.html",
"classpython_1_1uq_1_1operations_1_1forcePositivity_1_1localFullGridSearch_1_1LocalFullGridCandidates.html#ada5e0fbb38b88e5210335d4fcaed44a4",
"classpython_1_1uq_1_1parameters_1_1ParameterBuilder_1_1GeneralParameterBuilder.html#aa1c3a1cd0d1e3f85cfe36535d7d51535",
"classpython_1_1uq_1_1quadrature_1_1HashQuadrature_1_1HashQuadrature.html#abd85e184b7701fa0e2ab591e602c84e5",
"classpython_1_1uq_1_1refinement_1_1RefinementManagerDescriptor_1_1RefineCurrentNodesDescriptor.html#aa0b3503839da6becbb351c48514d5a48",
"classpython_1_1uq_1_1refinement_1_1RefinementStrategy_1_1VarianceOptRanking.html",
"classpython_1_1uq_1_1sampler_1_1asgc_1_1ASGCSamplerSpecification_1_1ASGCSamplerSpecification.html#a58198be4260db5860229952317368d45",
"classpython_1_1uq_1_1uq__setting_1_1UQSettingManager_1_1UQSettingManager.html",
"classpython_1_1uq_1_1uq__setting_1_1samplingresult_1_1Samplingresult.html#aa78b51f4fc7707c95f3238d9bc81c86f",
"classsgpp_1_1base_1_1BoundingBox.html#a55d2cdce67130e9e72cc96d065dcb347",
"classsgpp_1_1base_1_1BsplineModifiedClenshawCurtisBasis.html#a45817dff68f2d96f1ad010c64b0b795f",
"classsgpp_1_1base_1_1DataMatrix.html",
"classsgpp_1_1base_1_1DataMatrixSP.html#add8fc3bdbe00f42672677875baf092d6",
"classsgpp_1_1base_1_1DataVectorSP.html#af3cea273a4cc453e0f1e874db165ca3b",
"classsgpp_1_1base_1_1DehierarchisationPoly.html#ae820a702a69b7b46282b24ef26c2765a",
"classsgpp_1_1base_1_1DistributionsVector.html#adcf5b0b692fe9936d50366666c2b2d76",
"classsgpp_1_1base_1_1FundamentalSplineModifiedBasis.html#a06c738c5ed28c0705e7f7c1781fdb02a",
"classsgpp_1_1base_1_1Grid.html#aecc463b7ec3c0316b6b6dffec2a88a67",
"classsgpp_1_1base_1_1HashGridIterator.html#a9ac85483d60fe89501ab959f3eb542c7",
"classsgpp_1_1base_1_1HashRefinement.html#a0f5cba4e555787f1525aeae025885d85",
"classsgpp_1_1base_1_1HierarchisationModLinearClenshawCurtis.html#a74e29c4b1aed04506fab2f60ea4a5a2e",
"classsgpp_1_1base_1_1HierarchisationSLE.html#ad005e6173e91f9ee004cd7d46ab14c97",
"classsgpp_1_1base_1_1KernelSourceBuilderBase.html#ad88c9d7548306122b7ce30700c008e46",
"classsgpp_1_1base_1_1LinearModifiedClenshawCurtisBasis.html#a54cfb35b1a087f5e09d14ce38ac9d4ef",
"classsgpp_1_1base_1_1ModNakBsplineGrid.html#ae5025e7927fc2e00ddedc8b2e0399e7a",
"classsgpp_1_1base_1_1NakBsplineExtendedBasis.html#a336663b27d1a2824dbc61ad5db004c6b",
"classsgpp_1_1base_1_1OCLClonedBuffer.html#aa7592b1477376c1f458077cab4c4a2bb",
"classsgpp_1_1base_1_1OperationConvert.html#ae43d9059c2a68205ef5afbb68dba97e3",
"classsgpp_1_1base_1_1OperationEvalGradientModBsplineNaive.html#ad3f018c5056822eefb6ed1faebcb0aab",
"classsgpp_1_1base_1_1OperationEvalHessian.html#a162917293b5bc2a6eedfda99b99a6282",
"classsgpp_1_1base_1_1OperationEvalHessianWaveletBoundaryNaive.html#acddd45f3e24be35d19138e440a50fcd6",
"classsgpp_1_1base_1_1OperationEvalModLinearNaive.html#a9a45bcc08dd96169f59e4e96698df2dc",
"classsgpp_1_1base_1_1OperationEvalPartialDerivativeBsplineNaive.html#a833b339f773ad59463d7fc1ba4d3d156",
"classsgpp_1_1base_1_1OperationEvalPolyBoundary.html#a051f0ae22166c4292fe8e12543973de7",
"classsgpp_1_1base_1_1OperationFirstMomentModPoly.html#a8e5740e65ccec83c4f380e53b00809c1",
"classsgpp_1_1base_1_1OperationHierarchisationModPoly.html#a673dc08ac2000b5c7eaa889d44780b48",
"classsgpp_1_1base_1_1OperationMultipleEvalLinearBoundary.html#a7d95233e6058bbd014011d1ef6a074b3",
"classsgpp_1_1base_1_1OperationMultipleEvalNakBsplineModifiedNaive.html#a3bdf475fd348bddd3e4e1bcf92ef3915",
"classsgpp_1_1base_1_1OperationQuadratureMC.html#a0f40c1ed80dc017a00a4b8e8d0a57ad7",
"classsgpp_1_1base_1_1OperationSecondMomentLinearBoundary.html#a3b8dc91e2bc4e195800d2302a571edab",
"classsgpp_1_1base_1_1OperationWeightedQuadratureNakBsplineModified.html#aa7f1de11046f448848e1c2edb3b216e5",
"classsgpp_1_1base_1_1PolyClenshawCurtisGrid.html#af4b81fbcf4a472c0e5af8777fe37921c",
"classsgpp_1_1base_1_1Printer.html#a92f69c423f4035a16ae52d537c7b7812",
"classsgpp_1_1base_1_1ScaledScalarFunctionGradient.html#a4bb3e46dc5285474c3a950e783f82a33",
"classsgpp_1_1base_1_1StencilHierarchisationModLinear.html#a61d5c5c1e09a240f408a04d474e40c78",
"classsgpp_1_1base_1_1VectorFunctionGradient.html#a4a331ec3bde9cf57790a3fd41f89ff22",
"classsgpp_1_1base_1_1WeaklyFundamentalNakSplineModifiedBasisDeriv2.html#a367eae45a4dc29dc6864cb7bd1fa093b",
"classsgpp_1_1base_1_1not__implemented__exception.html#a1f9242274d25c4a14f74173325954886",
"classsgpp_1_1combigrid_1_1CombinationGrid.html#aa5fc35b885304646fc816d3f87ebd66f",
"classsgpp_1_1combigrid_1_1OperationPole.html#a3a7adf5d16233bbb5bea6ad3652c9261",
"classsgpp_1_1datadriven_1_1AlgorithmAdaBoostBase.html#abbb8c326984c2bc14f444e6128dad8f4",
"classsgpp_1_1datadriven_1_1ClassificationRefinementFunctor.html#a0a2d168fc118a5660c2c3ae087d95a70",
"classsgpp_1_1datadriven_1_1DBMatOffline.html#a5e86f9e29e937d0ba907ac8ba7cc5be2",
"classsgpp_1_1datadriven_1_1DBMatOnlineDE.html#a55779fe7691bf67acbfd4a6b7c8b3eb7",
"classsgpp_1_1datadriven_1_1DMSystemMatrixDRE.html#ad1edd9dd7631b01e5dbd112cb87fdfa4",
"classsgpp_1_1datadriven_1_1DataShufflingFunctor.html#aff3d626249d410c9595a139c781173bd",
"classsgpp_1_1datadriven_1_1DensityDifferenceEstimationMinerFactory.html#a7163492911c741ba01b7045f9b123ae4",
"classsgpp_1_1datadriven_1_1DiscreteParameter.html#a4445e51781344317c37f682c024dde86",
"classsgpp_1_1datadriven_1_1GridFactory.html#a396b9392d3ff23ea10d363502b4c5863",
"classsgpp_1_1datadriven_1_1KernelDensityEstimator.html#adfd48e29004a3ea72ff63b30c6f959d3",
"classsgpp_1_1datadriven_1_1LearnerSGDE.html",
"classsgpp_1_1datadriven_1_1LearnerSVM.html#a420d1b1504d81387a65f8d5ba3532eb4",
"classsgpp_1_1datadriven_1_1ModelFittingBase.html#a5a9d736cb616dd6e52064d0e7cbc7e76",
"classsgpp_1_1datadriven_1_1ModelFittingDensityEstimationCG.html#abbe652acdb7e8694040d3993056f80b0",
"classsgpp_1_1datadriven_1_1MultiSurplusRefinementFunctor.html#aeb78c7e20b34dd13e690c876caa5af09",
"classsgpp_1_1datadriven_1_1OperationInverseRosenblattTransformation1DLinear.html",
"classsgpp_1_1datadriven_1_1OperationInverseRosenblattTransformationPoly.html#a0453e6b887d571f8b6d5a582c57e91f5",
"classsgpp_1_1datadriven_1_1OperationMultiEvalCuda.html#a026ce575f39721d514818b82a85a12a0",
"classsgpp_1_1datadriven_1_1OperationMultiEvalStreamingModOCLFastMultiPlatform.html#a3ff89d6ba0f853b4b38f10a77a8c031c",
"classsgpp_1_1datadriven_1_1OperationMultipleEvalDistributed.html#a18ebfd292c0356f14a4a7cbc2abe20e8",
"classsgpp_1_1datadriven_1_1OperationRosenblattTransformation1DPolyClenshawCurtis.html#ab3758e25a0602e46ddeeb89282680aec",
"classsgpp_1_1datadriven_1_1OperationTest.html#ad9ac0ece67f07dbda4727ae1c26dcb4a",
"classsgpp_1_1datadriven_1_1PolynomialChaosExpansion.html#a3c84a2ff50e90a8e75574281cdbcf3bb",
"classsgpp_1_1datadriven_1_1Scorer.html#ad2e7685ae34496faf05ebacd011c56fe",
"classsgpp_1_1datadriven_1_1StreamingBSplineOCLKernelSourceBuilder.html#a11c53447099af0f991a397af65cdc070",
"classsgpp_1_1datadriven_1_1SubspaceNodeCombined.html#a14c8303c1fe020aa88888771fe9746fd",
"classsgpp_1_1datadriven_1_1VisualizerDensityEstimation.html#aa545a4b91799a00fe8a08d433d756c6b",
"classsgpp_1_1datadriven_1_1clusteringmpi_1_1MPIWorkerPackageBase.html#acb35630ef716f71eccf84ac4d0f41e87",
"classsgpp_1_1optimization_1_1FuzzyExtensionPrincipleViaTransformation.html#a841935740772d709073e2cb7fef56d44",
"classsgpp_1_1optimization_1_1IterativeGridGeneratorFullAdaptiveRitterNovak.html#aa7290fbd5d15f0ab9344d65ddee43584",
"classsgpp_1_1optimization_1_1OperationMultipleHierarchisationLinear.html#af5755c98a110b9105a379dda1defd8eb",
"classsgpp_1_1optimization_1_1OperationMultipleHierarchisationWaveletBoundary.html#a5a26522662ced678a90a6ea90435d2c8",
"classsgpp_1_1optimization_1_1TriangularFuzzyInterval.html#adf0d102b3ad15947783055e01fb4bb7b",
"classsgpp_1_1optimization_1_1optimizer_1_1DifferentialEvolution.html#a7ec9f99ab41450d2accac497d27c0bde",
"classsgpp_1_1optimization_1_1optimizer_1_1NLCG.html#a2e830cd8bb2818553241d5afd9717e75",
"classsgpp_1_1optimization_1_1optimizer_1_1UnconstrainedOptimizer.html#a9b3d9001e047087e3f953a70c0a16fda",
"classsgpp_1_1optimization_1_1test__problems_1_1Floudas.html#a40790bf8b050751e73605c751d614240",
"classsgpp_1_1optimization_1_1test__problems_1_1G05Objective.html#acac90ce1e8f6a942c5e2bbfa8a10ee45",
"classsgpp_1_1optimization_1_1test__problems_1_1G11Objective.html#ade6ecba939c514e8466925d6802dfc3c",
"classsgpp_1_1optimization_1_1test__problems_1_1IncreasingPowerObjective.html#a00d23c79a4b8a4b7f15153e6c53f88f3",
"classsgpp_1_1optimization_1_1test__problems_1_1SimionescuObjective.html#aa033a6075d943657b01675047aee22ea",
"classsgpp_1_1pde_1_1HeatEquationParabolicPDESolverSystemParallelOMP.html#adc5d253a0a56badd157d57a63bee86bb",
"classsgpp_1_1pde_1_1OperationEllipticPDESolverSystemDirichlet.html#a7d840dc4272115801585d4e1c41fdd2c",
"classsgpp_1_1pde_1_1OperationLaplaceModPoly.html",
"classsgpp_1_1pde_1_1OperationMatrixLTwoDotExplicitPolyClenshawCurtis.html#af5e0682ae03b0f558fddcfa86c9c32ef",
"classsgpp_1_1pde_1_1PhiPhiDownBBLinear.html#afcbd4a71ef55c45a33bf845a95326fd1",
"classsgpp_1_1pde_1_1UpDownFourOpDims.html#a8f2ebb69b88487dc8c04c2d5ab8594ad",
"classsgpp_1_1quadrature_1_1NaiveSampleGenerator.html#a1c5f63a969e3baa42ceab2f7f5481b8a",
"classsgpp_1_1solver_1_1OperationParabolicPDESolverSystem.html#aba561c5847117139e0ff984d355207b7",
"create__scripts_8py.html",
"dir_a36cb6e10a5fc103584463a67ba88be8.html",
"examples_cpp.html#examples_cpp_module_base",
"learner_8cpp.html",
"namespacepython_1_1classifier.html#a034d01c6ea70a746dcec6c971959657e",
"namespacepython_1_1uq_1_1analysis_1_1mc.html",
"namespacepython_1_1uq_1_1quadrature_1_1bilinearform_1_1BilinearGaussQuadratureStrategy.html",
"namespacesgpp_1_1base.html#a62c4aae76549a46278a02f2d19a091e9",
"namespacesgpp_1_1datadriven.html#ad622a0a106238178faf9e1477bbf2026a4c1a0de6eb38af2b9075c30b1ef2bfa4",
"parabola_8py.html#a751add509716188cd479c52c2b0bc8d1",
"structsgpp_1_1base_1_1AdaptivityConfiguration.html#ac34fe1e3913a826e645ee7315e8c7ea8",
"structsgpp_1_1datadriven_1_1HashGridPointCompare.html",
"structsgpp_1_1datadriven_1_1VisualizationParameters.html#ad97f12fb2f4eec7203a96bdc409c4e72"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';