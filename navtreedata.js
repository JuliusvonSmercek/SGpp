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
"LearnedKnowledgeFormatter_8py.html",
"NaturalBsplineBasis_8hpp.html#a7387c9b142ebc897c68f5835d17af91f",
"OperationEvalModPolyClenshawCurtisNaive_8hpp.html",
"OperationLaplacePolyClenshawCurtisBoundary_8hpp.html",
"OperationMultipleEvalSubspaceSimple__multImpl_8cpp.html",
"OperationUPCombinationGrid_8cpp.html",
"ScalarFunctionGradient_8hpp.html",
"WaveletBoundaryGrid_8cpp.html",
"classjson_1_1IDNode.html#ad55c2121e38542eb34363e87db0ee1f3",
"classpython_1_1controller_1_1CheckpointController_1_1CheckpointController.html#a0cf70cc24e0cf23c65ea31fc372cc0f7",
"classpython_1_1learner_1_1LearnerBuilder_1_1LearnerBuilder.html#a7134f8c4378b16b5a1da7a85d3812314",
"classpython_1_1learner_1_1Types_1_1BorderTypes.html",
"classpython_1_1uq_1_1analysis_1_1asgc_1_1ASGCAnalysisBuilder_1_1ASGCAnalysisBuilder.html#a7a91bf4abaa592a78a04ae881bd3d48a",
"classpython_1_1uq_1_1dists_1_1Beta_1_1Beta.html#a6fa88262e2f22d101a0eeb77fa62fced",
"classpython_1_1uq_1_1dists_1_1LibAGFDist_1_1LibAGFDist.html#a625f5eb26eed1ea3633e56e537b3dbe4",
"classpython_1_1uq_1_1dists_1_1TLognormal_1_1TLognormal.html#a78880888e5785d1aefbd7db636056336",
"classpython_1_1uq_1_1learner_1_1Learner_1_1Learner.html#a716345fa05538b76a6fa3e51a2d3cd7a",
"classpython_1_1uq_1_1learner_1_1builder_1_1RegressorSpecificationDescriptor_1_1FoldingDescriptor.html#a6277831a30fbb69d537845ee4bfdeba7",
"classpython_1_1uq_1_1models_1_1Model_1_1Model.html#a7b9a24473469f5e6a2e833b4ecd06145",
"classpython_1_1uq_1_1operations_1_1forcePositivity_1_1localFullGridSearch_1_1LocalFullGridCandidates.html#aedac724cf973e5b7786da6df8a668804",
"classpython_1_1uq_1_1parameters_1_1ParameterBuilder_1_1ParameterBuilder.html#a8af05df8bad14448e08631d411f576d7",
"classpython_1_1uq_1_1quadrature_1_1HashQuadrature_1_1HashQuadratureMap.html#a727dd5a1cdad4dd13067179d950eb458",
"classpython_1_1uq_1_1refinement_1_1RefinementManagerDescriptor_1_1RefineCurrentNodesDescriptor.html#ae2784f0537a49568c30dac0662aa3211",
"classpython_1_1uq_1_1refinement_1_1RefinementStrategy_1_1VarianceOptRanking.html#a404ecf560cfd481126a645e3d53a169a",
"classpython_1_1uq_1_1sampler_1_1asgc_1_1ASGCSamplerSpecification_1_1ASGCSamplerSpecification.html#a8f83fafcc3995be9ad814afe7614daee",
"classpython_1_1uq_1_1uq__setting_1_1UQSettingManager_1_1UQSettingManager.html#a42560920305154f756306c21e394af14",
"classpython_1_1uq_1_1uq__setting_1_1samplingresult_1_1Samplingresult.html#aeb88379cc9df60dbb7b46fe1cd12bc8d",
"classsgpp_1_1base_1_1BoundingBox.html#a71d79ffa3ad0a43b490552f9342d50d6",
"classsgpp_1_1base_1_1BsplineModifiedClenshawCurtisBasis.html#a648bd7536ca361cfb2d3905cbbae6281",
"classsgpp_1_1base_1_1DataMatrix.html#a04caec26c94ab5951b48927af32fcddd",
"classsgpp_1_1base_1_1DataMatrixSP.html#ae6452d6b88051bbb613c734330598705",
"classsgpp_1_1base_1_1DehierarchisationFundamentalNakSplineBoundary.html#a027257cebb95fd992f3061595bdef0a3",
"classsgpp_1_1base_1_1DehierarchisationPolyBoundary.html#a51feb221a354735a2a19a739657da3f2",
"classsgpp_1_1base_1_1EmptyVectorFunctionGradient.html",
"classsgpp_1_1base_1_1FundamentalSplineModifiedBasis.html#a3dea5818e1bbca8a80cb19c4750539c4",
"classsgpp_1_1base_1_1Grid.html#afb4bdca27a4de421add1e69f5c6ac4bd",
"classsgpp_1_1base_1_1HashGridIterator.html#ac8446db1f8eae7ede8c48471aab4554a",
"classsgpp_1_1base_1_1HashRefinement.html#a86048f1c85489a20315f9e00182c6aab",
"classsgpp_1_1base_1_1HierarchisationModLinearClenshawCurtis.html#adc1b507e3183cdb9ca98cac4b5bf4e48",
"classsgpp_1_1base_1_1HierarchisationSLE.html#ae0e93ef7384013078534d55b180da9e0",
"classsgpp_1_1base_1_1L0BoundaryGridGenerator.html#a4ab4ee3e2ac04a1167b796f22003eb42",
"classsgpp_1_1base_1_1LinearModifiedClenshawCurtisBasis.html#afe85424c10eeff440667794c8fc7ebc8",
"classsgpp_1_1base_1_1ModPolyClenshawCurtisGrid.html#a4dfac27275f719c97413dfda8c3aef9a",
"classsgpp_1_1base_1_1NakBsplineExtendedBasis.html#aa0999de3062359d95516ea92ebee26d8",
"classsgpp_1_1base_1_1OCLClonedBuffer.html#ad62d2d917287b5b35faacdb991d9fd33",
"classsgpp_1_1base_1_1OperationConvertPrewavelet.html#a6d54bffa9a96f5ec0d2da13012dbeadd",
"classsgpp_1_1base_1_1OperationEvalGradientModFundamentalSplineNaive.html#a472e39dcc4f310ceb3a35ea8fa3ae53a",
"classsgpp_1_1base_1_1OperationEvalHessianBsplineBoundaryNaive.html",
"classsgpp_1_1base_1_1OperationEvalHessianWaveletNaive.html#a692c7d8856594e1b6459f874ce3b525d",
"classsgpp_1_1base_1_1OperationEvalModLinearNaive.html#ada9c8b575c6139b4539313e72405b990",
"classsgpp_1_1base_1_1OperationEvalPartialDerivativeFundamentalNakSplineNaive.html",
"classsgpp_1_1base_1_1OperationEvalPolyBoundary.html#a1829956b36119b872947a0b263986c30",
"classsgpp_1_1base_1_1OperationFirstMomentModPolyClenshawCurtis.html#a3f2cf4c8b71d4e78ccbcf87ce9039129",
"classsgpp_1_1base_1_1OperationHierarchisationModPolyClenshawCurtis.html",
"classsgpp_1_1base_1_1OperationMultipleEvalLinearBoundaryNaive.html",
"classsgpp_1_1base_1_1OperationMultipleEvalNakBsplineModifiedNaive.html#ac7352652f77407f6b10dbc97bd7b9789",
"classsgpp_1_1base_1_1OperationQuadratureMC.html#a80b41bf3d4b80a40ac960edb6a3c5124",
"classsgpp_1_1base_1_1OperationSecondMomentModBspline.html#ac29db4f243204876cd5417425cd6b2ec",
"classsgpp_1_1base_1_1OperationWeightedQuadratureNakPBspline.html#a0b8a71739ac21745831963c0741c4384",
"classsgpp_1_1base_1_1PolyGrid.html#a367c171870ea6d564ce56aaf438479b7",
"classsgpp_1_1base_1_1Printer.html#aa6312b8709d430b1c42260d6a00df72d",
"classsgpp_1_1base_1_1ScaledScalarFunctionGradient.html#a815c004d0d89c2d147e5723c039421d7",
"classsgpp_1_1base_1_1StencilHierarchisationModLinear.html#a96f92f4d59a037dd7b43a68b12166090",
"classsgpp_1_1base_1_1VectorFunctionGradient.html#aaad7910329b62fc4114336cd3b2819e6",
"classsgpp_1_1base_1_1WeaklyFundamentalNakSplineModifiedBasisDeriv2.html#aa5c3946e4ec2fb04489e0d2bf9aa1159",
"classsgpp_1_1base_1_1not__implemented__exception.html#af2787e382b2afc57bf50708ecebbfebb",
"classsgpp_1_1combigrid_1_1FullGrid.html#a40d9dc501f57e6961e679ed744efee88",
"classsgpp_1_1combigrid_1_1OperationPoleHierarchisationGeneral.html#ad8f2b9732ab7e8578893ba4cef4208c1",
"classsgpp_1_1datadriven_1_1AlgorithmAdaBoostBase.html#aeb79456e7366468f55ae23ee21bedcb8",
"classsgpp_1_1datadriven_1_1ClassificationRefinementFunctor.html#a7a497c507bdf4bf9bf4725b0b7265cbf",
"classsgpp_1_1datadriven_1_1DBMatOffline.html#a93f1fa6d0120518a11bf210d55e13dd5",
"classsgpp_1_1datadriven_1_1DBMatOnlineDE.html#a9ddc723483252b25d442542b869ffb2d",
"classsgpp_1_1datadriven_1_1DataBasedRefinementFunctor.html#a31522b553656e5d2667ad1d38aa7e6c0",
"classsgpp_1_1datadriven_1_1DataShufflingFunctorRandom.html#a8ccb185671adfcac377716c1d3beedd6",
"classsgpp_1_1datadriven_1_1DensityEstimationMinerFactory.html",
"classsgpp_1_1datadriven_1_1EpanechnikovKernel.html#ac779b7aebfc1b72fd9749370b5426882",
"classsgpp_1_1datadriven_1_1GridPointBasedRefinementFunctor.html#a8060d88fb353f474735d1fbe1138dd76",
"classsgpp_1_1datadriven_1_1LearnerBase.html#a1f17ac75f2ff093caa9e06da513d2339",
"classsgpp_1_1datadriven_1_1LearnerSGDE.html#a47b35752c6812c42f60335b9d26e3cb6",
"classsgpp_1_1datadriven_1_1LearnerSVM.html#ae560cdebb9c88a0c3e84ef6566d2d5b0",
"classsgpp_1_1datadriven_1_1ModelFittingBase.html#acd8492f30eab4b7c7e64a05886ee7eeb",
"classsgpp_1_1datadriven_1_1ModelFittingDensityEstimationCombi.html#a4530a545997019d11b971d1cada36f08",
"classsgpp_1_1datadriven_1_1NearestNeighbors.html#ab7079113315ed451fb728754838c367c",
"classsgpp_1_1datadriven_1_1OperationInverseRosenblattTransformation1DModBsplineClenshawCurtis.html#ac93cce4cdc194bd6b8b4445e009c5f5e",
"classsgpp_1_1datadriven_1_1OperationInverseRosenblattTransformationPolyBoundary.html#acd252243b3b4c66f9c64b76bce6844fa",
"classsgpp_1_1datadriven_1_1OperationMultiEvalCuda.html#a7dc5582303596e86087fa2611e85ff82",
"classsgpp_1_1datadriven_1_1OperationMultiEvalStreamingModOCLFastMultiPlatform.html#ab1cb1aca878684c85cc835d7ac101489",
"classsgpp_1_1datadriven_1_1OperationMultipleEvalLinearDistributed.html#a7c1ad5932e5807c8f866a9bbfd9f0d05",
"classsgpp_1_1datadriven_1_1OperationRosenblattTransformationBspline.html#ae995ee7eb67a6420d932f05cb3637b63",
"classsgpp_1_1datadriven_1_1OperationTestLinearBoundary.html#a62deeeb6104a5c6bd078e8e13de74207",
"classsgpp_1_1datadriven_1_1PrimalDualSVM.html#a537fa76fee81f756d9f6d277c38ceaff",
"classsgpp_1_1datadriven_1_1SortedDataset.html#a0c72afaa837fad12f5027834bd55fb10",
"classsgpp_1_1datadriven_1_1StreamingModOCLFastMultiPlatform_1_1KernelMultTranspose.html#a7069c3cdc530572453c011b022dbf8a8",
"classsgpp_1_1datadriven_1_1SubspaceNodeCombined.html#a1fad7f511e02aaf95b7a9d3b642af42cacb4fb1757fb37c43cded35d3eb857c43",
"classsgpp_1_1datadriven_1_1VisualizerDummy.html",
"classsgpp_1_1datadriven_1_1clusteringmpi_1_1OperationDummy.html#a29e588edb09323aa7b2752085de2166c",
"classsgpp_1_1optimization_1_1FuzzyInterval.html",
"classsgpp_1_1optimization_1_1OperationMultipleHierarchisation.html#a6978aef706a396d94d7f3163c5ad36b3",
"classsgpp_1_1optimization_1_1OperationMultipleHierarchisationModNakBspline.html#acbce044330963258d842bda1a1aa6c60",
"classsgpp_1_1optimization_1_1ResponseSurfaceVector.html#ad76c72af526c851fa2cf61088b07b974",
"classsgpp_1_1optimization_1_1optimizer_1_1AugmentedLagrangian.html#aa49b59471c6a90a3118aa3865881130e",
"classsgpp_1_1optimization_1_1optimizer_1_1LevenbergMarquardt.html#a223e9c30204161e662e44c0618ad0acb",
"classsgpp_1_1optimization_1_1optimizer_1_1Rprop.html",
"classsgpp_1_1optimization_1_1test__problems_1_1Branin01Objective.html#add49d28bd8deb93c716c8b8e779d7967",
"classsgpp_1_1optimization_1_1test__problems_1_1G04EqualityConstraint.html#a9b6be095671bc6dd21f9cf9af320261d",
"classsgpp_1_1optimization_1_1test__problems_1_1G09EqualityConstraint.html#a5401e0b574884fe49539b8f60da30373",
"classsgpp_1_1optimization_1_1test__problems_1_1Griewank.html#a1fa98ae7d2a093bf1a28d3318aa71c90",
"classsgpp_1_1optimization_1_1test__problems_1_1SHCB.html#a41b2e0a36a07fdeaebff8c895158609b",
"classsgpp_1_1optimization_1_1test__problems_1_1TremblingParabola.html#afe0cf9a97a1eb006a6dca898fba83b21",
"classsgpp_1_1pde_1_1LaplaceEnhancedDownBBLinear.html#adeed72253d996e4308cffc2c3b860c43",
"classsgpp_1_1pde_1_1OperationLaplaceExplicitLinear.html",
"classsgpp_1_1pde_1_1OperationMatrixLTwoDotExplicitLinear.html#a5877cb80dcc284e84cd986af25a4ab24",
"classsgpp_1_1pde_1_1OperationParabolicPDESolverSystemDirichlet.html#a31c5761fedf78a273b838932789445a4",
"classsgpp_1_1pde_1_1PoissonEquationSolver.html",
"classsgpp_1_1pde_1_1UpDownOneOpDimWithShadow.html#a719cb1d4845e07923feff67b2ce2b17f",
"classsgpp_1_1solver_1_1ConjugateGradientsSP.html#ad8470fa3787ddb48bc04803fb46ac5ca",
"classsgpp_1_1solver_1_1StepsizeControl.html#ac818df3073d1611018beee02084f161f",
"dir_04eb28bab8d6590cfa0a645247f6bc17.html",
"example__constrainedOptimization__cpp_8doxy.html",
"functions_vars_g.html",
"namespacefriedman2__4d.html#aed36ecec9f471d4a858a8bf1f3de26fe",
"namespacepython_1_1plotDeltas3d.html#aeb8d7b3815ed7a3a79ca4c815202ab17",
"namespacepython_1_1uq_1_1operations_1_1sparse__grid.html#a0ba5a417a607bc753fc42ee677278d8a",
"namespacepython_1_1utils_1_1data__projections.html#ad7db37480293bee7b012dfcd87077feb",
"namespacesgpp_1_1datadriven.html#a12c64c71b70f5fc05f0c0f95ecdb1eef",
"namespacesgpp_1_1optimization_1_1file__io.html#a739f9f0a72bb241b0208ba1c763056d3",
"sg__projections_8py.html#a0f8eb46d09255e7e1f3fbb933984b6b0",
"structsgpp_1_1datadriven_1_1CrossvalidationConfiguration.html#a9ccba328a0ad54770d963bcb92877863",
"structsgpp_1_1datadriven_1_1RefinementResultNetworkMessage.html#a77eeca26f1780bd280409fb64b33eea0"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';