// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearModifiedClenshawCurtisBasis.hpp>
#include <sgpp/base/grid/generation/StandardGridGenerator.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * Grid with modified Clenshaw-Curtis poly basis functions
 */
class ModLinearClenshawCurtisGrid : public Grid {
 protected:
  /**
   * This constructor creates a new GridStorage out of the stream.
   *
   * @param istr inputstream that contains the grid information
   */
  explicit ModLinearClenshawCurtisGrid(std::istream& istr);

 public:
  /**
   * Constructor of grid with Clenshaw-Curtis poly basis functions
   *
   * @param dim the dimension of the grid
   */
  explicit ModLinearClenshawCurtisGrid(size_t dim);

  /**
   * Destructor.
   */
  ~ModLinearClenshawCurtisGrid() override;

  /**
   * @return string that identifies the grid type uniquely
   */
  sgpp::base::GridType getType() override;

  /**
   * @return B-spline basis
   */
  SBasis& getBasis() override;

  /**
   * @return pointer to a GridGenerator object
   */
  GridGenerator& getGenerator() override;

  /**
   * reads a grid out of a string
   *
   * @param istr string that contains the grid information
   * @return grid
   */
  static Grid* unserialize(std::istream& istr);

  /**
   * Serializes the grid.
   *
   * @param ostr stream to which the grid is written
   * @param version the serialization version of the file
   */
  void serialize(std::ostream& ostr, int version = SERIALIZATION_VERSION) override;

 protected:
  /// grid generator
  StandardGridGenerator generator;
  /// linear basis
  SLinearModifiedClenshawCurtisBase basis_;
};

}  // namespace base
}  // namespace sgpp
