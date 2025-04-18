// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/tools/json/IDNode.hpp>
#include <sgpp/base/tools/json/json_exception.hpp>

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace json {

IDNode::IDNode()
    : value(),
      //        internalType(InternalIDType::ID),
      isDouble(false),
      doubleValue(0.0),
      isUnsigned(false),
      unsignedValue(0),
      isSigned(false),
      signedValue(0),
      isBool(false),
      boolValue(false) {}

Node& IDNode::operator=(const Node& right) {
  const IDNode& idNode = dynamic_cast<const IDNode&>(right);
  this->operator=(idNode);
  return *this;
}

void IDNode::parse(std::vector<Token>& stream) {
  // create new text node
  if (stream[0].type == TokenType::ID) {
    this->value = stream[0].value;
    stream.erase(stream.begin());

    this->setupInternalType();
  } else {
    throw json_exception(stream[0], "expected id");
  }
}

void IDNode::setupInternalType() {
  if (this->value.compare("true") == 0) {
    this->isBool = true;
    this->boolValue = true;
  } else if (this->value.compare("false") == 0) {
    this->isBool = true;
    this->boolValue = false;
  }

  {
    std::stringstream conv_stream;
    conv_stream << this->value;
    uint64_t r;
    conv_stream >> r;
    if (!conv_stream.fail()) {
      this->isUnsigned = true;
      this->unsignedValue = r;
    }
  }
  {
    std::stringstream conv_stream;
    conv_stream << this->value;
    int64_t r;
    conv_stream >> r;
    if (!conv_stream.fail()) {
      this->isSigned = true;
      this->signedValue = r;
    }
  }
  {
    std::stringstream conv_stream;
    conv_stream << this->value;
    double r;
    conv_stream >> r;
    if (!conv_stream.fail()) {
      this->isDouble = true;
      this->doubleValue = r;
    }
  }
}

std::string& IDNode::get() { return this->value; }

void IDNode::set(const std::string& value) {
  this->value = value;

  this->setupInternalType();
}

double IDNode::getDouble() {
  //    if (this->internalType == InternalIDType::DOUBLE) {
  if (this->isDouble) {
    return this->doubleValue;
  } else {
    throw json_exception("node has not a numerical value");
  }
}

void IDNode::setDouble(double numericValue) {
  //    this->doubleValue = numericValue;
  //    this->internalType = InternalIDType::DOUBLE;
  std::stringstream stringstream;
  stringstream << numericValue;
  this->value = stringstream.str();
  this->setupInternalType();
}

uint64_t IDNode::getUInt() {
  //    if (this->internalType == InternalIDType::UINT) {
  if (this->isUnsigned) {
    return this->unsignedValue;
  } else {
    throw json_exception("node has not an unsigned integer value");
  }
}

void IDNode::setUInt(uint64_t uintValue) {
  //    this->unsignedValue = uintValue;
  //    this->internalType = InternalIDType::UINT;
  std::stringstream stringstream;
  stringstream << uintValue;
  this->value = stringstream.str();
  this->setupInternalType();
}

int64_t IDNode::getInt() {
  //    if (this->internalType == InternalIDType::INT) {
  if (this->isSigned) {
    return this->signedValue;
  } else {
    throw json_exception("node has not an integer value");
  }
}

void IDNode::setInt(int64_t intValue) {
  //    this->unsignedValue = intValue;
  //    this->internalType = InternalIDType::INT;
  std::stringstream stringstream;
  stringstream << intValue;
  this->value = stringstream.str();
  this->setupInternalType();
}

bool IDNode::getBool() {
  //    if (this->internalType == InternalIDType::BOOL) {
  if (this->isBool) {
    return this->boolValue;
  } else {
    throw json_exception("node has not a bool value");
  }
}

void IDNode::setBool(bool boolValue) {
  //    this->boolValue = boolValue;
  //    this->internalType = InternalIDType::BOOL;
  if (boolValue) {
    this->value = std::string("true");
  } else {
    this->value = std::string("false");
  }

  this->setupInternalType();
}

void IDNode::serialize(std::ostream& outFile, size_t indentWidth) { outFile << this->value; }

size_t IDNode::size() { return 1; }

Node* IDNode::clone() {
  IDNode* newNode = new IDNode(*this);
  return newNode;
}

}  // namespace json
