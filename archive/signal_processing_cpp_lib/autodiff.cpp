#include <iostream>
#include <cmath>
#include <vector>
#include <stdexcept>

// Forward Mode Automatic Differentiation using Dual Numbers
struct Dual {
  double value;
  double derivative;
};

Dual square(const Dual &x) {
  Dual result;
  result.value = x.value * x.value;
  result.derivative = 2 * x.value * x.derivative;
  return result;
}

Dual multiplyDual(const Dual &a, const Dual &b) {
  Dual result;
  result.value = a.value * b.value;
  result.derivative = a.value * b.derivative + a.derivative * b.value;
  return result;
}

Dual sinDual(const Dual &x) {
  Dual result;
  result.value = std::sin(x.value);
  result.derivative = std::cos(x.value) * x.derivative;
  return result;
}

double gradientDescentUpdate(double parameter, double gradient, double learningRate) {
  return parameter - learningRate * gradient;
}

std::vector<double> updateParameters(const std::vector<double>& parameters,
                                   const std::vector<double>& gradients,
                                   double learningRate) {
  if (parameters.size() != gradients.size()) {
      throw std::invalid_argument("Parameters and gradients vectors must have the same size.");
  }
  std::vector<double> updated;
  updated.reserve(parameters.size());
  for (std::size_t i = 0; i < parameters.size(); ++i) {
    updated.push_back(parameters[i] - learningRate * gradients[i]);
  }
  return updated;
}

// Reverse Mode Automatic Differentiation using a Computational Graph
struct Node {
  double value;
  double gradient; // adjoint value for reverse mode AD
  std::vector<Node*> children;
  std::vector<double> localDerivatives; // Store local derivatives w.r.t. children

  Node(double val) : value(val), gradient(0.0) {}
};

// Recursive backpropagation through the computational graph
void propagateGradient(Node* node) {
  for (std::size_t i = 0; i < node->children.size(); ++i) {
    node->children[i]->gradient += node->gradient * node->localDerivatives[i];
    propagateGradient(node->children[i]);
  }
}

int main() {
  // --- Forward Mode Automatic Differentiation Demonstration ---
  std::cout << "Forward Mode AD using Dual Numbers:" << std::endl;

  const double X_VAL = 3.0;
  Dual x;
  x.value = X_VAL;
  x.derivative = 1.0; // Seed derivative

  Dual xSquared = square(x);
  Dual f = sinDual(xSquared);

  std::cout << "f(x) = sin(x^2) evaluated at x = " << x.value << " is: " << f.value << std::endl;
  std::cout << "Derivative f'(x) computed with forward mode AD: " << f.derivative << std::endl;

  std::cout << "\nGradient Descent Parameter Update:" << std::endl;
  const double INITIAL_PARAM = 1.0;
  const double LEARNING_RATE = 0.01;
  double updatedParameter = gradientDescentUpdate(INITIAL_PARAM, f.derivative, LEARNING_RATE);
  std::cout << "Old parameter: " << INITIAL_PARAM << ", Updated parameter: " << updatedParameter << std::endl;

  // --- Reverse Mode Automatic Differentiation Demonstration ---
  std::cout << "\nReverse Mode AD using a Computational Graph:" << std::endl;

  Node* xNode = new Node(X_VAL);
  Node* squareNode = new Node(xNode->value * xNode->value);
  Node* sinNode = new Node(std::sin(squareNode->value));

  // Link nodes and set local derivatives
  squareNode->children.push_back(xNode);
  squareNode->localDerivatives.push_back(2 * xNode->value);

  sinNode->children.push_back(squareNode);
  sinNode->localDerivatives.push_back(std::cos(squareNode->value));

  sinNode->gradient = 1.0;

  propagateGradient(sinNode);

  std::cout << "f(x) = sin(x^2) evaluated at x = " << xNode->value << " is: " << sinNode->value << std::endl;
  std::cout << "Gradient computed via reverse mode AD (f'(x)): " << xNode->gradient << std::endl;

  delete sinNode;
  delete squareNode;
  delete xNode;

  return 0;
}
