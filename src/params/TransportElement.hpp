#ifndef __TRANSPORT_ELEMENT_HPP
#define __TRANSPORT_ELEMENT_HPP

#include <string>
#include <iostream>

// Carried per particle: Point, Vector (rho-hat, w, S) or Tensor (F)
enum class TransportElement { Point, Vector, Tensor };

inline TransportElement parse_transport_element(const std::string& s){
  if (s == "point") return TransportElement::Point;
  if (s == "vector") return TransportElement::Vector;
  if (s == "tensor") return TransportElement::Tensor;
  std::cerr << "transport must be point, vector or tensor, not " << s << std::endl;
  exit(1);
}

inline const char* to_string(const TransportElement e){
  switch (e){
    case TransportElement::Point: return "point";
    case TransportElement::Vector: return "vector";
    case TransportElement::Tensor: return "tensor";
  }
  return "?";
}

#endif
