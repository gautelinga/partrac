#ifndef __TRANSPORT_ELEMENT_HPP
#define __TRANSPORT_ELEMENT_HPP

#include <string>

// Carried per particle: Point, Vector (rho-hat, w, S) or Tensor (F)
enum class TransportElement { Point, Vector, Tensor };

inline const char* to_string(const TransportElement e){
  switch (e){
    case TransportElement::Point: return "point";
    case TransportElement::Vector: return "vector";
    case TransportElement::Tensor: return "tensor";
  }
  return "?";
}

#endif
