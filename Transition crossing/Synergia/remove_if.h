#ifndef CHEF_BASIC_TOOLKIT_REMOVE_IF_H
#define CHEF_BASIC_TOOLKIT_REMOVE_IF_H

#include <algorithm>

namespace chef {

  template <typename FI, typename PRED>
  FI remove_if(FI b, FI e, PRED p) {
    b = std::find_if(b, e, p);
    if (b != e) {
      for (FI i = b; ++i != e; /* nothing */) {
        if (!p(*i)) {
          *b++ = std::move(*i);
        }
      }
    }
    return b;
  }
}

#endif
