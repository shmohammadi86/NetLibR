#ifndef MIMUW_ADORATE_42_UTILS_H
#define MIMUW_ADORATE_42_UTILS_H

#include "constants.h"

#include <iostream>
#include <thread>

inline void log(const std::string msg)
{
    if (DEBUG_PARALLEL) {
        std::cout << "[" << std::this_thread::get_id() << "] " << msg << "\n";
    }
}

template<class C, class T>
inline auto contains(const C& v, const T& x)
-> decltype(end(v), true)
{
    return end(v) != std::find(begin(v), end(v), x);
}

#endif //MIMUW_ADORATE_42_UTILS_H
