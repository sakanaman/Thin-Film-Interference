#ifndef RAY_HPP
#define RAY_HPP

#include "api.hpp"

struct Ray
{
    Vec o, d;
    Ray(Vec o_, Vec d_) : o(o_), d(d_) {}
};

#endif