#pragma once
#ifndef STDAFX_H__INCLUDED
#define STDAFX_H__INCLUDED

#include <Halide.h>
using Halide::Image;

// objects/methods declared in the H will cause duplicate linker symbols, so just extern them here
// #include "../apps/support/image_io.h"
// image_io.h will be compiled/linked in stdafx.cpp 

// externs
template<typename T> Halide::Image<T> load(std::string filename);
template<typename T> void save(Halide::Image<T> im, std::string filename);

#endif STDAFX_H__INCLUDED