#include "stdafx.h"

#include "../apps/support/image_io.h"

// method to force the inclusion of the specified functions
void force_to_include_types() 
{
	Halide::Image<uint8_t> img = load<uint8_t>("foo");
	save(img, "bar");
}