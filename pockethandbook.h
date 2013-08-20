#pragma once
#ifndef POCKETHANDBOOK_H__INCLUDED
#define POCKETHANDBOOK_H__INCLUDED

#include <Halide.h>
#include <string>

namespace PocketHandbook {

// set of test suites for PocketHandbook algorithms
class TestSuite {
public:
	static void Run() {
		TestGraylevelHistogram();
	}
private:
	static void TestGraylevelHistogram();
};
	
// compute the "gray level histogram" of an image.. NOTE:
// in all cases, we'll assume RGB images instead of grayscale and just
// compute the "gray level" for each channel individually
class GraylevelHistogram
{
	bool _validate;
	
public:
	GraylevelHistogram()
		: _validate(false)
	{
	}
	
	// compute the histogram of the input 3 channel image and return a 256x3 output
	// results as the normalized histogram. 
	// Output Dimensions:
	// [0] = gray level
	// [1] = color channel
	// Values = sum of pixels at this gray level [0..255] divided by total number of pixels in the image
	Halide::Image<float> Run(Halide::Image<uint8_t> input);
	
	// convert the normalized histogram from Run into an image
	void SaveHistogramToImage(Halide::Image<float> histogram, std::string filename);
};

}

#endif //POCKETHANDBOOK_H__INCLUDED