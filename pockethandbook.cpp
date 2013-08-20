#include "pockethandbook.h"
#include <iostream>
#include <iomanip>
#include "stdafx.h"

namespace PocketHandbook { 
		
void TestSuite::TestGraylevelHistogram()
{
	GraylevelHistogram f;

	Halide::Image<uint8_t> input = load<uint8_t>("rgb.png");
	Halide::Image<float> histogram = f.Run(input);
	f.SaveHistogramToImage(histogram, "graylevelhistogram_x.png");
}

// GraylevelHistogram
Halide::Image<float> GraylevelHistogram::Run(Halide::Image<uint8_t> input)
{
	Halide::Var x, y, c;
	Halide::RDom r(0, input.width(), 0, input.height()); // reduction does not iterate over channels
	
	Halide::Func h; 					// histogram
	h(x,c) = 0; 						// initialize to zero
	h( input(r.x, r.y, c), c ) += 1; 	// increment each element
	
	// scale to 0..1 range
	int N = input.width() * input.height();
	Halide::Func hscale;
	hscale(x,c) = h(x,c)/(float)N; // 
	
	// compute the histogram
	Halide::Image<int> hist = h.realize( 256, 3 ); // every byte for every color channel
	Halide::Image<float> histscale = hscale.realize(256, 3); // realize as a scaled set of values
	
	if(_validate) {
		int N = 0; // check count!
		std::cout << "Graylevel Histogram (channel 0 only)" << std::endl;
		for(int i=0; i<hist.width(); i++) {
			int r = hist(i,0);
			int g = hist(i,1);
			int b = hist(i, 2);
			
			float rf = histscale(i,0); 
			float gf = histscale(i,1);
			float bf = histscale(i,2);
			
			std::cout << i << "   " << r << " " << g << " " << b
				<< std::setprecision(5)
				<< " " << rf << " " << gf << " " << bf
				<< std::endl;
		
			N += r + g + b;
		}
	
		if(N != input.width() * input.height() * 3) {
			throw "Incorrect histogram calculation";
		}
	}	
	
	return histscale;
}

void GraylevelHistogram::SaveHistogramToImage(Halide::Image<float> histogram, std::string filename)
{	
	std::cout << "Generating Image..." << std::endl;
	
	if(2 != histogram.dimensions()) 
		throw "histogram dimensions must be 2";
		
	Halide::Var x, y, c;
	
	// output an image based on the histogram
	Halide::Func gen;
	int GEN_W = 256;
	int GEN_H = 1024;
	 // middle gray if not set... otherwise it's the actual color represented by the histogram
	gen(x,y,c) = Halide::select( histogram(x,c) > y/(float)GEN_H, 255, 0);
	
	Halide::Func gen_convert;
	gen_convert(x,y,c) = Halide::cast<uint8_t>( clamp(gen(x,y,c), 0, 255));
	
	Halide::Image<uint8_t> output = gen_convert.realize(GEN_W, GEN_H, 3); // all channels
	
	save(output, filename); 
}

} // namespace PocketHandbook