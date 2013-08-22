#include "pockethandbook.h"
#include <iostream>
#include <iomanip>
#include "stdafx.h"
#include "test.h"

namespace PocketHandbook { 
		
void TestSuite::TestGraylevelHistogram()
{
	GraylevelHistogram f;

	Halide::Image<uint8_t> input = load<uint8_t>("rgb.png");
	Halide::Image<float> histogram = f.Run(input);
	f.SaveHistogramToImage(histogram, "graylevelhistogram_x.png");
}

void TestSuite::TestBrightness()
{
	Brightness f;
	
	Halide::Image<uint8_t> input = load<uint8_t>("rgb.png");
	Halide::Image<uint8_t> brightened = f.Run(input, 50);
	Halide::Image<uint8_t> darkened = f.Run(input, -50);
	
	save(brightened, "brightness_plus50.png");
	save(darkened, "brightness_minus50.png");
}

void TestSuite::TestBinaryDilationFilter()
{
	Halide::Image<uint8_t> input(5,5,3);
	input(2,2,0) = 10;
	input(0,0,0) = 20;
	
	BinaryDilationFilter bdf;
	Halide::Image<uint8_t> output = bdf.Run(input);
	
	// debug/diagnose the output
	for(int i=0; i<5; i++) {
		for (int j=0; j<5; j++) {
			
			char buf[256];
			sprintf(buf, "%.2d,%.2d", i, j);
			
			uint8_t r = output(i,j,0);
			std::cout << (int)r << " ";
			/*
			if(i==0 || i==4 || j==0 || j==4)
				assertEqual(0, output(i,j,0), buf);
			else
				assertEqual(255, output(i,j,0), buf);
			*/
		}
		std::cout << std::endl;
	}
	
	// the interesting point is [1,1] where both the 10 and 20 are available. binary dilation
	// should choose the larger value, 20
	assertEqual_u8(20, output(1,1,0), "competing value chooses larger");
	assertEqual_u8(0, output(4,4,0), "unneighbored cell doesn't change value");
	assertEqual_u8(20, output(0,0,0), "seed doesn't change");
	assertEqual_u8(20, output(0,1,0), "grow down");
	assertEqual_u8(20, output(1,0,0), "grow right");
}

void TestSuite::RunBinaryDilationOnRGBImage()
{
	Halide::Image<uint8_t> input = load<uint8_t>("binary_dilation_input.png");
	BinaryDilationFilter bdf;
	Halide::Image<uint8_t> output = bdf.Run(input);
	save( output, "binary_dilation_output.png" );
}

void TestSuite::TestFlip()
{
	std::cout << "TestFlip" << std::endl;
	
	Halide::Image<uint8_t> input(3,3,3);
	input(0,0,0) = 1;
	
	Flip f;
	Halide::Image<uint8_t> ud = f.Run(input, true);
	Halide::Image<uint8_t> lr = f.Run(input, false);
	
	assertEqual_u8(0, ud(0,0,0), "updown moves the source pixel");
	assertEqual_u8(1, ud(0,2,0), "pixel flips to bottom in up-down");
	assertEqual_u8(0, lr(0,0,0), "left-right moves source pixel");
	assertEqual_u8(1, lr(2,0,0), "pixel flips to top-right in left-right flip");
}

void TestSuite::RunFlipExample()
{
	Halide::Image<uint8_t> input = load<uint8_t>("rgb.png");
	Flip flip;
	
	save( flip.Run(input, true), "_flip_up_down.png");
	save( flip.Run(input, false), "_flip_left_right.png");
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

//
// Brightness
//
Halide::Image<uint8_t> Brightness::Run(Halide::Image<uint8_t> input, int brightness)
{
	Halide::Var x, y, c;
	
	Halide::Func f;
	f(x,y) = Halide::cast<uint8_t>(Halide::clamp(Halide::cast<int>(input(x,y))+brightness,0,255));
	
	return f.realize( input.width(), input.height(), input.channels() );
}

//
// BinaryDilationFilter
//
Halide::Image<uint8_t> BinaryDilationFilter::Run(Halide::Image<uint8_t> in)
{
	Halide::Var x, y, c;
	
	// bug: select with <uint8_t> arguments is failing, so everything here will be done in <int>

	// "circular" structuring function, used as a BINARY mask. Any mask element
	// that is non-zero will be used to select the maximum value from neighboring pixels
	Halide::Image<int> mask(3,3); 
	mask(0,0) = 1;
	mask(0,1) = 1;
	mask(0,2) = 1;
	
	mask(1,0) = 1;
	mask(1,1) = 1;
	mask(1,2) = 1;
	
	mask(2,0) = 1; 	
	mask(2,1) = 1;
	mask(2,2) = 1;
	Halide::RDom r(mask);
	
	// clamp input sampling
	Halide::Func input("input");
	int W = in.width();
	int H = in.height();
	input(x,y,c) = in(clamp(x,0,W-1), clamp(y,0,H-1), c);
	
	// main function 
	Halide::Func f("f");
	
	// works okay, but mask should really be a select(), not a multiply
	//f(x,y,c) = 0;
	//f(x,y,c) = Halide::max(mask(r.x,r.y) * input_clamp( x+r.x-1, y+r.y-1, c ), f(x,y,c));
	
	// NOTE:
	// select() with uint8_t doesn't work quite right, so cast to integer. final step g() will
	// cast this back to uint8_t
	f(x,y,c) = 0;
	f(x,y,c) = Halide::max(
		select(
			mask(r.x,r.y) > 0,
			Halide::cast<int>(input( x+r.x-1, y+r.y-1, c )), 
			0), 
		f(x,y,c)
	);

	// convert back to uint8_t
	Halide::Func g("g");
	g(x,y,c) = Halide::cast<uint8_t>( clamp(f(x,y,c), 0, 255) );
	
	// create the output image
	return g.realize(in.width(), in.height(), in.channels());
}

// FLIP
Halide::Image<uint8_t> Flip::Run(Halide::Image<uint8_t> input, bool upDown)
{
	Halide::Var x, y, c;
	Halide::Func f;
	int W = input.width();
	int H = input.height();
	
	if(upDown) {
		f(x,y,c) = input(x,H-y-1,c);
	}
	else {
		f(x,y,c) = input(W-x-1,y,c);
	}
	
	return f.realize(input.width(), input.height(), input.channels());
}

} // namespace PocketHandbook