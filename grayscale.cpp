// sample c++ gradient halide tutorial

// Demonstrates how to create create grayscale from a color image
// Techniques covered:
// 1. image load
// 2. image save
// 3. in-place editing of pixel values

// On os x, compile with:
// g++ grayscale.cpp -I ../src/include -L ../../bin -lHalide `libpng-config --cflags --ldflags` -o grayscale
// DYLD_LIBRARY_PATH=../../bin ./grayscale

#include <Halide.h>

using Halide::Image;
#include "../apps/support/image_io.h"

void simpleFunc();
void simple1DFunc();
void simple2DFunc();
void funcWithExplicitExpression();
void playWithExpressions();
void changeDimensions();
void official_argmax();
//void generateCheckerboard(); FAILING
void setAndReadImage();
void manualCheckerboard();
void imageMathXorCheckerboard();
void swizzleColors();
void compositionOfTwo1DFuncs();

int _errors = 0;

void runOtherTests()
{
	_errors = 0;
	printf("__Running Tests__\n");
	
	official_argmax();

	simpleFunc();
	simple1DFunc();
	simple2DFunc();
	funcWithExplicitExpression();
	playWithExpressions();
	changeDimensions();
	//generateCheckerboard();
	setAndReadImage();
	manualCheckerboard();
	imageMathXorCheckerboard();
	swizzleColors();
	compositionOfTwo1DFuncs();
	
	if(_errors == 0) {
		printf("Success! (0 errors)\n");
	}
	else {
		printf("FAIL! %d error(s)\n", _errors);
	}
}

void assertEqual(int expected, int actual, const char* message) 
{
	if(expected != actual) {
		printf("ERROR -- expected (%d) actual (%d) -- %s\n", expected, actual, message);
		++_errors;
	}
	else {
		//printf("Pass... %s\n", message);
	}
}

// simple 0-dimension function, and raw access to the results
void simpleFunc()
{
	printf("simpleFunc\n");
	
	Halide::Func f;
	f() = 3;
	
	Halide::Buffer buf = f.realize();
	int* pBuf = (int*) buf.host_ptr(); // use host_ptr, NOT raw_buffer
	int output = *pBuf; // there will be a buffer allocated for the result, even when there are 0-dimensions
	
	assertEqual(3, output, "simpleFunc output should be 3");
}

// Show a 1-dimensional function
void simple1DFunc()
{
	printf("simple1DFunc\n");
	
	Halide::Func f;
	Halide::Var x;
	f(x) = 2*x+1; // function of 1 variable requires the definition of a ::Var
	
	Halide::Buffer buf = f.realize(10);
	
	int* pBuf = (int*) buf.host_ptr(); // do NOT use ::raw_buffer
	
	for(int i=0; i<10; i++) {
		assertEqual(2*i+1, pBuf[i], "1dfunc");
	}
}

// Show a 2-dimensional function, wrapped in an ::Image for easy access to the
// results
void simple2DFunc()
{
	printf("simple2DFunc\n");
	
	Halide::Func f;
	Halide::Var x, y;
	f(x,y) = x+y;
	
	Halide::Buffer buf = f.realize(10,10);
	Halide::Image<int> img = buf; // cast to ::Image for easy access to data
	
	for(int i=0; i<10; i++) {
		for(int j=0; j<10; j++) {
			assertEqual(i+j, img(i,j), "2dfunc");
		}
	}
}

void funcWithExplicitExpression()
{
	printf("funcWithExplicitExpression\n");
	
	Halide::Func f;
	Halide::Var x;
	Halide::Expr e;
	
	// instead of defining the expression in-place, show that the expression
	// can be defined first
	e = 2*x + 1;
	f(x) = e;
	
	Halide::Buffer buf = f.realize(10);

	int* pBuf = (int*) buf.host_ptr(); // do NOT use ::raw_buffer
	for(int i=0; i<10; i++) {
		assertEqual(2*i+1, pBuf[i], "funcWithExplicitExpression");
	}
}

void playWithExpressions()
{
	printf("playWithExpressions\n");
	
	Halide::Func f;
	Halide::Var x;
	Halide::Expr a, b;
	
	a = 2*x + 1;
	b = x*x;
	f(x) = a + b;
	
	Halide::Image<int> results = f.realize(10);
	
	for(int i=0; i<10; i++) {
		assertEqual(i*i + 2*i + 1, results(i), "playWithExpressions");
	}
}

void changeDimensions()
{
	printf("changeDimensions...\n");
	Halide::Func f, g;
	Halide::Var x;
	
	f(x) = Halide::Tuple(x, x*x);
	//f.compute_root();
	
	//Assertion failed: (buffers.size() == 1 && "Can only cast single-element realizations to buffers or images"), function operator Halide::Buffer, file ../src/include/Halide.h, line 5479.
	//Halide::Buffer buf( f.realize(10) );
	//Halide::Buffer buf = f.realize(10) );
	
	// so the function f can be realized, but the data can't be extracted... have
	// to push it back through a second function to get the data back
	/*
	Halide::Tuple t = f(x);
	g(x) = t[0] + t[1];
	Halide::Image<int> results = g.realize(10);
	*/
	Halide::Tuple t = f(x);
	g(x) = t[0] + t[1];
	Halide::Image<int> results = g.realize(10);
	
	for(int i=0; i<10; i++) {
		//int a = results(0,i);
		//int aa = results(1,i);
		//assertEqual(i, a, "changeDimensions--UP");
		//assertEqual(i*i, aa, "changeDimensions--UP-squared");
		
		int a = results(i);
		assertEqual(i*i+i, a, "g(f(x))");
	}
}

void setAndReadImage()
{
	printf("setAndReadImage\n");
	Image<uint8_t> input(16);
	
	input(0) = 127;
	uint8_t v = input(0);
	
	assertEqual(127, v, "setAndReadImage");
}

void manualCheckerboard()
{
	printf("manualCheckerboard\n");
	Image<uint8_t> input(16,16,3);
	int block = 4;
	
	for(int i=0; i<16; i++) {
		for (int j=0; j<16; j++) {
			bool lit = ((i/block) & 1) ^ ((j/block)&1) == 1;
			input(i,j,0) = lit ? 255 : 0;
			input(i,j,1) = (uint8_t)(i/16.0f * 255);
			input(i,j,2) = (uint8_t)(j/16.0f * 255);
			
			//input(i,j) = 0; // only sets the 0 component
			//input(i,j) = Halide::Tuple(lit?255:0, 127, 0);
		}
	}
	
	save(input, "manual_checkerboard.png");
}

void imageMathXorCheckerboard()
{
	// create two images by hand, then xor them together to get 
	// a checkerboard pattern
	Image<uint8_t> horiz_bars(16,16,3);
	Image<uint8_t> vert_bars(16,16,3);
	Halide::Var x, y, c;
	int block = 4;
	
	for(int i=0; i<16; i++) {
		bool vlit = (i/block)&1 == 1;
		for(int j=0; j<16; j++) {
			bool hlit = (j/block)&1 == 1;
			horiz_bars(i,j,0) = hlit ? 255 : 0;
			vert_bars(i,j,1) = vlit ? 255 : 0;
		}
	}
	
	save(horiz_bars, "horiz_bars.png");
	save(vert_bars, "vert_bars.png");
	
	Halide::Func f;
	f(x,y,c) = horiz_bars(x,y,c) ^ vert_bars(x,y,c); // function must handle all 3 dimensions
	Image<uint8_t> output = f.realize(16,16,3);
	save(output, "xor_bars.png");
	
	// TODO: read back in and verify the colors
}

void swizzleColors()
{
	// swizzle the colors of an image
	Image<uint8_t> input = load<uint8_t>("rgb.png");
	
	Halide::Func f;
	Halide::Var x, y, c;
	
	f(x,y,c) = input(x,y,(c+1)%3); // simple shift swizzle
	
	Image<uint8_t> output = f.realize( input.width(), input.height(), input.channels() );
	
	save(output, "swizzle.png");
}

// from OFFICIAL /test/correctness folder
void official_argmax()
{
	printf("official_argmax\n");
	
	Halide::Func f, arg_max_f;
	Halide::Var x;
	Halide::RDom r(0,100);
	
	f(x) = x * (100-x);
	
	arg_max_f() = 0;
	Halide::Expr best_so_far = f(clamp(arg_max_f(),0,100));
	arg_max_f() = select( f(r) > best_so_far, r, arg_max_f() );
	
	f.compute_root();
	Halide::Image<int> result_f = arg_max_f.realize();
	
	printf ("%d\n", result_f(0));
}

void generateCheckerboard()
{
	printf("generateCheckerboard\n");
	
	Halide::Func f;
	Halide::Var x, y;
	//Halide::Expr is_odd_x, is_odd_y;
	//is_odd_x(x) = x/4 & 1 == 1;
	//is_odd_y(y) = y/4 & 1 == 1;
	
	f(x,y) = 127;
	
	Halide::Image<uint8_t> results = f.realize(64,64);
	//save(results, "checker.png");
}

// make a 2D function out of 2 1d functions
void compositionOfTwo1DFuncs()
{
	printf("compositionOfTwo1DFuncs\n");
	
	Halide::Func f,g,h;
	Halide::Var x, y;
	
	g(x) = x*x;
	h(x) = x+x;
	f(x,y) = g(x) + h(y); // x*x + y+y
	
	char buf[256];
	
	//NOTE: this MUST be <int>, not <uint8_t> or the bytes will
	//come back 
	Halide::Image<int> img = f.realize(10,10);
	for(int i=0; i<10; i++) {
		for (int j=0; j<10; j++) {
			
			int v = i*i + (j+j);
			
			sprintf(buf, "compositionOfTwo1DFuncs(%d,%d)", i,j);
			assertEqual(v, img(i,j), buf);
		}
	}
}

int main(int argc, char **argv) {
	// load is from image_io.h
	Halide::Image<uint8_t> input = load<uint8_t>("rgb.png");
	Halide::Var x, y, c;
	Halide::Expr to3;

	// long form of the program
	Halide::Expr value = input(x,y,c);
	value = Halide::cast<float>(value);
	
	printf("Image Dimension: %d\n", input.dimensions());
	for(int i=0; i<input.dimensions(); i++) {
		printf("   extent of dimension %d is %d\n", i, input.extent(i));
		printf("   stride of dimension %d is %d\n", i, input.stride(i));
	}
	
	// color values of the top row
	/*
	for(int i=0; i<input.extent(0); i++) {
		printf("(%d,0): ", i);
		for(int j=0; j<3; j++) {
			uint8_t p8 = input(i,0,j);
			printf("%.2x ", p8);
		}
		printf("\n");
	}
	*/
	
	Halide::Expr e_red = input(x,y,0);
	Halide::Expr e_green = input(x,y,1);
	Halide::Expr e_blue = input(x,y,2);
	
	// color of top left pixel
	uint8_t red = input(0,0,0);
	uint8_t green = input(0,0,1);
	uint8_t blue = input(0,0,2);
	printf("Color of top left pixel is %d %d %d\n", red, green, blue);
	
	// program details
	value = 255 * (1 - value / 255.0f); // invert
	//value = 1.5f * value; // brighten
	
	/*
	Halide::Func gs;
	gs() = Halide::Tuple(255.0f, 0.0f, 0.0f); //value(0); //(value(x,y,0) + value(x,y,1) + value(x,y,2))/3.0f; // grayscale
	value = gs();
	*/
	
	// convert back to 8-bit integer
	value = Halide::min(value,255.0f);
	value = Halide::cast<uint8_t>(value);
	
	Halide::Func brighter;
	
	brighter(x,y,c) = value;
	
	// realize the function
	Halide::Image<uint8_t> output = brighter.realize(input.width(), input.height(), input.channels());
	
	// from image_io.h
	save(output, "rgb_brighter.png");
	
	runOtherTests();
	
	return 0;
}
