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
void expressionsTakeParentVariables();
void expressionWithNestedScopeVariable();
void sobel();
void sobel_rdom();
void grayscale();
void rdom();
void rdom_simple();
void rdom_sum();
void fib();
//void histogram();
void dim2reduction();
void official_convolve();
void rdom_dimensionality();
void lambda_intro(); // first lambda test
void official_side_effects();
void extern_c();
void power_function(); // extern_c falloff image generator
void power_native(); // all-native power function
void checkerboard();
void colored_checkerboard();


int _errors = 0;

void runOtherTests()
{
	_errors = 0;
	printf("__Running Tests__\n");

	// newest first so they fail faster
	colored_checkerboard();
	
	#if 0 // disabled for development
	checkerboard();
	power_native();
	power_function();
	extern_c();
	official_side_effects();
	lambda_intro();
	rdom_dimensionality();
	official_convolve();
	sobel_rdom();
	dim2reduction();
	
	//histogram();
	fib();
	rdom_sum();
	rdom_simple();
	rdom();

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
	expressionsTakeParentVariables();
	expressionWithNestedScopeVariable();
	sobel();
	grayscale();
	#endif
	
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

// expressions should take variables from the parent function
void expressionsTakeParentVariables()
{
	printf("expressionsTakeParentVariables\n");
	
	Halide::Expr a, b;
	Halide::Var x, y;
	Halide::Func f;
	
	a = x*x;
	b = y+y;
	f(x,y) = a + b; // x*x + y+y = x^2 + 2y
	
	Halide::Image<int> img = f.realize(10,10);
	
	for(int i=0; i<img.width(); i++) {
		for (int j=0; j<img.height(); j++) {
			int v = i*i + j+j;
			assertEqual(v, img(i,j), "");
		}
	}
}

// expression of X in nested function takes local definition of X
void expressionWithNestedScopeVariable()
{
	printf("expressionWithNestedScopeVariable\n");
	
	Halide::Expr a, b;
	Halide::Var x, y;
	Halide::Func f, g,  h;
	
	a = x*x;
	b = x+x;
	
	g(x) = a;
	h(x) = b;
	f(x,y) = g(x) + h(y); // a(x) + b(y) = x*x + y+y = x^2 + 2y
	
	Halide::Image<int> img = f.realize(10,10);
	for(int i=0; i<img.width(); i++) {
		for (int j=0; j<img.height(); j++) {
			int v = i*i + 2*j;
			
			assertEqual(v, img(i,j), "");
		}
	}
} 

void grayscale()
{
	printf("grayscale\n");
	
	Halide::Image<uint8_t> input = load<uint8_t>("rgb.png");
	
	Halide::Func f;
	Halide::Var x, y, c;
	
	// each color component is the same
	// min() in case the values overflow
	// cast() to convert the floating point values back to U8 [0..255]
	f(x,y,c) = Halide::cast<uint8_t>(min( 
		0.299f * input(x,y,0) + // red -- converts to floating point
		0.587f * input(x,y,1) + // green
		0.114f * input(x,y,2) 	// blue
		, 255.0f)
	);
		
	Halide::Image<uint8_t> output = f.realize( input.width(), input.height(), input.channels());
	save(output, "grayscale.png");
}

// sobel filter 
void sobel()
{
	printf ("sobel\n");

	Halide::Image<uint8_t> input = load<uint8_t>("rgb.png");
	
	Halide::Var x, y, c;
	Halide::Func f;
	Halide::Expr vs, hs; // vertical sobel filter, horizontal sobel filter
	
	// this definition ensures that samples do not lie outside of the image domain
	// and will return a "clamped" value of the pixel
	Halide::Func in;
	int W = input.width();
	int H = input.height();
	in(x,y,c) = input(clamp(x,0,W-1), clamp(y,0,H-1), c);
	
	// horizontal sobel kernel
	hs = -1.0f * in(x-1,y-1,c) + /*0*/ 1.0f * in(x+1,y-1,c) +
	     -2.0f * in(x-1,y  ,c) + /*0*/ 2.0f * in(x+1,y  ,c) +
		 -1.0f * in(x-1,y+1,c) + /*0*/ 1.0f * in(x+1,y+1,c);
		
	// vertical sobel kernel
	vs = -1.0f * in(x-1,y-1,c) - 2.0f * in(x,y-1,c)   - 1.0f * in(x+1,y-1,c) +
		  /* 0 + 0 + 0 */
		  1.0f * in(x-1,y+1,c) + 2.0f * in(x,y+1,c)   + 1.0f * in(x+1,y+1,c);
	
	// compute the magnitude of the sobel filter (per component)
	f(x,y,c) = Halide::cast<uint8_t>( min( sqrt(vs*vs + hs*hs), 255.0f));
	
	// Compute the sobel filter
	Halide::Image<uint8_t> output = f.realize( W, H, input.channels() );
	save(output, "sobel.png");
}

// magnitude function that can be used by the pipeline
extern "C" float mag(float dx, float dy) {
	return (float)sqrt(dx*dx + dy*dy);
}
HalideExtern_2(float, mag, float, float);

void sobel_rdom()
{
	printf("sobel_rdom\n");
	
	// sobel filter using reduction domain
	Halide::Image<uint8_t> rgb = load<uint8_t>("rgb.png");
	int W = rgb.width();
	int H = rgb.height();
	int C = rgb.channels();
	
	// clamp edge values
	Halide::Func input("input");
	Halide::Var x("x"), y("y"), c("c");

	input(x,y,c) = rgb(clamp(x,0,W-1), clamp(y,0,H-1),c);
	
	// define the sobel kernels
	Halide::Image<float> hsi(3,3);
	hsi(0,0) = -1; hsi(1,0) = 0; hsi(2,0) = 1;
	hsi(0,1) = -2; hsi(1,1) = 0; hsi(2,1) = 2;
	hsi(0,2) = -1; hsi(1,2) = 0; hsi(2,2) = 1;
	
	Halide::Image<float> vsi(3,3);
	vsi(0,0) = -1; vsi(1,0) = -2; vsi(2,0) = -1;
	vsi(0,1) =  0; vsi(1,1) =  0; vsi(2,1) =  0;
	vsi(0,2) = +1; vsi(1,2) = +2; vsi(2,2) = +1;
	
	Halide::RDom hs( hsi );
	Halide::RDom vs( vsi );
	
	// Sample using reduction
	Halide::Func fh("fh"), fv("fv"), sobel("sobel"), f("f");
	
	// NOTE: Dimensionality of reduction  must match dimension of 'pure definition'
	//fh(x,y) += input(hs.x, hs.y) * input(x,y); // reduction over horizontal kernel
	fh(x,y,c) += hsi(hs.x, hs.y) * input(x + hs.x - 1, y + hs.y - 1, c); // offset to center the kernel
	fv(x,y,c) += vsi(vs.x, vs.y) * input(x + vs.x - 1, y + vs.y - 1, c); // reduction over vertical kernel
	
	// compute the magnitude.. actually this doesn't work
	/*
	Halide::Expr mag;
	mag(x,y) = Halide::sqrt(x*x + y*y);
	*/
	
	// compute sobel magnitude
	// NOTE: uses an "extern" function 'mag'
	sobel(x,y,c) = mag(fh(x,y,c), fv(x,y,c)); // sobel magnitude
	
	//sobel(x,y,c) = Halide::sqrt( fh(x,y,c)*fh(x,y,c) + fv(x,y,c)*fv(x,y,c) );
	
	// scale and convert back to u8
	f(x,y,c) = Halide::cast<uint8_t>( min(sobel(x,y,c), 255.0f)); 	// scale

	Halide::Image<uint8_t> output = f.realize( W, H, C );
	save(output, "sobel_rdom.png");
}

void rdom_dimensionality()
{
	printf("rdom_dimensionality\n");
	
	// originally couldn't convolve with a real image... so this test illustratest the dimensionality
	// problem of using a 2D reduction (RDom(x,y)) on a 3D image (x,y,c)
	Halide::Image<uint8_t> rgb = load<uint8_t>("rgb.png");
	int W = rgb.width();
	int H = rgb.height();
	
	Halide::Func input;
	Halide::Var x, y, c;
	input(x,y,c) = rgb(clamp(x,0,W-1), clamp(y,0,H-1), c);
	
	Halide::Image<uint8_t> tent(3,3);
	tent(0,0) = 1;
	tent(1,1) = 1;
	tent(2,2) = 1;
	
	for(int i=0; i<3; i++) {
		for (int  j=0; j<3; j++) {
			uint8_t r = tent(i,j);
			//uint8_t g = tent(i,j);
			//uint8_t b = tent(i,j);
			//printf("%d,%d,%d   ", r, g, b);
			printf("%d   ", r);
		}
		printf("\n");
	}
	
	Halide::RDom r(tent);
	Halide::Func f;
	
	f(x,y,c) += tent(r.x,r.y) * input(x+r.x-1,y+r.y-1,c);
	
	Halide::Image<uint8_t> output = f.realize( W, H, rgb.channels() );
	save(output, "rdom_dimensionality.png");
}

void official_convolve()
{
	printf("official_convolve\n");
	
    int W = 64*3, H = 64*3;

    Halide::Image<uint16_t> in(W, H);
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            in(x, y) = rand() & 0xff;
        }
    }


    Halide::Var x("x"), y("y");

    Halide::Image<uint16_t> tent(3, 3);
    tent(0, 0) = 1;
    tent(0, 1) = 2;
    tent(0, 2) = 1;
    tent(1, 0) = 2;
    tent(1, 1) = 4;
    tent(1, 2) = 2;
    tent(2, 0) = 1;
    tent(2, 1) = 2;
    tent(2, 2) = 1;
    
    Halide::Func input("input");
    input(x, y) = in(clamp(x, 0, W-1), clamp(y, 0, H-1));

    Halide::RDom r(tent);

    /* This iterates over r outermost. I.e. the for loop looks like:
     * for y:
     *   for x:
     *     blur1(x, y) = 0
     * for r.y: 
     *   for r.x: 
     *     for y: 
     *       for x: 
     *         blur1(x, y) += tent(r.x, r.y) * input(x + r.x - 1, y + r.y - 1)
     *     
     * In general, reductions iterate over the reduction domain outermost.
     */
    Halide::Func blur1("blur1");
    blur1(x, y) += tent(r.x, r.y) * input(x + r.x - 1, y + r.y - 1);
}

// reduction domains
void rdom()
{
	printf("rdom\n");
	
	//http://halide-lang.org/Halide/class_halide_1_1_r_dom.html#details
	Halide::Func f;
	Halide::RDom r(1,8);
	Halide::Var x;
	
	f(x) = x;
	f(r) = f(r-1) + f(r) + f(r+1);
	
	// NOTE !!!
	// sample from the first pass, f(x), and store back into the result. since
	// this is done from 1 to 8, any samples that have already been calculated such
	// as f(r-1) will be the previous 'reduced' values. any values that have not been
	// calculated such as f(r+1) will be the 'un-reduced' values (simple 'x'). this
	// looks like:
	/*
		0
		1  => 0 + 1+2 = 3
		2  => 3 + 2+3 = 8
		3  => 8 + 3+4 = 15
		4  => 15 + 4+5 = 24
		5  => 24 + 5+6 = 35
		6  => 35 + 6+7 = 48
		7  => 48 + 7+8 = 63
		8  => 63 + 8+9 = 80
		9
	*/
	
	Halide::Image<int> result = f.realize(10); 
	
	int prev = 0;
	assertEqual(0, result(0), "rdom[0]"); // edge condition, not affected by reduction
	for(int i=1; i<9; i++) {
		printf ("result[%d] = %d\n", i, result(i));
		
		int v = prev + i + (i+1);
		assertEqual(v, result(i), "rdom[i]");
		
		prev = v;
	}
	assertEqual(9, result(9), "rdom[9]"); // edge condition, not affected by reduction
}

// from the official Doxygen comments in Halide docs
void rdom_simple()
{
	printf("rdom_simple\n");
	
	Halide::Func f;
	Halide::Var x;
	Halide::RDom r(0,5); // 5 is the extent, max value is 0+5-1 = 4

	f(x) = x;
	f(r) = f(r) * 2;
	Halide::Image<int> result = f.realize(10);
	
	for(int i=0; i<5; i++) {
		printf("result[%d] = %d\n", i, result(i));
		assertEqual(i*2, result(i), "rdom_simple[...]");
	}
	for(int i=5; i<10; i++) {
		assertEqual(i, result(i), "rdom_simple[5..10]");
	}
}

void rdom_sum()
{
	printf("rdom_sum\n");
	
	// sum all numbers in the range
	Halide::Func f;
	Halide::Var x;
	Halide::RDom r(1,9); // start index, length (NOT start index, end index exclusive)
	
	f(x) = x;
	f(r) = f(r) + f(r-1); // apply over each element in place
	
	/*
	rdom_sum
	result[0] = 0
	result[1] = 1
	result[2] = 3
	result[3] = 6
	result[4] = 10
	result[5] = 15
	result[6] = 21
	result[7] = 28
	result[8] = 36
	result[9] = 45
	*/
		
	Halide::Image<int> result = f.realize(10);
	
	int sum = 0;
	for(int i=0; i<10; i++) {
		printf("result[%d] = %d\n", i, result(i));
		
		sum += i;
		int v = result(i);
		
		assertEqual(sum, v, "rdom_sum[...]");
	}
	assertEqual(45, sum, "rdom_sum final sum = (N-1) * N/2");
}

void fib()
{
	printf("fib\n");
	
	// from rdom doc (http://halide-lang.org/Halide/class_halide_1_1_r_dom.html#details)
	// compute fibonacci numbers
	Halide::Func f;
	Halide::RDom r(2,18); // first 20
	Halide::Var x;
	
	f(x) = 1;
	f(r) = f(r-1) + f(r-2);
	Halide::Image<int> result = f.realize(20);
	
	int a = 1, b = 1;
	for(int i=2; i<20; i++) {
		int c = a+b;
		a = b;
		b = c;
		assertEqual(c, result(i), "fib");
	}
}

#if 0
// can't get this to work...
void histogram()
{
	printf("histogram\n");
	
	// "scattering" from the official Doxygen comments

	// NOTE: UInt is a Halide constructor that returns a Type
	Halide::ImageParam input(Halide::UInt(8), 2); // 2D image of type 8-bit integer unsigned
	
	Halide::Func histogram;
	Halide::Var x;
	
	Halide::RDom r(input);
	histogram(x) = 0; // initialize
	histogram(input(r.x,r.y)) = histogram( input(r.x, r.y)) + 1; // scatter
	
	Halide::Image<int> result = histogram.realize(10,10);
	
}
#endif

// 2D reduction
void dim2reduction()
{
	printf("dim2reduction\n");
	
	Halide::ImageParam input(Halide::Int(8), 2); // 2D 8-bit signed integer
	
	input(0,0) = 1;
	input(0,1) = 2;
	input(1,0) = 3;
	input(1,1) = 4;
	
	Halide::Func f;
	Halide::Var x, y;
	
	f(x,y) = input(x,y);
	
	Halide::Image<int> result = f.realize(2,2);
	
	/*
	assertEqual(2, input.extent(0), "x_extent");
	assertEqual(2, input.extent(1), "y_extent");
	
	// grayscale.cpp:639: error: cannot convert ‘Halide::Expr’ to ‘int’ in initialization
	int v = input(0,0);
	
	for(int i=0; i<2; i++) {
		for (int j=0; j<2; j++) {
			//printf("input[%d,%d] = %d", i, j, input(i,j));
		}
	}
	*/
	
}

void lambda_intro()
{
	// first lambda test
	printf("lambda_intro\n");
	
	Halide::Var x,y;
	//Halide::Func f;
	//f(x,y) = x*y;
	
	// MUST BE <int>
	Halide::Image<int> output = lambda(x,y,x*y).realize(10,10);
	
	assertEqual(0, output(0,0), "0,0");
	assertEqual(1, output(1,1), "1,1");
	assertEqual(12, output(4,3), "4,3");
	
	for(int i=0; i<10; i++) {
		for(int j=0; j<10; j++) {
			assertEqual(i*j, output(i,j), "lambda_intro(i,j)");
		}
	}
}

// side_effects.cpp (Official from tests)
extern "C" int argmin(int x, float y) {
    static float minVal = 1e10;
    static int minIndex = 0;

    if (y < minVal) {
        minIndex = x;
        minVal = y;
    }

	//printf("extern_c_argmin(%d, %.3f) -> [%d] %.3f\n", x, y, minIndex, minVal);

    return minIndex;
}
HalideExtern_2(int, argmin, int, float);

void official_side_effects()
{
	printf("official_side_effects\n");
	
	    Halide::Var x, y;
	    Halide::Func f, g, h;

	    f(x) = Halide::sin(x/10.0f+17);

	    // Compute argmin of f over [-100..100]
	    Halide::RDom r(-100, 100);
	    g(x) = 0;
	    g(x) = argmin(r, f(r));

	    Halide::Image<int> result = g.realize(1);
	    int idx = result(0);
	    printf("sin(%d/10.0f+17) = %f\n", idx, sinf(idx/10.0f+17));

	    printf("Success!\n");
	
		assertEqual(-60, idx, "official_side_effects");
}

// extern_c method declaration
//
extern "C" int __extern_c_function(int x0, int y0, int x, int y)
{
	// compute distance from center
	int dx = x-x0;
	int dy = y-y0;
	return (int)sqrt( dx*dx + dy*dy );
}
HalideExtern_4(int, __extern_c_function, int, int, int, int);

void extern_c()
{
	printf("extern_c\n");
	
	Halide::Image<int> input(256,256);
	input(0,0) = 1;
	input(0,1) = 2;
	int x0 = 128;
	int y0 = 128;
	
	Halide::Func f; 
	Halide::Var x, y;
	f(x,y) = __extern_c_function(x0, y0, x, y);
	
	Halide::Image<int> output = f.realize(256,256);
	
	assertEqual(1, input(0,0), "0,0");
	assertEqual(0, output(x0,y0), "x0,y0");
	
	// check every pixel
	for(int i=0; i<256; i++) {
		for (int j=0; j<256; j++) {
			int expected = __extern_c_function(x0,y0,i,j); // call it directly, just to make sure Func is calling it
			assertEqual(expected, output(i,j), "i,j" );
		}
	}
}

extern "C" float __dist(int x0, int y0, int x, int y)
{
	float dx = (x-x0), 
		  dy = (y-y0);
	
	return (float)sqrt(dx*dx + dy*dy);
}

HalideExtern_4(float, __dist, int, int, int, int);

// Generate a power function falloff image (for exmaple for use in computer graphics)
void power_function() // extern_c falloff image generator
{
	printf("power_function\n");
	
	Halide::Var x, y, c;
	Halide::Func fpow, f;
	float k = 3.0f;
	
	fpow(x,y,c) = 255.0f * (1-Halide::pow(__dist(0, 0, x, y) / 256.0f, k ));
	
	f(x,y,c) = Halide::cast<uint8_t>( Halide::clamp(fpow(x,y,c), 0, 255.0f));
	
	Halide::Image<uint8_t> output = f.realize(256,256,3);
	
	//Halide::Image<uint8_t> output = Halide::cast<uint8_t>( min(255.0f, f.realize(256,256,3) ));
	save(output, "power_function.png");
}

// there actually *is* a distance function in halide, it's "hypot"
// offseting this, too, to create a power/falloff image centered at 0.5,0.5
void power_native()
{
	printf("power_native\n");
	
	Halide::Var x, y, c;
	Halide::Func fpow, f;
	float k = 3.0f;
	int x0 = 128, y0 = 128;
	
	fpow(x,y,c) = 255.0f * (1-Halide::pow(Halide::hypot(x-x0, y-y0) / 128.0f, k ));
	
	f(x,y,c) = Halide::cast<uint8_t>( Halide::clamp(fpow(x,y,c), 0, 255.0f));
	Halide::Image<uint8_t> output = f.realize(256,256,3);
	save(output, "power_native.png");
}

void checkerboard()
{
	printf("checkerboard\n");
	
	// generate a checkerboard pattern
	Halide::Var x, y, c;
	Halide::Expr is_lit;
	Halide::Func f;

	int shift = 5;
	
	is_lit = ((x>>shift&1) ^ (y>>shift&1))==1; // calculation of a checkerboard
	
	// NOTE: without uint8_t cast, interprets the <int> output as <u8> which makes a funny picture
	f(x,y,c) = Halide::cast<uint8_t>( Halide::select(is_lit, 255, 0) ); 
	
	
	Halide::Image<uint8_t> output = f.realize(256,256,3);
	
	/*Halide::cast<uint8_t>( 
		lambda(x, y, c, Halide::select(is_lit, 0, 255)).realize(256,256,3)
		);*/
	save(output, "checkerboard.png");
}

// displays the use of a select "case" statement for the color
void colored_checkerboard()
{
	printf("colored_checkerboard\n");
	
	// generate a colored checkerboard pattern
	Halide::Var x, y, c;
	Halide::Expr is_lit;
	Halide::Expr get_lit_color;
	Halide::Expr get_unlit_color;
	Halide::Func f;
	
	int lit_r = 255, lit_g = 0, lit_b = 0;
	int unlit_r = 0, unlit_g = 255, unlit_b = 0;

	int shift = 5;
	
	is_lit = ((x>>shift&1) ^ (y>>shift&1))==1; // calculation of a checkerboard
	
	get_lit_color = Halide::select(c==0, lit_r, select(c==1, lit_g, lit_b));
	get_unlit_color = Halide::select(c==0, unlit_r, select(c==1, unlit_g, unlit_b));
	
	// NOTE: without uint8_t cast, interprets the <int> output as <u8> which makes a funny picture
	f(x,y,c) = Halide::cast<uint8_t>( Halide::select(is_lit, get_lit_color, get_unlit_color) ); 
	
	Halide::Image<uint8_t> output = f.realize(256,256,3);
	
	save(output, "colored_checkerboard.png");
}
void run_original_test()
{
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
}

int main(int argc, char **argv) {
	
	//run_original_test();
	runOtherTests();
	
	return 0;
}
