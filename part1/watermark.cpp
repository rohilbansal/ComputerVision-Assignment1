//
// Watermark.cpp : Add watermark to an image, or inspect if a watermark is present.
//
// Based on skeleton code by D. Crandall, Spring 2017
//
// PUT YOUR NAMES HERE
//
//

//Link to the header file
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <SImage.h>
#include <SImageIO.h>
#include <fft.h>
#include <constants.h>
#include <vector>
#include <math.h>

using namespace std;

// This code requires that input be a *square* image, and that each dimension
//  is a power of 2; i.e. that input.width() == input.height() == 2^k, where k
//  is an integer. You'll need to pad out your image (with 0's) first if it's
//  not a square image to begin with. (Padding with 0's has no effect on the FT!)
//
// Forward FFT transform: take input image, and return real and imaginary parts.
//
void fft(const SDoublePlane &input, SDoublePlane &fft_real, SDoublePlane &fft_imag)
{
  fft_real = input;
  fft_imag = SDoublePlane(input.rows(), input.cols());

  FFT_2D(1, fft_real, fft_imag);
}

// Inverse FFT transform: take real and imaginary parts of fourier transform, and return
//  real-valued image.
//
void ifft(const SDoublePlane &input_real, const SDoublePlane &input_imag, SDoublePlane &output_real)
{
  output_real = input_real;
  SDoublePlane output_imag = input_imag;

  FFT_2D(0, output_real, output_imag);
}

// Write this in Part 1.1
SDoublePlane fft_magnitude(const SDoublePlane &fft_real, const SDoublePlane &fft_imag);

// Write this in Part 1.2
SDoublePlane remove_interference(const SDoublePlane &input);

// Write this in Part 1.3 -- add watermark N to image
SDoublePlane mark_image(const SDoublePlane &input, int N);

// Write this in Part 1.3 -- check if watermark N is in image
SDoublePlane check_image(const SDoublePlane &input, int N);


int main(int argc, char **argv)
{
  try {

    if(argc < 4)
      {
	cout << "Insufficent number of arguments; correct usage:" << endl;
	cout << "    p2 problemID inputfile outputfile" << endl;
	return -1;
      }

    string part = argv[1];
    string inputFile = argv[2];
    string outputFile = argv[3];
    cout << "In: " << inputFile <<"  Out: " << outputFile << endl;

    SDoublePlane input_image = SImageIO::read_png_file(inputFile.c_str());
    int N = atoi(argv[4]);
    if(part == "1.1")
      {
        // generate a binary vector of l length, we can have this l in a constant file

        // cout << l_param << endl;
        std::vector<int> v(l_param);
        srand(N);
        for(int i = 0; i < l_param; i++)
          v[i] = rand()&1;

         // print binary vector
        //for (std::vector<int>::const_iterator i = v.begin(); i != v.end(); ++i)
        //  std::cout << *i << ' ';
        // cout << endl;

        SDoublePlane fft_real, ft_imag;
        // convert image into frequency domain using FFT
        fft(input_image, fft_real, ft_imag);

        SDoublePlane input_real = SDoublePlane(fft_real.rows(), fft_real.cols());
        SDoublePlane input_imag = SDoublePlane(fft_real.rows(), fft_real.cols());

        // cout << fft_real.rows() << " " << fft_real.cols() <<  " " << fft_real[0][0] << endl;
        for(int i = 0; i < fft_real.rows(); i++){
          for(int j = 0; j < fft_real.cols(); j++){
            if(i >= 156 && j >= 156 && i <= 162 && j <= 162){
              fft_real[i][j] = fft_real[0][0];
              ft_imag[i][j] = ft_imag[0][0];
            }

            if(i >= 350 && j >= 350 && i <= 356 && j <= 356){
              fft_real[i][j] = fft_real[0][0];
              ft_imag[i][j] = ft_imag[0][0];
            }

            input_real[i][j] = log(sqrt((fft_real[i][j] * fft_real[i][j]) + (ft_imag[i][j] * ft_imag[i][j])));

          }
        }

        //SImageIO::
        SImageIO::write_png_file("cleaned_spectogram", input_real, input_real, input_real);

       SDoublePlane output_real;

       ifft(fft_real, ft_imag, output_real);
       SImageIO::write_png_file(outputFile.c_str(), output_real, output_real, output_real);

	// do something here!
      }
    else if(part == "1.2")
      {
	// do something here!
      }
    else if(part == "1.3")
      {
	if(argc < 6)
	  {
	    cout << "Need 6 parameters for watermark part:" << endl;
	    cout << "    p2 1.3 inputfile outputfile operation N" << endl;
	    return -1;
	  }
	string op(argv[4]);
	if(op == "add")
	  {
	    // add watermark
	  }
	else if(op == "check")
	  {
	    // check watermark
	  }
	else
	  throw string("Bad operation!");
	      }
    else
      throw string("Bad part!");

  }
  catch(const string &err) {
    cerr << "Error: " << err << endl;
  }
}
