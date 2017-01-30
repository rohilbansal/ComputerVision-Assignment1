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
//SDoublePlane fft_magnitude(const SDoublePlane &fft_real, const SDoublePlane &fft_imag);

// Write this in Part 1.2


// Write this in Part 1.3 -- add watermark N to image
SDoublePlane mark_image(const SDoublePlane &input, int N);

// Write this in Part 1.3 -- check if watermark N is in image
SDoublePlane check_image(const SDoublePlane &input, int N);



/*
gdhody 29-Jan-2017
*/
SDoublePlane fft_magnitude(const SDoublePlane &fft_real, const SDoublePlane &fft_imag){
  SDoublePlane spectrogram = SDoublePlane(fft_real.rows(),fft_real.cols());
  int rowLoop = 0, columnLoop = 0;
	double value , realComponent, imaginaryComponent;
	for (rowLoop = 0;rowLoop < fft_real.rows(); ++rowLoop){
		for(columnLoop = 0;columnLoop < fft_real.cols(); ++columnLoop){
			//spectrogram[rowLoop][columnLoop] = log(sqrt(fft_real[rowLoop][columnLoop]*fft_real[rowLoop][columnLoop]+fft_imag[rowLoop][columnLoop]*fft_imag[rowLoop][columnLoop]));
			realComponent = (double)fft_real[rowLoop][columnLoop];
			imaginaryComponent = (double)fft_imag[rowLoop][columnLoop];
			value = sqrt(realComponent + imaginaryComponent);
			if (value != 0)
				spectrogram[rowLoop][columnLoop] = (double)log(value);
			else
				spectrogram[rowLoop][columnLoop] = (double)0;
		}
	}
	printf("Magnitude Computed\n");
	return spectrogram;
}


SDoublePlane remove_interference(const SDoublePlane &input){
	int rowLoop = 0, columnLoop = 0;
	int interference [] = {156,160,352,256};
	int rowsInImage = input.rows();
	for (int _2moves = 0; _2moves < 2; ++_2moves){
		for (rowLoop = interference[0 + (2*_2moves)];rowLoop <= interference[1 + (2*_2moves)];++rowLoop){
			for (columnLoop = 0; columnLoop < rowsInImage; ++columnLoop){
				input[rowLoop][columnLoop] = (double)0;
			}
		}
	}
	printf("Interference Removed\n");
	return input;
}
/*
gdhody 29-Jan-2017
*/



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
    cout << "In: " << inputFile <<"  Out: " << outputFile;

    printf("%s\n",inputFile.c_str());
		SDoublePlane input_image = SImageIO::read_png_file(inputFile.c_str());
    
    if(strcmp(argv[1],"1.1")==0)
      {
				/*
				gdhody 29/1/2017
				*/
				SDoublePlane real,imagine,result;
				fft(input_image,real,imagine);
				result = fft_magnitude(real,imagine);
				SImageIO::write_png_file(argv[3],result,result,result);
				printf ("Magnitude 1.1 Module Finished\n");
				/*
				gdhody 29/1/2017
				*/	
      }
    else if(strcmp(argv[1],"1.2")==0)
      {
			/*
				gdhody 29/1/2017
				*/
				//Get the FFT
				SDoublePlane real,imagine,result;
				fft(input_image,real,imagine);
				result = fft_magnitude(real,imagine);
				//Remove the noise
				result = remove_interference(result);
				int interference [] = {155,161,350,256};
				int rowsInImage = real.rows();
				for (int _2moves = 0; _2moves < 2; ++_2moves){
					for (int rowLoop = interference[0 + (2*_2moves)];rowLoop <= interference[1 + (2*_2moves)];++rowLoop){
						for (int columnLoop = 0; columnLoop < rowsInImage; ++columnLoop){
							real[rowLoop][columnLoop] = (double)0;
							imagine[rowLoop][columnLoop] = (double)0;
						}
					}
				}
				for (int x = 0; x < 511; ++x){
					real[255][x]=0;
					imagine[255][x]=0;
					real[x][255]=0;
					imagine[x][255]=0;
					real[x][255]=0;
					imagine[x][255]=0;
					real[x][255]=0;
					imagine[x][255]=0;
				}
				//Do the IDFT transformation
				SDoublePlane Nresult;
				ifft(real,imagine,Nresult);
				SImageIO::write_png_file(argv[3],Nresult,Nresult,Nresult);
				/*
				gdhody 29/1/2017
				*/	
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
       
	int N = atoi(argv[5]);
      }
    else
      throw string("Bad part!");

  } 
  catch(const string &err) {
    cerr << "Error: " << err << endl;
  }
}








