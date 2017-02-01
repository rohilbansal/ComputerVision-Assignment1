//
// Watermark.cpp : Add watermark to an image, or inspect if a watermark is present.
//
// Based on skeleton code by D. Crandall, Spring 2017
//
// ctewani-gdhody-bansalro
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
#include <typeinfo>


#define PI 3.14159265


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


// Write this in Part 1.3 -- add watermark N to image


// Write this in Part 1.3 -- check if watermark N is in image
SDoublePlane check_image(const SDoublePlane &input, int N);



/*
gdhody 29-Jan-2017
*/
SDoublePlane mark_image(const SDoublePlane &input, int N){
	//Generating the v vector through N and l
	srand(N);
	SDoublePlane finalC = SDoublePlane(input.rows(),input.cols());
	//---L---
	int l = N;
	double *v;
	v = new double[l];
	//Assigning the v vector with key N
	for(int bitGen = 0;bitGen <= l; ++bitGen){
		v[bitGen] = (double) (rand() % 2);
		//printf("%f",v[bitGen]);	
	}
	//Constants Initialization
	SDoublePlane real,imagine;
	fft(input,real,imagine);
	int centerPoint = input.rows()/2 - 1;
	//---RADIUS---
	int radius = (int)(1 * (input.rows()/2))-50;
	//---ALPHA---
	double alpha = 25.00;
	//4 Quadrants
	int quadrantL = l / 4;
	//Angle, Real Value Holder, Absolute Constant
	double thetha = 0.0, Rvalue, negativeConstant = 1.0;
	int moveConstantOverV = 0;
	//Center of FFT
	int Xcenter = centerPoint, Ycenter = centerPoint;
	//To visualize circular watermark
	for(int i=0;i < input.rows();++i){
		for(int j=0;j < input.cols();++j){
			finalC[i][j] = 0;
		}
	}
	printf("Initialization For Watermark Done\n");
	//Updating Real Values to add watermark	
	for(int quad = 0; quad < 4; ++quad){
		for(int loop = 0; loop < quadrantL; ++loop){
			thetha = (PI* ((double)(quad*90.00) + ((loop+1)*(90.00/(double)quadrantL)) )) / 180.00;
			Rvalue = real[Xcenter + (int)ceil(radius * cos(thetha))][Ycenter + (int)ceil(radius * sin(thetha))];
			//printf("Previous%f\n", Rvalue);
			//printf("%f\n",(alpha * abs(Rvalue) * ((double)(v[moveConstantOverV]))));
			negativeConstant = 1.0;
			if (Rvalue < 0.0)
				negativeConstant = -1.0;
			//printf("%f,%f,%f\n",alpha,Rvalue*negativeConstant,v[moveConstantOverV]);
			real[Xcenter + (int)ceil(radius * cos(thetha))][Ycenter + (int)ceil(radius * sin(thetha))] = Rvalue + (alpha * Rvalue * negativeConstant * v[moveConstantOverV]);
			//printf("%d,%d,%d\n",quad,loop,moveConstantOverV);
			finalC[Xcenter + (int)ceil(radius * cos(thetha))][Ycenter + (int)ceil(radius * sin(thetha))] = v[moveConstantOverV]*255;	
			moveConstantOverV = moveConstantOverV + 1;
			//printf("Looped%f\n\n", real[Xcenter + (int)ceil(radius * cos(thetha))][Ycenter + (int)ceil(radius * sin(thetha))]);
		}
	}
	SDoublePlane watermarkFFT = fft_magnitude(real,imagine);
	SImageIO::write_png_file("Watermark_FFT.png",watermarkFFT,watermarkFFT,watermarkFFT);
	SDoublePlane outputReturnImage;
	ifft(real,imagine,outputReturnImage);
	SImageIO::write_png_file("Watermark_Circle_Illusration.png",finalC,finalC,finalC);
	printf("ADD-Watermark FFT and Circle Illustration Images Generated For Referral\n");
	return outputReturnImage;
}




SDoublePlane fft_magnitude(const SDoublePlane &fft_real, const SDoublePlane &fft_imag){
	
	//Declaring spectrogram that contains magnitude of FFT
	SDoublePlane spectrogram = SDoublePlane(fft_real.rows(),fft_real.cols());
  	//Declaring constants needed to run code
	int rowLoop = 0, columnLoop = 0;
	double minV = LONG_MAX,maxV = LONG_MIN, imageColorConstant = 255.00, value = 0.0;
	
	//Logic to calculate magnitude and find min and max value in our fft, min max used for normalization of image so we can display it correctly in output
	for (rowLoop = 0;rowLoop < fft_real.rows(); ++rowLoop){
		for(columnLoop = 0;columnLoop < fft_real.cols(); ++columnLoop){
			value = ( sqrt( (fft_real[rowLoop][columnLoop]*fft_real[rowLoop][columnLoop]) + (fft_imag[rowLoop][columnLoop]*fft_imag[rowLoop][columnLoop]) ) );
			
			//Log of 0 is -inf we don't want that
			if (value == 0)
				spectrogram[rowLoop][columnLoop] = LONG_MIN;
			else
				spectrogram[rowLoop][columnLoop] = log(value);

			//Find Min Max value except where value is LONG_MIN
			if (spectrogram[rowLoop][columnLoop] != LONG_MIN){
				if (spectrogram[rowLoop][columnLoop] < minV)
					minV = spectrogram[rowLoop][columnLoop];
				if (spectrogram[rowLoop][columnLoop] > maxV)
					maxV = spectrogram[rowLoop][columnLoop];
			}
		}
	}

	//updating the normalied values in spectrogram for displaying in png image
	for (rowLoop = 0;rowLoop < fft_real.rows(); ++rowLoop){
		for(columnLoop = 0;columnLoop < fft_real.cols(); ++columnLoop){
			//IF Long_Min means the magnitude was log of 0 so it is 0
			if (spectrogram[rowLoop][columnLoop] != LONG_MIN) 
				spectrogram[rowLoop][columnLoop] = ((spectrogram[rowLoop][columnLoop] - minV) * (imageColorConstant/(maxV-minV)));
			else
				spectrogram[rowLoop][columnLoop] = 0;	
		}
	}

	printf("FFT-Magnitude min max without normalization %f,%f; After normalization (0-255)\n",minV,maxV);
	printf("FFT-Magnitude Computed For Image FFT\n");
	return spectrogram;
}



SDoublePlane remove_interference(const SDoublePlane &input){	
	SDoublePlane real,imagine,result,Mresult;
	//Do The FFT
	fft(input,real,imagine);
	int rowLoop = 0, columnLoop = 0;
	//Box Coordinates To Remove Message Put Into Image As Noise - HI
	int interferenceBox [] = {155,350};
	int rowWidth = 5, columnWidth = 8;
	//Box Coordinates Defined
		//Logic To Remove Box Coordinates with 0 intensity
		for (rowLoop = 0;rowLoop <= rowWidth; ++rowLoop){
			for (columnLoop = 0; columnLoop <= columnWidth; ++columnLoop){
				real[interferenceBox[0] + rowLoop][interferenceBox[0] + columnLoop] = 0.0;
				real[interferenceBox[1] + rowLoop + 1][interferenceBox[1] + columnLoop] = 0.0;
				imagine[interferenceBox[0] + rowLoop][interferenceBox[0] + columnLoop] = 0.0;
				imagine[interferenceBox[1] + rowLoop + 1][interferenceBox[1] + columnLoop] = 0.0;
				/*Remove Complete || Horizontal Lines
				real[interference[0] + rowLoop][columnLoop] = (double)0;	
				real[interference[1] + rowLoop][columnLoop] = (double)0;
				imagine[interference[0] + rowLoop][columnLoop] = (double)0;    
			       	imagine[interference[1] + rowLoop][columnLoop] = (double)0;*/
			}
		}
	Mresult = fft_magnitude(real,imagine);
	printf("RMI-Interference Removed\n");
	SImageIO::write_png_file("Noise_Removed_fft.png",Mresult,Mresult,Mresult); 		
	printf("RMI-FFT Magnitude Image Generated After Noise Message - HI - Removal : %s\n","Noise_Removed_fft.png");
	//Parse The Image Back Through IFFT
	ifft(real,imagine,result);
	//Return The Grayscale Component To Be Written Back To PNG
	return result;
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
    cout << "In: " << inputFile <<"  Out: " << outputFile << endl;
    //printf("%s\n",inputFile.c_str());
    SDoublePlane input_image = SImageIO::read_png_file(inputFile.c_str());
    

    if(strcmp(argv[1],"1.1")==0)
      {
				/*
				gdhody 29/1/2017
				*/
				SDoublePlane real,imagine,result;
				//Do the fft
				fft(input_image,real,imagine);
				//Find the magnitude
				result = fft_magnitude(real,imagine);
				SImageIO::write_png_file(argv[3],result,result,result);
				printf ("Magnitude 1.1 Module Finished, with output Image : %s\n",argv[3]);
				/*
				gdhody 29/1/2017
				*/	
      }
    else if(strcmp(argv[1],"1.2")==0)
      {
				/*
				gdhody 29/1/2017
				*/
				SDoublePlane Nresult;
				//Remove the noise
				Nresult = remove_interference(input_image);
				//Save The Image
				SImageIO::write_png_file(argv[3],Nresult,Nresult,Nresult);
				printf ("Noise Removal 1.2 Module Finished, with output Image : %s\n",argv[3]);
				/*
				gdhody 29/1/2017
				*/	
      }
    else if(part == "1.3")
      {
	/*if(argc < 6)
	  {
	    cout << "Need 6 parameters for watermark part:" << endl;
	    cout << "    p2 1.3 inputfile outputfile operation N" << endl;
	    return -1;
	  }
	string op(argv[4]);*/
	if(strcmp(argv[4],"add")==0)
	  {
		SDoublePlane Aresult = mark_image(input_image,90);
		SImageIO::write_png_file(argv[3],Aresult,Aresult,Aresult);
		printf ("Watermark Add 1.3 Module Finished, with output Image : %s\n",argv[3]);
	  }
	else if(strcmp(argv[4],"check")==0)
	  {
	    // check watermark
	  }
	else
	  throw string("Bad operation!");
       
	//int N = atoi(argv[5]);
      }
    else
      throw string("Bad part!");

  } 
  catch(const string &err) {
    cerr << "Error: " << err << endl;
  }
}








