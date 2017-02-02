#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <SImage.h>
#include <SImageIO.h>
#include <fft.h>

using namespace std;

int main(int argc, char** argv){

	SDoublePlane gaussian = SDoublePlane(3,3);
	gaussian[0][0] = 0.075;
	gaussian[0][1] = 0.125;
	gaussian[0][2] = 0.075;
	gaussian[1][0] = 0.125;
	gaussian[1][1] = 0.200;
	gaussian[1][2] = 0.125;
	gaussian[2][0] = 0.075;
	gaussian[2][1] = 0.125;
	gaussian[2][2] = 0.075;
	
	string inputFile_1 = argv[1];
	string inputFile_2 = argv[2];

	SDoublePlane input_1 = SImageIO::read_png_file(inputFile_1.c_str());
	SDoublePlane input_2 = SImageIO::read_png_file(inputFile_2.c_str());
	SDoublePlane lowPass_1 = SDoublePlane(input_1.rows(), input_1.cols());
	SDoublePlane lowPass_2 = SDoublePlane(input_1.rows(), input_1.cols());
	SDoublePlane highPass = SDoublePlane(input_1.rows(), input_1.cols());
	SDoublePlane morphPass = SDoublePlane(input_1.rows(), input_1.cols());

	for(int i=1 ; i<input_1.rows()-1; i++){
		for(int j=1; j<input_1.cols()-1; j++){
			for(int k=0; k<gaussian.rows(); k++){
				for(int l=0; l<gaussian.cols(); l++){
					lowPass_1[i][j] += input_1[i][j] * gaussian[gaussian.rows()-k-1][gaussian.cols()-l-1];
				}
			}	
		}
	}

	for(int i=1 ; i<input_2.rows()-1; i++){
		for(int j=1; j<input_2.cols()-1; j++){
			for(int k=0; k<gaussian.rows(); k++){
				for(int l=0; l<gaussian.cols(); l++){
					lowPass_2[i][j] += input_2[i][j] * gaussian[gaussian.rows()-k-1][gaussian.cols()-l-1];
				}
			}	
		}
	}

	for(int i=1; i<input_2.rows()-1; i++){
		for(int j=1; j<input_2.cols()-1; j++){
			highPass[i][j] = input_2[i][j] - lowPass_2[i][j];
		}
	}

	for(int i=1; i<input_2.rows()-1; i++){
		for(int j=1; j<input_2.cols()-1; j++){
			morphPass[i][j] = (lowPass_1[i][j] + highPass[i][j]) / 2;
		}
	}
	
	SImageIO::write_png_file("output.png", morphPass, morphPass, morphPass);
	return 0;
}
