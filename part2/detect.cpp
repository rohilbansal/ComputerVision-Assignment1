//
// detect.cpp : Detect cars in satellite images.
//
// Based on skeleton code by D. Crandall, Spring 2017
//
// NAME
//
//


#include <SImage.h>
#include <SImageIO.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <limits.h>

#define PI 3.14159265


using namespace std;

// The simple image class is called SDoublePlane, with each pixel represented as
// a double (floating point) type. This means that an SDoublePlane can represent
// values outside the range 0-255, and thus can represent squared gradient magnitudes,
// harris corner scores, etc. 
//
// The SImageIO class supports reading and writing PNG files. It will read in
// a color PNG file, convert it to grayscale, and then return it to you in 
// an SDoublePlane. The values in this SDoublePlane will be in the range [0,255].
//
// To write out an image, call write_png_file(). It takes three separate planes,
// one for each primary color (red, green, blue). To write a grayscale image,
// just pass the same SDoublePlane for all 3 planes. In order to get sensible
// results, the values in the SDoublePlane should be in the range [0,255].
//

// Below is a helper functions that overlays rectangles
// on an image plane for visualization purpose. 

// Draws a rectangle on an image plane, using the specified gray level value and line width.
//
//COMMENTS - Final stop for drawing the box over the cars detected
void overlay_rectangle(SDoublePlane &input, int _top, int _left, int _bottom, int _right, double graylevel, int width)
{
  for(int w=-width/2; w<=width/2; w++) {
    int top = _top+w, left = _left+w, right=_right+w, bottom=_bottom+w;

    // if any of the coordinates are out-of-bounds, truncate them 
    top = min( max( top, 0 ), input.rows()-1);
    bottom = min( max( bottom, 0 ), input.rows()-1);
    left = min( max( left, 0 ), input.cols()-1);
    right = min( max( right, 0 ), input.cols()-1);
      
    // draw top and bottom lines
    for(int j=left; j<=right; j++)
	  input[top][j] = input[bottom][j] = graylevel;
    // draw left and right lines
    for(int i=top; i<=bottom; i++)
	  input[i][left] = input[i][right] = graylevel;
  }
}

// DetectedBox class may be helpful!
//  Feel free to modify.
//
//COMMENTS - Detected Box class details
class DetectedBox {
public:
  int row, col, width, height;
  double confidence;
};


//COMMENTS - write files code
// Function that outputs the ascii detection output file
void  write_detection_txt(const string &filename, const vector<DetectedBox> &cars)
{
  ofstream ofs(filename.c_str());

  for(vector<DetectedBox>::const_iterator s=cars.begin(); s != cars.end(); ++s)
    ofs << s->row << " " << s->col << " " << s->width << " " << s->height << " " << s->confidence << endl;
}

// Function that outputs a visualization of detected boxes
void  write_detection_image(const string &filename, const vector<DetectedBox> &cars, const SDoublePlane &input)
{
  SDoublePlane output_planes[3];

  for(int p=0; p<3; p++)
    {
      output_planes[p] = input;
      for(vector<DetectedBox>::const_iterator s=cars.begin(); s != cars.end(); ++s)
	overlay_rectangle(output_planes[p], s->row, s->col, s->row+s->height-1, s->col+s->width-1, p==2?255:0, 2);
    }

  SImageIO::write_png_file(filename.c_str(), output_planes[0], output_planes[1], output_planes[2]);
}
//End of write files


// The rest of these functions are incomplete. These are just suggestions to 
// get you started -- feel free to add extra functions, change function
// parameters, etc.
/*
gdhody Custom functions
*/
double gaussianValue_xy(const int x, const int y, const int sigma){
	return exp(-( ((double)(x * x + y * y)) / ((double)(2.00 * sigma * sigma)) )) / ((double)(2.00 * PI * sigma * sigma));
}


double gaussianValue_x(const int x, const int sigma){
	return exp(-( ((double)(x * x)) / ((double)(2.00 * sigma * sigma)) )) / sqrt((double)(2.00 * PI * sigma * sigma));
}



//Creates a Gaussian 1D kernel
SDoublePlane create_gaussian_kernel_1D(const int windowSize,const int sigma){
	printf("---***---gaussian kernel---***---\n");
	int middleSpot = (int)(windowSize / 2);
	SDoublePlane gaussianKernel = SDoublePlane(windowSize, windowSize);
	for (int rowloop = 0; rowloop < windowSize; ++rowloop){
		gaussianKernel[0][rowloop] = gaussianValue_x(rowloop - middleSpot,sigma);
		printf("%f  ",gaussianKernel[0][rowloop]);
		printf("\n");
	}	
	printf("---***---gaussian kernel---***---\n");
	printf("-----1D gaussian kernel with sigma %d and window size %d-----\n",sigma,windowSize);
	return gaussianKernel;		
}



//Creates a Gaussian 2D kernel
SDoublePlane create_gaussian_kernel_2D(const int windowSize,const int sigma){
	printf("---***---gaussian kernel---***---\n");
	int middleSpot = (int)(windowSize / 2);
	SDoublePlane gaussianKernel = SDoublePlane(windowSize, windowSize);
	for (int rowloop = 0; rowloop < windowSize; ++rowloop){
		for(int columnloop = 0; columnloop < windowSize; ++columnloop){
			gaussianKernel[rowloop][columnloop] = gaussianValue_xy(rowloop - middleSpot,columnloop - middleSpot,sigma);
			printf("%f  ",gaussianKernel[rowloop][columnloop]);
		}
		printf("\n");
	}	
	printf("---***---gaussian kernel---***---\n");
	printf("-----2D gaussian kernel with sigma %d and window size %d-----\n",sigma,windowSize);
	return gaussianKernel;
			
}

SDoublePlane extend_image_boundaries(const SDoublePlane input, int length_to_extend){

	SDoublePlane extended_image(input.rows() + (2 * length_to_extend), input.cols() + (2 * length_to_extend));
	
	//Fills the Tic Toe Positions of the 2D array with boundary values
	for (int outer_loop = length_to_extend; outer_loop > 0; --outer_loop){
		for(int inner_loop = length_to_extend; inner_loop < extended_image.cols() - length_to_extend; ++inner_loop){
			extended_image[outer_loop-1][inner_loop] = input[0][inner_loop - length_to_extend];
			extended_image[extended_image.rows() - outer_loop][inner_loop] = input[input.rows() - 1][inner_loop - length_to_extend];
		}
		for(int inner_loop = length_to_extend; inner_loop < extended_image.rows() - length_to_extend; ++inner_loop){ 
			extended_image[inner_loop][extended_image.cols() - outer_loop] = input[inner_loop - length_to_extend][input.cols() - 1];
			extended_image[inner_loop][outer_loop-1] = input[inner_loop - length_to_extend][0];
		}
	}
	

	//Fill The extended image with the original image
	for (int outer_loop = 0; outer_loop < input.rows(); ++outer_loop){
			for (int inner_loop = 0; inner_loop < input.cols(); ++inner_loop){
					extended_image[outer_loop + length_to_extend][inner_loop + length_to_extend] = input[outer_loop][inner_loop];
			}
	}


	//Fills the Corner Box of Tic Toe Positions With Diagonal Element
	for (int outer_loop = length_to_extend - 1; outer_loop >= 0; --outer_loop){
		for (int inner_loop = length_to_extend - 1; inner_loop >= 0; --inner_loop){
			extended_image[outer_loop][inner_loop] = extended_image[length_to_extend][length_to_extend];
			extended_image[outer_loop][extended_image.cols() - inner_loop - 1] = extended_image[length_to_extend][extended_image.cols() - length_to_extend];
			extended_image[extended_image.rows() - outer_loop - 1][inner_loop] = extended_image[extended_image.rows() - length_to_extend][length_to_extend];
			extended_image[extended_image.rows() - outer_loop - 1][extended_image.cols() - inner_loop - 1] = extended_image[extended_image.rows() - length_to_extend][extended_image.cols() - length_to_extend];	
		}
	}

	SImageIO::write_png_file("extended_image.png",extended_image,extended_image,extended_image);

return extended_image;
}



/*
Custom functions end
*/

void min_max(const SDoublePlane &magnitude,double &minV,double &maxV){
	//updating the normalied values in spectrogram for displaying in png image
	for (int row_loop = 0;row_loop < magnitude.rows(); ++row_loop){
	    for(int col_loop = 0;col_loop < magnitude.cols(); ++col_loop){
				if (magnitude[row_loop][col_loop] < minV)
						minV = magnitude[row_loop][col_loop];
				if (magnitude[row_loop][col_loop] > maxV)
					    maxV = magnitude[row_loop][col_loop];
		}
	}
}




void write_normalized_image(string file_name,SDoublePlane magnitude){
	double minV, maxV;
	min_max(magnitude, minV, maxV);
	//printf("%f,%f\n", minV, maxV);
	//updating the normalied values in spectrogram for displaying in png image
	for (int rowLoop = 0;rowLoop < magnitude.rows(); ++rowLoop){
	    for(int columnLoop = 0;columnLoop < magnitude.cols(); ++columnLoop){
	        if (magnitude[rowLoop][columnLoop] != LONG_MIN) 
	             magnitude[rowLoop][columnLoop] = ((magnitude[rowLoop][columnLoop] - minV) * (255.00/(maxV - minV)));
	    }
	}

	SImageIO::write_png_file(file_name.c_str(),magnitude,magnitude,magnitude);
}

// Convolve an image with a separable convolution kernel
//
SDoublePlane convolve_separable(const SDoublePlane &input, const double *row_filter, const double *column_filter, int filter_size)
{
	int center_point = ((int)(filter_size/2.00));
	SDoublePlane complete_image = extend_image_boundaries(input,center_point);
	SDoublePlane output = SDoublePlane(complete_image.rows(), complete_image.cols());
	SDoublePlane Routput = SDoublePlane(input.rows(), input.cols());

  for(int row_loop = 0; row_loop < complete_image.rows(); ++row_loop){
		  for(int col_loop = 0; col_loop < complete_image.cols(); ++col_loop){
				if (row_loop < center_point || row_loop >= input.rows() + center_point || col_loop < center_point || col_loop >= input.cols() + center_point)
						output[row_loop][col_loop] = complete_image[row_loop][col_loop];
		  }
  }




//Can be further improved with dynamic programming
  for(int row_loop = center_point; row_loop < complete_image.rows() - center_point; ++row_loop){
		  for(int col_loop = center_point; col_loop < complete_image.cols() - center_point; ++col_loop){
				  for(int filter_move = -center_point; filter_move <= center_point; ++filter_move){
							output[row_loop][col_loop] += complete_image[row_loop - filter_move][col_loop] * 
								  (*(column_filter + center_point + filter_move)); 
							//output[row_loop - center_point][col_loop - center_point] /= 8.00;
				  }
		  }
  }



  for(int row_loop = center_point; row_loop < complete_image.rows() - center_point; ++row_loop){
		  for(int col_loop = center_point; col_loop < complete_image.cols() - center_point; ++col_loop){
				  for(int filter_move = -center_point; filter_move <= center_point; ++filter_move){
							Routput[row_loop - center_point][col_loop - center_point] += output[row_loop][col_loop - filter_move] * 
								  (*(row_filter + center_point + filter_move));
							//output[row_loop - center_point][col_loop - center_point] /= 8.00;
				  }
		  }
  }

  // Convolution code here
  return Routput;
}

// Convolve an image with a  convolution kernel
//
SDoublePlane convolve_general(const SDoublePlane &input, const SDoublePlane &filter)
{
	int center_point = ((int)filter.rows() / 2);
	int filterW, filterH;
	filterW = filterH = filter.rows();
	SDoublePlane output = SDoublePlane(input.rows(), input.cols());        
	SDoublePlane complete_image = extend_image_boundaries(input,center_point);

	//Convolution Gaussian Procedure
        for(int i = center_point ; i < complete_image.rows() - center_point; i++){
                for(int j = center_point; j < complete_image.cols() - center_point; j++){
                        for(int k = 0; k < filterH; k++){
                                for(int l = 0; l < filterW; l++){
                                        output[i - center_point][j - center_point] = output[i - center_point][j - center_point] + (complete_image[k + i - center_point]
														[l + j - center_point] * filter[filterH - k - 1][filterW - l - 1]);
                                }
                        }
                }
        }

  
  return output;
}


// Apply a sobel operator to an image, returns the result
// 
SDoublePlane sobel_gradient_filter(const SDoublePlane &input, bool _gx)
{
  //SDoublePlane output(input.rows(), input.cols());
  double sobel_plane [] = {1,2,1};
  double sobel_spread [] = {-1,0,1};

  double *row_operation, *column_operation;
  if (_gx){
	//printf("dx\n");
	row_operation = sobel_spread;
	column_operation = sobel_plane;
  }
  else{
	//printf("dy\n");
	row_operation = sobel_plane;
	column_operation = sobel_spread;
  }

  SDoublePlane output = convolve_separable(input, row_operation, column_operation, 3);

  if (_gx)
		write_normalized_image("sobel_dx.png",output);
  else
		write_normalized_image("sobel_dy.png",output);
  // Implement a sobel gradient estimation filter with 1-d filters
  

  return output;
}

// Apply an edge detector to an image, returns the binary edge map
// 
SDoublePlane find_edges(const SDoublePlane &input, double thresh=0)
{
  SDoublePlane output(input.rows(), input.cols());

  // Implement an edge detector of your choice, e.g.
  // use your sobel gradient operator to compute the gradient magnitude and threshold
  
  return output;
}


//
// This min file just outputs a few test images. You'll want to change it to do 
//  something more interesting!
//
int main(int argc, char *argv[])
{
	if(!(argc == 2))
    	{
      		cerr << "usage: " << argv[0] << " input_image" << endl;
	        return 1;
    	}

  	string input_filename(argv[1]);
  	SDoublePlane input_image= SImageIO::read_png_file(input_filename.c_str());	
	SDoublePlane gaussian = create_gaussian_kernel_2D(5,1);
	SDoublePlane smoothed_image	= convolve_general(input_image,gaussian);
	SImageIO::write_png_file("output.png",smoothed_image,smoothed_image,smoothed_image);
	
	SDoublePlane sobel_dx = sobel_gradient_filter(input_image,true);
	SDoublePlane sobel_dy = sobel_gradient_filter(input_image,false);
	SDoublePlane magnitude = SDoublePlane(sobel_dx.rows(),sobel_dy.cols());
	for(int row_loop = 0; row_loop < magnitude.rows(); row_loop++){
		for(int col_loop = 0; col_loop < magnitude.cols(); col_loop++){
			magnitude[row_loop][col_loop] = sqrt((sobel_dx[row_loop][col_loop] * sobel_dx[row_loop][col_loop]) + (sobel_dy[row_loop][col_loop] * sobel_dy[row_loop][col_loop]));
		}
	}

	write_normalized_image("sobel_magnitude.png",magnitude);

/*
  // test step 2 by applying mean filters to the input image
  SDoublePlane mean_filter(3,3);
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      mean_filter[i][j] = 1/9.0;
  SDoublePlane output_image = convolve_general(input_image, mean_filter);

  
  // randomly generate some detected cars -- you'll want to replace this
  //  with your car detection code obviously!
  vector<DetectedBox> cars;
  for(int i=0; i<10; i++)
    {
      DetectedBox s;
      s.row = rand() % input_image.rows();
      s.col = rand() % input_image.cols();
      s.width = 20;
      s.height = 20;
      s.confidence = rand();
      cars.push_back(s);
    }

  write_detection_txt("detected.txt", cars);
  write_detection_image("detected.png", cars, input_image);
*/
}
