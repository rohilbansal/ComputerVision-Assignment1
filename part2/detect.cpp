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


#define HOG_CELL_SIZE 8
#define HOG_BIN_SIZE 9
#define HOG_BIN_SEPERATION 20
#define HOG_BLOCK_SIZE 2
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

SDoublePlane extend_image_boundaries(const SDoublePlane input, int length_to_extend, string file_name){

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

	SImageIO::write_png_file(file_name.c_str(),extended_image,extended_image,extended_image);

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
	SDoublePlane complete_image = extend_image_boundaries(input,center_point,"extended_convolve_image.png");
	SDoublePlane output = SDoublePlane(complete_image.rows(), complete_image.cols());
	SDoublePlane Routput = SDoublePlane(input.rows(), input.cols());


	//Fill the image content in the extended image
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
SDoublePlane convolve_general(SDoublePlane input, const SDoublePlane &filter)
{
	int center_point = ((int)filter.rows() / 2);
	int filterW, filterH;
	filterW = filterH = filter.rows();
	SDoublePlane output = SDoublePlane(input.rows(), input.cols());        
	SDoublePlane complete_image = extend_image_boundaries(input,center_point,"extended_convolve_image.png");

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


//Apply a threshold to binary Image
SDoublePlane binary_image(const char* image_name,const double threshold){
	SDoublePlane edge_image = SImageIO::read_png_file(image_name);
	for(int row_loop = 0;row_loop < edge_image.rows(); ++row_loop){
			for(int col_loop = 0;col_loop < edge_image.cols(); ++col_loop){
					if (edge_image[row_loop][col_loop] < threshold)
							edge_image[row_loop][col_loop] = 0.0;
					else
							edge_image[row_loop][col_loop] = 255.00;
		}
	}
	return edge_image;
}

SDoublePlane dilation(SDoublePlane edge_image,int filter_size){

	SDoublePlane morph_filter = SDoublePlane(filter_size,filter_size);
	for (int row = 0; row < morph_filter.rows(); ++row){
			for (int col = 0; col < morph_filter.cols(); ++col){
					morph_filter[row][col] = 1.00;
			}
	}
	
	int center_point = (int)(morph_filter.rows() / 2);

	SDoublePlane morph_dilation = SDoublePlane(edge_image.rows(),edge_image.cols());
	for (int row = center_point; row < edge_image.rows() - center_point; ++row){
			for (int col = center_point; col < edge_image.cols() - center_point; ++col){
					morph_dilation[row][col] = 0.0;
					for(int k = 0; k < morph_filter.rows(); k++){
                                for(int l = 0; l < morph_filter.cols(); l++){
                                        if (morph_filter[k][l] == 1.00 && edge_image[k + row - center_point][l + col - center_point] == 255.0){
												morph_dilation[row][col] = 255.0;
												break;
										}
                                }
								if (morph_dilation[row][col] == 255.00)
										break;
                    }
			}
	}

	return morph_dilation;
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


void hog_implementation(){
	SDoublePlane input_image = SImageIO::read_png_file("Car.png");	
	SDoublePlane sobel_dx = sobel_gradient_filter(input_image,true);
	SDoublePlane sobel_dy = sobel_gradient_filter(input_image,false);
	SDoublePlane magnitude = SDoublePlane(sobel_dx.rows(),sobel_dy.cols());
	SDoublePlane angle = SDoublePlane(sobel_dx.rows(),sobel_dy.cols());
	for(int row_loop = 0; row_loop < magnitude.rows(); row_loop++){
		for(int col_loop = 0; col_loop < magnitude.cols(); col_loop++){
			magnitude[row_loop][col_loop] = sqrt((sobel_dx[row_loop][col_loop] * sobel_dx[row_loop][col_loop]) + (sobel_dy[row_loop][col_loop] * sobel_dy[row_loop][col_loop]));
			if (sobel_dx[row_loop][col_loop] != 0.0)
				angle[row_loop][col_loop] = (double)(atan( ((double)sobel_dy[row_loop][col_loop]) / ((double)sobel_dx[row_loop][col_loop]) ) * 180.00) / PI;
			else if (sobel_dy[row_loop][col_loop] != 0.0)
				angle[row_loop][col_loop] = 0.0;				
			else
				angle[row_loop][col_loop] = -90.00;
			//since atan2 gives between -PI/2 to +PI/2 and histogram wants a range between 0 and 180
			angle[row_loop][col_loop] += 90.00;
		}
	}


	//angle max min
	double minA,maxA;
	min_max(angle,minA,maxA);
	printf ("Angle %f %f",minA,maxA);

	write_normalized_image("hog_magnitude.png",magnitude);
	SDoublePlane hog_image = SImageIO::read_png_file("hog_magnitude.png");
	printf("Hog Magnitude Computed\n");
	
	//HOG Constants
	int cell_size = HOG_CELL_SIZE;
	int bin_size = HOG_BIN_SIZE;
	int bin_seperation = HOG_BIN_SEPERATION;
	
	//HOG Histogram Initialization
	int size_histogram = (hog_image.rows()/cell_size) * (hog_image.cols()/cell_size);
	double ** histogram = new double*[size_histogram];
	for (int assign = 0;assign < size_histogram; ++assign){
		histogram[assign] = new double[bin_size];
		for (int init_zero = 0; init_zero < bin_size; ++init_zero)
			histogram[assign][init_zero] = 0.0;
	}
	printf("Hog Initialized\n");

	//HOG Histogram
	int index_x = 0, index_y = 0, angle_xy = 0, lower_bound = 0, upper_bound = 0;
	for (int row = 0; row < hog_image.rows()/cell_size; ++row){
		for (int col = 0; col < hog_image.cols()/cell_size; ++col){
			for (int cell_r = 0; cell_r < cell_size; ++cell_r){
				for (int cell_c = 0; cell_c < cell_size; ++cell_c){
					index_x = (cell_size * row) + cell_r;
					index_y = (cell_size * col) + cell_c;
					angle_xy = ((int)angle[index_x][index_y]);
					
					//finding bin number
					for (int move = -10; move <= 190; move += HOG_BIN_SEPERATION){
							if (angle_xy <= move){
									upper_bound = move;
									break;
							}
							if ((angle_xy - move) <= HOG_BIN_SEPERATION)
									lower_bound = move;		
					}
					if (upper_bound == 190)
							histogram[(row * ((int)(hog_image.cols()/cell_size))) + col][HOG_BIN_SIZE - 1] += magnitude[index_x][index_y]; 
					else if (lower_bound == -10)
							histogram[(row * ((int)(hog_image.cols()/cell_size))) + col][0] += magnitude[index_x][index_y];
					else{
							histogram[(row * ((int)(hog_image.cols()/cell_size))) + col][(int)(upper_bound / HOG_BIN_SEPERATION)] += magnitude[index_x][index_y] 
																			* ( (angle[index_x][index_y] - lower_bound) / ((float)bin_seperation) );
							//histogram[(row * (hog_image.rows()/cell_size)) + col][(upper_boundary - 10) / bin_seperation] += 0;
							histogram[(row * ((int)(hog_image.cols()/cell_size))) + col][(int)(lower_bound / HOG_BIN_SEPERATION)] += magnitude[index_x][index_y] 
																			* ( (upper_bound - angle[index_x][index_y]) / ((float)bin_seperation) );
					}
					//if (lower_bound >= 90.00)
						//printf("index_x%dindex_y%dLB%dUB%d\n",index_x,index_y,lower_bound,upper_bound);
				}
			}
		}
	}

	
	
	//HOG Blocks
	int block_bin_size = HOG_BIN_SIZE * 4;
	int size_block = (hog_image.rows()/cell_size - 1) * (hog_image.cols()/cell_size - 1);


	double **blocks = new double*[size_block];
	for (int assign = 0;assign < size_block; ++assign){
		blocks[assign] = new double[block_bin_size];
		for (int init_zero = 0; init_zero < block_bin_size; ++init_zero)
			blocks[assign][init_zero] = 0.0;
	}

	double **blocks_normalized = new double*[size_block];
	for (int assign = 0;assign < size_block; ++assign){
		blocks_normalized[assign] = new double[block_bin_size];
		for (int init_zero = 0; init_zero < block_bin_size; ++init_zero)
			blocks_normalized[assign][init_zero] = 0.0;
	}

	//Finding the blocks
	int blocks_per_row = (hog_image.cols()/cell_size);
	for (int row = 0; row < (hog_image.rows()/cell_size) - 1; ++row){
		for (int col = 0; col < (hog_image.cols()/cell_size) - 1; ++col){
			for (int block_r = 0; block_r < HOG_BLOCK_SIZE; ++block_r){
				for (int block_c = 0; block_c < HOG_BLOCK_SIZE; ++block_c){
					for (int cell = 0; cell < HOG_BIN_SIZE; ++cell){
							//printf("%d , %d\n",(row * (blocks_per_row - 1)) + col,(((block_r * HOG_BLOCK_SIZE) + block_c) * HOG_BIN_SIZE) + cell);
							blocks[(row * (blocks_per_row - 1)) + col][(((block_r * HOG_BLOCK_SIZE) + block_c) * HOG_BIN_SIZE) + cell] 
									= histogram[(row * blocks_per_row) + col + (block_r * blocks_per_row) + block_c][cell];
						
					}
				}
			}
		}
	}
	printf("HOG Block Size (C,R) : %d,%d\n",((hog_image.cols()/cell_size)-1) , ((hog_image.rows()/cell_size)-1));

	//Normalizing the blocks
	double magnitude_block = 0.0;
	for (int block_head = 0; block_head < size_block; ++block_head){
			magnitude_block = 0.0;
			for (int block_val = 0; block_val < block_bin_size; ++block_val)
				magnitude_block += (blocks[block_head][block_val] * blocks[block_head][block_val]);
			magnitude_block = sqrt(magnitude_block);
			for (int block_val = 0; block_val < block_bin_size; ++block_val)
				blocks_normalized[block_head][block_val] = blocks[block_head][block_val] / magnitude_block;
	}


	ofstream angle_value;
	angle_value.open("angle.txt");
	for (int row = 0; row < hog_image.rows(); ++row){
			for (int col = 0; col < hog_image.cols(); ++col){
					angle_value << angle[row][col] << ",";
			}
			angle_value << "\n";
	}
	angle_value.close();


	//output to txt file
	ofstream hog_vector_values;
	hog_vector_values.open("hog_vectors.txt");
	for (int block_head = 0; block_head < size_block; ++block_head){
			for (int block_val = 0; block_val < block_bin_size; ++block_val)
				hog_vector_values << blocks_normalized[block_head][block_val] << ",";
			hog_vector_values << "\n";
	}
	hog_vector_values.close();

	ofstream histogram_vector_values;
	histogram_vector_values.open("histogram_vectors.txt");
	for (int block_head = 0; block_head < size_histogram; ++block_head){
			for (int block_val = 0; block_val < HOG_BIN_SIZE; ++block_val)
				histogram_vector_values << histogram[block_head][block_val] << ",";
			histogram_vector_values << "\n";
	}
	histogram_vector_values.close();



	printf("Hog Out\n");
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
  	SDoublePlane input_image = SImageIO::read_png_file(input_filename.c_str());	
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
	SDoublePlane edge_image = binary_image("sobel_magnitude.png",40.00);
	SImageIO::write_png_file("binary_image.png",edge_image,edge_image,edge_image);
	SDoublePlane morph_dilation = dilation(edge_image,11);
	SDoublePlane pass_image = SDoublePlane(edge_image.rows(),edge_image.cols());
	for (int row = 0;row < edge_image.rows(); ++row){
			for (int col = 0;col < edge_image.cols(); ++col){
					if (morph_dilation[row][col] == 255.00)
							pass_image[row][col] = input_image[row][col];
			}
	}
	//morph_dilation = dilation(morph_dilation,9);
	/*DoublePlane morph_erosion = SDoublePlane(edge_image.rows(),edge_image.cols());
	for (int row = 1; row < edge_image.rows() - 1; ++row){
			for (int col = 1; col < edge_image.cols() - 1; ++col){
					morph_erosion[row][col] = 255.0;
					for(int k = 0; k < morph_erode_filter.rows(); k++){
                                for(int l = 0; l < morph_erode_filter.cols(); l++){
                                        if (morph_erode_filter[k][l] == 1.00 && edge_image[k + row - 1][l + col - 1] == 0.0)
												morph_erosion[row][col] = 0.0;
                                }
                    }
			}
	}*/

//	SImageIO::write_png_file("morph_erosion.png",morph_erosion,morph_erosion,morph_erosion);
	
	SImageIO::write_png_file("morph_dilation.png",morph_dilation,morph_dilation,morph_dilation);
	SImageIO::write_png_file("pass_image.png",pass_image,pass_image,pass_image);
	hog_implementation();




/*

 SDoublePlane car_template = SImageIO::read_png_file("template-car.png");

  // calculate car template's mean
  int car_template_x = car_template.rows();
  int car_template_y = car_template.cols();
  cout << car_template_x << " " <<  car_template_y << endl;
  double mean_car_template;
  for(int i = 0; i < car_template_x; i++){
    for(int j = 0; j < car_template_y; j++){
        mean_car_template += car_template[i][j];
    }
  }
  cout << "sum " << mean_car_template << endl;
  mean_car_template = mean_car_template / (car_template_x * car_template_y);
  cout << mean_car_template << endl;

  // calculate sample window's mean

  int window_x_size[2] = {32, 74};
  int window_y_size[2] = {89, 108};

  double mean_window;

  for(int i = window_x_size[0]; i < window_x_size[1]; i++){
    for(int j = window_y_size[0]; j < window_y_size[1]; j++){
      mean_window += input_image[i][j];
    }
  }

  cout << "sum " << mean_window << endl;
  mean_window = mean_window / ((window_y_size[1] - window_y_size[0]) * (window_x_size[1] - window_x_size[0]));
  cout << mean_window << endl;

  // variance of car template
  double S_ff;

  for(int i = 0; i < car_template_x; i++){
    for(int j = 0; j < car_template_y; j++){
      S_ff += (car_template[i][j] - mean_car_template) * (car_template[i][j] - mean_car_template);
    }
  }

  cout << "sum S_ff " << S_ff << endl;
  S_ff = S_ff / (car_template_x * car_template_y);

  cout << "S_ff " << S_ff << endl;
  // variance of window
  double S_gg;

  for(int i = window_x_size[0]; i < window_x_size[1]; i++){
    for(int j = window_y_size[0]; j < window_y_size[1]; j++){
      S_gg += (input_image[i][j] - mean_window) * (input_image[i][j] - mean_window);
    }
  }
  cout << "sum S_gg " << S_gg << endl;
  S_gg = S_gg / (car_template_x * car_template_y);
  cout << "S_gg " << S_gg << endl;

  double S_fg;
  for(int i_template = 0, i_window = window_x_size[0]; i_template < car_template_x && i_window < window_x_size[1]; i_template++, i_window++){
    for(int j_template = 0, j_window = window_y_size[0]; j_template < car_template_y && j_window < window_y_size[1]; j_template++, j_window++){
      //cout << "template " << i_template << " " << j_template << " window " << i_window << " " << j_window << endl;
      // NOTE -  not sure if I should use abs
      S_fg += abs(car_template[i_template][j_template] - mean_car_template) * abs(input_image[i_window][j_window] - mean_window);
    }
  }

  cout << "sum S_fg " << S_fg << endl;

  S_fg = S_fg / (car_template_x * car_template_y);

  cout << "S_fg " << S_fg << endl;

  // correlation

  double corr;
  corr = S_fg / (sqrt(S_ff) * sqrt(S_gg));
  cout << corr << endl;

  // run a sliding window on the pass image of the same dimension as template image


*/









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
