

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <opencv2/core/core.hpp>

// cov
#include <opencv2/ml/ml.hpp>

//fft
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/imgcodecs/imgcodecs.hpp"
#include "opencv2/highgui/highgui.hpp"

#include <iostream>


// x86
// #define FILEPATH "/home/echocare/Downloads/for_trego/input_COV.csv"
//ARM
#ifdef __ARM
#define COV_FILEPATH "/home/debian/inputs/input_COV.csv"
#define FFT_FILEPATH "/home/debian/inputs/input_FIR_FFT.csv"
#define OUTPUT_COV_FILEPATH "/home/debian/outputs/output_COV.csv"
#define OUTPUT_FFT_FILEPATH "/home/debian/outputs/output_FFT.csv"
#else
/* x86 */
#define COV_FILEPATH "/home/echocare/Desktop/winshare/input_COV.csv"
#define FFT_FILEPATH "/home/echocare/Desktop/winshare/input_FIR_FFT.csv"
#define OUTPUT_COV_FILEPATH "/home/echocare/Desktop/winshare/output_COV.csv"
#define OUTPUT_FFT_FILEPATH "/home/echocare/Desktop/winshare/output_FFT.csv"
#endif


using namespace cv;
using namespace std;


#define COLS 288
#define ROWS 200

#include <fstream>

void writeCSV(string filename, Mat m)
{
   ofstream myfile;
   myfile.open(filename.c_str());
   myfile<< cv::format(m, cv::Formatter::FMT_CSV) << std::endl;
   myfile.close();
}

void covar_opencv()
{
	// read data from csv
#if 0 // cv 2.x
	CvMLData mlData;
	mlData.read_csv(COV_FILEPATH);
#endif

	Ptr<ml::TrainData> tdata = ml::TrainData::loadFromCSV(COV_FILEPATH,0,0,-1);
	Mat_<float> samples   = tdata->getTrainSamples();
#if 0
	const CvMat* tmp = mlData.get_values();
	cv::Mat samples(tmp, true);
	tmp->CvMat::~CvMat();
#endif
//	samples = cv::Mat(ROWS, COLS, CV_32FC1, InputArray);
	// process data

	//cv::Mat_<uchar> samples(2,9);  samples << 1,3,2,5,8,7,12,2,4,8,6,9,4,3,3,2,7,7;
	cv::Mat_<float> covar, mean;
	std::cout << "\nsamples\n" << samples;
	cv::calcCovarMatrix( samples, covar, mean, cv::COVAR_NORMAL|cv::COVAR_COLS|cv::COVAR_SCALE, CV_32FC1);
	std::cout << "\nMean\n" << mean << "\nCovar\n" << covar << std::endl;
	 writeCSV(OUTPUT_COV_FILEPATH, covar);
}

void fft_opencv()
{
	// Ptr<ml::TrainData> tdata = ml::TrainData::loadFromCSV(FFT_FILEPATH,0,1,1);
	Ptr<ml::TrainData> tdata = ml::TrainData::loadFromCSV(FFT_FILEPATH,-1,-1,-1);
	Mat_<float> Full   = tdata->getTrainSamples();
	Mat I = Full.row(0);

	// cv::hconcat(I,Last,I);
	std::cout << "Input:\n" << I;
	std::cout << "\nlast:\n" << Full.at<float>(0,287);
	std::cout << "\nlast-1:\n" << Full.at<float>(0,286);
	std::cout << "\nI last:\n" << I.at<float>(0,287);
	std::cout << "\nI last-1:\n" << I.at<float>(0,286);

    Mat padded;                            //expand input image to optimal size
 //   int m = getOptimalDFTSize( I.rows );
 //   int n = getOptimalDFTSize( I.cols ); // on the border add zero values
    int m = getOptimalDFTSize( I.rows ); // each time one raw
    int n = getOptimalDFTSize( I.cols ); // each time one raw

    copyMakeBorder(I, padded, 0, m - I.rows, 0, n - I.cols, BORDER_CONSTANT, Scalar::all(0));
    Mat planes[] = {Mat_<float>(padded), Mat::zeros(padded.size(), CV_32F)};
    Mat complexI;
    merge(planes, 2, complexI);         // Add to the expanded another plane with zeros
    dft(complexI, complexI);           // this way the result may fit in the source matrix
    // dft(padded, complexI);
    // dft(I,complexI);
    // compute the magnitude and switch to logarithmic scale
    // => log(1 + sqrt(Re(DFT(I))^2 + Im(DFT(I))^2))

    writeCSV(OUTPUT_FFT_FILEPATH, complexI);

	// std::cout << "dft:\n" << complexI;
}
extern "C"
{
	void opencv_tests()
	{
		//	ReadFileToInputArray("/home/echocare/Downloads/for_trego/input_COV.csv");
		// fft_opencv();
		covar_opencv();

	}
}
