#include "/reg/neh/home/khegazy/baseScripts/imageProcessing.h"
#include "/reg/neh/home/khegazy/baseScripts/plotClass.h"
#include "/reg/neh/home/khegazy/baseScripts/tools.h"
#include <fftw3.h>
#include <valarray>

using namespace std;


struct imgInfoStruct {

  string date;
  string path;
  string fileName;
  string scan;
  int run;
  float stagePos;

  imgInfoStruct () {
    run = -999;
  }

};


namespace ppFunct {

  bool getScanRunInfo(std::vector<imgInfoStruct> &imgINFO, std::string runListName, bool verbose);
  void makeRunLists(std::vector<imgInfoStruct> &imgINFO);
}


/*
  class centerfnctr {

    public:
	vector< vector<double> >* img;
	int count;
	//centerfnctr(vector< vector<double> > &img_i) {img = &img_i;}
	//centerfnctr(); 
	double operator() (vector<double> vect) { 
	  count++;
	  //cout<<"calling fxn          ";
	  return imgProc::centerSymXsqr(img, vect[0], vect[1], 8, 2, 80, 20);
	  //return fabs(1-imgProc::centerSymXsqr(img, vect[0], vect[1], 8, 3, 100, 20));
	}

  };

*/
