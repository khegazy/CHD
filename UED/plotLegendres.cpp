#include "chd.h" 


using namespace std;


int main(int argc, char* argv[]) {

  if (argc<2) {
    cerr<<"ERROR: Missing input arguments, must run code ./analysis.exe 'fileList.txt' !!!"<<endl;
    cerr<<"         Can also run ./analysis.exe 'fileList.txt' 'treeName'"<<endl;
    exit(0);
  }

  string fileList(argv[1]);
  string treeName("physics");
  if (argc==3) string treeName(argv[2]);
  analysisClass analysis(fileList, treeName);  // Use this when specifying treeName


  ///// Load environment and get the number of events /////
  uint64_t Nentries;
  Nentries = analysis.setupEnvironment();  // Alter this function as needed for specific setup

  cout.setf(ios::scientific);


  /////  Flags  /////
  bool findVarRegions = true;

  
  /////  Setting up variables  /////
  int NradBins = 50;
  int Nlg = 6;

  string fileName = "alignment.txt";
  string outputDir = "output/data/";
  for (int iarg=2; iarg<argc; iarg+=2) {
    if (strcmp(argv[iarg],"-Ofile")==0) {
      string str(argv[iarg+1]); 
      fileName = str + ".txt";
    }
    else if (strcmp(argv[iarg],"-Odir")==0) {
      string str(argv[iarg+1]); 
      outputDir = str;
    }
    else {
      cerr<<"ERROR!!! Option "<<argv[iarg]<<" does not exist!"<<endl;
      exit(0);
    }
  }

  double dshift = 4e-3;
  double seed = clock();

  ////////////////////////////////////////////
  /////  Finding Regions that Fluctuate  /////
  ////////////////////////////////////////////

  double qMax = 11.3;

  int ir, ic, itr = 0;
  int curRun = -1;
  string runInd = "";
  string curDate = "";
  string curScan = "";
  double curPosition = 0;
  std::vector<PLOToptions> opts(6);
  std::vector<std::string> vals(6);
  std::vector<PLOToptions> oppts(4);
  std::vector<std::string> vaals(4);
  oppts[0] = xLabel;   vaals[0] = "Time [ps]";
  oppts[1] = yLabel;   vaals[1] = "Scattering Q [inv Angs]";
  oppts[2] = ySpan;    vaals[2] = "0,"+to_string(qMax);
  oppts[3] = draw;     vaals[3] = "CONT4Z";
  //oppts[3] = logz;     vaals[3] = true;
  opts[0] = xLabel;   vals[0] = "Time [ps]";
  opts[1] = yLabel;   vals[1] = "Scattering Q [inv Angs]";
  opts[2] = ySpan;    vals[2] = "0,"+to_string(qMax);
  opts[3] = draw;     vals[3] = "CONT4Z";
  opts[4] = minimum;
  opts[5] = maximum;

  std::vector<double> cbar;
  cbar.push_back(10e18);
  cbar.push_back(25e18);
  cbar.push_back(10e18);
  cbar.push_back(30);
  cbar.push_back(20);
  cbar.push_back(10);


  int arrSize = 188448;
  bool filledRun = false;
  double* refImg = NULL;
  double prevStagePos, stageDiff;
  double stageCut = 1.8e10; //0.0041;
  double compare;

  bool newEntry;
  double rad;
  std::map< double, std::vector<double> > diffPs, legCoeff_map;
  std::map< double, double > counts;

  std::vector<double> atmDiff(NradBins);
  FILE* atmFile = fopen("../simulation/diffractionPattern/output/references/atomicScattering_CHD.dat", "rb");
  fread(&atmDiff[0], sizeof(double), atmDiff.size(), atmFile);
  fclose(atmFile);

  ///////////////////////////
  /////  Aligning Runs  /////
  ///////////////////////////

  std::vector<string> runInds;
  std::map< string, std::map< double, double* > >  diffP_arrays;
  std::map< string, int > runShifts;

  ///// Loop through events in the file /////
  for (uint64_t ievt=0; ievt<Nentries; ievt++) {
    analysis.loadEvent(ievt);

    //if ((*date != "20161105")) {
    //if ((*date != "20161104") || (*scan != "LongScan3")) {
    //  continue;
    //}
    //cout << imgNum << "  " << runNum << "  "<< ievt << " "<<Nentries<<endl;
    // Expect the reference image to be the first image
    if (imgIsRef) {
      continue;
    }

    // Ignore reference images taken before scan
    //if (stagePos < (t0StagePos - 0.021)) {
    //  continue;
    //}

    /////  Make new root file and initialize variables  /////
    if ((curRun != runNum) || (curDate != *date) || (curScan != *scan) || (ievt == Nentries-1)) {

      if  ((curDate != *date) || (curScan != *scan) || (ievt == Nentries-1)) {
        cout<<"start plot "<<curRun<<"  "<<legCoeff_map.size()<<endl;
        if ((curRun != -1) && (legCoeff_map.size() > 1)) {
          cout<<"Plotting"<<endl;
          legCoeff_map.erase(legCoeff_map.begin());
          legCoeff_map.erase(--legCoeff_map.end());
          //legCoeff_map.erase(--legCoeff_map.end());
          //legCoeff_map.erase(--legCoeff_map.end());
          for (auto itr = legCoeff_map.begin(); itr!=legCoeff_map.end(); itr++) {
            cout<<"test pos: "<<itr->first<<endl;
          }
          int ind = 0;
          double *delays = new double[legCoeff_map.size() + 1];
          std::vector< std::vector< std::vector<double> > > img(Nlg); //legCoeff_map.size());
          for (uint i=0; i<Nlg; i++) {
            cout<<"i: "<<i<<endl;
            ind = 0;
            img[i].resize(legCoeff_map.size());
            for (auto coeffs : legCoeff_map) {
              cout<<"pos: "<<coeffs.first<<endl;
              delays[ind] = coeffs.first/(3e-1);
              img[i][ind].resize(NradBins, 0);
              for (uint ir=0; ir<NradBins; ir++) {
                img[i][ind][ir] = coeffs.second[i*NradBins + ir]/counts[coeffs.first];
              }
              ind++;
            }

            for (uint iy=0; iy<img[i][0].size(); iy++) {
              double sum = 0;
              for (uint ix=0; ix<img[i].size(); ix++) {
                sum += img[i][ix][iy];
              }
              sum /= (double)img[i].size();
              for (uint ix=0; ix<img[i].size(); ix++) {
                // Subtract average
                //img[i][ix][iy] -= sum;
                // Scaling
                img[i][ix][iy] /= atmDiff[iy]*(qMax*(iy+0.5)/NradBins);
                if (iy < 2) {
                  img[i][ix][iy] = 0;
                }
              }
            }

            cout<<"last entry"<<endl;
            delays[ind] = 2*delays[ind-1] - delays[ind-2];
            cout<<"plotting"<<endl;

              vals[4] = to_string(-1*cbar[i]);
              vals[5] = to_string(cbar[i]);
              plt.printXY(img[i], delays, curDate + "_" + curScan 
                    + "_Leg" + to_string(i), opts, vals);
              /*
            if ((i==0 || i==2) && (cbar.find(curDate)!=cbar.end()) && (cbar[curDate].find(curScan)!=cbar[curDate].end())) {
              vals[2] = to_string(-1*cbar[curDate][curScan][i/2]);
              vals[3] = to_string(cbar[curDate][curScan][i/2]);
              plt.printXY(img[i], delays, curDate + "_" + curScan 
                    + "_Leg" + to_string(i), opts, vals);
            }
            else {
              plt.printXY(img[i], delays, curDate + "_" + curScan 
                    + "_Leg" + to_string(i), oppts, vaals);
            }
            */
            cout<<"plotted"<<endl;
          }
          cout<<"del del"<<endl;
          delete[] delays;
        }
        
        cout<<"clearing"<<endl;
        legCoeff_map.clear();
        counts.clear();
      }

          cout<<"reset"<<endl;
      filledRun = false;
      curPosition = 0;
      curRun  = runNum;
      curDate = (*date);
      curScan = (*scan);

      runInd = curDate + "_" + curScan + "_" + to_string(curRun);
      runInds.push_back(runInd);
      cout<<runInd<<endl;
      runShifts[runInd] = 0;
      prevStagePos = -1;


    }

    // Don't fill run with large time steps
    if (filledRun) {
      //continue;
    }

    // Need initial stage position (skips first image)
    if (prevStagePos == -1) {
      prevStagePos = stagePos;
      continue;
    }

    // Round the stage difference to the correct decimal place
    compare = 1;
    itr = -1;
    while ((fabs(compare) > 0.01) && (itr < 4)) {
      itr++;
      compare = (tools::round((stagePos - prevStagePos)*pow(10, itr))/pow(10, itr) 
          - (stagePos - prevStagePos))
          /(stagePos - prevStagePos);
    }
    stageDiff = ((double)tools::round((stagePos - prevStagePos)*pow(10, itr)))/pow(10, itr);
    //cout<<"stage diff: "<<stagePos<<"  "<<stagePos - prevStagePos<<"  "<<stageDiff<<endl;
    prevStagePos = stagePos;

    stageDiff = ((double)((int)(1000*stageDiff)))/1000;

    curPosition += stageDiff;
    //cout<<endl;

    // Make sure we are still in regime of small step sizes
    if (stageDiff > stageCut) {
      filledRun = true;
      //continue;
    }


    //curPosition = ((double)((int)(1000*curPosition)))/1000.;
    //if (legCoeff_map.find(curPosition) == legCoeff_map.end()) {
    newEntry = true;
      //cout<<"adding : "<<curPosition<<"   "<<endl;
          for (auto itr = legCoeff_map.begin(); itr!=legCoeff_map.end(); itr++) {
            //cout<<"search pos: "<<itr->first<<endl;
            if (fabs((curPosition - itr->first)/curPosition) < 0.01) {
              curPosition = itr->first;
              newEntry = false;
            }
          }
    if (newEntry) {
      //cout<<"ADDED"<<endl;
      legCoeff_map[curPosition].resize((*legCoeffs).size(), 0);
      counts[curPosition] = 0;
    }
    //cout<<"added"<<endl;
    for (uint i=0; i<(*legCoeffs).size(); i++) {
      legCoeff_map[curPosition][i] += imgNorm*(*legCoeffs)[i];
    }
    counts[curPosition] += imgNorm;
  }

  // Clean up
  for (auto &runItr : diffP_arrays) {
    for (auto &imgItr : runItr.second) {
      delete[] imgItr.second;
    }
  }

  return 1;
}
