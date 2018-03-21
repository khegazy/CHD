#include "preProcessing.h"


////////////////////////////////
////  Retrieving File Info  ////
////////////////////////////////

bool ppFunct::getScanRunInfo(std::vector<imgInfoStruct> &imgINFO, std::string runListName, bool verbose) {

  ifstream fileNames(runListName.c_str());
  if (!fileNames.is_open()) {
    cerr << "ERROR: Cannot open file " << runListName << endl;
    exit(0);
  }

  size_t ipos, fpos;
  string fileName;
  std::map<float, imgInfoStruct> imgInfoMap;

  std::string curScan="";
  std::string curDate="";
  int curRun=-1;

  while (getline(fileNames, fileName)) {

    cout<<"Now looking at file "<<fileName<<endl;

    imgInfoStruct imgInfo;

    // Finding date (2016)
    ipos = fileName.find("2016",0);
    imgInfo.date = fileName.substr(ipos, 8);

    // Finding Scan name
    ipos = fileName.find("scan",0);
    if (ipos == string::npos) {
      ipos = fileName.find("Scan",0);
      if (ipos == string::npos) {
        ipos = fileName.find("Range",0);
        if (ipos == string::npos) {
          cerr << "ERROR: Cannot find 'scan' or 'Scan' in image name!!!\n\n";
          exit(0);
        }
      }
    }
    while (fileName.at(ipos-1) != '/') ipos--;
    fpos = fileName.find("/", ipos);
    imgInfo.scan = fileName.substr(ipos, fpos-ipos);
    if (verbose) cout << "\tSCAN: " << imgInfo.scan;

    // Finding run number
    ipos = fileName.find("run",0);
    if (ipos == string::npos) {
      ipos = fileName.find("Run",0);
    }
    if (ipos != string::npos) {
      fpos = fileName.find("/", ipos);
      imgInfo.run = stoi(fileName.substr(ipos+3, fpos-ipos-3));
    }
    if (verbose) cout << "\tRUN: " << imgInfo.run;

    // Finding file name
    ipos = fileName.rfind("/");
    imgInfo.fileName = fileName.substr(ipos, fileName.length()-ipos);
    imgInfo.path = fileName.substr(0, ipos);
    if (verbose) cout << "\tFILEPATH: " << imgInfo.path;
    if (verbose) cout << "\tFILENAME: " << imgInfo.fileName;

    // Finding stage position
    imgInfo.stagePos = stof(imgInfo.fileName.substr(imgInfo.fileName.length()-17, 8));
    if (verbose) cout << "\tSTAGE POSITION: " << imgInfo.stagePos;

    // Fill imgINFO with ordered events
    if ((curScan != imgInfo.scan) || (curRun != imgInfo.run) || (curDate != imgInfo.date)) {

      if (verbose) cout << "\n\nFilling imgINFO\n";

      for (auto itr: imgInfoMap) {
        imgINFO.push_back(itr.second);
      }
      imgInfoMap.clear();

      curScan = imgInfo.scan;
      curRun = imgInfo.run;
      curDate = imgInfo.date;
    }

    if (verbose) cout << endl;
    imgInfoMap[imgInfo.stagePos] = imgInfo;
  }

  for (auto itr: imgInfoMap) {
    imgINFO.push_back(itr.second);
  }

  fileNames.close();
  if (verbose) cout << "\nINFO: Finished retrieving info from runList!\n\n";

  // Check stagePos is ordered
  if (verbose) {
    for (uint i=0; i<imgINFO.size(); i++) {
      cout << imgINFO[i].scan << "  " << imgINFO[i].run
        << "  " << imgINFO[i].stagePos << endl;
    }
  }
 
  return true;
}




//////////////////////////////////////////////////////
/////  Creating smaller runList files if needed  /////
//////////////////////////////////////////////////////

void ppFunct::makeRunLists(std::vector<imgInfoStruct> &imgINFO) {
  ofstream outList;
  std::string curScan="";
  std::string curDate="";
  int curRun=-1;
 
  for (int ifl=0; ifl<imgINFO.size(); ifl++) {
    if ((curScan != imgINFO[ifl].scan) || (curRun != imgINFO[ifl].run)
          || (curDate != imgINFO[ifl].date)) {

      curScan = imgINFO[ifl].scan;
      curRun = imgINFO[ifl].run;
      curDate = imgINFO[ifl].date;

      if (outList.is_open()) {
        outList.close();
      }

      outList.open(("runLists/runList_Date-" + imgINFO[ifl].date
          + "_Scan-" + imgINFO[ifl].scan
          + "_Run-" + to_string(imgINFO[ifl].run)
          + ".txt").c_str());
    }
    cout<<"input for: "<<"runLists/runList_Date-" + imgINFO[ifl].date
              + "_Scan-" + imgINFO[ifl].scan
              + "_Run-" + to_string(imgINFO[ifl].run)
              + ".txt"<<endl;
    cout<<imgINFO[ifl].path + imgINFO[ifl].fileName<<endl;
    outList << imgINFO[ifl].path + imgINFO[ifl].fileName << endl;
  }
  outList.close();

  // Exiting program so you can use correct runLists
  exit(1);
}
       
