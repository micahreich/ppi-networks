#include "dataExplorer.h"
#include <fstream>
#include <dirent.h>
#include <errno.h>
#include <string>
#include <iostream>
#include <vector>

using namespace std;

string moduleName = "[DataExplorer] ";
int STATUS_CODE = 0;

DataExplorer::DataExplorer(string _dataDirectory) {
    string* requiredFiles[4] = {&proteinLinksPath, &proteinAnnotationsPath, &proteinSequencesPath, &clusterSizesPath};
    dataDirectory = _dataDirectory;

    DIR* dp = opendir(_dataDirectory.c_str());
    struct dirent* dirp;

    if(dp == NULL) {
        cout << moduleName << "Error(" << errno << ") opening " << dataDirectory << endl;
        STATUS_CODE = 1;
    } else {
        cout << moduleName << "Retrieving data from " << dataDirectory << endl;
    }

    while((dirp = readdir(dp)) != NULL) {
        string currFileName = string(dirp->d_name);
        
        // iterate through data directory, assign file paths to appropriate data depending on names
        if(substring("protein.links", currFileName)) {
            proteinLinksPath = dataDirectory + "/" + currFileName;
            cout << moduleName << "Found: " <<  proteinLinksPath << " !" << endl;
        } else if(substring("clusters.proteins", currFileName)) {
            proteinAnnotationsPath = dataDirectory + "/" + currFileName;
            cout << moduleName << "Found: " <<  proteinAnnotationsPath << " !" << endl;
        } else if(substring("protein.sequences", currFileName)) {
            proteinSequencesPath = dataDirectory + "/" + currFileName;
            cout << moduleName << "Found: " <<  proteinSequencesPath << " !" << endl;
        } else if(substring("clusters.info", currFileName)) {
            clusterSizesPath = dataDirectory + "/" + currFileName;
            cout << moduleName << "Found: " <<  clusterSizesPath << " !" << endl;
        }
    }

    for(string* fname : requiredFiles) {
        string fnameString = *fname;
        if(!substring(dataDirectory, fnameString)) {
            cout << moduleName << "Error! No file path found for " << fnameString << " !" << endl;
            STATUS_CODE = 1;
        }
    }
}

bool DataExplorer::substring(string s1, string s2) {
    return s2.find(s1) != string::npos;
}

int main() {
    DataExplorer d("data/Ecoli");
}