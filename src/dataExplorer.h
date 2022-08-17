#ifndef DataExplorer_H
    #define DataExplorer_H

#include <fstream>
#include <dirent.h>
#include <errno.h>
#include <string>
#include <iostream>
#include <vector>
#include <map>

using namespace std;

class DataExplorer {
    string dataDirectory;

    string proteinLinksPath = "proteinLinksPath";
    map<string, vector<string> > proteinLinks;

    string proteinAnnotationsPath = "proteinAnnotationsPath";
    map<string, vector<string> > proteinAnnotations;

    string proteinSequencesPath = "proteinSequencesPath";
    map<string, string> proteinSequences;

    string clusterSizesPath = "clusterSizesPath";
    map<string, int> clusterSizes;

    bool substring(string s1, string s2);

    public:
        DataExplorer(string _dataDirectory);

        // define static getter methods
        map<string, vector<string> > getProteinLinks() {return proteinLinks;}
        map<string, vector<string> > getProteinAnnotations() {return proteinAnnotations;}
        map<string, string> getProteinSequences() {return proteinSequences;}
        map<string, int> getClusterSizes() {return clusterSizes;}
};

#endif