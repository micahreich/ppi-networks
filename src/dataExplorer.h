#ifndef DataExplorer_H
    #define DataExplorer_H

#include <fstream>
#include <dirent.h>
#include <errno.h>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <tuple>
#include <regex>

using namespace std;

class DataExplorer {
    string dataDirectory;

    const regex proteinNameRegex = regex("([0-9]+)(.)([a-zA-Z0-9]+)(_)[0-9]+");
    const regex proteinScoreRegex = regex("[0-9][0-9][0-9]([0-9]?)");
    const regex clusterIDRegex = regex("(CL):([0-9])+");
    const regex clusterSizeRegex = regex("[0-9]+");
    const regex aminoAcidRegex = regex("[A-Z]+");
    

    string proteinLinksPath = "proteinLinksPath";
    map<string, vector< tuple<string, int> > > proteinLinks;
    void populateProteinLinks();

    string proteinAnnotationsPath = "proteinAnnotationsPath";
    map<string, vector<string> > proteinAnnotations;
    void populateProteinAnnotations();

    string proteinSequencesPath = "proteinSequencesPath";
    map<string, string> proteinSequences;
    void populateProteinSequences();

    string clusterSizesPath = "clusterSizesPath";
    map<string, int> clusterSizes;
    void populateClusterSizes();

    bool substring(string s1, string s2);
    void timedFunction(DataExplorer* ob, void (DataExplorer::*fn)());

    public:
        DataExplorer(string _dataDirectory);

        // define static getter methods
        map<string, vector< tuple<string,int> > > getProteinLinks() {return proteinLinks;}
        map<string, vector<string> > getProteinAnnotations() {return proteinAnnotations;}
        map<string, string> getProteinSequences() {return proteinSequences;}
        map<string, int> getClusterSizes() {return clusterSizes;}
};

#endif