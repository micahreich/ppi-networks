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
    float EDGE_WEIGHT_THRESHOLD; // add protein links if their edge weight is or exceeds this threshold
    string dataDirectory; // directory (relative file path starting with `data/`) where network files are stored

    const regex proteinNameRegex = regex("([0-9]+)(.)([a-zA-Z0-9]+)(_)[0-9]+");
    const regex proteinScoreRegex = regex("[0-9][0-9][0-9]([0-9]?)");
    const regex clusterIDRegex = regex("(CL):([0-9])+");
    const regex clusterSizeRegex = regex("[0-9]+");
    const regex aminoAcidRegex = regex("[A-Z]+");
    
    string proteinLinksPath = "proteinLinksPath";
    map<string, map<string, float> > proteinLinks;
    void populateProteinLinks(); // read in data from protein links file and populate a map of tuples

    string proteinAnnotationsPath = "proteinAnnotationsPath";
    map<string, vector<string> > proteinAnnotations;
    void populateProteinAnnotations(); // read in data from protein annotations file and populate a map of vectors

    string proteinSequencesPath = "proteinSequencesPath";
    map<string, string> proteinSequences; 
    void populateProteinSequences(); // read in data from protein sequences file and populate a map of strings

    string clusterSizesPath = "clusterSizesPath";
    map<string, int> clusterSizes;
    vector<int> clusterSizesSorted;
    void populateClusterSizes(); // read in data from clusters info file and populate a map of ints
    void populateClusterSizesSorted(); // use cluster sizes map to create a decreasing-sorted vector of cluster sizes

    bool substring(string s1, string s2); // helper function returns true iff s2 is substr of s1

    public:
        DataExplorer(string _dataDirectory, float _EDGE_WEIGHT_THRESHOLD);

        // define static getter methods
        map<string, map<string, float> > getProteinLinks() {return proteinLinks;}
        map<string, vector<string> > getProteinAnnotations() {return proteinAnnotations;}
        map<string, string> getProteinSequences() {return proteinSequences;}
        map<string, int> getClusterSizes() {return clusterSizes;}
        vector<int> getClusterSizesSorted() {return clusterSizesSorted;}
};