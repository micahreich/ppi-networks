#include "dataExplorer.h"
#include <string>
#include <vector>
#include <map>
#include <tuple>
#include <unordered_set>

using namespace std;

class FunctionPrediction {
    DataExplorer* de;
    
    // perform breadth first search on the network to get all proteins within distance `radius` of `protein`
    unordered_set<string> getNeighborhood(string protein, int radius, bool lastLevel);

    // majority approach helper functions
    int expectedNumOfClusters(string clusterID, int n);
    float edgeCapacity(string proteinA, string proteinB, string clusterID, int t);

    public:
        FunctionPrediction(DataExplorer* _de); // class takes in pointer to DataExplorer object made using `new` keyword

        // perform majority appraoch at radius `radius` for all proteins in `unknownProteinList`
        map<string, vector<string> > majorityApproach(vector<string> unknownProteinList, int radius);

        // perform functional flow appraoch at radius `radius` for all proteins in `unknownProteinList`
        map<string, vector<string> > functionalFlowApproach(vector<string> unknownProteinList, int radius);

        // perform alignment appraoch at radius `radius` for all proteins in `unknownProteinList`
        map<string, vector<string> > alignmentApproach(vector<string> unknownProteinList, int radius);
};