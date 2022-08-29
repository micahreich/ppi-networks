#include "functionPrediction.h"
#include <string>
#include <vector>
#include <map>
#include <tuple>
#include <unordered_set>
#include <queue>
#include <algorithm>

using namespace std;

string moduleName = "[FunctionPrediction] ";

// map<string, vector<string> > t = {make_pair("1", vector<string>{"2"}),
//                                   make_pair("2", vector<string>{"3", "A"}),
//                                   make_pair("3", vector<string>{"4"}),
//                                   make_pair("4", vector<string>{"5", "3"}),
//                                   make_pair("5", vector<string>{"6"}),
//                                   make_pair("6", vector<string>{"7"}),
//                                   make_pair("A", vector<string>{"B"})};

FunctionPrediction::FunctionPrediction(DataExplorer* _de) {
    de = _de;
}

unordered_set<string> FunctionPrediction::getNeighborhood(string protein, int radius, bool lastLevel=false) {
    DataExplorer dataExplorer = *de;
    auto proteinLinks = dataExplorer.getProteinLinks();

    queue<string> bfsQueue;
    bfsQueue.push(protein);

    unordered_set<string> neighbors;

    unordered_set<string> visited; 
    visited.insert(protein);

    int currDepth = 1;

    while (!bfsQueue.empty() && currDepth <= radius) {
        int levelSize = bfsQueue.size();

        while (levelSize > 0) {
            string currProtein = bfsQueue.front(); bfsQueue.pop(); // return and remove the top element
            if(!lastLevel) {
                neighbors.insert(currProtein);
            } else if(lastLevel && currDepth == radius) {
                neighbors.insert(currProtein);
            }

            for(auto& neighborPair : proteinLinks[currProtein]) {
                string neighbor = get<0>(neighborPair);

                if(!visited.count(neighbor)) {
                    visited.insert(neighbor);
                    bfsQueue.push(neighbor);
                }
            }

            // for(string neighbor : ::t[currProtein]) {
            //     if(visited.find(neighbor) == visited.end()) {
            //         visited.insert(neighbor);
            //         bfsQueue.push(neighbor);
            //     }
            // }

            levelSize--;
        }

        currDepth++;
    }

    return neighbors;
}

int FunctionPrediction::expectedNumOfClusters(string clusterID, int n) {
    return (int)(n * ((float)(de->getClusterSizes()[clusterID]) / (float)(de->getProteinLinks().size())));
}

map<string, vector<string> > FunctionPrediction::majorityApproach(vector<string> unknownProteinList, int radius) {
    string approachName = "[Majority Approach] ";
    map<string, vector<string> > unknownProteinAssignments;

    cout << string(80, '.') << endl;

    int proteinCount = 1;

    for(string protein : unknownProteinList) {
        unordered_set<string> neighborhood = getNeighborhood(protein, radius);

        for(string& unknownProtein : unknownProteinList) { // removes all unknown proteins from the neighborhood
            if(neighborhood.count(unknownProtein)) {
                neighborhood.erase(unknownProtein);
            }
        }

        cout << approachName << "Protein: " << protein << "(" << proteinCount << "/" << unknownProteinList.size() << ")" << endl;
        cout << string(approachName.size(), ' ') << "Neighborhood size: " << neighborhood.size() << endl;

        vector<string> observedClustersList;
        unordered_set<string> observedClustersSet;

        vector<tuple<string, float> > observedClustersScores;

        for(string neighbor : neighborhood) { // list all observed clusters
            vector<string> neighborAnnotations = de->getProteinAnnotations()[neighbor];
            observedClustersList.insert(observedClustersList.end(), neighborAnnotations.begin(), neighborAnnotations.end()); // add to observed clusters

            for(string& annotation : neighborAnnotations) {
                observedClustersSet.insert(annotation);
            }
        }

        for(string cluster : observedClustersSet) { // populate observedClustersFreq with frequency of all observed clusters
            int observedClusterFreq = count(observedClustersList.begin(),
                                            observedClustersList.end(), cluster);
            int expectedClusterFrequency = expectedNumOfClusters(cluster, neighborhood.size());

            float clusterScore = (float)(observedClusterFreq - expectedClusterFrequency) / expectedClusterFrequency;
            observedClustersScores.push_back(make_pair(cluster, clusterScore));
        }

        // sort vector of observed cluster scores 
        sort(observedClustersScores.begin(), observedClustersScores.end(), [](const tuple<string, float>& t1, // lambda function to compare tuples
                                                                              const tuple<string, float>& t2) {
                                                                                return (get<1>(t1) > get<1>(t2)); 
                                                                              });

        // extract winners from the descending-sorted vector of tuples
        vector<string> winners = {};

        for(auto& clusterPair : observedClustersScores) {
            winners.push_back(get<0>(clusterPair));
        }

        // assign current protein this winners vector
        unknownProteinAssignments[protein] = winners;

        proteinCount++;
    }

    return unknownProteinAssignments;
}

// reservoirs dictionaries contains keys with protein, values as dictionaries with clusters as keys & values as water amounts
map<string, map<string, float>> reservoirsPrev{};
map<string, map<string, float>> reservoirsCurr{};

map<string, vector<string> > FunctionPrediction::functionalFlowApproach(vector<string> unknownProteinList, int radius) {
    string approachName = "[Functional Flow Approach] ";
    map<string, vector<string> > unknownProteinAssignments;

    map<string, map<string, float>> totalEnteredAmounts{};

    cout << string(80, '.') << endl;

    int proteinCount = 1;

    unordered_set<string> knownNeighborhood;
    unordered_set<string> knownNeighborhoodClusters;

    // make totalNeighborhood the union of all radius `radius` neighborhoods for each unknown protein
    for(string& unknownProtein : unknownProteinList) {
        unordered_set<string> neighborhood = getNeighborhood(unknownProtein, radius);
        knownNeighborhood.insert(neighborhood.begin(), neighborhood.end());
    }

    for(string& unknownProtein : unknownProteinList) { // removes all unknown proteins from the known neighborhood
        if(knownNeighborhood.count(unknownProtein)) {
            knownNeighborhood.erase(unknownProtein);
        }
    }

    for(string knownProtein : knownNeighborhood) {
        vector<string> knownProteinAnnotation = de->getProteinAnnotations()[knownProtein];
        knownNeighborhoodClusters.insert(knownProteinAnnotation.begin(), knownProteinAnnotation.end());
    }

    // populate iteration 0 starting water amounts for known proteins
    for(string knownProtein : knownNeighborhood) {
        map<string, float> proteinFunctionMap;
        unordered_set<string> negativeClusters = knownNeighborhoodClusters;

        vector<string> knownProteinAnnotation = de->getProteinAnnotations()[knownProtein];

        // fill clusters of knownProtein with MAX water
        for(string cluster : knownProteinAnnotation) {
            negativeClusters.erase(cluster);
            proteinFunctionMap[cluster] = (float) INT_MAX;
        }

        // fill non-clusters of knownProtein with zero water
        for(string cluster : negativeClusters) {
            proteinFunctionMap[cluster] = (float) 0;
        }

        reservoirsPrev[knownProtein] = proteinFunctionMap;
    }
    
    // populate iteration 0 starting water amounts for unknown proteins
    for(string unknownProtein : unknownProteinList) {
        map<string, float> proteinFunctionMap;

        for(string cluster : knownNeighborhoodClusters) {
            proteinFunctionMap[cluster] = (float) 0;
        }

        reservoirsPrev[unknownProtein] = proteinFunctionMap;
    }

    // begin main algorithm
    unordered_set<string> iteratedProteins = knownNeighborhood;
    iteratedProteins.insert(unknownProteinList.begin(), unknownProteinList.end());

    map<string, map<string, float> > proteinLinks = de->getProteinLinks();
    
    cout << approachName << " Working neighborhood size: " << iteratedProteins.size() << endl;
    
    for(int t = 1; t <= radius; t++) {
        cout << approachName << " Iteration (" << t << "/" << radius << ")" << endl;

        // iterate through all considered proteins
        for(string protein : iteratedProteins) {
            // iterate through all observed clusters
            for(string cluster : knownNeighborhoodClusters) {
                // determine the amount entered/left for this cluster for this protein at iteration `i`
                float amtEntered = 0;

                for(auto& proteinLink : proteinLinks[protein]) {
                    string neighborProtein = proteinLink.first;

                    if(iteratedProteins.count(neighborProtein)) {
                        amtEntered += edgeCapacity(neighborProtein, protein, cluster, t);
                    }
                }

                float amtLeft = 0;

                for(auto& proteinLink : proteinLinks[protein]) {
                    string neighborProtein = proteinLink.first;

                    if(iteratedProteins.count(neighborProtein)) {
                        amtLeft += edgeCapacity(neighborProtein, protein, cluster, t);
                    }
                }

                // the amount of water in this reservoir = water previously + (amount entered - amount left)
                reservoirsCurr[protein][cluster] = reservoirsPrev[protein][cluster] + (amtEntered - amtLeft);

                // update the score of this cluster by adding the amount entered for this cluster
                if(count(unknownProteinList.begin(), unknownProteinList.end(), protein)) {
                    totalEnteredAmounts[protein][cluster] += amtEntered;
                }
            }
        }

        reservoirsPrev = reservoirsCurr;
        reservoirsCurr.clear();
    }

    // convert dictionary of dictionaries into dictionary of vector of tuples
    map<string, vector<tuple<string, float> > > totalEnteredAmountsVectorized;

    for(string& protein : unknownProteinList) {
        map<string, float> proteinTotalEnteredAmounts = totalEnteredAmounts[protein];
        vector<tuple<string, float> > proteinTotalEnteredAmountsVectorized;

        for(auto& clusterPair : proteinTotalEnteredAmounts) {
            proteinTotalEnteredAmountsVectorized.push_back(
                make_tuple(clusterPair.first, clusterPair.second)
            );
        }

        totalEnteredAmountsVectorized[protein] = proteinTotalEnteredAmountsVectorized;
    }

    // sort vector of tuples and put winning proteins into assignments in sorted order without scores
    for(auto& proteinClusterPair : totalEnteredAmountsVectorized) {
        string protein = proteinClusterPair.first;
        vector<tuple<string, float> >& clusterScores = totalEnteredAmountsVectorized[protein];

        sort(clusterScores.begin(), clusterScores.end(), [](const tuple<string, float>& t1, // lambda function to compare tuples
                                                            const tuple<string, float>& t2) {
                                                                return (get<1>(t1) > get<1>(t2)); 
                                                            });
        
        vector<string> winners{};

        for(auto& clusterScorePair : clusterScores) {
            string clusterID = get<0>(clusterScorePair);
            winners.push_back(clusterID);
        }

        unknownProteinAssignments[protein] = winners;
    }

    return unknownProteinAssignments;
} 

float FunctionPrediction::edgeCapacity(string proteinA, string proteinB, string clusterID, int t) {
    if(t == 0) { return 0; }

    if(reservoirsPrev[proteinA][clusterID] <= reservoirsPrev[proteinB][clusterID]) {
        return 0;
    } else {
        map<string, map<string, float> > proteinLinks = de->getProteinLinks();

        float edgeCapacity = proteinLinks[proteinA][proteinB];
        return edgeCapacity;

        /*
        // update proteinA proteinB edge capacity
        if(proteinLinks[proteinA].find(proteinB) != proteinLinks[proteinA].end()) {
            edgeCapacity = proteinLinks[proteinA][proteinB];
        }
        */
    }
}

int main() {
    // unordered_set<string> res = getNeighborhood("1", 5, true);

    // for(auto& s : res) {
    //     cout << s << " ";
    // }
    // cout << endl;
    map<string, float> t;

    int x = 1;
    int y = 2;
    int z = 5;

    t["a"] = (float)(x - y) / z;

    // cout << t["a"] << endl;

    unordered_set<int> set1 = {1, 2};
    unordered_set<int> set2 = set1;
    set2.erase(2);

    for(int i : set2) {
        cout << i << " ";
    } cout << endl;

    return 0;
}