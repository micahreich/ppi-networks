#include "dataExplorer.h"
#include <fstream>
#include <dirent.h>
#include <errno.h>
#include <iostream>
#include <vector>
#include <regex>
#include <tuple>
#include <string>
#include <chrono>

using namespace std;

string moduleName = "[DataExplorer] ";
int STATUS_CODE = 0;
double MS_TO_S = 1 / 1e3; 

bool DataExplorer::substring(string s1, string s2) {
    return s2.find(s1) != string::npos;
}

DataExplorer::DataExplorer(string _dataDirectory, float _EDGE_WEIGHT_THRESHOLD) {
    string* requiredFiles[4] = {&proteinLinksPath, &proteinAnnotationsPath, &proteinSequencesPath, &clusterSizesPath};
    dataDirectory = _dataDirectory;
    EDGE_WEIGHT_THRESHOLD = _EDGE_WEIGHT_THRESHOLD;

    DIR* dp = opendir(_dataDirectory.c_str());
    struct dirent* dirp;

    if(dp == NULL) {
        cout << moduleName << "Error(" << errno << ") opening " << dataDirectory << endl;
        STATUS_CODE = 1;
    } else {
        cout << moduleName << "Retrieving data from " << dataDirectory << "..." << endl;
    }

    cout << string(80, '.') << endl;

    // iterate through data directory, assign file paths to appropriate data depending on names
    while((dirp = readdir(dp)) != NULL) {
        string currFileName = string(dirp->d_name);
        
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

    // populate protein links
    cout << string(80, '.') << endl;

    auto start = chrono::high_resolution_clock::now();
    cout << moduleName << "Populating protein links..." << endl;

    // populateProteinLinks();

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start).count() * MS_TO_S;
    cout << moduleName << "Finished! Took " << duration << " seconds " << endl;

    // populate protein annotations
    start = chrono::high_resolution_clock::now();
    cout << moduleName << "Populating protein annotations..." << endl;
    
    // populateProteinAnnotations();

    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start).count() * MS_TO_S;
    cout << moduleName << "Finished! Took " << duration << " seconds " << endl;

    // populate protein sequences
    start = chrono::high_resolution_clock::now();
    cout << moduleName << "Populating protein AA sequences..." << endl;
    
    // populateProteinSequences();

    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start).count() * MS_TO_S;
    cout << moduleName << "Finished! Took " << duration << " seconds " << endl;

    // populate cluster sizes
    start = chrono::high_resolution_clock::now();
    cout << moduleName << "Populating cluster sizes..." << endl;
    
    populateClusterSizes();

    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start).count() * MS_TO_S;
    cout << moduleName << "Finished! Took " << duration << " seconds " << endl;

    // populate sorted cluster sizes
    start = chrono::high_resolution_clock::now();
    cout << moduleName << "Populating sorted cluster sizes array..." << endl;

    populateClusterSizesSorted();

    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start).count() * MS_TO_S;
    cout << moduleName << "Finished! Took " << duration << " seconds " << endl;

    // fin
    cout << string(80, '.') << endl;
}

void DataExplorer::populateProteinLinks() {
    ifstream file(proteinLinksPath);

    if(file.is_open()) {
        int cnt = 0;

        string line;
        const string delim = " ";

        vector<string> tokenizedLine(3);

        while(getline(file, line)) {
            size_t pos = 0;
            size_t numTokens = 0;

            tokenizedLine = vector<string>{"", "", ""};
            
            while((pos = line.find(delim)) != string::npos) {
                string token = line.substr(0, pos);
                tokenizedLine[numTokens] = token;
                ++numTokens;

                line.erase(0, pos + delim.length());
            }
            tokenizedLine[2] = line; // set final token to the link score

            // verify that the tokenized line is a correct representation of a protein link
            string currProtein1, currProtein2;
            float currScore;

            if(regex_match(tokenizedLine[0], proteinNameRegex) &&
               regex_match(tokenizedLine[1], proteinNameRegex) &&
               regex_match(tokenizedLine[2], proteinScoreRegex)) {
                currProtein1 = tokenizedLine[0];
                currProtein2 = tokenizedLine[1];
                
                int currScoreRaw =  stoi(tokenizedLine[2]);
                currScore = (currScoreRaw / 1e2);

                if(proteinLinks.find(currProtein1) != proteinLinks.end()) {
                    // extend entry in the hashmap
                    if(currScore >= EDGE_WEIGHT_THRESHOLD) {
                        proteinLinks[currProtein1].push_back(make_tuple(currProtein2, currScore));
                    }
                } else {
                    // make a new entry in the hashmap
                    if(currScore >= EDGE_WEIGHT_THRESHOLD) {
                        proteinLinks[currProtein1] = vector<tuple <string,float> >{
                            make_tuple(currProtein2, currScore)
                        };
                    } else {
                        proteinLinks[currProtein1] = vector<tuple <string,float> >{};
                    }
                }
            }            
        }
        file.close();
    }
}

void DataExplorer::populateProteinAnnotations() {

    ifstream file(proteinAnnotationsPath);

    if(file.is_open()) {
        string line;
        const string delim = " ";

        vector<string> tokenizedLine(3);

        while(getline(file, line)) {
            size_t pos = 0;
            size_t numTokens = 0;

            tokenizedLine = vector<string>{"", "", ""};
            
            while((pos = line.find(delim)) != string::npos) {
                string token = line.substr(0, pos);
                tokenizedLine[numTokens] = token;
                ++numTokens;

                line.erase(0, pos + delim.length());
            }
            tokenizedLine[2] = line; // set final token to the link score

            // verify that the tokenized line is a correct representation of a protein link
            string currProtein, currCluster;

            if(regex_match(tokenizedLine[1], clusterIDRegex) &&
               regex_match(tokenizedLine[2], proteinNameRegex)) {

                currProtein = tokenizedLine[2];

                size_t colonDelimIdx = tokenizedLine[1].find(":");
                currCluster = tokenizedLine[1].substr(colonDelimIdx + 1, 
                                                      tokenizedLine[1].length() - (colonDelimIdx + 1));

                cout << currProtein << " " << currCluster << endl;
                if(proteinAnnotations.find(currProtein) != proteinAnnotations.end()) {
                    // extend entry in the hashmap
                    proteinAnnotations[currProtein].push_back(currCluster);
                } else {
                    // make a new entry in the hashmap
                    proteinAnnotations[currProtein] = vector<string>{ currCluster };
                }
            }
        }
        file.close();
    }
}

void DataExplorer::populateProteinSequences() {
    ifstream file(proteinSequencesPath);

    if(file.is_open()) {
        string line;
        const string delim = " ";

        vector<string> tokenizedLine(2);

        while(getline(file, line)) {
            size_t pos = 0;
            size_t numTokens = 0;

            tokenizedLine = vector<string>{"", ""};
            
            while((pos = line.find(delim)) != string::npos) {
                string token = line.substr(0, pos);
                tokenizedLine[numTokens] = token;
                ++numTokens;

                line.erase(0, pos + delim.length());
            }
            tokenizedLine[1] = line; // set final token to the link score

            // verify that the tokenized line is a correct representation of a protein link
            string currProtein, currSequence;

            if(regex_match(tokenizedLine[0], proteinNameRegex) &&
               regex_match(tokenizedLine[1], aminoAcidRegex)) {
                currProtein = tokenizedLine[0];
                currSequence = tokenizedLine[1];

                // make a new entry in the hashmap
                proteinSequences[currProtein] = currSequence;
            }
        }
        file.close();
    }
}

void DataExplorer::populateClusterSizes() {
    ifstream file(clusterSizesPath);

    if(file.is_open()) {
        string line;
        const string delim = " ";

        vector<string> tokenizedLine{};

        while(getline(file, line)) {
            size_t pos = 0;
            size_t numTokens = 0;

            tokenizedLine = vector<string>{};
            
            while((pos = line.find(delim)) != string::npos) {
                string token = line.substr(0, pos);
                tokenizedLine.push_back(token);
                ++numTokens;

                line.erase(0, pos + delim.length());
            }
            tokenizedLine.push_back(line); // set final token to the link score

            // verify that the tokenized line is a correct representation of a protein link
            string currCluster;
            int currClusterSize;

            if(regex_match(tokenizedLine[1], clusterIDRegex) &&
               regex_match(tokenizedLine[2], clusterSizeRegex)) {
                currClusterSize = stoi(tokenizedLine[2]);

                size_t colonDelimIdx = tokenizedLine[1].find(":");
                currCluster = tokenizedLine[1].substr(colonDelimIdx + 1, 
                                                      tokenizedLine[1].length() - (colonDelimIdx + 1));

                // make a new entry in the hashmap
                clusterSizes[currCluster] = currClusterSize;
            }
        }
        file.close();
    }
}

void DataExplorer::populateClusterSizesSorted() {
    DataExplorer::clusterSizesSorted = vector<int>{};

    for (const auto &clusterPair : DataExplorer::clusterSizes) {
        clusterSizesSorted.push_back(clusterPair.second);
    }

    sort(clusterSizesSorted.begin(), clusterSizesSorted.end(), greater<int>());
}

int main() {
    DataExplorer d("data/ecoli", 0);
    // vector<int> s = d.getClusterSizesSorted();
    // cout << s.size() << endl;
    // cout << s[0] << " " << s[s.size() - 1] << endl;
}