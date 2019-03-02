/*
Code to compute persistent path homology of a network.

Author: Samir Chowdhury
Date: ()


*/

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<algorithm>
#include<tuple>
#include<sstream>


// compile as 
// clang++ -std=c++11 computePPH.cpp -o computePPH
// ./computePPH

using namespace std;

// FUNCTION DECLARATIONS HERE

// main functions for pph
tuple<int, vector<vector<double> > > readFile(string fileName);
vector<int> getp0(int numNodes);
tuple<vector<double>, vector<vector<int>>> getp1(vector<vector<double>> &edgeData);
tuple<vector<double>, vector<vector<int>>> getp2(vector<vector<double>> &edgeData);

// helper functions
vector<vector <double> > sortByCol(int sortBy, vector<vector<double>> &unsrtMat);
void print2DVector(vector<vector<double> > &myvec);
void print2DIntVector(vector<vector<int> > &myvec);
void print1DIntVector(vector<int> &myvec);
void print1DVector(vector<double> &myvec);
tuple<bool, int> findPathInList(vector<int> &path, vector<vector<int>> &list);


// Start here
int main()
{
    int numNodes;
    vector<vector<double> > edgeData;
    string fileName;
    fileName = "cyc3.txt";

    tie(numNodes,edgeData) = readFile(fileName);
    vector<int> p0 = getp0(numNodes);
    vector<vector<int>> p1;
    vector<double> p1time;
    tie(p1time,p1) = getp1(edgeData);
    
    vector<vector<int>> p2;
    vector<double> p2time;
    tie(p2time,p2) = getp2(edgeData);
    
    print2DIntVector(p1);
    vector<int> tmp;
    tmp.push_back(2);
    tmp.push_back(3);
    print1DIntVector(tmp);

    bool found_tmp;
    int index_found_tmp;
    tie(found_tmp,index_found_tmp)=findPathInList(tmp,p1);
    cout << "found_tmp " << found_tmp << endl;
    cout << "index " << index_found_tmp << endl;
    //isPathInList(tmp,p1);
    //cout << "here goes" << endl;
    //print2DIntVector(p2);
    //print1DVector(p2time);
    return 0;
}

tuple<bool, int> findPathInList(vector<int> &path, vector<vector<int>> &list){
/* function that searches for p-path in list of p-paths
If found, return (found, index of path in list).
Otherwise return (!found, list.end())
*/
    // initialize return objects
    bool found;
    int index;
    // create iterator and perform find
    vector<vector<int>>::iterator it;
    it = find(list.begin(), list.end(), path);
    if (it != list.end()){
        found = true;
        index = distance(list.begin(),it);
        //cout << "Found item at position " << index << endl;
    }
    else {
        found = false;
        index = distance(list.begin(),it);
        //cout << "item not found. distance" << index << endl;
    }
    return make_tuple(found,index);
}




tuple<vector<double>, vector<vector<int>>> getp2(vector<vector<double> > &edgeData){
    
    vector<vector<double>> p2temp;

    vector<double> temp(4);
    //double tmpTime;

    vector<double> helpTemp(2);

    for (int i = 0; i < edgeData.size(); i++){
        helpTemp[0] = edgeData[i][0];
        helpTemp[1] = edgeData[i][1];
        for (int j = 0; j < edgeData.size();j++){
            if (edgeData[j][0] == helpTemp[1]){
                temp[0] = helpTemp[0];//first vertex
                temp[1] = helpTemp[1];//second vertex
                temp[2] = edgeData[j][1];//third vertex
                //max of first/second pair
                temp[3] = max(edgeData[j][2],edgeData[i][2]);
                p2temp.push_back(temp);
                //p2time.push_back();
                //cout << helpTemp[0] << helpTemp[1] << edgeData[j][1] << endl;
            }
        }
    }
    //print2DVector(p2temp);
    p2temp = sortByCol(3,p2temp);
    //cout << "try again" << endl;
    //print2DVector(p2temp);
    vector<vector<int>> p2;
    vector<double> p2time(p2temp.size());
    vector<int> tmp(3); //tmp(3) for 2-paths
    
    for (int i = 0; i < p2temp.size() ; ++i){
        p2time[i] = p2temp[i].back();
        for (int j = 0; j < 3; ++j){
            // j< 3 above for 2-paths
            tmp[j] = p2temp[i][j]; 
        }
        p2.push_back(tmp);
        
    }

    return make_tuple(p2time,p2);
}



tuple<vector<double>, vector<vector<int>>> getp1(vector<vector<double> > &edgeData){
    // get 1 paths and entry times
    vector<vector<int>> p1;
    vector<double> p1time;
    vector<int> temp(2);
    double tmpTime;
    for (int i =0; i < edgeData.size() ; i++){
        temp[0] = edgeData[i][0];
        temp[1] = edgeData[i][1];
        tmpTime = edgeData[i][2];
        p1.push_back(temp);//append
        p1time.push_back(tmpTime);
    }
    return make_tuple(p1time,p1);
}

vector<int> getp0(int numNodes){
// get 0 paths
    vector<int> p0;
    for (int i=1; i<= numNodes; i++){
        p0.push_back(i);
    }
    return p0;
}



tuple<int, vector<vector<double> > > readFile(string fileName){
/*
Function that reads a provided text file.
Top row of file contains number of nodes
Remaining rows contain directed edge pairs and 
entry time.
Important: node labels are 1, 2, 3, ... (do not start from 0!)
Multiple values returned, hence using tuple class 
edgeData returned sorted by time
*/
    // create ifstream object
    ifstream file;
    string currentLine;

    // initialize objects that will be returned
    int numNodes;
    vector<vector<double>> edgeData;
    double a, b,t ;

    file.open(fileName);
    getline(file,currentLine);//first row has number of nodes
    numNodes = stoi(currentLine);//convert to integer
    
    // now append each row to vector of vectors
    vector<double> tempVect(3); //to store each (edge, time) row
    while (getline(file,currentLine)){
        stringstream lineSS(currentLine); //build stringstream object
        lineSS >> a >> b >> t;//pass stringstream into variables
        tempVect[0] = a;
        tempVect[1] = b;
        tempVect[2] = t;
        edgeData.push_back(tempVect);//append
    } 
    file.close();
    edgeData = sortByCol(2,edgeData);//sort by time
    return make_tuple(numNodes, edgeData);
}


void print2DVector(vector<vector<double> > &myvec){
/* Little helper function for printing 
out all the values of a 2d vector of doubles
*/

    for (int i = 0; i < myvec.size(); i++){
        for (int j = 0; j < myvec[i].size(); j++){
            cout << myvec[i][j] << " ";
        }
        cout << endl;
    }
}

void print2DIntVector(vector<vector<int> > &myvec){
/* Little helper function for printing 
out all the values of a 2d vector of ints
*/

    for (int i = 0; i < myvec.size(); i++){
        for (int j = 0; j < myvec[i].size(); j++){
            cout << myvec[i][j];
        }
        cout << endl;
    }
}

void print1DIntVector(vector<int> &myvec){
/* Print all values of int vector
*/
    for (int i = 0; i < myvec.size(); i++){
        cout << myvec[i] << endl;
    }
}

void print1DVector(vector<double> &myvec){
/* Print all values of int vector
*/
    for (int i = 0; i < myvec.size(); i++){
        cout << myvec[i] << endl;
    }
}



vector<vector <double> > sortByCol(int sortBy, vector<vector <double> > &unsrtMat){
/*
    Function that sorts a matrix of doubles by a particular column. E.g.
        121
        232
        311 becomes
        121
        311
        232
    Uses sort() from <algorithm> library.
    Third argument in sort() is a comparator; in this case
    it is a lambda expression.
*/
    
    // [sortBy] puts variable inside capture-list
    sort(unsrtMat.begin(),unsrtMat.end(), 
    [sortBy] (const std::vector<double> &a, const std::vector<double> &b)
          {
              return a[sortBy] < b[sortBy];
          });

    return unsrtMat;
}
