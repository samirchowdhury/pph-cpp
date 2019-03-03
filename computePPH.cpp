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
vector<vector<int>> getp0(int numNodes);
tuple<vector<double>, vector<vector<int>>> getp1(vector<vector<double>> &edgeData);
tuple<vector<double>, vector<vector<int>>> getp2(vector<vector<double>> &edgeData);
vector<vector<int>> simpleBoundary(vector<int> &path);
void addInfAPaths(vector<vector<int>> &p1, vector<double> &p1time, vector<vector<int>> &p2);

tuple<vector<vector<int>>, int>  basisChange(vector<int> &kpath, 
vector<bool> &pjmarked, vector<vector<int>> &pj, vector<vector<vector<int>>> &tj);

vector<vector<int>> removeUnmarked(vector<vector<int>> &pjsummands, vector<bool> &pjmarked, vector<vector<int>> &pj);

tuple<int, int> getMaxIdx(vector<vector<int>> &pjsummands, vector<vector<int>> &pj);


/*void buildSlots(vector<vector<int>> &p0, vector<vector<int>> &p1, vector<vector<int>> &p2,
vector<vector<vector<int>>> &t0, vector<vector<vector<int>>> &t1, vector<vector<vector<int>>> &t2);
*/

// helper functions
vector<vector <double> > sortByCol(int sortBy, vector<vector<double>> &unsrtMat);
void print2DVector(vector<vector<double> > &myvec);
void print2DIntVector(vector<vector<int> > &myvec);
void print1DIntVector(vector<int> &myvec);
void print1DVector(vector<double> &myvec);
void print1DBoolVector(vector<bool> &myvec);
tuple<bool, int> findPathInList(vector<int> &path, vector<vector<int>> &list);
bool hasRepeats(vector<int> &myvec);
vector<vector<int>> symDiff(vector<vector<int>> &arr1, vector<vector<int>> &arr2);

// Global variable
const int maxTime = 10000;

// Start here
int main()
{
    
    int numNodes;
    vector<vector<double> > edgeData;
    string fileName;
    fileName = "cyc3.txt";

    tie(numNodes,edgeData) = readFile(fileName);
    vector<vector<int>> p0 = getp0(numNodes);
    vector<vector<int>> p1;
    vector<double> p1time;
    tie(p1time,p1) = getp1(edgeData);
    
    vector<vector<int>> p2;
    vector<double> p2time;
    tie(p2time,p2) = getp2(edgeData);
    
    //print2DIntVector(p1);
    /*
    vector<int> tmp;
    tmp.push_back(1);
    tmp.push_back(2);
    tmp.push_back(1);
    */
    //print1DIntVector(tmp);

    //bool found_tmp;
    //int index_found_tmp;
    //tie(found_tmp,index_found_tmp)=findPathInList(tmp,p1);

    //cout << "found_tmp " << found_tmp << endl;
    //cout << "index " << index_found_tmp << endl;
    //isPathInList(tmp,p1);


    //vector<vector<int>> summands = simpleBoundary(tmp);
    //print2DIntVector(summands);

    addInfAPaths(p1,p1time,p2);
    //print2DIntVector(p1);
    //print1DVector(p1time);
    //print2DIntVector(p0);



/*
Marked paths
Different from python implementation.
Instead of checking if a summand is in 
the pool of marked paths, we actually find
the index of the summand in the corresponding list.
then we set that index to true (to denote marked)
in the vector of bools labeled pXmarked.
Indices of pXmarked agree with pX.
*/
    // all p0 are marked.
    vector<bool> p0marked(p0.size(),true);
    vector<bool> p1marked(p1.size(),false);
    vector<bool> p2marked(p2.size(),false);

// Build slots
    vector<vector<vector<int>>> t0(p0.size());
    vector<vector<vector<int>>> t1(p1.size());
    vector<vector<vector<int>>> t2(p2.size());
    
// this goes inside loop
    print2DIntVector(p1);

    // adding 0 at end bc all are 0-paths
    vector<vector<int>> summands0;
    int maxIdx0;

    int j = 0;

    tie(summands0, maxIdx0) = basisChange(p1[0],p0marked, p0,t0);
    t0[maxIdx0] = summands0;
    cout << "at t0 index " << maxIdx0 << " we have " << endl;
    print2DIntVector(t0[maxIdx0]);
    cout << "next"<<endl;

    j = 1;

    tie(summands0, maxIdx0) = basisChange(p1[1],p0marked, p0,t0);
    t0[maxIdx0] = summands0;
    cout << "at t0 index " << maxIdx0 << " we have " << endl;
    print2DIntVector(t0[maxIdx0]);
    cout << "next"<<endl;

    j = 2;

    tie(summands0, maxIdx0) = basisChange(p1[2],p0marked, p0,t0);
    if (summands0.empty()){
        cout << "empty sum" << endl;
        p1marked[j] = true;

        t0[maxIdx0] = summands0;
        cout << "at t0 index " << maxIdx0 << " we have " << endl;
        print2DIntVector(t0[maxIdx0]);
    }
    
    cout << "marked 0" << endl;
    print1DBoolVector(p0marked);
    cout << "marked 1" << endl;
    print1DBoolVector(p1marked);





    return 0;
}


tuple<vector<vector<int>>, int> basisChange(vector<int> &kpath, 
vector<bool> &pjmarked, vector<vector<int>> &pj, vector<vector<vector<int>>> &tj){
/*
Column reduction operation. Here j = k-1.
Take a k-path, compute its boundary summands.

Differs from PPH algorithm 2, only returns u and i.
u = "good" summands
i = max index of a summand

Thus need an extra function to compute et, the entry time.
*/
    int maxIdx;
    vector<vector<int>> summands;
    summands = simpleBoundary(kpath);
/*
Remove unmarked summands. All 0-paths are marked, so don't
need this if k = 1;
*/

    if (kpath.size()>1){
        summands = removeUnmarked(summands, pjmarked, pj);
    }
    print1DIntVector(kpath);
    cout << "boundary" << endl;
    print2DIntVector(summands);

    // Loop will start here
    while (!summands.empty()){ 
    /*
    Find term with highest index
    */
        int summandIdx;
        tie(maxIdx,summandIdx)=getMaxIdx(summands,pj);
        cout << "maxidx is " << maxIdx << endl;
        cout << "summandidx is " << summandIdx << endl;


        cout << "is spot empty?" << tj[maxIdx].empty() << endl;

        if (tj[maxIdx].empty()){
            // if empty, break out of loop
            break;
        }
        // column reduction step
        // 
        summands = symDiff(summands, tj[maxIdx]);
        cout << "summands are " << endl;
        print2DIntVector(summands);
    }

    if (summands.empty()){
        cout << "ok got empty summand" << endl;
    }

    return make_tuple(summands, maxIdx);
}

vector<vector<int>> symDiff(vector<vector<int>> &arr1, vector<vector<int>> &arr2){
/* builds the symmetric difference of two lists of j-paths
*/
    vector<vector<int>> res;
    vector<int> tmp;

    bool found_tmp;
    int index_found_tmp;

    // for each row of arr1, if not in arr2, add to res
    for (int i = 0; i<arr1.size(); ++i){
        tmp = arr1[i];
        tie(found_tmp,index_found_tmp)=findPathInList(tmp,arr2);
        if (!found_tmp){
            res.push_back(tmp);
        }
    }
    // for each row of arr2, if not in arr1, add to res
    for (int i = 0; i<arr2.size(); ++i){
        tmp = arr2[i];
        tie(found_tmp,index_found_tmp)=findPathInList(tmp,arr1);
        if (!found_tmp){
            res.push_back(tmp);
        }
    }
    return res;
}



tuple<int, int> getMaxIdx(vector<vector<int>> &pjsummands, vector<vector<int>> &pj){
/* given a list of j-paths,
pick out the j-path which has the 
greatest index in pj. return 
this index as maxIdx. also return the 
index of the summand achieving maxIdx
as summandIdx */
    int maxIdx, summandIdx;
    vector<int> tmp;
    bool found_tmp;
    int index_found_tmp;

    // initialize using 0th summand provided
    summandIdx = 0;
    tmp = pjsummands[summandIdx];
    tie(found_tmp,index_found_tmp)=findPathInList(tmp,pj);
    maxIdx = index_found_tmp;

    // see if other summands have higher index
    for (int i =1; i < pjsummands.size(); ++i){
        tmp = pjsummands[i];
        tie(found_tmp,index_found_tmp)=findPathInList(tmp,pj);
        if (index_found_tmp > maxIdx){
            // if bigger index found, update
            summandIdx = i;
            maxIdx = index_found_tmp;
        }
    }

    return make_tuple(maxIdx,summandIdx);
}

vector<vector<int>> removeUnmarked(vector<vector<int>> &pjsummands, vector<bool> &pjmarked, vector<vector<int>> &pj){
    vector<vector<int>> newSummands;
    vector<int> tmp;
    bool found_tmp;
    int index_found_tmp;

    for (int i = 0; i < pjsummands.size(); ++i){
        tmp = pjsummands[i];
        tie(found_tmp,index_found_tmp)=findPathInList(tmp,pj);
        if (pjmarked[index_found_tmp]==true){
            newSummands.push_back(tmp);
        }
    }

    return newSummands;
}

void addInfAPaths(vector<vector<int>> &p1, vector<double> &p1time, vector<vector<int>> &p2){
    /* go through 2-paths, see which "missing"
    paths need to be added at infinity*/
    vector<int> tmp;
    bool hasRep;
    bool inp1;
    int index_found_tmp;
    int lp = p2[0].size();

    for (int i = 0; i<p2.size(); ++i){
        for (int j = 0; j < lp; ++j){
            for (int k = 0; k< lp; ++k){
                if (k!=j){
                    tmp.push_back(p2[i][k]);
                }         
            }
            //print1DIntVector(tmp);
            hasRep = hasRepeats(tmp);
            tie(inp1,index_found_tmp)=findPathInList(tmp,p1);
                if (hasRep==0 && inp1==0){
                    p1.push_back(tmp);
                    p1time.push_back(maxTime);
                }
                tmp.clear();
            /*
            cout << "has rep" << hasRepeats(tmp) << endl;
            cout << "inp1 " << inp1 << endl;
            cout << "next" << endl;
            */
            tmp.clear();          
        }     
    }
}



vector<vector<int>> simpleBoundary(vector<int> &path){
/*
Compute the boundary summands of a p-path, p>= 1.
*/
    vector<vector<int>> summands;
    vector<int> tmp;
    /* 
    outer loop is the hat operator.
    it says which index to exclude.
    */
    for (int i = 0; i< path.size();++i){
        for (int j = 0; j < path.size(); ++j){
            if (i!=j){
                tmp.push_back(path[j]);
            }
        }
        // summand built, now test for repeats
        if (hasRepeats(tmp)==false){
            summands.push_back(tmp);
        }
        tmp.clear();
    }
    return summands;
}


bool hasRepeats(vector<int> &myvec){
    bool rep = false;
    for (int i =0; i< myvec.size()-1; ++i){
        if (myvec[i] == myvec[i+1]){
            rep = true;
        }
    }
    return rep;
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

vector<vector <int>> getp0(int numNodes){
// get 0 paths
    vector<vector<int>> p0;
    vector<int> tmp;
    for (int i=1; i<= numNodes; i++){
        tmp.push_back(i);
        p0.push_back(tmp);
        tmp.clear();
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

void print1DBoolVector(vector<bool> &myvec){
/* Print all values of bool vector
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
