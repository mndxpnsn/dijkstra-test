/*
 * fib_heap_comp.cpp
 *
 *  Created on: 9 Oct 2021
 *      Author: derekharrison
 */

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <chrono>
#include <map>
#include <set>

using namespace std;

vector<int> shortestReach(int n, vector< vector<int> > edges, int s) {
    vector<int> minimumPath(n, INT_MAX);
    set<int> pikedNodes;

    /* it could be done with 2D array two, zero no edge, non-zero weigh of excited edge */
    vector<map< int , int> > minWeightEdgesFrom(n);
    map<int, int> * srcEdges;
    map<int, int> * dstEdges;
    int edgeWeight, edgeSrcNodeIndex, edgeDstNodeIndex;

    for(auto edge: edges)
    {
        if(edge[0] == 9 && edge[1] == 1)
            cout << "";
        edgeSrcNodeIndex = edge[0] - 1;
        edgeDstNodeIndex = edge[1] - 1;
        edgeWeight = edge[2];

        srcEdges = &minWeightEdgesFrom[edgeSrcNodeIndex];
        dstEdges = &minWeightEdgesFrom[edgeDstNodeIndex];

        auto foundEdge = srcEdges->find(edgeDstNodeIndex);
        if(foundEdge == srcEdges->end() || foundEdge->second > edgeWeight) {
            (*srcEdges)[edgeDstNodeIndex] = edgeWeight;
            (*dstEdges)[edgeSrcNodeIndex] = edgeWeight;
        }
    }

    minimumPath[s - 1] = 0;

    int minDistNodeIndex, minValue;
    while (pikedNodes.size() != n) {
        /* it could be a heap to improve time complexity from (n^2) to (nlog(n)) */
        minDistNodeIndex = -1;
        minValue = INT_MAX;
        for(int i = 0; i < n; i ++)
        {
            if(pikedNodes.find(i) == pikedNodes.end() && minValue > minimumPath[i])
            {
                minValue = minimumPath[i];
                minDistNodeIndex = i;
            }
        }
        /* no connected node anymore */
        if(minValue == INT_MAX)
            break;

        pikedNodes.insert(minDistNodeIndex);
        /* update picked node adjacent weights */
        for(auto edgeFrom: minWeightEdgesFrom[minDistNodeIndex]) {
            if(minimumPath[minDistNodeIndex] + edgeFrom.second < minimumPath[edgeFrom.first])
                minimumPath[edgeFrom.first] = minimumPath[minDistNodeIndex] + edgeFrom.second;
        }
    }
    for(int i = s - 1; i < n - 1; i++)
        minimumPath[i] = minimumPath[i + 1];
    minimumPath.resize(n - 1);

    for(int i = 0; i < n - 1; i++)
        if(minimumPath[i] == INT_MAX)
            minimumPath[i] = -1;

    return minimumPath;
}

