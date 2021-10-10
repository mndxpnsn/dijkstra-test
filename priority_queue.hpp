/*
 * bin_heap_comp.hpp
 *
 *  Created on: 9 Oct 2021
 *      Author: mndx
 */

#ifndef PRIORITY_QUEUE_HPP_
#define PRIORITY_QUEUE_HPP_

#include <iostream>
#include <vector>
#include <queue>
#include <list>

using namespace std;

typedef pair<int,int> myPair;
class Graph{
    int V;
    list<myPair> *adj;
public:
    Graph(int V);
    void addEdge(int u,int v,int w);
    vector<int> shortestPath(int src);
};

#endif /* PRIORITY_QUEUE_HPP_ */
