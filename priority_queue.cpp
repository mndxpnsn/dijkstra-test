/*
 * Implementation of Dijkstra's algorithm
 * using priority queues. The code below
 * was gotten from github.
 */

#include <iostream>
#include <vector>
#include <queue>
#include <list>

#include "priority_queue.hpp"

#define INF 0x3f3f3f3f
using namespace std;

Graph::Graph(int V){
    this->V = V;
    adj = new list<myPair>[this->V];
}

void Graph::addEdge(int u,int v,int w){
    adj[u].push_back({v,w});
    adj[v].push_back({u,w});
}

vector<int> Graph::shortestPath(int src){
    //Dijkstra's algorithm
    priority_queue<myPair,vector<myPair>,greater<myPair> > pq;
    vector<int> dist(this->V,INF);
    dist[src] = 0;
    list<myPair>::iterator it;

    pq.push({0,src});
    while(!pq.empty()){
        int u = pq.top().second;
        pq.pop();

        for(it = adj[u].begin();it!=adj[u].end();++it){
            int v = it->first;
            int w = it->second;
            if(dist[v] > dist[u] + w){
                dist[v] = dist[u] + w;
                pq.push({dist[v],v});
            }
        }
    }

    //Set results
    vector<int> result;
    for(int i=0;i<this->V;i++){
        if(i != src) {
            int distance = dist[i];
            if(distance == INF) {
                result.push_back(-1);
            }
            else {
                result.push_back(distance);
            }
        }
    }

    return result;
}
