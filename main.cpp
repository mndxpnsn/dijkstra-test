/*
 * main.cpp
 *
 *  Created on: 30 Sep 2021
 *      Author: derekharrison
 */


#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include <stdlib.h>
#include <time.h>
#include <vector>

#include "bin_heap.hpp"
#include "fib_heap.hpp"
#include "user_types.hpp"

int main(int argc, char* argv[]) {

    //Declarations
    int s = 2; //Start vertex must be greater or equal to 1
    int n = 2499; //Number of vertices
    int num_edges = 3125; //Number of edges

    //Create edges
//    srand(time(NULL));
    std::vector< std::vector<int> > edges;
    for(int i = 0; i < num_edges; ++i) {
        int start_vert = rand() % n + 1;
        int end_vert = rand() % n + 1;
        int weight = rand() % 200 + 1;

        std::vector<int> edge_elem;
        edge_elem.push_back(start_vert);
        edge_elem.push_back(end_vert);
        edge_elem.push_back(weight);
        edges.push_back(edge_elem);
    }

    clock_t tv1, tv2;
    double time;
    tv1 = clock();

    //Compute distances to nodes from start vertex
    std::vector<int> results = shortest_reach2(n, edges, s);

    tv2 = clock();
    time = (tv2 - tv1)/(CLOCKS_PER_SEC / (double) 1000.0);

    //Print results
    float tot_num_ops_est = 7*n + 3*num_edges + 6.4*n*log(n)/log(2);
    float tot_num_ops_est_bin_min = 10*n + num_edges + 0.9*n*log(n)/log(2) + 0.18*num_edges*log(n)/log(2);
    float complexity_ratio = tot_num_ops / tot_num_ops_est;
    int size_results = (int) results.size();
    for(int i = 0; i < size_results; ++i) {
        std::cout << results[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "timing execution: " << time << std::endl;
    std::cout << "number of operations estimated fib heap 7V + 3E + 6.4VlgV: " << tot_num_ops_est << std::endl;
    std::cout << "number of operations measured: " << tot_num_ops << std::endl;
    std::cout << "number of operations estimated bin heap 10V + E + 0.9VlgV + 0.18ElgV: " << tot_num_ops_est_bin_min << std::endl;
    std::cout << "complexity ratio fib heap: " << complexity_ratio << std::endl;
    std::cout << "complexity ratio bin min heap: " << (tot_num_ops/tot_num_ops_est_bin_min) << std::endl;
    std::cout << "done" << std::endl;

    return 0;
}
