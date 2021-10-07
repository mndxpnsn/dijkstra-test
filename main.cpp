/*
 * main.cpp
 *
 *  Created on: 30 Sep 2021
 *      Author: derekharrison
 *
 *      Compare complexity of Dijkstra's algorithm
 *      using a Fibonacci heap and a binary min
 *      heap.
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
    std::vector<int> results = shortest_reach(n, edges, s);

    tv2 = clock();
    time = (tv2 - tv1)/(CLOCKS_PER_SEC / (double) 1000.0);

    int tot_num_ops1 = tot_num_ops;

    clock_t tv12, tv22;
    double time2;
    tv12 = clock();

    //Reset total number of operations counter
    tot_num_ops = 0;

    //Compute distances to nodes from start vertex
    std::vector<int> results2 = shortest_reach2(n, edges, s);

    tv22 = clock();
    time2 = (tv22 - tv12)/(CLOCKS_PER_SEC / (double) 1000.0);

    //Print results
    float tot_num_ops_est = 5*n + 4*num_edges + 6.4*n*log(n)/log(2);
    float complexity_ratio = tot_num_ops1 / tot_num_ops_est;
    float tot_num_ops_est_bin_min = 9*n + 2*num_edges + 0.9*n*log(n)/log(2) + 0.18*num_edges*log(n)/log(2);
    float complexity_ratio2 = tot_num_ops / tot_num_ops_est_bin_min;
    std::cout << "timing execution fib heap: " << time << std::endl;
    std::cout << "timing execution bin heap: " << time2 << std::endl;
    std::cout << "number of operations estimated fib heap 5V + 4E + 6.4VlgV: " << tot_num_ops_est << std::endl;
    std::cout << "number of operations estimated bin heap 9V + 2E + 0.9VlgV + 0.18ElgV: " << tot_num_ops_est_bin_min << std::endl;
    std::cout << "number of operations measured fib heap: " << tot_num_ops1 << std::endl;
    std::cout << "number of operations measured bin heap: " << tot_num_ops << std::endl;
    std::cout << "complexity ratio fib heap: " << complexity_ratio << std::endl;
    std::cout << "complexity ratio bin min heap: " << complexity_ratio2 << std::endl;

    int size_results = (int) results.size();
    for(int i = 0; i < size_results; ++i) {
    	std::cout << results[i] << " ";
    }
    std::cout << std::endl;

    int size_results2 = (int) results2.size();
    for(int i = 0; i < size_results2; ++i) {
    	std::cout << results2[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "done" << std::endl;

    return 0;
}
