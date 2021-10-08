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
#include <chrono>

#include "bin_heap.hpp"
#include "fib_heap.hpp"
#include "user_types.hpp"

int main(int argc, char* argv[]) {

    //Declarations
    int s = 2; //Start vertex must be greater or equal to 1
    int n = 2499; //Number of vertices
    int num_edges = 3125; //Number of edges

    //Create edges
    srand(time(NULL));
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

    //Time results based on Fibonacci heap
    clock_t start_time_fib_heap, end_time_fib_heap;
    double time;
    start_time_fib_heap = clock();

    //Compute distances to nodes from start vertex using a fibonacci heap
    std::vector<int> results = shortest_reach(n, edges, s);

    end_time_fib_heap = clock();
    time = (double) (end_time_fib_heap - start_time_fib_heap) / CLOCKS_PER_SEC * 1000.0;

    //Store results computation based on Fibonacci heap
    int tot_num_ops1 = tot_num_ops;

    //Reset total number of operations counter
    tot_num_ops = 0;

    //Time results based on binary heap
    clock_t start_time_bin_heap, end_time_bin_heap;
    double time2;
    start_time_bin_heap = clock();

    //Compute distances to nodes from start vertex using a binary heap
    std::vector<int> results2 = shortest_reach2(n, edges, s);

    end_time_bin_heap = clock();
    time2 = (double) (end_time_bin_heap - start_time_bin_heap) / CLOCKS_PER_SEC * 1000.0;

    //Check if results of both methods are the same
    bool results_match = true;
    int size_res1 = (int) results.size();
    int size_res2 = (int) results2.size();
    if(size_res1 != size_res2) { results_match = false; }
    for(int i = 0; i < size_res1; ++i) {
        if(results[i] != results2[i]) {
            results_match = false;
        }
    }

    //Print results
    float tot_num_ops_est = 5*n + 4*num_edges + 6.4*n*log(n)/log(2);
    float complexity_ratio = tot_num_ops1 / tot_num_ops_est;
    float tot_num_ops_est_bin_min = 6*n + 2*num_edges + 0.9*n*log(n)/log(2) + 0.18*num_edges*log(n)/log(2);
    float complexity_ratio2 = tot_num_ops / tot_num_ops_est_bin_min;
    std::cout << "results obtained from binary heap and fibonacci heap match: " << results_match << std::endl;
    std::cout << "timing execution fibonacci heap: " << time << std::endl;
    std::cout << "timing execution binary heap: " << time2 << std::endl;
    std::cout << "number of operations estimated fibonacci heap 5V + 4E + 6.4VlgV: " << tot_num_ops_est << std::endl;
    std::cout << "number of operations estimated binary heap 6V + 2E + 0.9VlgV + 0.18ElgV: " << tot_num_ops_est_bin_min << std::endl;
    std::cout << "number of operations measured fibonacci heap: " << tot_num_ops1 << std::endl;
    std::cout << "number of operations measured binary heap: " << tot_num_ops << std::endl;
    std::cout << "complexity ratio fibonacci heap: " << complexity_ratio << std::endl;
    std::cout << "complexity ratio binary min heap: " << complexity_ratio2 << std::endl;
    std::cout << std::endl;

    //Print results based on Fibonacci heap
    std::cout << "Shortest distances from start vertex:" << std::endl;
    int size_results = (int) results.size();
    for(int i = 0; i < size_results; ++i) {
        std::cout << results[i] << " ";
    }
    std::cout << std::endl;

    //Print results based on binary heap
    int size_results2 = (int) results2.size();
    for(int i = 0; i < size_results2; ++i) {
        std::cout << results2[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "done" << std::endl;

    return 0;
}
