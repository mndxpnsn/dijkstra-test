/*
 * main.cpp
 *
 *  Created on: 30 Sep 2021
 *      Author: derekharrison
 *
 *      Compare complexity of Dijkstra's algorithm
 *      using a Fibonacci heap and a binary min
 *      heap.
 *
 *      Also compare running time using priority queues
 *      and arrays.
 */

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <chrono>

#include "array.hpp"
#include "bin_heap.hpp"
#include "fib_heap.hpp"
#include "user_types.hpp"
#include "priority_queue.hpp"

int main(int argc, char* argv[]) {

    //Declarations
    int s = 1; //Start vertex. The minimum index for vertices is 1
    int n = 2449; //Number of vertices
    int num_edges = 100*n; //Number of edges

    //Initialize graph for priority queue implementation
    Graph graph(n);

    //Create edges
    srand(time(NULL));
    std::vector< std::vector<int> > edges;
    for(int i = 0; i < num_edges; ++i) {
        int start_vert = rand() % n + 1;
        int end_vert = rand() % n + 1;
        int weight = rand() % n + 1;

        std::vector<int> edge_elem;
        edge_elem.push_back(start_vert);
        edge_elem.push_back(end_vert);
        edge_elem.push_back(weight);
        edges.push_back(edge_elem);

        graph.addEdge(start_vert-1, end_vert-1, weight);
    }

    //Time results based on arrays
    clock_t start_time_fib_heap_ver, end_time_fib_heap_ver;
    double time_ver;
    start_time_fib_heap_ver = clock();
    //Compute distances to nodes from start vertex using arrays
    std::vector<int> results_ver = shortestReach(n, edges, s);
    end_time_fib_heap_ver = clock();
    time_ver = (double) (end_time_fib_heap_ver - start_time_fib_heap_ver) / CLOCKS_PER_SEC * 1000.0;

    //Time results based on priority queue
    clock_t start_time_bin_heap_ver, end_time_bin_heap_ver;
    double time_ver2;
    start_time_bin_heap_ver = clock();
    //Compute distances to nodes from start vertex using priority queues
    vector<int> results_ver2 = graph.shortestPath(s - 1);
    end_time_bin_heap_ver = clock();
    time_ver2 = (double) (end_time_bin_heap_ver - start_time_bin_heap_ver) / CLOCKS_PER_SEC * 1000.0;

    //Reset operation counters
    tot_num_ops = 0;
    num_ops_relax = 0;
    num_ops_decrease_key = 0;
    num_ops_extract_min = 0;
    num_ops_v_overhead = 0;
    num_ops_e_overhead = 0;

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
    int num_ops_relax1 = num_ops_relax;
    int num_ops_decrease_key1 = num_ops_decrease_key;
    int num_ops_extract_min1 = num_ops_extract_min;
    int num_ops_v_overhead1 = num_ops_v_overhead;
    int num_ops_e_overhead1 = num_ops_e_overhead;

    //Reset operation counters
    tot_num_ops = 0;
    num_ops_relax = 0;
    num_ops_decrease_key = 0;
    num_ops_extract_min = 0;
    num_ops_v_overhead = 0;
    num_ops_e_overhead = 0;

    //Time results based on binary heap
    clock_t start_time_bin_heap, end_time_bin_heap;
    double time2;
    start_time_bin_heap = clock();
    //Compute distances to nodes from start vertex using a binary heap
    std::vector<int> results2 = shortest_reach2(n, edges, s);
    end_time_bin_heap = clock();
    time2 = (double) (end_time_bin_heap - start_time_bin_heap) / CLOCKS_PER_SEC * 1000.0;

    //Store results computation based on binary heap
    int tot_num_ops2 = tot_num_ops;
    int num_ops_relax2 = num_ops_relax;
    int num_ops_decrease_key2 = num_ops_decrease_key;
    int num_ops_extract_min2 = num_ops_extract_min;
    int num_ops_v_overhead2 = num_ops_v_overhead;
    int num_ops_e_overhead2 = num_ops_e_overhead;

    //Check if results of the various methods are the same
    bool results_match = true;
    int size_res1 = (int) results.size();
    int size_res2 = (int) results2.size();
    int size_ver = (int) results_ver.size();
    int size_ver2 = (int) results_ver2.size();
    if(size_res1 != size_res2) { results_match = false; }
    if(size_res2 != size_ver) { results_match = false; }
    if(size_res1 != size_ver2) { results_match = false; }
    for(int i = 0; i < size_res1; ++i) {
        if(results[i] != results2[i]) { results_match = false; }
        if(results2[i] != results_ver[i]) { results_match = false; }
        if(results[i] != results_ver2[i]) { results_match = false; }
    }

    //Print results
    float tot_num_ops_est = 5*n + 4*num_edges + 6.4*n*log(n)/log(2);
    float complexity_ratio = tot_num_ops1 / tot_num_ops_est;
    float tot_num_ops_est_bin_min = 9*n + 2*num_edges + 0.9*n*log(n)/log(2) + 0.18*num_edges*log(n)/log(2);
    float complexity_ratio2 = tot_num_ops2 / tot_num_ops_est_bin_min;

    std::cout << "Execution results:" << std::endl;
    std::cout << "results obtained from the various methods match: " << results_match << std::endl;
    std::cout << "timing execution arrays: " << time_ver << std::endl;
    std::cout << "timing execution priority queue: " << time_ver2 << std::endl;
    std::cout << "timing execution fibonacci heap: " << time << std::endl;
    std::cout << "timing execution binary heap: " << time2 << std::endl;
    std::cout << std::endl;

    std::cout << "Measurements number of operations:" << std::endl;
    std::cout << "total number of operations measured fibonacci heap: " << tot_num_ops1 << std::endl;
    std::cout << "total number of operations measured binary heap: " << tot_num_ops2 << std::endl;
    std::cout << "number of operations relax fibonacci heap: " << num_ops_relax1 << std::endl;
    std::cout << "number of operations relax binary heap: " << num_ops_relax2 << std::endl;
    std::cout << "number of operations decrease key fibonacci heap: " << num_ops_decrease_key1 << std::endl;
    std::cout << "number of operations decrease key binary heap: " << num_ops_decrease_key2 << std::endl;
    std::cout << "number of operations extract min fibonacci heap: " << num_ops_extract_min1 << std::endl;
    std::cout << "number of operations extract min binary heap: " << num_ops_extract_min2 << std::endl;
    std::cout << "number of overhead operations fibonacci heap: " << num_ops_v_overhead1 + num_ops_e_overhead1 << std::endl;
    std::cout << "number of overhead operations binary heap: " << num_ops_v_overhead2 + num_ops_e_overhead2 << std::endl;
    std::cout << "number of operations extract min / VlgV, fibonacci heap: " << (float) num_ops_extract_min1/(n * log(n)/log(2)) << std::endl;
    std::cout << "number of operations extract min / VlgV, binary heap: " << (float) num_ops_extract_min2/(n * log(n)/log(2)) << std::endl;
    std::cout << std::endl;

    std::cout << "Complexity and total number of operations ratios:" << std::endl;
    std::cout << "complexity ratio fibonacci heap: " << complexity_ratio << std::endl;
    std::cout << "complexity ratio binary min heap: " << complexity_ratio2 << std::endl;
    std::cout << "fibonacci heap tot num ops ratio: " << ((float) tot_num_ops1) / (num_ops_relax1 + num_ops_decrease_key1 + num_ops_extract_min1 + num_ops_v_overhead1 + num_ops_e_overhead1) << std::endl;
    std::cout << "binary heap tot num ops ratio: " << ((float) tot_num_ops2) / (num_ops_relax2 + num_ops_decrease_key2 + num_ops_extract_min2 + num_ops_v_overhead2 + num_ops_e_overhead2) << std::endl;
    std::cout << "total number of ops fib heap / total number of ops binary heap: " << (float) tot_num_ops1 / tot_num_ops2 << std::endl;
    std::cout << std::endl;

//    //Print results based on Fibonacci heap
//    std::cout << "Shortest distances from start vertex:" << std::endl;
//    int size_results = (int) results.size();
//    for(int i = 0; i < size_results; ++i) {
//        std::cout << results[i] << " ";
//    }
//    std::cout << std::endl;
//
//    //Print results based on binary heap
//    int size_results2 = (int) results2.size();
//    for(int i = 0; i < size_results2; ++i) {
//        std::cout << results2[i] << " ";
//    }
//    std::cout << std::endl;
//
//    //Print results based on arrays
//    int size_results_ver = (int) results_ver.size();
//    for(int i = 0; i < size_results_ver; ++i) {
//        std::cout << results_ver[i] << " ";
//    }
//    std::cout << std::endl;
//
//    //Print results based on priority queues
//    int size_results_ver2 = (int) results_ver2.size();
//    for(int i = 0; i < size_results_ver2; i++){
//        std::cout << results_ver2[i] << " ";
//    }
//    std::cout << std::endl;

    std::cout << "done" << std::endl;

    return 0;
}
