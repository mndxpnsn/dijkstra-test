/*
 * main.cpp
 *
 *  Created on: 30 Sep 2021
 *      Author: mndx
 *
 *      Compare number of operations and running time
 *      of Dijkstra's algorithm using a Fibonacci heap
 *      and a binary min heap.
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

    //Input parameters
    int s = 31; //Start vertex. The minimum index for vertices is 1
    int n = 2499; //Number of vertices
    int num_edges = 100*n; //Number of edges

    //Initialize graph for priority queue implementation
    Graph graph(n);

    //Create edges
    srand((unsigned) time(NULL));
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

        //Add edge to priority queue
        graph.addEdge(start_vert-1, end_vert-1, weight);
    }

    //Time results based on arrays
    clock_t start_time, end_time;
    double time_arrays;
    start_time = clock();

    //Compute distances to nodes from start vertex using arrays
    std::vector<int> results_arrays = shortestReach(n, edges, s);

    end_time = clock();
    time_arrays = (double) (end_time - start_time) / CLOCKS_PER_SEC * 1000.0;

    //Time results based on priority queue
    double time_prio_queue;
    start_time = clock();

    //Compute distances to nodes from start vertex using priority queues
    std::vector<int> results_prio_queue = graph.shortestPath(s - 1);

    end_time = clock();
    time_prio_queue = (double) (end_time - start_time) / CLOCKS_PER_SEC * 1000.0;

    //Reset operation counters
    tot_num_ops = 0;
    num_ops_relax = 0;
    num_ops_decrease_key = 0;
    num_ops_extract_min = 0;
    num_ops_v_overhead = 0;
    num_ops_e_overhead = 0;

    //Time results based on Fibonacci heap
    double time_fib;
    start_time = clock();

    //Compute distances to nodes from start vertex using a Fibonacci heap
    std::vector<int> results_fib = shortest_reach(n, edges, s);

    end_time = clock();
    time_fib = (double) (end_time - start_time) / CLOCKS_PER_SEC * 1000.0;

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
    double time_bin;
    start_time = clock();

    //Compute distances to nodes from start vertex using a binary heap
    std::vector<int> results_bin = shortest_reach2(n, edges, s);

    end_time = clock();
    time_bin = (double) (end_time - start_time) / CLOCKS_PER_SEC * 1000.0;

    //Store results computation based on binary heap
    int tot_num_ops2 = tot_num_ops;
    int num_ops_relax2 = num_ops_relax;
    int num_ops_decrease_key2 = num_ops_decrease_key;
    int num_ops_extract_min2 = num_ops_extract_min;
    int num_ops_v_overhead2 = num_ops_v_overhead;
    int num_ops_e_overhead2 = num_ops_e_overhead;

    //Check if results of the various methods are the same
    bool results_match = true;
    int size_fib = (int) results_fib.size();
    int size_bin = (int) results_bin.size();
    int size_arrays = (int) results_arrays.size();
    int size_prio = (int) results_prio_queue.size();
    if(size_fib != size_bin) { results_match = false; }
    if(size_bin != size_arrays) { results_match = false; }
    if(size_fib != size_prio) { results_match = false; }
    for(int i = 0; i < size_fib; ++i) {
        if(results_fib[i] != results_bin[i]) { results_match = false; }
        if(results_bin[i] != results_arrays[i]) { results_match = false; }
        if(results_fib[i] != results_prio_queue[i]) { results_match = false; }
    }

    //Print results
    std::cout << "Execution results:" << std::endl;
    std::cout << "results obtained from the various methods match: " << results_match << std::endl;
    std::cout << "timing execution arrays: " << time_arrays << std::endl;
    std::cout << "timing execution priority queue: " << time_prio_queue << std::endl;
    std::cout << "timing execution fibonacci heap: " << time_fib << std::endl;
    std::cout << "timing execution binary heap: " << time_bin << std::endl;
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

    std::cout << "Total number of operations ratios:" << std::endl;
    std::cout << "fibonacci heap tot num ops ratio: " << ((float) tot_num_ops1) / (num_ops_relax1 + num_ops_decrease_key1 + num_ops_extract_min1 + num_ops_v_overhead1 + num_ops_e_overhead1) << std::endl;
    std::cout << "binary heap tot num ops ratio: " << ((float) tot_num_ops2) / (num_ops_relax2 + num_ops_decrease_key2 + num_ops_extract_min2 + num_ops_v_overhead2 + num_ops_e_overhead2) << std::endl;
    std::cout << "total number of ops fib heap / total number of ops binary heap: " << (float) tot_num_ops1 / tot_num_ops2 << std::endl;
    std::cout << std::endl;

    //Print results based on Fibonacci heap
    std::cout << "Shortest distances from start vertex:" << std::endl;
    for(int i = 0; i < size_fib; ++i) {
        std::cout << results_fib[i] << " ";
    }
    std::cout << std::endl;

    //Print results based on binary heap
    for(int i = 0; i < size_bin; ++i) {
        std::cout << results_bin[i] << " ";
    }
    std::cout << std::endl;

    //Print results based on arrays
    for(int i = 0; i < size_arrays; ++i) {
        std::cout << results_arrays[i] << " ";
    }
    std::cout << std::endl;

    //Print results based on priority queues
    for(int i = 0; i < size_prio; i++){
        std::cout << results_prio_queue[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "done" << std::endl;

    return 0;
}
