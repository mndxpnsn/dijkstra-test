/*
 * bin_heap.cpp
 *
 *  Created on: 7 Oct 2021
 *      Author: derekharrison
 */

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include <stdlib.h>
#include <time.h>
#include <vector>

#include "fib_heap.hpp"
#include "user_types.hpp"

Heap::Heap(int size) {
    heap_size = size;
    A = new node*[size+1];
    element_map = new int[size+1];
    size_array = size + 1;
    element_map[0] = 0;
    for(int i = 1; i <= heap_size; ++i) {
    	num_ops_v_overhead++;
        tot_num_ops++;
        element_map[i] = i;
        A[i] = new node;
        A[i]->pi = NULL;

    }
}

Heap::~Heap() {
//    for(int i = 0; i < size_array; ++i) {
//        delete A[i];
//    }
//
//    delete [] A;
//    delete [] element_map;
}

int Heap::parent(int i) {
    return i/2;
}

int Heap::left(int i) {
    return 2*i;
}

int Heap::right(int i) {
    return 2*i + 1;
}

node* Heap::get_heap_element(int node_index) {
    int index_in_heap = element_map[node_index];
    return A[index_in_heap];
}

int Heap::get_heap_index(int node_index) {
    int index_in_heap = element_map[node_index];
    return index_in_heap;
}

int Heap::get_root_index() {
    return A[1]->index;
}

void Heap::min_heapify(node* A[], int i) {
    num_ops_extract_min++;
    tot_num_ops++;
    int l, r, smallest;
    l = Heap::left(i);
    r = Heap::right(i);
    if(l < heap_size + 1 && A[l]->key < A[i]->key) {
        smallest = l;
    }
    else {
        smallest = i;
    }
    if(r < heap_size + 1 && A[r]->key < A[smallest]->key) {
        smallest = r;
    }
    if(smallest != i) {
        node* dummy;
        dummy = A[i];

        element_map[A[smallest]->index] = i;
        element_map[A[i]->index] = smallest;

        A[i] = A[smallest];
        A[smallest] = dummy;

        Heap::min_heapify(A, smallest);
    }
}

void Heap::build_min_heap() {
    for(int i = heap_size/2; i > 0; --i) {
        Heap::min_heapify(A, i);
    }
}

void Heap::set_heap(node* B[]) {
    for(int i = 1; i < heap_size + 1; ++i) {
    	num_ops_v_overhead++;
        tot_num_ops++;
        A[i] = B[i];
    }
}

void Heap::get_heap(node B[]) {
    for(int i = 1; i < heap_size + 1; ++i) {
        B[i] = *A[i];
    }
}

int Heap::get_heap_size() {
    return heap_size;
}

bool Heap::min_heap_verify() {
    bool is_min_heap = true;
    for(int i = (heap_size - 1)/2; i > 0; --i) {
        int l, r;
        l = Heap::left(i);
        r = Heap::right(i);
        if(A[i]->key > A[l]->key || A[i]->key > A[r]->key) {
            is_min_heap = false;
        }
    }

    return is_min_heap;
}

void Heap::print_heap() {
    for(int i = heap_size/2; i > 0; --i) {
        int l, r;
        l = Heap::left(i);
        r = Heap::right(i);
        if(l < heap_size + 1 && r < heap_size + 1) {
            printf("node: %i, key: %i, key left child: %i, key right child: %i\n", i, A[i]->key,  A[l]->key,  A[r]->key);
        }
    }
}

void Heap::print_element_map() {
    for(int i = 1; i <= heap_size; ++i) {
        int index_loc = element_map[i];
        std::cout << "index: " << i << ", A[index]->index: " << A[index_loc]->index
                  << ", key: " << A[index_loc]->key << ", i: " << i << std::endl;
    }
}

node* Heap::heap_extract_min() {

    if(heap_size < 1) {
        std::cout << "heap size is less than 1" << std::endl;
    }
    node* min = A[1];

    element_map[A[heap_size]->index] = 1;
    A[1] = A[heap_size];
    heap_size = heap_size - 1;

    Heap::min_heapify(A, 1);

    return min;
}

void Heap::heap_decrease_key(int index, double key) {
    if(key > A[index]->key) {
        printf("new key is larger than current key\n");
    }
    else {
        A[index]->key = key;
        while(index > 1 && A[parent(index)]->key > A[index]->key) {
            tot_num_ops++;
            num_ops_decrease_key++;
            element_map[A[index]->index] = parent(index);
            element_map[A[parent(index)]->index] = index;

            node* dummy = A[index];
            A[index] = A[parent(index)];
            A[parent(index)] = dummy;

            index = parent(index);
        }
    }
}

void relax(node* u, node* v, int** w, Heap* heap) {
    int index_in_heap = heap->get_heap_index(v->index);
    if(v->key > u->key + w[u->index][v->index]) {
        int weight = u->key + w[u->index][v->index];
        heap->heap_decrease_key(index_in_heap, weight);
        v->pi = u;
    }
}

int map_index2(int n, int index, int s) {
    int r;

    if(index >= s) { r = index - s + 1; }
    else { r = n - s + index + 1; }

    return r;
}

int map_inverse(int n, int index, int s) {
    int r;

    r = s + index - 1;
    if(r > n) {
        r = r - n;
    }

    return r;
}

void populate_adj_and_weight_hr(int** adj_mat,
                                int** weight_mat,
                                int size_graph,
                                std::vector< std::vector<int> > &edges,
                                node** heap,
                                int s) {

    for(int i = 1; i < size_graph + 1; ++i) {
    	num_ops_v_overhead++;
        tot_num_ops++;
        heap[i] = new node;
        heap[i]->key = INF;
        heap[i]->pi = NULL;
        heap[i]->index = i;
        heap[i]->index_og = map_inverse(size_graph, i, s);
    }

    heap[1]->key = 0;

    int num_edges = edges.size();
    for(int i = 0; i < num_edges; ++i) {
    	num_ops_e_overhead++;
        tot_num_ops++;
        int start_index = edges[i][0];
        int end_index = edges[i][1];
        int weight = edges[i][2];

        int start = map_index2(size_graph, start_index, s);
        int end = map_index2(size_graph, end_index, s);

        heap[start]->adj_nodes.push_back(end);
        heap[end]->adj_nodes.push_back(start);

        weight_mat[start][end] = weight;
        weight_mat[end][start] = weight;

        adj_mat[start][end] = SETVAR;
        adj_mat[end][start] = SETVAR;
    }

    for(int i = 0; i < num_edges; ++i) {
    	num_ops_e_overhead++;
        tot_num_ops++;
        int start_index = edges[i][0];
        int end_index = edges[i][1];
        int weight = edges[i][2];

        int start = map_index2(size_graph, start_index, s);
        int end = map_index2(size_graph, end_index, s);

        if(weight_mat[start][end] >= weight) {
            weight_mat[start][end] = weight;
            weight_mat[end][start] = weight;
        }
    }
}

std::vector<int> shortest_reach2(int n, std::vector< std::vector<int> > &edges, int s) {

    //Time results
    clock_t start_time_fib_heap, end_time_fib_heap;
    double time;
    start_time_fib_heap = clock();

    std::vector<node*> rs_S;

    //Set index maps
    int* index_map_end = new int[n+1];

    //Initialize weight and adjacency matrices and binary min heap
    node** heap = new node*[n + 1];
    int** adj_mat = int2D(n + 1);
    int** weight_mat = int2D(n + 1);

    //Populate weight and adjacency matrices and initialize heap
    populate_adj_and_weight_hr(adj_mat, weight_mat, n, edges, heap, s);

    //End time measurements
    end_time_fib_heap = clock();
    time = (double) (end_time_fib_heap - start_time_fib_heap) / CLOCKS_PER_SEC * 1000.0;
    std::cout << "execution time binary heap: " << time << std::endl;

    //Set heap and build heap
    Heap min_heap(n);
    min_heap.set_heap(heap);
    min_heap.build_min_heap();

    //Perform Dijkstra's algorithm
    int heap_size = min_heap.get_heap_size();
    int rs_elem_counter = 0;
    while(heap_size > 0) {

        node* u = min_heap.heap_extract_min();
        heap_size = min_heap.get_heap_size();

        int u_index = u->index;
        int num_adj_nodes = u->adj_nodes.size();

        for(int i = 0; i < num_adj_nodes; ++i) {
            tot_num_ops++;
            num_ops_relax++;
            int it = u->adj_nodes[i];
            node* v = min_heap.get_heap_element(it);
            int v_index = v->index;

            //Extracted nodes always point to node 1 in the heap,
            //and the node at 1 may not be an adjacent node
            //Therefore adjacency must be verified with adj_mat
            if(adj_mat[u_index][v_index] == SETVAR) {
                relax(u, v, weight_mat, &min_heap);
            }
        }

        rs_S.push_back(u);
        index_map_end[u->index_og] = rs_elem_counter;
        rs_elem_counter++;
    }

    //Reorder results
    int size_results = rs_S.size();
    std::vector<int> rs_S_reordered;

    for(int i = 1; i <= size_results; ++i) {
    	num_ops_v_overhead++;
        tot_num_ops++;
        int j = index_map_end[i];
        if(rs_S[j]->index_og != s) {
            int key = rs_S[j]->key;
            if(key == INF) { key = -1; }
            rs_S_reordered.push_back(key);
        }
    }

    //Deallocate memory
    free_int2D(adj_mat, n + 1);
    free_int2D(weight_mat, n + 1);
    delete [] heap;
    delete [] index_map_end;

    return rs_S_reordered;
}
