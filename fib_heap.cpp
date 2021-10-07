/*
 * fib_heap.cpp
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

#include "user_types.hpp"

int tot_num_ops = 0;

bool** bool2D(const int size) {
    bool** p = new bool*[size];

    for(int i = 0; i < size; ++i) {
        tot_num_ops++;
        p[i] = new bool[size];
    }

    return p;
}

int** int2D(const int size) {
    int** p = new int*[size];

    for(int i = 0; i < size; ++i) {
        tot_num_ops++;
        p[i] = new int[size];
    }

    return p;
}

void free_bool2D(bool** p, int size) {
    for(int i = 0; i < size; ++i) {
        tot_num_ops++;
        delete [] p[i];
    }

    delete [] p;
}

void free_int2D(int** p, int size) {
    for(int i = 0; i < size; ++i) {
        tot_num_ops++;
        delete [] p[i];
    }

    delete [] p;
}

void free_node_ref(node** v_ref, int size) {
    for(int i = 0; i < size; ++i) {
        tot_num_ops++;
        delete v_ref[i];
    }

    delete [] v_ref;
}

void fib_heap_insert(FibHeap* H, node* x) {
    x->degree = 0;
    x->p = NULL;
    x->child = NULL;
    x->mark = false;

    if(H->min == NULL) {
        x->left = x;
        x->right = x;
        H->min = x;
        H->n = 0;
    }
    else {
        x->left = H->min;
        x->right = H->min->right;
        H->min->right->left = x;
        H->min->right = x;
        if(x->key < H->min->key) {
            H->min = x;
        }
    }

    H->n = H->n + 1;
}

void print_root_list(node* z) {
    node* xt = z;
    if(xt != NULL) {
        do {
            std::cout << "xt->key: " << xt->key;
            std::cout << ", xt->degree: " << xt->degree << std::endl;
            xt = xt->right;
        } while(xt != z);
    }
}

void make_child_of(FibHeap* H, node* y, node* x) {

    //Remove node from root list
    y->left->right = y->right;
    y->right->left = y->left;

    if(x->child == NULL) {
        x->child = y;
        y->p = x;
        y->left = y;
        y->right = y;
    }
    else {
        y->left = x->child;
        y->right = x->child->right;
        y->p = x;
        x->child->right->left = y;
        x->child->right = y;
    }

    //Set mark
    y->mark = false;

    x->degree = x->degree + 1;
}

void link(FibHeap* H, node** A, node* x, node* y) {

    int d = x->degree;

    //Make y child of x;
    make_child_of(H, y, x);

    A[d] = NULL;
    A[d+1] = x;

    if(y == H->min) {
        H->min = x;
    }
}

void link_dup_deg(FibHeap* H, node** A, node*& x, bool& there_is_dup) {
    int d = x->degree;
    //There is a node with the same degree and node is not node x
    if(A[d] != NULL && A[d] != x) {
        there_is_dup = true;
        node* y = A[d];
        //Link x and y
        if(y->key > x->key) {
            //Make y child of x
            link(H, A, x, y);
        }
        else {
            //Make x child of y
            link(H, A, y, x);
            x = y;
        }
    }
    //There is no node with the same degree or node is node x
    else {
        A[d] = x;
    }
}

void consolidate(FibHeap* H) {

    //Compute upper bound root list
    double golden = (1.0 + sqrt(5.0)) / 2.0;
    double f = log(H->n) / log(golden);
    int D = floor(f + 0.01) + 1;

    //Allocate memory for root list construction
    node** A = new node*[D + 1];
    for(int i = 0; i < D + 1; ++i) {
        tot_num_ops++;
        A[i] = NULL;
    }

    //Ensure all root nodes have unique degrees
    node* x = H->min;
    if(x != NULL) {
        bool there_is_dup = true;
        while(there_is_dup) {
            there_is_dup = false;
            x = H->min;
            do {
                tot_num_ops++;
                link_dup_deg(H, A, x, there_is_dup);
                x = x->right;
            } while(x != H->min);
        }
    }

    //Reconstruct root list
    H->min = NULL;
    for(int i = 0; i < D + 1; ++i) {
        tot_num_ops++;
        if(A[i] != NULL) {
            if(H->min == NULL) {
                A[i]->left = A[i];
                A[i]->right = A[i];
                A[i]->p = NULL;
                H->min = A[i];
            }
            else {
                A[i]->left = H->min;
                A[i]->right = H->min->right;
                H->min->right->left = A[i];
                H->min->right = A[i];
                A[i]->p = NULL;
                if(A[i]->key < H->min->key) {
                    H->min = A[i];
                }
            }
        }
    }

    //Free root list reference
    delete [] A;
}

void print_child_list(node* child) {
    node* xt = child;
    if(xt != NULL) {
        do {
            std::cout << "xt->child->key: " << xt->key;
            std::cout << ", xt->child->degree: " << xt->degree << std::endl;
            if(xt->child != NULL) {
                std::cout << "xt->child->child->key: " << xt->child->key << std::endl;
            }
            xt = xt->right;
        } while(xt != child);
    }
}

void print_list(node* z) {
    node* xt = z;
    if(xt != NULL) {
        do {
            std::cout << "xt->key: " << xt->key;
            std::cout << ", xt->degree: " << xt->degree << std::endl;
            if(xt->child != NULL) {
                print_child_list(xt->child);
            }
            xt = xt->right;
        } while(xt != z);
    }
}

bool numbers_children_match(node* z, int& num_nodes) {
    bool nums_match = true;
    int num_of_nodes = 0;

    node* xt = z->child;
    if(xt != NULL) {
        do {
            num_of_nodes++;
            if(xt->child != NULL) {
                nums_match = numbers_children_match(xt, num_nodes);
                if(!nums_match) { return false; }
            }
            xt = xt->right;
        } while(xt != z->child);

        num_nodes = num_nodes + num_of_nodes;

        if(num_of_nodes == z->degree) { nums_match = true; }
        else { nums_match = false; }
    }

    return nums_match;
}

fib_props numbers_match(node* z) {

    bool nums_match = true;
    int num_nodes = 0;
    fib_props fib_heap_props = { nums_match, num_nodes };

    node* xt = z;
    if(xt != NULL) {
        do {
            num_nodes++;
            nums_match = numbers_children_match(xt, num_nodes);
            fib_heap_props.deg_is_num_child = nums_match;
            fib_heap_props.num_nodes = num_nodes;
            if(!nums_match) { return fib_heap_props; }
            xt = xt->right;
        } while(xt != z);
    }

    fib_heap_props.deg_is_num_child = nums_match;
    fib_heap_props.num_nodes = num_nodes;

    return fib_heap_props;
}

bool is_fib_heap_children(node* z) {
    bool is_fibheap = true;

    node* xt = z->child;
    if(xt != NULL) {
        do {
            if(xt->p->key > xt->key) {
                return is_fibheap = false;
            }
            if(xt->child != NULL) {
                is_fibheap = is_fib_heap_children(xt);
                if(!is_fibheap) { return false; }
            }
            xt = xt->right;
        } while(xt != z->child);
    }

    return is_fibheap;
}

void nullify_children_parent_node(node* z) {
    node* xt = z->child;
    if(xt != NULL) {
        do {
            xt->p = NULL;
            xt = xt->right;
        } while(xt != z->child);
    }
}

bool is_fib_heap(node* z) {
    bool is_fibheap = true;

    node* xt = z;
    if(xt != NULL) {
        do {
            is_fibheap = is_fib_heap_children(xt);
            if(!is_fibheap) { return false; }
            xt = xt->right;
        } while(xt != z);
    }

    return is_fibheap;
}

node* fib_heap_extract_min(FibHeap* H) {

    node* z = H->min;

    if(z != NULL) {
        //Add each child of z to root list
        node* y = z->child;
        if(y != NULL) {
            //Set children's parent node to NULL
            nullify_children_parent_node(z);

            y->left->right = z->right;
            z->right->left = y->left;
            y->left = z;
            z->right = y;
            z->degree = 0;

            z->child = NULL;
        }

        //Remove z from root list
        z->left->right = z->right;
        z->right->left = z->left;

        if(z == z->right) {
            H->min = NULL;
        }
        else {
            H->min = z->right;
            consolidate(H);
        }

        H->n = H->n - 1;

    }

    return z;

}

void cut(FibHeap* H, node* x, node* y) {

    //If x is only child set child of parent to null
    if(x == x->right) {
        y->child = NULL;
        y->degree = 0;
    }
    else {
        y->child = x->right;
        y->degree = y->degree - 1;
    }

    //Remove x from child list of y and add x to root list of H
    x->left->right = x->right;
    x->right->left = x->left;

    x->right = H->min->right;
    x->left = H->min;

    H->min->right->left = x;
    H->min->right = x;

    x->p = NULL;
    x->mark = false;
}

void cascading_cut(FibHeap* H, node* y) {
    node* z = y->p;
    if(z != NULL) {
        if(y->mark == false) {
            y->mark = true;
        }
        else {
            cut(H, y, z);
            cascading_cut(H, z);
        }
    }
}

void fib_heap_decrease_key(FibHeap* H, node* x, int k) {
    if(k > x->key) {
        const char* s = "new key is greater than current key";
        std::cout << s << std::endl;
        throw s;
    }

    x->key = k;
    node* y = x->p;
    if(y != NULL && x->key < y->key) {
        cut(H, x, y);
        cascading_cut(H, y);
    }

    if(x->key < H->min->key) {
        H->min = x;
    }
}

void relax(node* u, node* v, int** w, FibHeap* H) {

    if(v->key > u->key + w[u->index][v->index]) {
        int weight = u->key + w[u->index][v->index];
        fib_heap_decrease_key(H, v, weight);
        v->key = weight;
    }
}

int map_index(int n, int index, int s) {
    int r;

    if(index >= s) { r = index - s; }
    else { r = n - s + index; }

    return r;
}

void set_weight_mat_and_ref(int size_graph,
                            std::vector< std::vector<int> >& edges,
                            int start_vertex,
                            FibHeap* H,
                            int** weight_mat,
                            node** node_refs) {


    //Initialize and construct heap
    for(int i = 0; i < size_graph; ++i) {
        tot_num_ops++;
        node_refs[i] = new node;
        node_refs[i]->key = inf;
        node_refs[i]->index = i;
        if(i == 0) {
            node_refs[i]->key = 0;
        }
        fib_heap_insert(H, node_refs[i]);
    }

    //Set weight  matrix and adjacent nodes
    int num_edges = (int) edges.size();
    for(int i = 0; i < num_edges; ++i) {
        tot_num_ops++;
        int start_index = edges[i][0] - 1;
        int end_index = edges[i][1] - 1;
        int weight = edges[i][2];

        int start = map_index(size_graph, start_index, start_vertex);
        int end = map_index(size_graph, end_index, start_vertex);

        node_refs[start]->adj_nodes.push_back(end);
        node_refs[end]->adj_nodes.push_back(start);

        weight_mat[start][end] = weight;
        weight_mat[end][start] = weight;
    }

    //Traverse edges again to pick minimum weights
    for(int i = 0; i < num_edges; ++i) {
        tot_num_ops++;
        int start_index = edges[i][0] - 1;
        int end_index = edges[i][1] - 1;
        int weight = edges[i][2];

        int start = map_index(size_graph, start_index, start_vertex);
        int end = map_index(size_graph, end_index, start_vertex);

        bool is_greater = weight_mat[start][end] >= weight;
        if(is_greater) {
            weight_mat[start][end] = weight;
            weight_mat[end][start] = weight;
        }
    }
}

bool check_fib_heap(FibHeap* H) {
    //This is the general test for the fibonacci heap.
    //The function returns true if the heap satisfies
    //the fibonacci heap properties

    //Compute heap properties
    fib_props fh_props = numbers_match(H->min);
    bool heap_is_fibheap = is_fib_heap(H->min);

    //Check if number of children equal node degrees
    bool deg_is_num_childs = fh_props.deg_is_num_child;

    //Check if number of nodes counted in heap equals H.n
    int num_nodes = fh_props.num_nodes;
    bool num_nodes_match = (num_nodes == H->n);

    //Check to see if heap is properly structured
    bool heap_is_ok = num_nodes_match && deg_is_num_childs && heap_is_fibheap;

    return heap_is_ok;
}

void dijkstra(FibHeap* H, int** w, node** node_refs) {

    //Perform Dijkstra's algorithm
    while(H->n > 0) {
        node* u = fib_heap_extract_min(H);

        int num_adj_nodes = (int) u->adj_nodes.size();
        for(int i = 0; i < num_adj_nodes; ++i) {
            tot_num_ops++;
            int index_ref = u->adj_nodes[i];
            node* v = node_refs[index_ref];
            relax(u, v, w, H);
        }
    }
}

std::vector<int> reorder_results(FibHeap* H, int n, int s, node** node_refs) {

    std::vector<int> results;
    for(int i = 0; i < n; ++i) {
        tot_num_ops++;
        if(i != s) {
            int index = map_index(n, i, s);
            if(node_refs[index]->key == inf) {
                results.push_back(-1);
            }
            else {
                results.push_back(node_refs[index]->key);
            }
        }
    }

    return results;
}

std::vector<int> shortest_reach(int n, std::vector< std::vector<int> >& edges, int s) {

    //Declarations
    FibHeap H;
    std::vector<int> results;

    //Map start vertex index
    s = s - 1;

    //Initialize weight matrix and node references
    int** weight_mat = int2D(n);
    node** node_refs = new node*[n];

    //Set weight matrix and create heap
    set_weight_mat_and_ref(n, edges, s, &H, weight_mat, node_refs);

    //Perform Dijkstra's algorithm
    dijkstra(&H, weight_mat, node_refs);

    //Reorder results
    results = reorder_results(&H, n, s, node_refs);

    //Deallocate memory
    free_int2D(weight_mat, n);
    free_node_ref(node_refs, n);

    return results;
}

