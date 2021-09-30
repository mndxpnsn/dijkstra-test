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

const int SETVAR = 314159;
const int INF = 3e+8;
const int inf = 3e+8;
int tot_num_ops = 0;

typedef struct FibHeapProperties {
    bool deg_is_num_child;
    int num_nodes;
} fib_props;


typedef struct Node {
    Node* left;
    Node* right;
    Node* p;
    Node* child;

    std::vector<int> adj_nodes;

    int key;
    int degree;
    int index;
    bool mark;

    int index_og;
    Node* pi;
} node;


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

//Dijkstra's algorithm using binary min heaps below

class Heap {
private:
    int heap_size;
    int size_array;
    node** A;
    node** n_ref;
    int* element_map;

    int parent(int i);
    int left(int i);
    int right(int i);
    void min_heapify(node* A[], int i);

public:
    Heap(int size);
    ~Heap();

    void set_heap(node B[]);
    void get_heap(node B[]);
    int get_heap_size();
    node* heap_extract_min();
    void heap_decrease_key(int index, double key);

    node* get_heap_element(int index);
    int get_root_index();
    void print_element_map();
    int get_heap_index(int index);

    void build_min_heap();
    bool min_heap_verify();
    void print_heap();
};

Heap::Heap(int size) {
    heap_size = size;
    A = new node*[size+1];
    element_map = new int[size+1];
    size_array = size + 1;
    n_ref = new node*[size+1];
    element_map[0] = 0;
    for(int i = 1; i <= heap_size; ++i) {
    	tot_num_ops++;
        element_map[i] = i;
        A[i] = new node;
        A[i]->pi = NULL;

    }
}

Heap::~Heap() {
//    for(int i = 0; i < size_array; ++i) {
//    	delete A[i];
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

void Heap::set_heap(node B[]) {
    for(int i = 1; i < heap_size + 1; ++i) {
    	tot_num_ops++;
        A[i] = &B[i];
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

void set_index_map(int size_graph, int* index_map, int* index_map_inverse, int* index_map_end, int s) {
    index_map[0] = index_map_inverse[0] = 0; //Point to zero for unused element

    int index_track = 1;
    for(int i = s; i <= size_graph; ++i) {
    	tot_num_ops++;
        index_map[i] = index_track;
        index_map_inverse[index_track] = i;
        index_map_end[i] = 0;
        index_track++;
    }
    for(int i = 1; i <= s - 1; ++i) {
    	tot_num_ops++;
        index_map[i] = index_track;
        index_map_inverse[index_track] = i;
        index_map_end[i] = 0;
        index_track++;
    }
}

void populate_adj_and_weight_hr(int* index_map,
		                        int* index_map_inverse,
		                        int** adj_mat,
		                        int** weight_mat,
								int size_graph,
								std::vector< std::vector<int> > &edges,
								node* heap,
								int s) {

    int** elem_is_set = int2D(size_graph + 1);

    for(int i = 1; i < size_graph + 1; ++i) {
    	tot_num_ops++;
        heap[i].key = INF;
        heap[i].pi = NULL;
        heap[i].index = i;
        heap[i].index_og = index_map_inverse[i];
    }
    heap[1].key = 0;

    int num_edges = edges.size();
    for(int i = 0; i < num_edges; ++i) {
    	tot_num_ops++;
        int start = index_map[edges[i][0]];
        int end = index_map[edges[i][1]];
        int weight = edges[i][2];

        heap[start].adj_nodes.push_back(end);
        heap[end].adj_nodes.push_back(start);

        if(elem_is_set[start][end] != SETVAR) {
            weight_mat[start][end] = weight_mat[end][start] = weight;
            elem_is_set[start][end] = elem_is_set[end][start] = SETVAR;
        }
        else if(elem_is_set[start][end] == SETVAR && weight_mat[start][end] >= weight) {
            weight_mat[start][end] = weight_mat[end][start] = weight;
        }
        adj_mat[start][end] = adj_mat[end][start] = SETVAR;
    }

    //Deallocate set flags
    free_int2D(elem_is_set, size_graph + 1);
}

std::vector<int> shortest_reach2(int n, std::vector< std::vector<int> > &edges, int s) {

    std::vector<node> rs_S;

    //Set index maps
    int* index_map = new int[n+1];
    int* index_map_inverse = new int[n+1];
    int* index_map_end = new int[n+1];
    set_index_map(n, index_map, index_map_inverse, index_map_end, s);

    //Initialize weight and adjacency matrices and binary min heap
    node* heap = new node[n + 1];
    int** adj_mat = int2D(n + 1);
    int** weight_mat = int2D(n + 1);

    populate_adj_and_weight_hr(index_map, index_map_inverse, adj_mat, weight_mat, n, edges, heap, s);

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

        rs_S.push_back(*u);
        index_map_end[u->index_og] = rs_elem_counter;
        rs_elem_counter++;
    }

    //Reorder results
    int size_results = rs_S.size();
    std::vector<int> rs_S_reordered;

    for(int i = 1; i <= size_results; ++i) {
    	tot_num_ops++;
        int j = index_map_end[i];
        if(rs_S[j].index_og != s) {
            int key = rs_S[j].key;
            if(key == INF) { key = -1; }
            rs_S_reordered.push_back(key);
        }
    }

    //Deallocate memory
    free_int2D(adj_mat, n + 1);
    free_int2D(weight_mat, n + 1);
    delete [] heap;
    delete [] index_map;
    delete [] index_map_inverse;
    delete [] index_map_end;

    return rs_S_reordered;
}


//Dijkstra's algorithm using fibonacci heaps below

class FibHeap {
public:
    int n;
    node* min;
    FibHeap() { min = NULL; n = 0; }
};



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
        if(xt->right != z) {
            while(xt->right != z) {
                std::cout << "xt->key: " << xt->key;
                std::cout << ", xt->degree: " << xt->degree << std::endl;
                xt = xt->right;
            }
            if(xt->right == z) {
                std::cout << "xt->key: " << xt->key;
                std::cout << ", xt->degree: " << xt->degree << std::endl;
            }
        }
        else {
            std::cout << "X == X->RIGHT" << std::endl;
            std::cout << "xt->key: " << xt->key;
            std::cout << ", xt->degree: " << xt->degree << std::endl;
        }
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

void consolidate(FibHeap* H) {

    double golden = (1.0 + sqrt(5.0)) / 2.0;
    double f = log(H->n) / log(golden);
    int D = floor(f + 0.01) + 1;

    node** A = new node*[D + 2];
    for(int i = 0; i < D + 2; ++i) {
        tot_num_ops++;
        A[i] = NULL;
    }

    node* x = H->min;
    if(x != NULL) {
        //Root list has more than one node
        if(x->right != H->min) {

            //Ensure all root nodes have unique degrees
            bool there_is_dup = true;
            while(there_is_dup) {
                there_is_dup = false;
                x = H->min;
                while(x->right != H->min) {
                    tot_num_ops++;
                    int d = x->degree;
                    if(A[d] != NULL && A[d] != x) {
                        there_is_dup = true;
                        node* y = A[d];
                        if(y->key > x->key) {
                            //Make y child of x;
                             make_child_of(H, y, x);

                             A[d] = NULL;
                             A[d+1] = x;

                            if(y == H->min) {
                                H->min = x;
                            }
                        }
                        else {
                            //Make x child of y;
                            make_child_of(H, x, y);

                            A[d] = NULL;
                            A[d+1] = y;

                            if(x == H->min) {
                                H->min = y;
                            }

                            x = y;
                        }
                    }
                    else {
                        A[d] = x;
                    }
                    x = x->right;
                }

                if(x->right == H->min) {
                    tot_num_ops++;
                    int d = x->degree;
                    if(A[d] != NULL && A[d] != x) {
                        there_is_dup = true;
                        node* y = A[d];
                        if(y->key > x->key) {
                            //Make y child of x;
                            make_child_of(H, y, x);

                            A[d] = NULL;
                            A[d+1] = x;

                            if(y == H->min) {
                                H->min = x;
                            }
                        }
                        else {
                            //Make x child of y;
                            make_child_of(H, x, y);

                            A[d] = NULL;
                            A[d+1] = y;

                            if(x == H->min) {
                                H->min = y;
                            }

                            x = y;
                        }
                    }
                    else {
                        A[d] = x;
                    }
                }
            }
        }
        //Root list has only one node
        else {
            int d = x->degree;
            A[d] = x;
        }
    }

    //Reconstruct root list
    H->min = NULL;
    for(int i = 0; i < D + 2; ++i) {
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
}

void print_child_list(node* child) {
    node* xt = child;
    if(xt != NULL) {
        if(xt->right != child) {
            while(xt->right != child) {
                std::cout << "xt->child->key: " << xt->key;
                std::cout << ", xt->child->degree: " << xt->degree << std::endl;
                if(xt->child != NULL) {
                    std::cout << "xt->child->child->key: " << xt->child->key << std::endl;
                }
                xt = xt->right;
            }
            if(xt->right == child) {
                std::cout << "xt->child->key: " << xt->key;
                std::cout << ", xt->child->degree: " << xt->degree << std::endl;
                if(xt->child != NULL) {
                    if(xt->child != NULL) {
                        std::cout << "xt->child->child->key: " << xt->child->key << std::endl;
                    }
                }
            }
        }
        else {
            std::cout << "X->CHILD == X->CHILD->RIGHT" << std::endl;
            std::cout << "xt->child->key: " << xt->key;
            std::cout << ", xt->child->degree: " << xt->degree << std::endl;
        }
    }
}

void print_list(node* z) {
    node* xt = z;
    if(xt != NULL) {
        if(xt->right != z) {
            while(xt->right != z) {
                std::cout << "xt->key: " << xt->key;
                std::cout << ", xt->degree: " << xt->degree << std::endl;
                if(xt->child != NULL) {
                    print_child_list(xt->child);
                }
                xt = xt->right;
            }
            if(xt->right == z) {
                std::cout << "xt->key: " << xt->key;
                std::cout << ", xt->degree: " << xt->degree << std::endl;
                if(xt->child != NULL) {
                    if(xt->child != NULL) {
                        print_child_list(xt->child);
                    }
                }
            }
        }
        else {
            std::cout << "X == X->RIGHT" << std::endl;
            std::cout << "xt->key: " << xt->key;
            std::cout << ", xt->degree: " << xt->degree << std::endl;
            if(xt->child != NULL) {
                print_child_list(xt->child);
            }
        }
    }
}

bool numbers_children_match(node* z, int& num_nodes) {
    bool nums_match = true;
    int num_of_nodes = 0;

    node* xt = z->child;
    if(xt != NULL) {
        while(xt->right != z->child) {
            num_of_nodes++;
            if(xt->child != NULL) {
                nums_match = numbers_children_match(xt, num_nodes);
                if(!nums_match) { return false; }
            }
            xt = xt->right;
        }
        if(xt->right == z->child) {
            num_of_nodes++;
            if(xt->child != NULL) {
                nums_match = numbers_children_match(xt, num_nodes);
                if(!nums_match) { return false; }
            }
        }

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
        while(xt->right != z) {
            num_nodes++;
            nums_match = numbers_children_match(xt, num_nodes);
            fib_heap_props.deg_is_num_child = nums_match;
            fib_heap_props.num_nodes = num_nodes;
            if(!nums_match) { return fib_heap_props; }
            xt = xt->right;
        }
        if(xt->right == z) {
            num_nodes++;
            nums_match = numbers_children_match(xt, num_nodes);
            fib_heap_props.deg_is_num_child = nums_match;
            fib_heap_props.num_nodes = num_nodes;
            if(!nums_match) { return fib_heap_props; }
        }
    }

    fib_heap_props.deg_is_num_child = nums_match;
    fib_heap_props.num_nodes = num_nodes;

    return fib_heap_props;
}

bool is_fib_heap_children(node* z) {
    bool is_fibheap = true;

    node* xt = z->child;
    if(xt != NULL) {
        while(xt->right != z->child) {
            if(xt->p->key > xt->key) {
                return is_fibheap = false;
            }
            if(xt->child != NULL) {
                is_fibheap = is_fib_heap_children(xt);
                if(!is_fibheap) { return false; }
            }
            xt = xt->right;
        }
        if(xt->right == z->child) {
            if(xt->p->key > xt->key) {
                return is_fibheap = false;
            }
            if(xt->child != NULL) {
                is_fibheap = is_fib_heap_children(xt);
                if(!is_fibheap) { return false; }
            }
        }
    }

    return is_fibheap;
}

void nullify_children_parent_node(node* z) {
    node* xt = z->child;
    if(xt != NULL) {
        while(xt->right != z->child) {
            xt->p = NULL;
            xt = xt->right;
        }
        if(xt->right == z->child) {
            xt->p = NULL;
        }
    }
}

bool is_fib_heap(node* z) {
    bool is_fibheap = true;

    node* xt = z;
    if(xt != NULL) {
        while(xt->right != z) {
            is_fibheap = is_fib_heap_children(xt);
            if(!is_fibheap) { return false; }
            xt = xt->right;
        }
        if(xt->right == z) {
            is_fibheap = is_fib_heap_children(xt);
            if(!is_fibheap) { return false; }
        }
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

void set_index_map(int size_graph, int* index_map, int s) {

    int index_track = 0;
    for(int i = s; i < size_graph; ++i) {
        tot_num_ops++;
        index_map[i] = index_track;
        index_track++;
    }
    for(int i = 0; i < s; ++i) {
        tot_num_ops++;
        index_map[i] = index_track;
        index_track++;
    }
}

void set_adj_and_weight_mat_and_ref(FibHeap* H,
                                    int* index_map,
                                    int** adj_mat,
                                    int** weight_mat,
                                    int size_graph,
                                    std::vector< std::vector<int> >& edges,
                                    node** v_ref) {


    for(int i = 0; i < size_graph; ++i) {
        tot_num_ops++;
        v_ref[i] = new node;
        v_ref[i]->key = inf;
        v_ref[i]->index = i;
        if(i == 0) {
            v_ref[i]->key = 0;
        }
        fib_heap_insert(H, v_ref[i]);
    }

    //Set weight and adjacency matrices and node references
    int num_edges = (int) edges.size();
    int** elem_is_set = int2D(size_graph);
    for(int i = 0; i < num_edges; ++i) {
        tot_num_ops++;
        int start_index = edges[i][0] - 1;
        int end_index = edges[i][1] - 1;
        int weight = edges[i][2];

        int start = index_map[start_index];
        int end = index_map[end_index];

        v_ref[start]->adj_nodes.push_back(end);
        v_ref[end]->adj_nodes.push_back(start);

        if(elem_is_set[start][end] != SETVAR) {
            weight_mat[start][end] = weight_mat[end][start] = weight;
            elem_is_set[start][end] = elem_is_set[end][start] = SETVAR;
        }
        else if(elem_is_set[start][end] == SETVAR && weight_mat[start][end] >= weight) {
            weight_mat[start][end] = weight_mat[end][start] = weight;
        }
        adj_mat[start][end] = adj_mat[end][start] = SETVAR;
    }

    //Deallocate node flags
    free_int2D(elem_is_set, size_graph);
}

bool check_fib_heap(FibHeap* H) {
    /*This is the general test for the fibonacci heap.
      The function returns true if the heap satisfies
      the fibonacci heap properties
     */

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

void dijkstra(FibHeap* H, int** w, node** v_ref) {

    //Perform Dijkstra's algorithm
    while(H->n > 0) {
        node* u = fib_heap_extract_min(H);

        int num_adj_nodes = (int) u->adj_nodes.size();
        for(int i = 0; i < num_adj_nodes; ++i) {
            tot_num_ops++;
            int index_ref = u->adj_nodes[i];
            node* v = v_ref[index_ref];
            relax(u, v, w, H);
        }
    }
}

std::vector<int> shortest_reach(int n, std::vector< std::vector<int> >& edges, int s) {

    //Declarations
    FibHeap H;
    const int inf = 3e+8;

    //Set index map
    s = s - 1; //Substract 1 from start index
    int* index_map = new int[n];
    set_index_map(n, index_map, s);

    //Initialize heap
    int num_nodes = n;
    node** v_ref = new node*[num_nodes];

    //Initialize weight and adjacency matrices
    int** adj_mat = int2D(n);
    int** weight_mat = int2D(n);

    set_adj_and_weight_mat_and_ref(&H, index_map, adj_mat, weight_mat, n, edges, v_ref);

    //Perform Dijkstra's algorithm
    dijkstra(&H, weight_mat, v_ref);

    //Reorder results
    std::vector<int> rs_S_reordered;
    for(int i = 0; i < n; ++i) {
        tot_num_ops++;
        if(i != s) {
            int index = index_map[i];
            if(v_ref[index]->key == inf) {
                rs_S_reordered.push_back(-1);
            }
            else {
                rs_S_reordered.push_back(v_ref[index]->key);
            }
        }
    }

    //Deallocate memory
    free_int2D(adj_mat, n);
    free_int2D(weight_mat, n);
    free_node_ref(v_ref, n);
    delete [] index_map;

    return rs_S_reordered;
}

int main(int argc, char* argv[]) {

    //Declarations
    int s = 2; //Start vertex must be greater or equal to 1
    int n = 2499; //Number of vertices
    int num_edges = 3125250; //Number of edges

    //Create edges
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

    //Print results
    float tot_num_ops_est = 10*n + 3*num_edges + 6.4*n*log(n)/log(2);
    float complexity_ratio = tot_num_ops / tot_num_ops_est;
//    int size_results = (int) results.size();
//    for(int i = 0; i < size_results; ++i) {
//        std::cout << results[i] << " ";
//    }
    std::cout << std::endl;
    std::cout << "timing execution: " << time << std::endl;
    std::cout << "number of operations estimated 10V + 3E + 6.4VlgV: " << tot_num_ops_est << std::endl;
    std::cout << "number of operations measured: " << tot_num_ops << std::endl;
    std::cout << "complexity ratio: " << complexity_ratio << std::endl;
    std::cout << "done" << std::endl;

    return 0;
}
