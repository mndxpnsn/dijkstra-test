/*
 * fib_heap.hpp
 *
 *  Created on: 7 Oct 2021
 *      Author: mndx
 */

#ifndef FIB_HEAP_HPP_
#define FIB_HEAP_HPP_

#include <vector>
#include "user_types.hpp"

bool** bool2D(const int size);
int** int2D(const int size);
void free_bool2D(bool** p, int size);
void free_int2D(int** p, int size);
void free_node_ref(node** v_ref, int size);
int map_index(int n, int index, int s);
std::vector<int> shortest_reach(int n, std::vector< std::vector<int> >& edges, int s);

#endif /* FIB_HEAP_HPP_ */
