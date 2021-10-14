/*
 * fib_heap_support.hpp
 *
 *  Created on: 11 Oct 2021
 *      Author: mndx
 */

#ifndef FIB_HEAP_SUPPORT_HPP_
#define FIB_HEAP_SUPPORT_HPP_

#include <vector>

#include "user_types.hpp"

int size_root_list(node* z);
void print_root_list(node* z);
void print_child_list(node* child);
void print_list(node* z);
bool numbers_children_match(node* z, int& num_nodes);
fib_props numbers_match(node* z);
bool is_fib_heap_children(node* z);
bool is_fib_heap(node* z);
bool check_fib_heap(FibHeap* H);

#endif /* FIB_HEAP_SUPPORT_HPP_ */
