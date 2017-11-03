/*
 *  node.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 1/2/13.
 *  Copyright (C) 2016 Mathieu Fourment. All rights reserved.
 *
 *  This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along with this program; if not,
 *  write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#include "node.h"

#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "parameters.h"


//#define HEIGHT_MAX 1000000

#pragma mark -

Node * new_Node( Node *parent, const char *nodename, const int counter ){
	Node *n = (Node*)malloc(sizeof(Node));
	assert(n);
	n->parent = parent;
	n->right = NULL;
	n->left = NULL;
	
	n->depth = -1;
	n->id = counter;
	
	n->class_id = -1;
	
	char name[50];
	if( nodename == NULL ){
		strcpy(name, "node");
		sprintf(name+4, "%d", counter);
		n->name = String_clone(name);
	}
	else{
		n->name = String_clone(nodename);
	}
	
	n->distance = new_Parameter_with_postfix(n->name, POSTFIX_DISTANCE, BL_DEFAULT, new_Constraint(BL_MIN, BL_MAX));
	n->height   = new_Parameter_with_postfix(n->name, POSTFIX_HEIGHT, 0, new_Constraint(0,0));
	
	n->postorder_idx = 0;
	n->preorder_idx  = 0;
	
	n->time = 0;
	n->info = NULL;
    n->annotation = NULL;
	n->poly = false;
    
#ifdef TIMETEST
    n->timeParameter = new_Parameter_with_postfix(n->name, POSTFIX_DISTANCE, 0.1, new_Constraint(0.0001, 0.9999));
#endif
	
	return n;
}

void free_Node( Node *node ){
	free(node->name);
	if( node->height != NULL){
		free_Parameter( node->height );
	}
	if( node->distance != NULL){
		free_Parameter( node->distance );
	}
	if( node->info != NULL ) free(node->info);
    if ( node->annotation != NULL ) {
        free_Hashtable(node->annotation);
    }
#ifdef TIMETEST
    free_Parameter(node->timeParameter);
#endif
	free(node);
	node = NULL;
}

// TODO: clone Hashtable
Node * clone_Node( const Node *node){
	Node *n = (Node*)malloc(sizeof(Node));
	assert(n);
	
	n->id = node->id;
	n->name = String_clone( node->name);
	
	n->parent = NULL;
	n->right  = NULL;
	n->left   = NULL;
	
	n->depth = node->depth;
	n->time  = node->time;
	
	n->postorder_idx = node->postorder_idx;
	n->preorder_idx  = node->preorder_idx;
	
	n->class_id = node->class_id;
	
	n->height   = ( node->height == NULL ? NULL : clone_Parameter( node->height, true));
	n->distance = ( node->distance == NULL ? NULL : clone_Parameter( node->distance, true));
	
	n->info = NULL;
	if( node->info != NULL ) n->info = String_clone(node->info);
    n->annotation = NULL;
    if ( node->annotation != NULL ) {
        // TODO
    }
	n->poly = node->poly;
    
    
#ifdef TIMETEST
    n->timeParameter = clone_Parameter( node->timeParameter, true);
#endif
	return n;
}

#pragma mark -

void Node_set_height( Node *node, const double value ){
	Parameter_set_value(node->height, value);
}

double Node_height( const Node *node ){
	return Parameter_value(node->height);
}

double Node_time_elapsed( Node *node ){
	if( Node_isroot(node) ){
		return -1;
	}
//#ifdef TIMETEST
//    return Node_t(node);
//#endif
	return (Node_height( Node_parent(node) ) - Node_height(node));
}

#ifdef TIMETEST

void height_to_time(Node *node){
    if( !Node_isroot(node) ){
        node->timeParameter->value = (Node_height( Node_parent(node) ) - Node_height(node));
    }
    
    if( !Node_isleaf(node) ){
        height_to_time(Node_left(node));
        height_to_time(Node_right(node));
    }
}

void Node_set_t( Node *node, const double value ){
    Parameter_set_value(node->timeParameter, value);
}


double Node_t( const Node *node ){
    return Parameter_value(node->timeParameter);
}
#endif

void Node_set_distance( Node *node, const double value ){
    Parameter_set_value(node->distance, value);
}


double Node_distance( const Node *node ){
    return Parameter_value(node->distance);
}

void Node_set_parent( Node *node, Node *parent ){
	node->parent = parent;
}

Node * Node_parent( Node *node ){
	return node->parent;
}

bool Node_isroot( const Node *node ){
	return ( node->parent == NULL ? true : false );
}

bool Node_isleaf( const Node *node ){
	return ( node->left == NULL ? true : false );
}

Node * Node_left( Node *node ){
	return node->left;
}

Node * Node_right( Node *node ){
	return node->right;
}

Node * Node_sibling( Node *node ){
    if ( Node_isroot(node)) {
        return NULL;
    }
    if ( node->parent->left == node ) {
        return node->parent->right;
    }
	return node->parent->left;
}

bool Node_isancestor(  Node  * const node, const Node *putative_ancestor ){
	if ( node == putative_ancestor ) {
		return false;
	}
	Node *n = node;
	while ( n != NULL && n != putative_ancestor ) {
		n = n->parent;
	}
	if( n == putative_ancestor ) return true;
	else return false;
}

bool Node_isrelated( Node *aNode1, Node *aNode2 ){
	if ( aNode1 == aNode2 ) {
		return false;
	}
	
	Node *a1 = aNode1;
	Node *a2 = aNode2;
	
	bool flag1 = false;
	bool flag2 = false;
	
	while( a1->depth != a2->depth ){
		if( a1->depth > a2->depth){
			a1 = a1->parent;
			flag1 = true;
		}
		else{
			a2 = a2->parent;
			flag2 = true;
            
		}
	}
	if( flag1 && flag2 ) return 0;
	else if( flag1 && a1 == a2 ) return 2; // node1 is the ancestor
	else if( flag2 && a1 == a2 ) return 1;
	//fprintf(stderr, "%d %d %s %s\n",flag1,flag2,a1->name,a2->name);
	return 0;
}

// Number of edges in the the shortest path between 2 nodes
int Node_graph_distance( Node *node1, Node *node2 ){
    Node *a = node1;
    Node *b = node2;
    int d = 0;
    while ( a != b ) {
        d++;
        if ( a->depth < b->depth ) {
            b = Node_parent(b);
        }
        else if( a->depth > b->depth ){
            a = Node_parent(a);
        }
        else {
            d++;
            a = Node_parent(a);
            b = Node_parent(b);
        }
    }
    return d;
}

static void _count_tips( Node *node, int *count ){
    if( !Node_isleaf(node)){
        _count_tips(node->left, count);
        _count_tips(node->right, count);
    }
    else{
        (*count)++;
    }
}
int Node_tip_count( Node *node ){
    int count = 0;
    _count_tips(node, &count);
	return count;
}

void Node_rotate( Node *node ){
    Node *temp = node->left;
    node->left = node->right;
    node->right = temp;
}

void Node_swap_parents( Node *node1, Node *node2 ){
    Node *parent1 = Node_parent(node1);
    Node *parent2 = Node_parent(node2);
    
    Node_set_parent(node2, parent1);
    Node_set_parent(node1, parent2);
    
    Node_removeChild(parent1, node1);
    Node_addChild(parent1, node2);
    
    Node_removeChild(parent2, node2);
    Node_addChild(parent2, node1);
}

double Node_time( Node *node ){
	return node->time;
}

void Node_set_time( Node *node, const double time ){
	node->time = time;
}

int Node_id( const Node *node ){
	return node->id;
}

void Node_set_id( Node *node, const int id ){
	node->id = id;
}

int Node_class_id( Node *node ){
	return node->class_id;
}

void Node_set_class_id( Node *node, const int class_id ){
	node->class_id = class_id;
}

int Node_depth( Node *node ){
	return node->depth;
}

void Node_set_depth( Node *node, const int depth ){
	node->depth = depth;
}

char * Node_name( Node *node ){
	return node->name;
}

void Node_set_name( Node *node, const char *name ){
	if(node->name == NULL ){
		node->name = String_clone( name);
	}
	else {
		if( strlen(node->name) != strlen(name) ){
			node->name = realloc( node->name, (strlen(name)+1) * sizeof(char));
			assert(node->name);
		}
		strcpy(node->name, name);
	}
	
}

void Node_rename( Node *n, const int newPos ){
	char name[50] = "node";
	sprintf(name+4, "%d", newPos);
	if( strlen(n->name ) != strlen(name)  ){
		n->name = (char *)realloc(n->name, (strlen(name)+1)*sizeof(char));
		assert(n->name);
	}
	strcpy(n->name, name);
}

void Node_set_annotation( Node *node, const char *key, const char *value ){
    if( node->annotation == NULL ){
        node->annotation = new_Hashtable_string(5);
    }
    Hashtable_add(node->annotation, String_clone(key), String_clone(value));
}

void Node_empty_annotation( Node *node ){
    if( node->annotation != NULL ){
        Hashtable_empty( node->annotation );
    }
}

char * Node_annotation( Node *node, const char *key ){
    if( node->annotation != NULL ){
        if ( Hashtable_exists(node->annotation, key) ) {
            return Hashtable_get(node->annotation, key);
        }
    }
    return NULL;
}



double Node_get_double_from_info( const Node *node, const char *str ){
	double value = NAN;
	if ( node->info != NULL ) {
		char *ptr = node->info;
		if ( String_contains_str(ptr, str) ) {
			StringBuffer *buffer = new_StringBuffer(20);
			ptr += String_index_of_str(ptr, str) + strlen(str);
			while ( *ptr == ' ' || *ptr == '\t' ) ptr++;
            if( *ptr == '-' ){
                StringBuffer_append_char(buffer, *ptr);
                ptr++;
                if ( *ptr == '\0') { // str cannot be "-"
                    return NAN;
                }
            }
            int npoint = 0;
			while ( (*ptr >= 48 && *ptr <= 57) || *ptr == '.' ) {
                if ( *ptr == '.' ) npoint++;
                if (npoint > 1 ) {
                    
                    free_StringBuffer(buffer);
                    return NAN;
                }
				StringBuffer_append_char(buffer, *ptr);
				ptr++;
			}
			if ( *ptr == 'e' || *ptr == 'E' ) {
                StringBuffer_append_char(buffer, *ptr);
                ptr++;
                if( *ptr == '-' || *ptr == '+' ){
                    StringBuffer_append_char(buffer, *ptr);
                    ptr++;
                }
                if( *ptr >= 48 && *ptr <= 57){
                    while ( (*ptr >= 48 && *ptr <= 57) ) {
                        StringBuffer_append_char(buffer, *ptr);
                        ptr++;
                    }
                }
                else {
                    free_StringBuffer(buffer);
                    return NAN;
                }
            }
			value = atof(buffer->c);
			free_StringBuffer(buffer);
		}
	}
	return value;
}

int Node_get_int_from_info( const Node *node, const char *str ){
	int value = -1;
	if ( node->info != NULL ) {
		char *ptr = node->info;
		if ( String_contains_str(ptr, str) ) {
			StringBuffer *buffer = new_StringBuffer(20);
			ptr += String_index_of_str(ptr, str) + strlen(str);
			while ( *ptr == ' ' || *ptr == '\t' ) ptr++;
			while ( (*ptr >= 48 && *ptr <= 57) || *ptr == '.' || *ptr == '-' ) {
				StringBuffer_append_char(buffer, *ptr);
				ptr++;
			}
			
			value = atoi(buffer->c);
			free_StringBuffer(buffer);
		}
	}
	return value;
}

char * Node_get_string_from_info( const Node *node, const char *str ){
	char *value = NULL;
	if ( node->info != NULL ) {
		char *ptr = node->info;
		if ( String_contains_str(ptr, str) ) {
			StringBuffer *buffer = new_StringBuffer(20);
			ptr += String_index_of_str(ptr, str) + strlen(str);
			
			if( *ptr == '\"' || *ptr == '\'' ){
				ptr++;
				while( *ptr != '\'' && *ptr != '\"' ){
					StringBuffer_append_char(buffer, *ptr);
					ptr++;					
				}
				if (buffer->length !=0 ) {
					value = String_clone(buffer->c);
				}
			}
			else {
				
			}
			
			free_StringBuffer(buffer);
		}
	}
	return value;
}

void Node_removeChild( Node *parent, Node *node){
	if( parent->left == node ){
		parent->left = NULL;
	}
	else if( parent->right == node ) {
		parent->right = NULL;
	}
}

bool Node_addChild( Node *parent, Node *node){
	if( parent->left == NULL ){
		parent->left = node;
		return true;
	}
	else if( parent->right == NULL ) {
		parent->right = node;
		return true;
	}
	return false;
}

void Node_add_listener( Node *node, Model *model ){
	
}

void Node_print( Node *node ){
    printf("Name %s\n", Node_name(node));
    printf("Left  %s\n", Node_name(Node_left(node)));
    printf("Right %s\n", Node_name(Node_right(node)));
    printf("Parent %s\n", Node_name(Node_parent(node)));
    Parameter_print(node->distance);
    printf("\n");
}
