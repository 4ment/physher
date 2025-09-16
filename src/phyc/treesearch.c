/*
 *  treesearch.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 7/08/13.
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

#include <stdio.h>
#include <assert.h>

#include "treesearch.h"
#include "treeio.h"

void printall(Tree* tree){
	for (int i = 0; i < Tree_node_count(tree); i++) {
		Node* node = Tree_node(tree, i);
		printf("%s %s %s %s\n", node->name, (node->parent!=NULL?node->parent->name:"NULL"),
			   (node->left!=NULL?node->left->name:"NULL"),
			   (node->right!=NULL?node->right->name:"NULL"));
	}
}

// root is the root node => its parent must be NULL
void _reroot(Node* root, Node* node){
	Node newroot = {0, NULL, NULL, node->parent, node, 0,0.0,0,0,NULL, NULL, 0, NULL, BL_DEFAULT, NULL, false};
	
	double branchLength = Node_distance(node->parent); // save branch length
	double midpoint = Node_distance(node)/2;
	
	Node_set_distance(node, midpoint);
	Node_set_distance(node->parent, midpoint);
	
	Node *n = node->parent;     // the node to which we need to add its parent as a child recurssively
	Node *nparent = n->parent;  // the parent which is a ref to the rest of the old tree
	
	Node_removeChild(nparent, n);
	Node_removeChild(n, node);
	
	node->parent = &newroot;
	n->parent = &newroot;
	
	while ( !Node_isroot(nparent) ) {
		Node_addChild(n, nparent);
		
		Node *temp = nparent->parent;
		double bl = Node_distance(nparent);
		Node_set_distance(nparent, branchLength);
		branchLength = bl;
		
		nparent->parent = n;
		
		n = nparent;
		nparent = temp;
		
		Node_removeChild(nparent, n);
	}
	
	if( nparent->left == NULL){
		nparent = nparent->right;
	}
	else {
		nparent = nparent->left;
	}
	Node_set_distance(nparent, Node_distance(nparent) + branchLength );
	Node_addChild(n, nparent);
	nparent->parent = n;
	
	root->left = newroot.left;
	root->right = newroot.right;
	Node_set_parent(root->left, root);
	Node_set_parent(root->right, root);
}

// Should not be called on root and its right node
Node* SPR_move( Tree *tree, Node *prune, Node *graft ){
    Node *oldroot = Tree_root(tree);
	Node* right = Node_right(oldroot);
    assert(oldroot != prune);
    
    if( oldroot == prune ){
        fprintf(stderr, "No SPR move on the root node!!!\n");
        exit(1);
    }
	Node* tripod = NULL;
	
    // prune is an ancestor of graft
	if( Node_isancestor(graft, prune) ){
        Node *node = Node_parent(prune);
		prune->parent = NULL;
		
		Node *pnode = Node_parent(graft);
		
		_reroot(prune, graft);
		Node_addChild(pnode, Node_parent(graft)); //parent->graft is actually prune
		prune->parent = node;
		tripod = prune;
    }
    else {
        Node *parent1 = Node_parent(prune);
        Node *sibling = Node_sibling(prune);
        Node *parent2 = Node_parent(graft);
		
		Node* left = Node_left(oldroot);
		
        // When the child of the root is grafted on a branch => parent1 is the root
        if( Node_isroot(parent1) ){
			Node* pgraft = Node_parent(graft);
			
			Node* n = graft;
			while(Node_parent(n) != sibling) n = Node_parent(n);
			Node* right = Node_parent(n); // this node should be fixed and never moved
			Node* tomove = Node_sibling(n);
			Node_set_distance(tomove, Node_distance(tomove) + Node_distance(n));
			
			// remove prune from root and replace it with tomove
			Node_removeChild(parent1, prune);
			Node_addChild(parent1, tomove);
			Node_set_parent(tomove, parent1);
			
			// remove kids from right
			Node_removeChild(right, tomove);
			Node_removeChild(right, n);
			// add kids of the node n below
			Node_addChild(right, Node_left(n));
			Node_addChild(right, Node_right(n));
			Node_set_parent(Node_left(n), right);
			Node_set_parent(Node_right(n), right);
			
			// remove kids from n
			Node_removeChild(n, Node_left(n));
			Node_removeChild(n, Node_right(n));
			
			// add prune and graft
			Node_addChild(n, prune);
			Node_addChild(n, graft);
			Node_set_parent(prune, n);
			Node_set_parent(graft, n);
			
			// next insert n above graft
			Node_removeChild(pgraft, graft);
			Node_addChild(pgraft, n);
			Node_set_parent(n, pgraft);
			Node_set_distance(graft, Node_distance(graft)/2.0);
			Node_set_distance(n, Node_distance(graft));
			tripod = prune->parent;
        }
		// this is an NNI involving the root
		else if((parent1 == right && parent2 == left) || (parent1 == left && parent2 == right)){
			Node_swap_parents(prune, graft);
			Parameter_fire(prune->distance, Node_id(prune));
			Parameter_fire(graft->distance, Node_id(graft));
			tripod = left;
		}
		// NNI involving the right node
		else if(Node_parent(parent2) == right && parent1 == right){
			Node_swap_parents(prune, graft);
			Parameter_fire(prune->distance, Node_id(prune));
			Parameter_fire(graft->distance, Node_id(graft));
			tripod = Node_sibling(prune);
		}
		// prune is a child of the right node
		// We do not want to detach the right node to reconnect it somewhere else
		// as we do in the simple case since this node should never change or move
		else if(parent1 == right){
			bool sameSide = false;
			
			// check if graft is also on the right side
			Node* n = graft;
			while (n->parent != oldroot) {
				n = n->parent;
			}
			sameSide = (n == right);
			
			// graft is no a descendent of the right node
			if(!sameSide){
				Node_set_distance(sibling, Node_distance(sibling)+Node_distance(left));
				// disconnect left node form its parent and its children
				Node_removeChild(oldroot, left);
				Node* n2 = graft;
				// n2 node to be attached to the root
				while (n2->parent != left) {
					n2 = n2->parent;
				}
				Node* n = Node_sibling(n2); // node to be moved to the right side
				Node_removeChild(left, n);
				Node_removeChild(left, n2);
				
				// replace prune with moving node n
				Node_removeChild(right, prune);
				Node_addChild(right, n);
				Node_set_parent(n, right);
				// attach child of left node to root node
				Node_addChild(oldroot, n2);
				Node_set_parent(n2, oldroot);
				// attach prune and graft to left node
				Node_addChild(left, prune);
				Node_set_parent(prune, left);
				Node_addChild(left, graft);
				Node_set_parent(graft, left);
				// insert left between graft and its parent
				Node_removeChild(parent2, graft);
				Node_addChild(parent2, left);
				Node_set_parent(left, parent2);
				double d = Node_distance(graft)/2.0;
				Node_set_distance(graft, d);
				Node_set_distance(left, Node_distance(graft)-d);
				tripod = left;
			}
			else{
				Node* siblingPrune = Node_sibling(prune);
				Node_set_distance(left, Node_distance(siblingPrune)+Node_distance(left));
				// disconnect sibling of prune form its parent and its children
				Node_removeChild(right, siblingPrune);
				Node* n2 = graft;
				// n2 node to be attached to the root
				while (n2->parent != siblingPrune) {
					n2 = n2->parent;
				}
				// if n2===graft then it is simply an NNI and it is taken care of above
				Node* n = Node_sibling(n2); // node to be moved to the right side
				Node_removeChild(siblingPrune, n);
				Node_removeChild(siblingPrune, n2);
				
				// replace prune with moving node n
				Node_removeChild(right, prune);
				Node_addChild(right, n);
				Node_set_parent(n, right);
				// attach child of left node to right node
				Node_addChild(right, n2);
				Node_set_parent(n2, right);
				// attach prune and graft to left node
				Node_addChild(siblingPrune, prune);
				Node_set_parent(prune, siblingPrune);
				Node_addChild(siblingPrune, graft);
				Node_set_parent(graft, siblingPrune);
				// insert left between graft and its parent
				Node_removeChild(parent2, graft);
				Node_addChild(parent2, siblingPrune);
				Node_set_parent(siblingPrune, parent2);
				double d = Node_distance(graft)/2.0;
				Node_set_distance(graft, d);
				Node_set_distance(siblingPrune, Node_distance(graft)-d);
				tripod = siblingPrune;
			}
		}
		// Simple case
		else {
			Node *grandpa = Node_parent(parent1);
			
			if( grandpa->left == parent1 ){
				grandpa->left = sibling;
			}
			else {
				grandpa->right = sibling;
			}
			
			Node_set_distance(sibling, Node_distance(parent1) + Node_distance(sibling));
			
			Node_removeChild(parent1, sibling);
			sibling->parent = grandpa;
			
			// destination
			Node_removeChild(parent2, graft);

			double d = Node_distance(graft);
			double d2 = d/2.0;
			
			Node_set_distance(graft,d2);
			
			Node_addChild(parent1, graft);
			Node_set_parent(graft, parent1);
			
			Node_addChild(parent2, parent1);
			Node_set_parent(parent1, parent2);
			
			Node_set_distance(parent1,d-d2);
			tripod = prune->parent;
		}
	}
	
	if (Node_distance(right) != 0.0) {
		Node* left = Node_left(Tree_root(tree));
		Node_set_distance(left, Node_distance(left)+Node_distance(right));
		Node_set_distance(right, 0);
	}
    Tree_set_topology_changed(tree);

    return tripod;
}



// They can be switched back since node1 and node2 conserve their branch length
void NNI_move( Tree *tree, Node *node1, Node *node2 ){
    
    Node_swap_parents(node1,node2);
	Parameter_fire(node1->distance, Node_id(node1));
	Parameter_fire(node2->distance, Node_id(node2));
	Tree_set_topology_changed(tree);
}

