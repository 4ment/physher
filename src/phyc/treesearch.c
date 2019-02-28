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

// Should not be called where prune is the root. graft can be the root
bool SPR_move( Tree *tree, Node *prune, Node *graft ){
    Node *oldroot = Tree_root(tree);
    bool simple = false;
    
    assert(oldroot != prune);
    
    if( oldroot == prune ){
        fprintf(stderr, "No SPR move on the root node!!!\n");
        exit(1);
    }
    
    // prune is an ancestor of graft
    // prune becomes the root
    if( Node_isancestor(graft, prune) ){
        //printf("case 1\n");
        Node *node = Node_parent(graft);
        Node *n = node;// node to graft using the root node
        Node *sibling = Node_sibling(node);
        
        // Disconnect below prune
        Node_removeChild(node, graft);
        graft->parent = NULL;
        
        double d = Node_distance(prune);
        
        
        Node *parent = NULL;
		
		parent = Node_parent(node);
		sibling = Node_sibling(node);
		
		double dd = Node_distance(node);
		
		while ( parent != prune ) {
			dd = Node_distance(parent);
			
			sibling = Node_sibling(parent);
//			printf("%s %s %s\n", node->name, parent->name, sibling->name);
			
			Node_addChild(node, parent);
			Node_removeChild(parent, node);
			
			Node *temp = Node_parent(parent);
			
			Node_set_distance(parent, Node_distance(node));
			
			Node_set_parent(parent, node);
			
			node   = parent;
			parent = temp;
			
			Node_removeChild(parent, node);
		}
		// parent is the root
		Node_removeChild(parent, sibling);
		Node_addChild(node, sibling);
		Node_set_parent(sibling, node);
		Node_addChild(parent, graft);
		Node_set_parent(graft, parent);
		
		 // went through the while loop above
		if( Node_parent(node) != parent ){
			Node_addChild(parent, n);
			Node_set_parent(n, parent);
		}
		Node_set_distance(n, d);
		Node_set_distance(sibling, Node_distance(sibling)+dd);
//		printall(tree);
		Tree_print_newick_subtree(stdout, oldroot, true);
		printf("\n");
//
//        // Go up to the root and "rotate" each node so the parent becomes the child
//        if( !Node_isroot( node ) ){
//            parent = Node_parent(node);
//            sibling = Node_sibling(node);
//
//            double dd = Node_distance(node);
//
//            while ( !Node_isroot(parent) ) {
//                dd = Node_distance(parent);
//
//                sibling = Node_sibling(parent);
//
//                Node_addChild(node, parent);
//                Node_removeChild(parent, node);
//
//                Node *temp = Node_parent(parent);
//
//                Node_set_distance(parent, Node_distance(node));
//
//                Node_set_parent(parent, node);
//
//                node   = parent;
//                parent = temp;
//
//                Node_removeChild(parent, node);
//            }
//
//            // parent is the root
//            Node_removeChild(parent, sibling);
//            Node_addChild(node, sibling);
//            Node_set_parent(sibling, node);
//
//            // went through the while loop above
//            if( Node_parent(node) != parent ){
//                Node_addChild(parent, n);
//                Node_set_parent(n, parent);
//            }
//
//            Node_set_distance(n, d);
//            Node_set_distance(sibling, Node_distance(sibling)+dd);
//
//        }
//        else {
//            Node_set_distance(sibling, Node_distance(sibling)+d);
//        }
//
//
//        // graft
//        parent = Node_parent(graft);
//
//        double dgraft = Node_distance(graft);
//
//        Node_set_distance(oldroot, dgraft/2.0);
//        Node_set_distance(graft, dgraft/2.0);
//
//        Node_addChild(oldroot, graft);
//        Node_set_parent(graft, oldroot);
//
//        Node_removeChild(parent, graft);
//        Node_addChild(parent, oldroot);
//        Node_set_parent(oldroot, parent);
//
//        Node_set_distance(prune, -1);
//
//        Tree_set_root(tree, prune);
    }
    else {
        Node *parent1 = Node_parent(prune);
        Node *sibling = Node_sibling(prune);
        Node *parent2 = Node_parent(graft);
        
        // graft is root
        // parent1 becomes the root
        // Undo the case below and when prune is the ancestor of graft
        // WARNING: the orignal rooting is lost
		if( Node_isroot(graft) ){
			Node *grandpa = Node_parent(parent1);
			Node* right = prune;
			while (!Node_isroot(Node_parent(right))) right = Node_parent(right);
			Node* tomove = Node_sibling(right);
			// remove the node from root
			Node_removeChild(graft, tomove);
		
			// disconnect parent1, will be used to reconnect tomove to root
			Node_removeChild(parent1, sibling);
			Node_removeChild(parent1, prune);
			Node_removeChild(grandpa, parent1);
			// grandpa takes care of sibling of prune
			Node_addChild(grandpa, sibling);
			Node_set_parent(sibling, grandpa);
			// prune is reconnected under the root
			Node_addChild(graft, prune);
			Node_set_parent(prune, graft);
			
			Node_addChild(parent1, Node_left(right));
			Node_addChild(parent1, Node_right(right));
			Node_set_parent(Node_left(right), parent1);
			Node_set_parent(Node_right(right), parent1);
			Node_removeChild(right, Node_left(right));
			Node_removeChild(right, Node_right(right));
			Node_addChild(right, tomove);
			Node_addChild(right, parent1);
			Node_set_parent(tomove, right);
			Node_set_parent(parent1, right);
			
//            Node *grandpa = Node_parent(parent1);
//            Node_set_distance(sibling, Node_distance(sibling)+Node_distance(parent1));
//
//            Node_set_distance(prune, Node_distance(prune)/2);
//            Node_set_distance(graft, Node_distance(prune));
//
//            if( grandpa->left == parent1 ){
//                grandpa->left = sibling;
//            }
//            else {
//                grandpa->right = sibling;
//            }
//            sibling->parent = grandpa;
//
//            Node_removeChild(parent1, sibling);
//
//            Node_addChild(parent1, graft);
//            Node_set_parent(graft, parent1);
//
//            Tree_set_root(tree, parent1);
//            parent1->parent = NULL;
        }
        // When the child of the root is grafted on a branch => parent1 is the root
        // Sibling becomes the root
        else if( Node_isroot(parent1) ){
			Node* sibgraft = Node_sibling(graft);
			Node* pgraft = Node_parent(graft);
			printf("%s\n", sibgraft->name);
			// For NNIs
//			Node_removeChild(parent1, prune);
//			Node_removeChild(parent2, sibgraft);
//			Node_addChild(parent1, sibgraft);
//			Node_addChild(parent2, prune);
//			Node_set_parent(prune, parent2);
//			Node_set_parent(sibgraft, parent1);
			
			Node* n = graft;
			while(Node_parent(n) != sibling) n = Node_parent(n);
			Node* right = Node_parent(n); // this node should be fixed and never moved
			Node* tomove = Node_sibling(n);
			
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
//			printall(tree);
//            Tree_set_root(tree, sibling);
//
//            Node_set_distance(prune, Node_distance(prune) + Node_distance(sibling));// not meaningful, sibling might be 0
//
//            // disconnect subtree
//            Node_removeChild(parent1, sibling);
//            sibling->parent = NULL;
//
//            // destination
//            Node_removeChild(parent2, graft);
//
//            double d = Node_distance(graft);
//            double d2 = d/2.0;
//
//            Node_set_distance(graft,d2);
//
//            Node_addChild(parent1, graft);
//            Node_set_parent(graft, parent1);
//
//            Node_addChild(parent2, parent1);
//            Node_set_parent(parent1, parent2);
//
//            Node_set_distance(parent1,d-d2);
			
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
            sibling->parent = grandpa;
            
            Node_set_distance(sibling, Node_distance(parent1) + Node_distance(sibling));
            
            Node_removeChild(parent1, sibling);
            
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
            simple = true;
        }
    }
    Tree_set_topology_changed(tree);
    
    return simple;
}



// They can be switched back since node1 and node2 conserve their branch length
void NNI_move( Tree *tree, Node *node1, Node *node2 ){
    
    Node_swap_parents(node1,node2);
    
	Tree_set_topology_changed(tree);
}

