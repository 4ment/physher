///*
// *  model.c
// *  PhyC
// *
// *  Created by Mathieu Fourment on 31/1/12.
// *  Copyright (C) 2016 Mathieu Fourment. All rights reserved.
// *
// *  This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
// *  as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
// *
// *  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
// *
// *  You should have received a copy of the GNU General Public License along with this program; if not,
// *  write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
// */
//
//
//#include "model.h"
//
//#include <stdlib.h>
//#include <assert.h>
//
//#include "mstring.h"
//
//static void dummyUpdate( Model *self, Model *model, int index ){}
//static double _logP(Model *model){return 0;}
//static double _dlogP(Model *model, const Parameter* p){return 0;}
//
//#pragma mark -
//
//
//Model * new_Model( const char *name, void *obj ){
//	Model *model = (Model*)malloc(sizeof(Model));
//    assert(model);
//	model->name = String_clone(name);
//	model->obj = obj;
//	model->logP = _logP;
//	model->logP = _dlogP;
//	model->update  = dummyUpdate;
//	model->free = free_Model;
//	
//	model->listeners = new_ListenerList(1);
//	return model;
//}
//
//void free_Model( Model *model ){
//	free(model->name);
//	model->listeners->free(model->listeners);
//	free(model);
//}
//
//#pragma mark -
//
//static void _free_ListenerList( ListenerList *listeners ){
//	free(listeners->models);
//	free(listeners);
//}
//
//static void _ListenerList_fire( ListenerList *listeners, Model* model, int index){
//	for ( int i = 0; i < listeners->count; i++ ) {
//		listeners->models[i]->update( listeners->models[i], model, index );
//	}
//}
//
//static void _ListenerList_remove( ListenerList *listeners, Model* model ){
//	int i = 0;
//	for ( ; i < listeners->count; i++ ) {
//		if ( listeners->models[i] == model ) {
//			break;
//		}
//	}
//	if ( i == listeners->count) {
//		return;
//	}
//	i++;
//	for ( ; i < listeners->count; i++ ) {
//		listeners->models[i-1] = listeners->models[i];
//	}
//	listeners->models[listeners->count-1] = NULL;
//	listeners->count--;
//}
//
//static void _ListenerList_remove_all( ListenerList *listeners ){
//	int i = 0;
//	for ( ; i < listeners->count; i++ ) {
//		listeners->models[i] = NULL;
//	}
//	listeners->count = 0;
//}
//
//static void _ListenerList_add( ListenerList *listeners, Model *model ){
//	if ( listeners->count == listeners->capacity) {
//		listeners->capacity++;
//		listeners->models = realloc(listeners->models, listeners->capacity*sizeof(Model*));
//	}
//	listeners->models[listeners->count] = model;
//	listeners->count++;
//}
//ListenerList * new_ListenerList( const unsigned capacity ){
//	ListenerList *listeners = (ListenerList*)malloc( sizeof(ListenerList));
//    assert(listeners);
//	listeners->capacity = (capacity == 0 ? 1 : capacity);
//	listeners->count = 0;
//	listeners->models = (Model**)malloc( listeners->capacity * sizeof(Model*));
//	assert(listeners->models);
//	listeners->free = _free_ListenerList;
//	listeners->add = _ListenerList_add;
//	listeners->remove = _ListenerList_remove;
//	listeners->removeAll = _ListenerList_remove_all;
//	listeners->fire = _ListenerList_fire;
//	return listeners;
//}
//
//
//
