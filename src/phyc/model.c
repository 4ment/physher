/*
 *  model.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 31/1/12.
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


#include "model.h"

#include <stdlib.h>
#include <assert.h>

#include "mstring.h"


static void dummy( Model *model ){}
static void dummyUpdate( Model *self, Model *model, int index ){}

//static void noAction( Listeners *model ){}

#pragma mark -

Model * new_Model( const char *name, void *obj, int paramCount ){
	Model *model = (Model*)malloc(sizeof(Model));
    assert(model);
	model->name = String_clone(name);
	model->obj = obj;
	model->update  = dummyUpdate;
	model->store   = dummy;
	model->restore = dummy;
	
	model->listeners = new_ListenerList(paramCount);
	return model;
}

void free_Model( Model *model ){
	free(model->name);
	free_ListenerList(model->listeners);
	free(model);
}

#pragma mark -

ListenerList * new_ListenerList( const unsigned capacity ){
	ListenerList *listeners = (ListenerList*)malloc( sizeof(ListenerList));
    assert(listeners);
	listeners->capacity = (capacity == 0 ? 1 : capacity);
	listeners->count = 0;
	listeners->list = (Listener**)malloc( listeners->capacity * sizeof(Listener*));
    assert(listeners->list);
	return listeners;
}

void free_ListenerList( ListenerList *listeners ){
	for ( int i = 0; i < listeners->count; i++ ) {
		free_Listener(listeners->list[i]);
	}
	free(listeners->list);
	free(listeners);
}

void ListenerList_add( ListenerList *listeners, Model *model ){
	if ( listeners->count == listeners->capacity) {
		listeners->capacity++;
		listeners->list = realloc(listeners->list, listeners->capacity*sizeof(Listener*));
	}
	listeners->list[listeners->count] = new_Listener(model);
	listeners->count++;
}

void ListenerList_remove( ListenerList *listeners, Listener *listener ){
	int i = 0;
	for ( ; i < listeners->count; i++ ) {
		if ( listeners->list[i] == listener ) {
			break;
		}
	}
	if ( i == listeners->count) {
		return;
	}
	free_Listener(listeners->list[i]);
	i++;
	for ( ; i < listeners->count; i++ ) {
		listeners->list[i-1] = listeners->list[i];
	}
	listeners->list[listeners->count-1] = NULL;
	listeners->count--;
}

void ListenerList_remove_all( ListenerList *listeners ){
	int i = 0;
	for ( ; i < listeners->count; i++ ) {
		free_Listener(listeners->list[i]);
		listeners->list[i] = NULL;
	}
	listeners->count = 0;
}


void ListenerList_fire( ListenerList *listeners, Model *model, int index ){
	for ( int i = 0; i < listeners->count; i++ ) {
		listeners->list[i]->fire( listeners->list[i]->model, model, index );
	}
}

#pragma mark -

Listener * new_Listener( Model *model ){
	Listener *listener = (Listener*)malloc( sizeof(Listener));
    assert(listener);
	listener->model = model;
	listener->fire = NULL;
	return listener;
}


void free_Listener( Listener *listener ){
	free(listener);
}



