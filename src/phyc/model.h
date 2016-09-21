/*
 *  model.h
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

#ifndef Math_model_h
#define Math_model_h



struct _ListenerList;
typedef struct _ListenerList ListenerList;

struct _Listeners;
typedef struct _Listener Listener;

struct _Model;
typedef struct _Model Model;

struct _Event;
typedef struct _Event Event;

typedef void (*pAction)( Model *, Model *, int index );

struct _ListenerList {
	Listener **list;
	int count;
	int capacity;
};

struct _Listener {
	Model *model;
	pAction fire;
};




typedef void (*pUpdate)( Model *, Model *, int );

typedef void (*pStore)( Model * );
typedef void (*pRestore)( Model * );

struct _Model {
	void *obj; // pointer to model
	char *name;
	pUpdate update;
	pStore store;
	pRestore restore;
	
	//Parameters *parameters; // pointers to the parameters of the model
	ListenerList *listeners;
};


#pragma mark -

Model * new_Model( const char *name, void *obj, int paramCount );

void free_Model( Model *model );

#pragma mark -

ListenerList * new_ListenerList( const unsigned capacity );

void free_ListenerList( ListenerList *listeners );

void ListenerList_add( ListenerList *listeners, Model *model );

void ListenerList_remove( ListenerList *listeners, Listener *listener );

void ListenerList_remove_all( ListenerList *listeners );

void ListenerList_fire( ListenerList *listeners, Model *model, int index );


#pragma mark -

Listener * new_Listener( Model *model );

void free_Listener( Listener *listener );


#endif
