///*
// *  model.h
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
//#ifndef Math_model_h
//#define Math_model_h
//
//#include "parameters.h"
//
//struct _ListenerList;
//typedef struct _ListenerList ListenerList;
//
//struct _Listeners;
//typedef struct _Listener Listener;
//
//struct _Model;
//typedef struct _Model Model;
//
//struct Parameter;
//
//struct _ListenerList {
//	Model** models;
//	int count;
//	int capacity;
//	void (*free)( ListenerList*);
//	void (*fire)( ListenerList*, Model*, int );
//	void (*add)( ListenerList*, Model* );
//	void (*remove)( ListenerList*, Model* );
//	void (*removeAll)( ListenerList*);
//};
//
//struct _Model {
//	void *obj; // pointer to model
//	char *name;
//	double (*logP)( Model * );
//	double (*dlogP)( Model *, const struct Parameter* );
//	void (*free)( Model * );
//	void (*update)( Model *, Model *, int );
//	
//	//Parameters *parameters; // pointers to the parameters of the model
//	ListenerList *listeners;
//};
//
//
//#pragma mark -
//
//Model * new_Model( const char *name, void *obj );
//
//void free_Model( Model *model );
//
//#pragma mark -
//
//ListenerList * new_ListenerList( const unsigned capacity );
//
//
//#endif
