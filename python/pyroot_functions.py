#!/usr/bin/env python
# -*- coding: utf-8 -*-

def loopdir_new(keys, all_names):  # loop through all subdirectories of the root file and add all fitting object to a list
	object_list = []
	for key_object in keys:
		if ('TDirectory' in key_object.GetClassName()):
			object_list= object_list+loopdir_new(key_object.ReadObj().GetListOfKeys(), all_names)
		else:
			if (all(name in key_object.GetName() for name in all_names)):
				object_list.append(key_object)
	return object_list

def plot_tree(keys, conditions, details, sizerange, drawOption):
	for key_object in keys:
		if ('TTree' not in key_object.GetClassName()):
			print 'Object ' , key_object, ' is not a TTree object'
		else:
			
