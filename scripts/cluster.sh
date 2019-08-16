#!/bin/bash

pedestal_extension=.tree_pedestal.root

./bin/pedestal_tree $1 $2
./bin/external_trig_debug $1 $2 $1$pedestal_extension
