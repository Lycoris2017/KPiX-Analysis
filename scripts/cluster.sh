#!/bin/bash

pedestal_extension=.tree_pedestal.root

./bin/pedestal_tree $1 $2
./bin/cluster_analysis $1 $2 $1$pedestal_extension
