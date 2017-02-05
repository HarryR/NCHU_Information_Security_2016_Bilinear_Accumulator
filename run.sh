#!/bin/bash

if [ $# -gt 1 ]; then
  gcc -I$HOME/.local/include/pbc -L$HOME/.local/lib -Wl,-rpath $HOME/.local/lib -I./pbc -lpbc -lgmp -std=c99 -o $1 $1.c && clear && ./$1 < ./pbc/param/a.param
  exit 0
fi

gcc -I$HOME/.local/include/pbc -L$HOME/.local/lib -Wl,-rpath $HOME/.local/lib -I./pbc -lpbc -lgmp -std=c99 -o bilinear_acc bilinear_acc.c && clear && ./bilinear_acc < ./pbc/param/a.param
