#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
cd "$DIR/src"
g++ -o ../euler_fluid_simulation main.cpp Fluid.cpp Field.cpp -lsfml-graphics -lsfml-window -lsfml-system
../euler_fluid_simulation