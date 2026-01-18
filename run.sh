#!/bin/bash
clang++ Src/*.cpp -o myapp -limgui -lSDL3_image -lSDL3 -lGLEW -lGL && SDL_VIDEODRIVER=x11 ./myapp 
