#!/bin/sh
g++ L1JetEmulator.C $(root-config --cflags --libs) -Werror -Wall -Wextra -O2 -o L1JetEmulator.exe || exit 1

./L1JetEmulator.exe 0
./L1JetEmulator.exe 1

