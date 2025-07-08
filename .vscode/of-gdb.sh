#!/bin/bash
. /home/newuser/OpenFOAM/OpenFOAM-9/etc/bashrc WM_NCOMPROCS=2; export WM_COMPILE_OPTION=Debug
/usr/bin/gdb "$@"