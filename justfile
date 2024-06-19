#!/usr/bin/env -S just --justfile
# just reference  : https://just.systems/man/en/

set quiet := true

alias default := run

list:
    just --list

# Set up the environment with Rye
[unix]
build:
    mise run --timings setup:unix

[windows]
build:
    mise run --timings setup:windows

# Run all plotting scripts
run: && copy
    mise run --timings run

[private]
copy:
    mise run --timings copy
