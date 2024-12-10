#!/usr/bin/env -S just --justfile
# just reference  : https://just.systems/man/en/

set quiet := true

alias default := run

list:
    just --list

# Set up the environment
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

sync:
    mise run --timings sync

plot: sync
    mise run --timings run
