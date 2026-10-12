CXX ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Wpedantic
CPPFLAGS += -Iinclude
PYTHON ?= python3

.PHONY: all demo checks reproduce paper paper-affine papers
all: build/demo build/test_matching build/matching_experiments

build:
	mkdir -p build

build/demo: examples/demo.cpp include/recursive_tree_matching.hpp | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

build/test_matching: tests/test_matching.cpp include/recursive_tree_matching.hpp | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

build/matching_experiments: verification/matching_experiments.cpp include/recursive_tree_matching.hpp | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

demo: build/demo
	./build/demo

checks: build/test_matching
	./build/test_matching
	$(PYTHON) verification/verify_extensions.py

reproduce:
	$(PYTHON) reproduce.py

paper:
	cd paper && latexmk -pdf -interaction=nonstopmode -halt-on-error recursive_matching.tex

paper-affine:
	cd paper_affine && latexmk -pdf -interaction=nonstopmode -halt-on-error affine_matching.tex

papers: paper paper-affine
