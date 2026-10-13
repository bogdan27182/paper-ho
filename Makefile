PYTHON ?= python3

.PHONY: all check

all: check

check:
	$(PYTHON) reproduce.py
