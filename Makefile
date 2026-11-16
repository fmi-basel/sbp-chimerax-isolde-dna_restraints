# Targets:
#   wheel (default)        build a wheel in dist/
#   install / app-install  build and install into ChimeraX
#   uninstall              remove the installed bundle
#   test                   run the headless smoke test (install first)
#   test-pure              run the pure-Python tests (no ChimeraX needed)
#   run                    launch ChimeraX
#   clean                  remove build artifacts
#
# Uses the stable ChimeraX by default. "make DAILY=1 <target>" uses the daily
# build, and "make CHIMERAX_EXE=/path/to/ChimeraX <target>" uses any other.
# On Windows use make_win.bat.

BUNDLE_NAME = ChimeraX-ISOLDE-DNA_RESTRAINTS
BUNDLE_BASE_NAME = $(subst ChimeraX-,,$(BUNDLE_NAME))

OS = $(shell uname -s)
ifeq ($(OS),Darwin)
  ifdef DAILY
    CHIMERAX_EXE ?= /Applications/ChimeraX_Daily.app/Contents/bin/ChimeraX
  else
    CHIMERAX_EXE ?= /Applications/ChimeraX.app/Contents/bin/ChimeraX
  endif
else
  ifdef DAILY
    CHIMERAX_EXE ?= chimerax-daily
  else
    CHIMERAX_EXE ?= chimerax
  endif
endif

.DEFAULT_GOAL := wheel

wheel:
	$(CHIMERAX_EXE) --nogui --safemode --exit --cmd "devel build ."

install app-install:
	$(CHIMERAX_EXE) --nogui --safemode --exit --cmd "devel install ."

uninstall app-uninstall:
	$(CHIMERAX_EXE) --nogui --safemode --exit --cmd "toolshed uninstall $(BUNDLE_BASE_NAME)"

test:
	$(CHIMERAX_EXE) --nogui --exit --script tests/test_commands.py

# Pure-Python tests (no ChimeraX needed)
test-pure:
	python3 tests/test_pure.py

run:
	$(CHIMERAX_EXE)

clean:
	$(CHIMERAX_EXE) --nogui --safemode --exit --cmd "devel clean ."

.PHONY: wheel install app-install uninstall app-uninstall test test-pure run clean
