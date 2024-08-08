#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Convenience wrapper for running orthohmm directly from source tree."""
import sys

from orthohmm.orthohmm import main

if __name__ == "__main__":
    main(sys.argv[1:])
