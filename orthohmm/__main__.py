"""orthohmm.__main__: executed when orthohmm is called as script"""
import sys

from .orthohmm import main

main(sys.argv[1:])
