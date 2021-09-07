Python Wrapper around SIXTE: The SIXTEsoft package
==================================================
v1.1

This package is a wrapper around the tools of the SIXTE software
package. You can add the python wrapper to your packaging manager by
executing (in this folder) `pip3 install -e .` and then load the
package in Python3 via `import sixtesoft`.


# CHANGELOG

- 2021-09-07: v1.1
  Change test default to boolean to ensure that
  SIXTE_USE_PSEUDO_RNG is not initialized. Add possibility
  to use the python wrapper as module.


# How to add new functions

* Go into the sixtesoft/ folder
* Add new functions in the sixte.py script, or create new scripts
* Add functions and variables you want the user to be able to access
  to the __init__.py
* Increase the version number in setup.py and add a changelog entry in
  this README.md file
* The functions can be accessed from within python3 by `import sixtesoft`

