#!/bin/bash
#######################################
# Modeller
# Add the path of the required parameter file or script to the environment variable
SOFTWARE_DIR=$PWD
echo "#################### Modeller ########################" >> ~/.bashrc
echo "export PATH=\$PATH:$SOFTWARE_DIR/scripts" >> ~/.bashrc
echo "export modeller_prog=$SOFTWARE_DIR/program" >> ~/.bashrc
echo "source ~/.bashrc"
