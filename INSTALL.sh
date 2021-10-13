#!/bin/bash
#######################################
# Modeller								  #
# 2019-4-18 by LY					  #
#######################################

#将所需参数文件或脚本的路径加入环境变量
SOFTWARE_DIR=$PWD
echo "#################### Modeller ########################" >> ~/.bashrc
echo "export PATH=\$PATH:$SOFTWARE_DIR/scripts" >> ~/.bashrc
echo "export modeller_prog=$SOFTWARE_DIR/program" >> ~/.bashrc
echo "source ~/.bashrc"
