#this project uses trained data of H2O from deepmd documentation
#to run mc simulations

##########################################
To manually perform each step of computations
1. python launch_one_run.py ./path/to/mc.conf
2. make run_mc
3. ./run_mc ./path/to/cppIn.txt
4. python check_after_one_run.py ./path/to/mc.conf lastFileNum
5. go to 1, until no more data points are needed

#########################################