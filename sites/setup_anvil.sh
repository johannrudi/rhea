# load system defaults
module reset

# load modules
module load gcc/11.2.0
module load openblas/0.3.17
module load boost/1.74.0
module load valgrind/3.15.0
module load python/3.9.5

# load rhea-kit
module load anaconda/2021.05-py38
module use /home/x-johann/privatemodules
module load conda-env/rhea-kit-py3.8.8
