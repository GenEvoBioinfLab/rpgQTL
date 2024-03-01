### Install from scratch
This is a tutorial for people without any experience about conda-based python environment.  
This is mainly for Yale HPC user, but also applicable to general linux system with nvidia GPU.  
The package version is based on those on Yale HPC Farnam cluster.

#### 1. Install miniconda
* Go to https://docs.conda.io/en/latest/miniconda.html#linux-installers, download the Linux installers Python 3.8 file.  
(e.g. [Miniconda3-py38_4.12.0-Linux-x86_64.sh](https://repo.anaconda.com/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh))  
* Install the file by `bash Miniconda3-py38_4.12.0-Linux-x86_64.sh`.
    * Follow the prompts on the installer screens. Just use all default option.
    * On Yale HPC Farnam, with default option, your miniconda core directory will be in `home` and all python package will be in `project`.
* Close and re-open the terminal. Check if conda is successfully installed by `which python`. This should be under your miniconda core directory.
* (Highly recommand) Install Mamba by `conda install mamba -n base -c conda-forge`. The conda package manager will freeze if your conda environment gets too big.

#### 2. Create conda environment and install basic packages
*If you did not install mamba, replace all `mamba` by `conda` in all following parts.*
* `mamba create -n rpgQTL python=3.8`.
   * You can change the environment name `rpgQTL` to anything you want. If you do, remember to change for all the following parts.
* `conda activate rpgQTL`. Confirm the activation by `which python`. This should be under the rpgQTL environment folder.
* `mamba install numpy=1.20.1 scipy=1.7.0 pandas=1.2.3 matplotlib=3.3.4 jupyter ipython rpy2=3.4.2 -c conda-forge`
   * The version of these packages should not matter, but had not been tested.

#### 3. Install pytorch
* Install CUDA and cuDNN (skip if you are on Yale HPC). See [here](https://developer.nvidia.com/cuda-downloads).
* `mamba install cudatoolkit=11.1.1 pytorch=1.8.0 torchaudio=0.8.0 torchvision=0.2.2 -c pytorch -c conda-forge`
  * Unlike the previoius section, the package versions DO matters. You should use the pytorch version that matched to the CUDA version you have on your system. See more details in [Installing Pytorch](https://pytorch.org/get-started/locally/)

#### 4. Install tensorQTL
* `mamba install xarray=0.17.0 -c conda-forge`
   * It seems the newer version of xarray will not work (at least with the numpy version showed above).
* `which pip`. Make sure we are using the pip under rpgQTL environment folder.
* `pip install tensorqtl==1.0.5`
* Compatability with later version of tensorqtl has NOT beed tested.
* See details in [tensorQTL github](https://github.com/broadinstitute/tensorqtl)

#### 5. Install rpgQTL.py
* (Easy way) Simply copy `rpgQTL.py` file into your working directory. 
* Install rpgQTL as package.

#### 6. Running the scripts in general
* To run any scripts in the rpgQTL environment and make used of GPU:
```
source /home/jg2447/miniconda3/etc/profile.d/conda.sh (change to your miniconda core directory)
conda activate tensorQTL
module load CUDA/11.1.1-GCC-10.2.0 (skip if NOT on Yale HPC)
module load cuDNN/8.0.5.39-CUDA-11.1.1 (skip if NOT on Yale HPC)
python your_file.py
```
* For Yale HPC Farnam, add `#SBATCH --gpus-per-node titanv:1` to the sbatch jobfile to indicate using of one titanv gpu.
  * You can change the name titanv to other gpu names (see below).
  * Or you can change it to `#SBATCH --gpus-per-node 1` to require any type of gpu.
* `sinfo -N -O NodeHost:.9,Partition:.19,AllocMem:.11,FreeMem:.11,Memory:.11,CPUsState:.15,Gres:.22,GresUsed:.30 | grep "pi_gerstein_gpu\|HOSTNAMES" | (sed -u 1q; sort -k1,1)`
  * You can use this command to check all pi_gersetin_gpu node current states. 
  * The available gpu names are also shown.
