### Install from scratch

#### 1. Install miniconda
* Go to https://docs.conda.io/en/latest/miniconda.html#linux-installers, download the Linux installers Python 3.8 file.  
(e.g. [Miniconda3-py38_4.12.0-Linux-x86_64.sh](https://repo.anaconda.com/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh))  
* Install the file by `bash Miniconda3-py38_4.12.0-Linux-x86_64.sh`.
* Follow the prompts on the installer screens.
* Close and re-open the terminal. By default conda auto activate the base environment. Check if conda is successfully installed by `which python`. This should be under your miniconda core directory.
* (Highly recommand) Install Mamba by `conda install mamba -n base -c conda-forge`. The conda package manager will freeze if your conda environment gets too big.

#### 2. Create conda environment and install basic packages
*If you did not install mamba, replace all `mamba` by `conda` in all following parts.*
* `mamba create -n rpgQTL python=3.8`.
   * You can change the environment name `rpgQTL` to anything you want. If you do, remember to change for all the following parts.
* `conda activate rpgQTL`. Confirm the activation by `which python`. This should be under the rpgQTL environment folder.
* `mamba install numpy=1.20.1 scipy=1.7.0 pandas=1.2.3 matplotlib=3.3.4 jupyter ipython rpy2=3.4.2 -c conda-forge`
   * The version of these packages should not matter, but had not been tested.

#### 3. Install pytorch
* Make sure you have CUDA and cuDNN in your system. See [here](https://developer.nvidia.com/cuda-downloads).
* Install pytorch. See more details in [Installing Pytorch](https://pytorch.org/get-started/locally/). Unlike the previoius section, the package versions DO matter. You should use the pytorch version that matched to the CUDA version you have on your system.
   * Example: `mamba install cudatoolkit=11.1.1 pytorch=1.8.0 torchaudio=0.8.0 torchvision=0.2.2 -c pytorch -c conda-forge`

#### 4. Install tensorQTL
* `mamba install xarray=0.17.0 -c conda-forge`
   * It seems the newer version of xarray will not work (at least with the numpy version showed above, probably fixed with newer ones).
* `which pip`. Make sure we are using the pip under rpgQTL environment folder.
* `pip install tensorqtl==1.0.5`
* Compatability with later version of tensorqtl has NOT beed tested.
* See details in [tensorQTL github](https://github.com/broadinstitute/tensorqtl)

#### 5. Install rpgQTL.py
* Easy way: simply copy `rpgQTL.py` file into your working directory. 
* Alternatively, you can install rpgQTL as a package.
