### Install from scratch

#### 1. Install miniconda
* Go to [Conda](https://docs.conda.io/en/latest/miniconda.html#linux-installers), download the Linux installers Python 3.8 file.  
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
* `mamba install numpy scipy pandas matplotlib jupyter xarray rpy2 r-essentials -c conda-forge`
   * The version of these packages should not matter, but had not been tested.

#### 3. Install pytorch
* Make sure you have CUDA in your system.
* Install pytorch. See more details in [Installing Pytorch](https://pytorch.org/get-started/locally/). You should use the pytorch version that matched to the CUDA version you have on your system.

#### 4. Install tensorQTL
* `which pip`. Make sure we are using the pip under rpgQTL environment folder.
* `pip install tensorqtl==1.0.5`
* Compatability with later version of tensorqtl has NOT beed tested.
* See details in [tensorQTL github](https://github.com/broadinstitute/tensorqtl)

#### 5. Install R packages
* `which R`. Make sure we are using the R under rpgQTL environment folder.
* In R:
 ```
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("qvalue")
  ```

#### 5. Install rpgQTL
* `which python`. Make sure we are using the python under rpgQTL environment folder.
* Install rpgQTL:
```
git clone https://github.com/GenEvoBioinfLab/rpgQTL
cd rpgQTL
python3 setup.py install
```
