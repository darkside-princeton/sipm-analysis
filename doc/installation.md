# Installation
This python package is used for the SiPM data analysis of the DarkSide Princeton group and can be easily install with pip from within your conda environment. 

We need to install a new version Miniconda via the following two commands

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Make sure to run the `conda init` command to populate the necessary commands to your `.bashrc` so that upon a new login you will be able to use `conda` right away. 

After installation let's make sure we have the most up-to-date version
```
conda update conda
```

Now we need to create a new conda environement via 
```
conda env create -f conda.yml
```

and then you can activate it via 
```
conda activate sipm
```

Before we get the actual analysis code we can make sure that the conda environment also shows up as an available Jupyter kernel by running
```
ipython kernel install --name "sipm" --user
```

Next, we need to download the code repository
```
git clone https://github.com/darkside-princeton/sipm-analysis.git
```

Lets change to our new directory
```
cd sipm-analysis/
```

and then finally run
```
pip install -e .
```

Now you should be able to run one of the example Jupyter notebooks to test that your installation worked.