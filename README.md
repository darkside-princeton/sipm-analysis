# Installation
This python package is used for the SiPM data analysis of the DarkSide Princeton group and can be easily install with pip from within your conda environment. 

First create a new conda environement via 

```
conda env create -f conda.yml
```
and then you can activate it via 
```
conda activate ds-pu
```
Before we get the actual analysis code we can make sure that the conda environment also shows up as an available Jupyter kernel by running
```
ipython kernel install --name "ds-pu" --user
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

# Problems 
For any question reach out to Ako Jamil ([ako.jamil@princeton.edu](mailto:ako.jamil@princeton.edu)).