[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7892408.svg)](https://doi.org/10.5281/zenodo.7892408)
# BRaVDA - a solar wind data assimilation scheme


## Introduction

This repository provides a python implementation of the BRaVDA solar wind data assimilation scheme described in:

Lang, M., & Owens, M. J. (2019). A variational approach to data assimilation in the solar wind. Space Weather, 17, 59–83. https://doi.org/10.1029/2018SW001857

It combines the output of a coronal model, such as WSA or MAS, with the in situ observations, typically near 1 AU. It returns the optimum reconstruction of the solar wind, acounting for errors in both models and observations.

## Installation
 `BRaVDA` is written in Python 3.9.13 and has a range of dependencies, which are listed in `bravda_env.yml` files. Because of these dependencies, the simplest way to work with `BRaVDA` in `conda` is to create its own environment. With the anaconda prompt, in the root directory of `BRaVDA`, this can be done as:
```
>>conda env create -f bravda_env.yml
>>conda activate bravda_env
``` 

## Contact
Please contact either [Matthew Lang](https://github.com/SOJC6) or [Mathew Owens](https://github.com/mathewjowens). 

## Citations

If you use BRaVDA in a publication or presentation, please cite the software using the Zenodo reference with DOI:[10.5281/zenodo.7892408](https://doi.org/10.5281/zenodo.7892408) 

To cite this project, including the scientific basis and functionality of BRaVDA, please use: 

Lang, M., & Owens, M. J. (2019). A variational approach to data assimilation in the solar wind. Space Weather, 17, 59–83. https://doi.org/10.1029/2018SW001857

and

Lang, M., Witherington, J., Turner, H., Owens, M. J., & Riley, P. (2021). Improving solar wind forecasting using data assimilation. Space Weather, 19, e2020SW002698. https://doi.org/10.1029/2020SW002698
