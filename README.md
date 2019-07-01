Citation
--------
Please cite one or more of the following relevant papers if you use this code.

Code for Square Root Graphical Model (SQR) is based on:

David I. Inouye, Pradeep Ravikumar, Inderjit S. Dhillon.
Square Root Graphical Models: Multivariate Generalizations of Univariate Exponential Families that Permit Positive Dependencies
International Conference on Machine Learning (ICML), 2016.
https://www.davidinouye.com/publication/inouye-2016-square/inouye-2016-square.pdf

Code and data are also provided for the following review paper.

David I. Inouye, Eunho Yang, Genevera I. Allen, Pradeep Ravikumar.  
A review of multivariate distributions for count data derived from the Poisson distribution.   
*Wiley Interdisciplinary Reviews (WIREs): Computational Statistics*, 9:3, 2017. doi: 10.1002/wics.1398   
arXiv preprint: https://arxiv.org/pdf/1609.00066.pdf

The implementation of Square Root Graphical Models for the Poisson distribution is based 
on the following arXiv paper:

David I. Inouye, Pradeep Ravikumar, Inderjit S. Dhillon
Generalized Root Models: Beyond Pairwise Graphical Models for Univariate Exponential Families
*arXiv preprint arXiv:1606.00813*, 2016.  
arXiv preprint: https://arxiv.org/pdf/1606.00813.pdf

Installation
------------
You must install the R packages `VineCopula` and `XMRF` for the vine copula and TPGM models to work respectively.

Demo
----
The main demo file is `demo_comparison.m` but the `demo_comparison_check.m` file checks that all the methods run to completion for a really small dataset.

Data
----
The 6 datasets used in the paper are included as simple MAT files in the data folder.
