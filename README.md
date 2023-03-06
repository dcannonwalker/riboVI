# riboVI

Variational inference for detecting differential translation in ribosome profiling studies

### Authors: 
> David Walker and Peng Liu

### Contact:
> [dcwalker@iastate.edu] (David Walker)

### Citation: 
> Walker, D.C., Lozier, Z.R., Bi, R., Kanodia, P., Miller, W.A., Liu, P., 2023.
Variational inference for detecting differential translation in ribosome profiling studies. Submitted. 

### Installation
Using the `remotes` package: 
```
remotes::install_github("dcannonwalker/riboVI")
```
### Getting started

`ribovi()` contains the core functionality of the package. Below is a short demonstration using the `ribo_example` data set: 

```r
library(riboVI)
data(ribo_example)
head(ribo_example)

# Console out:
  geneid RiboCtrl1 RiboCtrl2 RiboTrt1 RiboTrt2 RNACtrl1 RNACtrl2 RNATrt1 RNATrt2
1  gene1        97       122      146      171      103      160     231     137
2  gene2         3         5        2        4        4        3       2       2
3  gene3       309      1025      708      598      118      497     359     295
4  gene4      2976      7218     7923     5408     3304     8327    8954    5303
5  gene5       149       155      141      200      207      214      85     180
6  gene6        97        56      151       83      114       59     156     122
```

`ribovi()` requires as arguments:
* `counts`, data frame with a gene ID column and columns of counts for each sample
* `X`, the fixed effects design matrix (with intercept)
* `Z`, the sample effects design matrix (with no intercept)
* Other arguments are optional, and we'll leave them set to their default values in this example

```r
treatment <- factor(rep(rep(c('ctrl', 'trt'), each = 2), 2))
preparation <- factor(rep(c('ribo', 'rna'), each = 4))
X <- model.matrix(~treatment + preparation + treatment:preparation)

sample <- rep(1:4, 2)
Z <- model.matrix(~0 + sample)

ribovi_out <- ribovi(counts = ribo_example, X = X, Z = Z)
```



