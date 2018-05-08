-   [**MOVords**](#MOVords)
    -   [Installation](#installation)
    -   [Example Usage](#example-usage)
    -   [**MOVpca**](#MOVpca)
    -   [**MOVpnmds**](#MOVpnmds)
    -   [License](#license)

<!-- README.md is generated from README.Rmd. Please edit that file -->
**MOVords**
==========

The package provides six functions: `MOVpca()`, `MOVca()`, `MOVpcoa()`, `MOVnmds()`, `MOVrda()`, `MOVcca()`. You can customize the display labels and themes. Labels can contain formulas. Image themes can be set by themes in ggplots or by other packages(ggthemr).

Usage
------------
Copy these R files into your working dictionary

Example Usage
-------------

``` r
library(MOVords)
require(vegan)
#> Loading required package: vegan
#> Loading required package: permute
#> Loading required package: lattice
#> This is vegan 2.4-4
require(ggplot2)
#> Loading required package: ggplot2
data(Cheminds)
```

**MOVpca**
---------

### 1. get group factor

``` r

Chemin.w <- hclust(dist(scale(Cheminds)), "ward.D")
gr <- cutree(Chemin.w , k=4)
grl <- factor(gr)
```

### 2. Compute PCA

``` r
Chemin.pca <- rda(Cheminds,scale = TRUE)
#head(summary(Chemin.pca))
```

### 3. Produce a plot

``` r
MOVpca(Chemin.pca)
```

<img src="README-unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

### 4. Add a group

``` r
MOVpca(Chemin.pca, group = grl) 
```

<img src="README-unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

### 5. Set a theme

``` r
MOVpca(Chemin.pca, group = grl, spacol = "white") + theme_dark()
```

<img src="README-unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

### 6. Set a theme (use MOVthemr packages, more themes in MOVthemer)

[](https://github.com/cttobin/ggthemr)

``` r
library(MOVthemr)

chalk_theme <- ggthemr('chalk', set_theme = FALSE)
fd_theme <- ggthemr('flat dark', set_theme = FALSE)

p <- MOVpca(Chemin.pca, group = grl, spacol = "white")
p + chalk_theme$theme
```

<img src="README-unnamed-chunk-9-1.png" style="display: block; margin: auto;" />

``` r
p + fd_theme$theme
```

<img src="README-unnamed-chunk-9-2.png" style="display: block; margin: auto;" />

### 7. Remove the arrow

``` r
MOVpca(Chemin.pca, group = grl, spearrow = NULL)
```

<img src="README-unnamed-chunk-10-1.png" style="display: block; margin: auto;" />

### 8. Modify legend title, group color and point shape

``` r
MOVpca(Chemin.pca, group = grl, spearrow = NULL) + 
  scale_color_manual(name = "Groups",values = c("red2", "purple1", "grey20","cyan")) +
  scale_shape_manual(name = "Groups",values = c(8,15,16,17))
```

<img src="README-unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

### 9. Add confidence ellipses

``` r
MOVpca(Chemin.pca, group = grl, spearrow = NULL, ellipse = TRUE) +
  scale_colour_hue(l = 70, c = 300)
```

<img src="README-unnamed-chunk-12-1.png" style="display: block; margin: auto;" />

------------------------------------------------------------------------

**MOVpnmds**
-----------

### 1. Compute NMDS

``` r
Chemin.nmds <- metaMDS(Cheminds, distance="bray")
#> Square root transformation
#> Wisconsin double standardization
#> Run 0 stress 0.04321381 
#> Run 1 stress 0.04321545 
#> ... Procrustes: rmse 0.000332486  max resid 0.001544699 
#> ... Similar to previous best
#> Run 2 stress 0.04321384 
#> ... Procrustes: rmse 0.0005743842  max resid 0.002679384 
#> ... Similar to previous best
#> Run 3 stress 0.06230603 
#> Run 4 stress 0.04321479 
#> ... Procrustes: rmse 0.00074061  max resid 0.003457211 
#> ... Similar to previous best
#> Run 5 stress 0.07681417 
#> Run 6 stress 0.0623054 
#> Run 7 stress 0.06972094 
#> Run 8 stress 0.0432159 
#> ... Procrustes: rmse 0.0004013882  max resid 0.001869575 
#> ... Similar to previous best
#> Run 9 stress 0.04321421 
#> ... Procrustes: rmse 0.0001015489  max resid 0.0004638241 
#> ... Similar to previous best
#> Run 10 stress 0.04321421 
#> ... Procrustes: rmse 0.000647809  max resid 0.003023385 
#> ... Similar to previous best
#> Run 11 stress 0.04321383 
#> ... Procrustes: rmse 5.693312e-05  max resid 0.0001668148 
#> ... Similar to previous best
#> Run 12 stress 0.04321469 
#> ... Procrustes: rmse 0.0007662745  max resid 0.003576096 
#> ... Similar to previous best
#> Run 13 stress 0.04321539 
#> ... Procrustes: rmse 0.0003309146  max resid 0.001540436 
#> ... Similar to previous best
#> Run 14 stress 0.06230623 
#> Run 15 stress 0.04321394 
#> ... Procrustes: rmse 0.0005394657  max resid 0.002445097 
#> ... Similar to previous best
#> Run 16 stress 0.04321387 
#> ... Procrustes: rmse 3.773925e-05  max resid 0.0001624311 
#> ... Similar to previous best
#> Run 17 stress 0.04321547 
#> ... Procrustes: rmse 0.0003366437  max resid 0.001564668 
#> ... Similar to previous best
#> Run 18 stress 0.07548226 
#> Run 19 stress 0.07779553 
#> Run 20 stress 0.2242599 
#> *** Solution reached
#head(summary(Chemin.nmds))
```

### 2. Produce a plot

``` r
MOVnmds(Chemin.nmds)
```

<img src="README-unnamed-chunk-14-1.png" style="display: block; margin: auto;" />

### 3. Add a group

``` r
MOVnmds(Chemin.nmds, group = grl) 
```

<img src="README-unnamed-chunk-15-1.png" style="display: block; margin: auto;" />

### 4. Set a theme

``` r
MOVnmds(Chemin.nmds, group = grl, spacol = "white") + theme_dark()
```

<img src="README-unnamed-chunk-16-1.png" style="display: block; margin: auto;" />

### 5. Set a theme (use ggthemr packages, more themes in ggthemer)

[](https://github.com/cttobin/ggthemr)

``` r

library(ggthemr)

chalk_theme <- ggthemr('chalk', set_theme = FALSE)
fd_theme <- ggthemr('flat dark', set_theme = FALSE)

p <- MOVnmds(Chemin.nmds, group = grl, spacol = "white")
p + chalk_theme$theme
```

<img src="README-unnamed-chunk-17-1.png" style="display: block; margin: auto;" />

``` r
p + fd_theme$theme
```

<img src="README-unnamed-chunk-17-2.png" style="display: block; margin: auto;" />

### 6. Remove the arrow

``` r
MOVnmds(Chemin.nmds, group = grl, spearrow = NULL)
```

<img src="README-unnamed-chunk-18-1.png" style="display: block; margin: auto;" />

### 7. Set labels

``` r
mlabs<-c("NH[4]^{`+`}" , "NO[3]^{`-`}" ,"delta^13*C","A[1]","sqrt(2*pi)","frac(x^2,2)",
         "sin(x)","hat(x)","bar(xy)","90*degree","x^{y+z}")

MOVnmds(Chemin.nmds, group = grl, spearrow = NULL, msplabs = mlabs)
```

<img src="README-unnamed-chunk-19-1.png" style="display: block; margin: auto;" />

### 8. Modify legend title, group color and point shape

``` r
MOVnmds(Chemin.nmds, group = grl, spearrow = NULL) + 
  scale_color_manual(name = "Groups",values = c("red2", "purple1", "grey20","cyan")) +
  scale_shape_manual(name = "Groups",values = c(8,15,16,17))
```

<img src="README-unnamed-chunk-20-1.png" style="display: block; margin: auto;" />

### 9. Add confidence ellipses

``` r
MOVnmds(Chemin.nmds, group = grl, spearrow = NULL, ellipse = TRUE) +
  scale_colour_hue(l = 70, c = 300)
```

<img src="README-unnamed-chunk-21-1.png" style="display: block; margin: auto;" />

License
-------

Released under GPL-3.
