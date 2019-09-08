
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mincountr

The goal of mincountr is to provide a simple point-counting mechanism
for microscope thin section images of rocks

## Installation

You can install mincountr from github with:

``` r
# install.packages("devtools")
devtools::install_github("Ignimbrit/mincountr")
```

## Example

Say you have collected an image from a mineral-bearing thin section
e.g. in an electron microprobe session and wonder about the relative
share of the different phases in your sample. Your image might look
something like this:

``` r
library(mincountr)

myimage <- mcr_load_image(
  system.file("extdata", "testim.png", package = "mincountr")
  )

plot(myimage)
```

![](README-load%20image-1.png)<!-- -->

As you can see, the different levels of brightness in the image allows
the observer to distinguish between several distinct phases. There is a
large, very bright mineral, some lightgrey minerals, a darkgrey, glassy
matrix, some black holes and so on. With `mincountr` we are able to
translate this qualitative optical assessment into practical numbers.
First let’s have a look at the density-distribution (like a continuous
histogram) of the images brightness.

``` r
mcr_inspect_phases(myimage)
```

![](README-plot%20brightness-1.png)<!-- -->

Here we can see that the four phases we distinguished in the image above
show up in the density distribution as distinct peaks of specific
brightness values. We can use this to assign a brightness-range to
certain phases. The peak on the far left (most dark) corresponds likely
to the hole we’ve seen in our thin section image. The peak’s value range
lays between value \~0-0.05. Other peaks (from left to right) range from
0.3-0.45 (glassy matrix?), from 0.5-0.65 (light gray minerals?) and from
0.92-1 (bright minerals?). Let me illustrate what I mean:

``` r
mcr_inspect_phases(myimage) +
  ggplot2::geom_vline(
    xintercept = c(0, 0.05, 0.3, 0.45, 0.5, 0.65, 0.92, 1),
    color = "red"
      )
```

![](README-illustrate%20peakborders-1.png)<!-- -->

Now in this example we just chose the borders of the peak by hand. This
is probably the safest method of constraining your brightness-levels, as
it allows you to chip in your personal mineralogical expertise. However,
sometimes you have a large stack of images you want to work with and
probably not the time to constrain peak-ranges by hand every single
time. This is why mincountr comes with an automatic mechanism to
generate those numbers for you.

``` r
myconstrains <- mcr_autoconstrain(myimage)
#> Warning in min(.): no non-missing arguments to min; returning Inf

print(myconstrains)
#> # A tibble: 4 x 4
#>      x1 peakpos     x2    ID
#>   <dbl>   <dbl>  <dbl> <int>
#> 1 0       0     0.0283     1
#> 2 0.338   0.382 0.437      2
#> 3 0.501   0.559 0.612      3
#> 4 0.963   0.989 1          4
```

The `mcr_autoconstrain` function automatically detects peaks and notes
their position (`peakpos`) and then goes on and calculates their borders
both on the left-hand-side (`x1`) and on the right hand side (`x2`).
Under the hood, `mcr_autoconstrain` identifies turning points in the
brightness-“spectra”, cuts it into pieces, one piece per peak, and then
loops over the single-peak spectra-pieces to calculate the
half-height-width.

Now we have all the information we need to assign certain areas of our
original image to distinct phases. The mincountr-package comes with a
function that lets you inspect what this assignement looks like.

``` r
mcr_inspect_assignement(
  myimage,
  lhs = myconstrains$x1,
  rhs = myconstrains$x2
)
```

![](README-check%20assignement-1.png)<!-- -->

As every pixel in the original image was now assigned to one of the 7
levels shown in the picture above, we can go ahead and just count the
pixels and then calculate the relative share of each group.

``` r
myresult <- mcr_herd_minerals(
  myimage,
  lhs = myconstrains$x1,
  rhs = myconstrains$x2
)

print(myresult)
#> # A tibble: 7 x 3
#>   Phase_ID pixels proportion_percentage
#>      <int>  <int>                 <dbl>
#> 1        1  34862                  9.88
#> 2        2  23583                  6.68
#> 3        3  97086                 27.5 
#> 4        4  18786                  5.32
#> 5        5  70327                 19.9 
#> 6        6  23485                  6.66
#> 7        7  84749                 24.0
```
