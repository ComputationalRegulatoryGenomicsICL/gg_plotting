# ggplot-based plotting for genomic data

*Warning: very unstable, function names and functionality subject to change*

This code is based on the idea of applying the 'grammar of graphics' to genomic 
data. Currently, it takes as basic input a ChIPprofile object from the [soGGi](http://www.bioconductor.org/packages/release/bioc/html/soGGi.html) 
package, and provides functions to transform this data into a long form dataframe
suitable for plotting with ggplot (with and without data summarisation), as well 
as some functions for simple plots. Plots should be easily modified with the usual 
ggplot syntax.

The basic workflow is:

ChIPprofile object with one or more assays 
  --> reshape with `reshape_chipprofile()`
          --> plot a heatmap with `gg_heatmap()`
      --> summarise signal by position and optionally other metadata columns using `summarise_signal()`
          --> plot a graph of average signal using ggplot!
          
Unfortunately it's not possible to make both these types of graphs from the same 
data.frame... You could however approximate the second type of graph using 
`geom_smooth()` on the un-summarised data.frame.

I've implemented a few ways of summarising data by position, but you can also 
summarise data in any way you like using dplyr functions.