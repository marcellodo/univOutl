# univOutl 0.5.0

## New features
* Added Walker et al. (2018) modified boxplot to the `boxB()` function.


# univOutl 0.4

## Improvements
* Updated the `plot4sizes()` function and its output.
* Updated compatibility with `robustbase`.
* Updated the `plot4sizes()` function and its output.
* Updated code depending on the `robustbase` package. Starting from `robustbase` v0.95-0 (2022-04-02), `mc()` now uses `doScale = FALSE` by default to guarantee convergence in extreme cases when used with `c.huberize = 1e11`. These changes are not backward compatible.

# univOutl 0.3

## Improvements
* Modified `boxB()` and `LocScaleB()` to return outliers by tail of the distribution (thanks to Eliot McGinnis for suggesting the improvement and providing modified code).

# univOutl 0.2

## New features
* Added the `skew.misc()` function to compute skewness measures.
* Added the possibility in `HBmethod()` to detect outliers in Escores using a boxplot adjusted for skewness.
* Added a new method in `LocScaleB()` to account for slightly skewed distributions.

# univOutl 0.1-5

## Improvements
* Modified `HBmethod()` according to a note from Hidiroglou and Emond (2018).
* Modified the data frame returned by `ratioSize()` when `return.data = TRUE`: the column `"size"` is now named `"sizeU"` for consistency with `HBmethod()`.
* Added the `plot4sizes()` function to plot output from `HBmethod()` or `ratioSize()`.

# univOutl 0.1-4

## Bug fixes
* Corrected a bug in `LocScaleB()` when computing the score.
* Added the possibility in `LocScaleB()` to estimate the scale using the Inter-Decile Range (`method = "IDR"`) or based on Gini's mean difference.