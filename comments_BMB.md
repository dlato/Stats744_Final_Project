## BMB comments

* Not sure I agree with your discussion of positive/negative/purifying selection. I thought positive selection referred to loci that are under selection for *variation* (i.e., novel variants are likely to be selected), while negative (=purifying) selection refers to loci where novel variancts are selected against.  You say that genes that are under negative selection are likely to be lost from the genome; that doesn't make sense to me. Is that what you actually meant?

Figure 1 (vio_str_box):

* y-axis-labels are actually too large (the '100' got omitted), and/or too closely spaced (although it's admittedly an awkward range; you could use {0.001,1,1000} [too sparse?] or {0.01,1,100} [not wide enough?], but the intermediate ticks aren't pretty)

* think I said this before: do you really need the right axis? Both scales are unitless (I think: "numbers per site" doesn't *really* have units?). Right-hand scale seems distracting/unnecessary.

* colours aren't distinguishable in grayscale (doesn't really matter though, as they're redundant/decorative here, unless you go on to use them in non-redundant way in later figures)

Figure 2

not bad to colour values with ratio>1, although this does contribute a little bit to a possibly unnecessary dichotomization (i.e. two points could be very close, one with omega=0.999 and the other with omega=1.001, and they'd be drawn in different colours)

maybe make trend lines slightly darker than the points (or use alpha<1 for points), to make them easier to distinguish?

Don't overuse "interestingly"

Figures 3-x

nice looking. What about the idea of binomial CIs on the means of bins?
Could these be condensed/combined?  Hard to compare across bacteria.

Agree that you need to do something to compare genomes across different scales (since they differ so much in size).  Might be a better way to do this (but can't think of it at the moment).

Not sure that blue/purple works as a pair, but whatever.

It's nice that ggplot provides a convenient way to show zero values on a log scale (i.e., stuck at the bottom axis). Might be necessary to emphasize this convention more for a general readership (unfortunately there's not a super-easy way to add an axis break in ggplot)

Good thoughts about ordering

Reasonable argumenat bout multiple scales even though I don't agree

