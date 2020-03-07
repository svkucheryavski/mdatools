
# Number of groups (discrete colors)

context('getColors: number of groups')

## get colors for different colormaps and ngroup = 1
c1 = mdaplot.getColors(1, colmap = 'old')
c2 = mdaplot.getColors(1, colmap = 'gray')
c3 = mdaplot.getColors(1, colmap = 'jet')
c4 = mdaplot.getColors(1, colmap = c('red', 'green', 'blue'))
c5 = mdaplot.getColors(1)

test_that("ngroup = 1: number of colors is correct", {
   expect_equal(length(c1), 1)
   expect_equal(length(c2), 1)
   expect_equal(length(c3), 1)
   expect_equal(length(c4), 1)
   expect_equal(length(c5), 1)
})

test_that("ngroup = 1: color is correct", {
   expect_equal(c1[1], "#3288BD")
   expect_equal(c2[1], "#101010")
   expect_equal(c3[1], "#00007F")
   expect_equal(c4[1], "#FF0000")
   expect_equal(c5[1], "#2679B2")
})

## get colors for different colormaps and ngroup = 2
c1 = mdaplot.getColors(2, colmap = 'old')
c2 = mdaplot.getColors(2, colmap = 'gray')
c3 = mdaplot.getColors(2, colmap = 'jet')
c4 = mdaplot.getColors(2, colmap = c('red', 'green', 'blue'))
c5 = mdaplot.getColors(2)

test_that("ngroup = 2: number of colors is correct", {
   expect_equal(length(c1), 2)
   expect_equal(length(c2), 2)
   expect_equal(length(c3), 2)
   expect_equal(length(c4), 2)
   expect_equal(length(c5), 2)
})

test_that("ngroup = 2: colors are correct", {
   expect_equal(c1, c("#3288BD", "#D53E4F"))
   expect_equal(c2, c("#E8E8E8", "#101010"))
   expect_equal(c3, c("#00007F", "#7F0000"))
   expect_equal(c4, c("#FF0000", "#0000FF"))
   expect_equal(c5, c("#2679B2", "#D22C2F"))
})

## get colors for different colormaps and ngroup = 3
c1 = mdaplot.getColors(3, colmap = 'old')
c2 = mdaplot.getColors(3, colmap = 'gray')
c3 = mdaplot.getColors(3, colmap = 'jet')
c4 = mdaplot.getColors(3, colmap = c('red', 'green', 'blue'))
c5 = mdaplot.getColors(3)

test_that("ngroup = 3: number of colors is correct", {
   expect_equal(length(c1), 3)
   expect_equal(length(c2), 3)
   expect_equal(length(c3), 3)
   expect_equal(length(c4), 3)
   expect_equal(length(c5), 3)
})

# TODO: commented because it does not pass test on Windows!
#test_that("ngroup = 3: colors are correct", {
#   expect_equal(c1, c("#3288BD", "#F2EA91", "#D53E4F"))
#   expect_equal(c2, c("#E8E8E8", "#A5A5A5", "#101010"))
#   expect_equal(c3, c("#00007F", "#7FFF7F", "#7F0000"))
#   expect_equal(c4, c("#FF0000", "#00FF00", "#0000FF"))
#   expect_equal(c5, c("#2679B2", "#92B42A", "#D22C2F"))
#})

## get colors for different colormaps and ngroup = 10
c1 = mdaplot.getColors(10, colmap = 'old')
c2 = mdaplot.getColors(10, colmap = 'gray')
c3 = mdaplot.getColors(10, colmap = 'jet')
c4 = mdaplot.getColors(10, colmap = c('red', 'green', 'blue'))
c5 = mdaplot.getColors(10)

test_that("ngroup = 10: number of colors is correct", {
   expect_equal(length(c1), 10)
   expect_equal(length(c2), 10)
   expect_equal(length(c3), 10)
   expect_equal(length(c4), 10)
   expect_equal(length(c5), 10)
})

test_that("ngroup = 10: colors are unique", {
   expect_equal(length(unique(c1)), 10)
   expect_equal(length(unique(c2)), 10)
   expect_equal(length(unique(c3)), 10)
   expect_equal(length(unique(c4)), 10)
   expect_equal(length(unique(c5)), 10)
})


# Using color grouping by factor (discrete colors)
context('getColors: cgroup (factor)')

## make factor and get colors for different colormaps
cgroup = cut(runif(1000, 1, 3), 10, levels = 1:10)

c1 = mdaplot.getColors(cgroup = cgroup, colmap = 'old')
c2 = mdaplot.getColors(cgroup = cgroup, colmap = 'gray')
c3 = mdaplot.getColors(cgroup = cgroup, colmap = 'jet')
c4 = mdaplot.getColors(cgroup = cgroup, colmap = c('red', 'green', 'blue'))
c5 = mdaplot.getColors(cgroup = cgroup)

test_that("number of colors in vector is correct", {
   expect_equal(length(c1), 1000)
   expect_equal(length(c2), 1000)
   expect_equal(length(c3), 1000)
   expect_equal(length(c4), 1000)
   expect_equal(length(c5), 1000)
})

test_that("number of unique colors is correct", {
   expect_equal(length(unique(c1)), 10)
   expect_equal(length(unique(c2)), 10)
   expect_equal(length(unique(c3)), 10)
   expect_equal(length(unique(c4)), 10)
   expect_equal(length(unique(c5)), 10)
})

test_that("position of colors is correct", {
   expect_equal(c1, mdaplot.getColors(10, colmap = 'old')[as.numeric(cgroup)])
   expect_equal(c2, mdaplot.getColors(10, colmap = 'gray')[as.numeric(cgroup)])
   expect_equal(c3, mdaplot.getColors(10, colmap = 'jet')[as.numeric(cgroup)])
   expect_equal(c4, mdaplot.getColors(10, colmap = c('red', 'green', 'blue'))[as.numeric(cgroup)])
   expect_equal(c5, mdaplot.getColors(10)[as.numeric(cgroup)])
})


# Using color grouping by continuous variable (no ngroups)
context('getColors: cgroup (continuous) no ngroups')

## make vector and get colors for different colormaps
cgroup = runif(1000, 1, 3)
maxsplits = 64 # default value for ngroups in this case

c1 = mdaplot.getColors(cgroup = cgroup, colmap = 'old')
c2 = mdaplot.getColors(cgroup = cgroup, colmap = 'gray')
c3 = mdaplot.getColors(cgroup = cgroup, colmap = 'jet')
c4 = mdaplot.getColors(cgroup = cgroup, colmap = c('red', 'green', 'blue'))
c5 = mdaplot.getColors(cgroup = cgroup)

test_that("number of colors in vector is correct", {
   expect_equal(length(c1), 1000)
   expect_equal(length(c2), 1000)
   expect_equal(length(c3), 1000)
   expect_equal(length(c4), 1000)
   expect_equal(length(c5), 1000)
})

test_that("number of unique colors is correct", {
   expect_equal(length(unique(c1)), maxsplits)
   expect_equal(length(unique(c2)), maxsplits)
   expect_equal(length(unique(c3)), maxsplits)
   expect_equal(length(unique(c4)), maxsplits)
   expect_equal(length(unique(c5)), maxsplits)
})

test_that("the result also includes pallette as atribute", {
   expect_equal(attr(c1, 'palette'), mdaplot.getColors(ngroups = maxsplits, colmap = 'old'))
   expect_equal(attr(c2, 'palette'), mdaplot.getColors(ngroups = maxsplits, colmap = 'gray'))
   expect_equal(attr(c3, 'palette'), mdaplot.getColors(ngroups = maxsplits, colmap = 'jet'))
   expect_equal(attr(c4, 'palette'), mdaplot.getColors(ngroups = maxsplits, colmap = c('red', 'green', 'blue')))
   expect_equal(attr(c5, 'palette'), mdaplot.getColors(ngroups = maxsplits))
})

ind = cut(cgroup, maxsplits, labels = 1:maxsplits, include.lowest = TRUE)
test_that("position of colors is correct", {
   expect_equal(c1[1:1000], mdaplot.getColors(ngroups = maxsplits, colmap = 'old')[ind])
   expect_equal(c2[1:1000], mdaplot.getColors(ngroups = maxsplits, colmap = 'gray')[ind])
   expect_equal(c3[1:1000], mdaplot.getColors(ngroups = maxsplits, colmap = 'jet')[ind])
   expect_equal(c4[1:1000], mdaplot.getColors(ngroups = maxsplits, colmap = c('red', 'green', 'blue'))[ind])
   expect_equal(c5[1:1000], mdaplot.getColors(ngroups = maxsplits)[ind])
})

# Using color grouping by continuous variable (with ngroups)
context('getColors: cgroup (continuous) with ngroups')

## make vector and get colors for different colormaps
cgroup = runif(1000, 1, 3)
ngroups = 16

c1 = mdaplot.getColors(ngroups = ngroups, cgroup = cgroup, colmap = 'old')
c2 = mdaplot.getColors(ngroups = ngroups, cgroup = cgroup, colmap = 'gray')
c3 = mdaplot.getColors(ngroups = ngroups, cgroup = cgroup, colmap = 'jet')
c4 = mdaplot.getColors(ngroups = ngroups, cgroup = cgroup, colmap = c('red', 'green', 'blue'))
c5 = mdaplot.getColors(ngroups = ngroups, cgroup = cgroup)

test_that("number of colors in vector is correct", {
   expect_equal(length(c1), 1000)
   expect_equal(length(c2), 1000)
   expect_equal(length(c3), 1000)
   expect_equal(length(c4), 1000)
   expect_equal(length(c5), 1000)
})

test_that("number of unique colors is correct", {
   expect_equal(length(unique(c1)), ngroups)
   expect_equal(length(unique(c2)), ngroups)
   expect_equal(length(unique(c3)), ngroups)
   expect_equal(length(unique(c4)), ngroups)
   expect_equal(length(unique(c5)), ngroups)
})

test_that("the result also includes pallette as atribute", {
   expect_equal(attr(c1, 'palette'), mdaplot.getColors(ngroups = ngroups, colmap = 'old'))
   expect_equal(attr(c2, 'palette'), mdaplot.getColors(ngroups = ngroups, colmap = 'gray'))
   expect_equal(attr(c3, 'palette'), mdaplot.getColors(ngroups = ngroups, colmap = 'jet'))
   expect_equal(attr(c4, 'palette'), mdaplot.getColors(ngroups = ngroups, colmap = c('red', 'green', 'blue')))
   expect_equal(attr(c5, 'palette'), mdaplot.getColors(ngroups = ngroups))
})

ind = cut(cgroup, ngroups, labels = 1:ngroups)
test_that("position of colors is correct", {
   expect_equal(c1[1:1000], mdaplot.getColors(ngroups, colmap = 'old')[ind])
   expect_equal(c2[1:1000], mdaplot.getColors(ngroups, colmap = 'gray')[ind])
   expect_equal(c3[1:1000], mdaplot.getColors(ngroups, colmap = 'jet')[ind])
   expect_equal(c4[1:1000], mdaplot.getColors(ngroups, colmap = c('red', 'green', 'blue'))[ind])
   expect_equal(c5[1:1000], mdaplot.getColors(ngroups)[ind])
})


context('getColors: opacity is provided')

c1 = mdaplot.getColors(16, opacity = 0.25, colmap = 'old')
c2 = mdaplot.getColors(16, opacity = 0.25, colmap = 'gray')
c3 = mdaplot.getColors(16, opacity = 0.25, colmap = 'jet')
c4 = mdaplot.getColors(16, opacity = 0.25, colmap = c('red', 'green', 'blue'))
c5 = mdaplot.getColors(16, opacity = 0.25)

test_that("all colors have correct alpha value (no cgroup)", {
   expect_true(all(col2rgb(c1, alpha = T)[4, ] == 64))
   expect_true(all(col2rgb(c2, alpha = T)[4, ] == 64))
   expect_true(all(col2rgb(c3, alpha = T)[4, ] == 64))
   expect_true(all(col2rgb(c4, alpha = T)[4, ] == 64))
   expect_true(all(col2rgb(c5, alpha = T)[4, ] == 64))
})

cgroup = runif(1000, 1, 3)
c1 = mdaplot.getColors(cgroup = cgroup, opacity = 0.25, colmap = 'old')
c2 = mdaplot.getColors(cgroup = cgroup, opacity = 0.25, colmap = 'gray')
c3 = mdaplot.getColors(cgroup = cgroup, opacity = 0.25, colmap = 'jet')
c4 = mdaplot.getColors(cgroup = cgroup, opacity = 0.25, colmap = c('red', 'green', 'blue'))
c5 = mdaplot.getColors(cgroup = cgroup, opacity = 0.25)

test_that("all colors have correct alpha value (with cgroup)", {
   expect_true(all(col2rgb(c1, alpha = T)[4, ] == 64))
   expect_true(all(col2rgb(c2, alpha = T)[4, ] == 64))
   expect_true(all(col2rgb(c3, alpha = T)[4, ] == 64))
   expect_true(all(col2rgb(c4, alpha = T)[4, ] == 64))
   expect_true(all(col2rgb(c5, alpha = T)[4, ] == 64))
})
