context("mdaplotg: plots with matrix as data and groupby parameter")

par(mfrow = c(2, 2))
data(people)
data <- people[, -6]
groupby <- as.factor(people[, "Sex"])
groupby_wrong <- people[, "Sex"]

groupby_df <- as.data.frame(people[, c("Sex", "Region")])
groupby_df$Sex <- factor(groupby_df$Sex, labels = c("M", "F"))
groupby_df$Region <- factor(groupby_df$Region, labels = c("S", "M"))

## shortcut for plotting function
tf <- function(groupby, ...) mdaplotg(data, groupby = groupby, ...)

# simple plots with groupby parameter

test_thatf("providing correct 'groupy' parameter works for any plot", {
   expect_silentf(tf(groupby, type = "p"))
   expect_silentf(tf(groupby, type = "l"))
   expect_silentf(tf(groupby, type = "b"))
   expect_silentf(tf(groupby, type = "h"))
})

test_thatf("providing 'groupy' as data frame parameter works for any plot", {
   expect_silentf(tf(groupby_df, type = "p"))
   expect_silentf(tf(groupby_df, type = "l"))
   expect_silentf(tf(groupby_df, type = "b"))
   expect_silentf(tf(groupby_df, type = "h"))
})

test_thatf("providing wrong 'groupy' parameter raises error for any plot", {
   expect_error(tf(groupby_wrong, type = "p"))
   expect_error(tf(groupby_wrong, type = "l"))
   expect_error(tf(groupby_wrong, type = "b"))
   expect_error(tf(groupby_wrong, type = "h"))
})

## checking labels

tfl <- function(groupby, ...) {
   mdaplotg(data, groupby = groupby, showl.labels = TRUE, ...)
}

test_thatf("parameter 'show.labels' is accepted", {
   expect_silentf(tfl(groupby, type = "p"))
   expect_silentf(tfl(groupby, type = "l"))
   expect_silentf(tfl(groupby, type = "b"))
   expect_silentf(tfl(groupby, type = "h"))
})

test_thatf("values can be used as labels", {
   expect_silentf(tfl(groupby, type = "p", labels = "values"))
   expect_silentf(tfl(groupby, type = "l", labels = "values"))
   expect_silentf(tfl(groupby, type = "b", labels = "values"))
   expect_silentf(tfl(groupby, type = "h", labels = "values"))
})

test_thatf("indices can be used as labels", {
   expect_silentf(tfl(groupby, type = "p", labels = "indices"))
   expect_silentf(tfl(groupby, type = "l", labels = "indices"))
   expect_silentf(tfl(groupby, type = "b", labels = "indices"))
   expect_silentf(tfl(groupby, type = "h", labels = "indices"))
})

## checking legend

test_thatf("legend can be hidden or have another position", {
   expect_silentf(tf(groupby, type = "p", show.legend = FALSE))
   expect_silentf(tf(groupby, type = "l", legend.position = "topleft"))
   expect_silentf(tf(groupby, type = "b", legend.position = "top"))
   expect_silentf(tf(groupby, type = "h", legend.position = "bottomleft"))
})

legend <- c("a", "b")
test_thatf("legend values can be modified", {
   expect_silentf(tf(groupby, type = "p", legend = legend))
   expect_silentf(tf(groupby, type = "l", legend = legend))
   expect_silentf(tf(groupby, type = "b", legend = legend))
   expect_silentf(tf(groupby, type = "h", legend = legend))
})

wrong_legend <- c("a", "b", "c", "d")
test_thatf("if number of provided legend values is wrong - error raises", {
   expect_error(tfgroupby, type = "p", legend = wrong_legend)
   expect_error(tfgroupby, type = "l", legend = wrong_legend)
   expect_error(tfgroupby, type = "b", legend = wrong_legend)
   expect_error(tfgroupby, type = "h", legend = wrong_legend)
})

## checking paremeters for plot series

test_thatf('"cex" can be specified as a single value or individual for each group"', {
   expect_silentf(mdaplotg(data, groupby = groupby, cex = 2, type = 'p'))
   expect_silentf(mdaplotg(data, groupby = groupby, cex = 2, type = 'b'))
   expect_silentf(mdaplotg(data, groupby = groupby, cex = 1:2, type = 'p'))
   expect_silentf(mdaplotg(data, groupby = groupby, cex = 1:2, type = 'b'))
   expect_error( mdaplotg(data, groupby = groupby, cex = 1:3, type = 'b'))
})

test_thatf('"pch" can be specified as a single value or individual for each group"', {
   expect_silentf(mdaplotg(data, groupby = groupby, pch = 2, type = 'p'))
   expect_silentf(mdaplotg(data, groupby = groupby, pch = 2, type = 'b'))
   expect_silentf(mdaplotg(data, groupby = groupby, pch = 1:2, type = 'p'))
   expect_silentf(mdaplotg(data, groupby = groupby, pch = 1:2, type = 'b'))
   expect_error( mdaplotg(data, groupby = groupby, pch = 1:3, type = 'b'))
})

test_thatf('"lty" can be specified as a single value or individual for each group"', {
   expect_silentf(mdaplotg(data, groupby = groupby, lty = 1, type = 'l'))
   expect_silentf(mdaplotg(data, groupby = groupby, lty = 1, type = 'b'))
   expect_silentf(mdaplotg(data, groupby = groupby, lty = 1:2, type = 'l'))
   expect_silentf(mdaplotg(data, groupby = groupby, lty = 1:2, type = 'b'))
   expect_error( mdaplotg(data, groupby = groupby, lty = 1:3, type = 'l'))
})

test_thatf('"lwd" can be specified as a single value or individual for each group"', {
   expect_silentf(mdaplotg(data, groupby = groupby, lwd = 1, type = 'l'))
   expect_silentf(mdaplotg(data, groupby = groupby, lwd = 1, type = 'b'))
   expect_silentf(mdaplotg(data, groupby = groupby, lwd = 1:2, type = 'l'))
   expect_silentf(mdaplotg(data, groupby = groupby, lwd = 1:2, type = 'b'))
   expect_error( mdaplotg(data, groupby = groupby, lwd = 1:3, type = 'l'))
})

test_thatf('"col" can be specified individual for each group"', {
   expect_silentf(mdaplotg(data, groupby = groupby, col = c('red', 'green'), type = 'p'))
   expect_silentf(mdaplotg(data, groupby = groupby, col = c('red', 'green'), type = 'h'))
   expect_silentf(mdaplotg(data, groupby = groupby, col = c('red', 'green'), type = 'l'))
   expect_silentf(mdaplotg(data, groupby = groupby, col = c('red', 'green'), type = 'b'))
   expect_error( mdaplotg(data, groupby = groupby, col = c('red', 'green', 'blue'), type = 'p'))
})

test_thatf('"col" can be specified as a single value"', {
   expect_silentf(mdaplotg(data, groupby = groupby, col = c('red'), type = 'p', pch = 1:2))
   expect_silentf(mdaplotg(data, groupby = groupby, col = c('red'), type = 'h', border = c('green', 'blue')))
   expect_silentf(mdaplotg(data, groupby = groupby, col = c('red'), type = 'l', lty = 1:2))
   expect_silentf(mdaplotg(data, groupby = groupby, col = c('red'), type = 'b', lty = 1:2, pch = 1:2))   
})


## checking ticks and ticklabels for x-axis

xticks = c(2, 4, 6, 8, 10, 12)
xticks_p = c(150, 170, 190)
xticklabels = colnames(people)[xticks]
xticklabels_p = c('Small', 'Normal', 'Large')

test_thatf('"xticks" can be user defined"', {
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'l', xticks = xticks))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'b', xticks = xticks))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'h', xticks = xticks))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'p', xticks = xticks_p))
})

test_thatf('"xticklabels" can be user defined"', {
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'l', xticks = xticks, xticklabels = xticklabels))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'b', xticks = xticks, xticklabels = xticklabels))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'h', xticks = xticks, xticklabels = xticklabels))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'p', xticks = xticks_p, xticklabels = xticklabels_p))
})

test_thatf('"xticklabels" can not be used without "xticks" or if they have different length"', {
   expect_error(mdaplotg(data, groupby = groupby, type = 'l', xticklabels = xticklabels))
   expect_error(mdaplotg(data, groupby = groupby, type = 'p', xticklabels = xticklabels_p))
   expect_error(mdaplotg(data, groupby = groupby, type = 'l', xticks = 1:10, xticklabels = xticklabels))
   expect_error(mdaplotg(data, groupby = groupby, type = 'p', xticks = 1:10, xticklabels = xticklabels_p))
})

## checking ticks and ticklabels for y-axis

yticks = c(10000, 30000)
yticks_p = c(50, 70, 90)
yticklabels = c('Mini', 'Mega')
yticklabels_p = c('Small', 'Normal', 'Large')

test_thatf('"yticks" can be user defined"', {
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'l', yticks = yticks))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'b', yticks = yticks))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'h', yticks = yticks))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'p', yticks = yticks_p))
})

test_thatf('"yticklabels" can be user defined"', {
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'b', yticks = yticks,   yticklabels = yticklabels))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'l', yticks = yticks,   yticklabels = yticklabels))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'h', yticks = yticks,   yticklabels = yticklabels))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'p', yticks = yticks_p, yticklabels = yticklabels_p))
})

test_thatf('"yticklabels" can not be used without "yticks" or if they have different length"', {
   expect_error(mdaplotg(data, groupby = groupby, type = 'l', yticklabels = yticklabels))
   expect_error(mdaplotg(data, groupby = groupby, type = 'p', yticklabels = yticklabels_p))
   expect_error(mdaplotg(data, groupby = groupby, type = 'l', yticks = 1:10, yticklabels = yticklabels))
   expect_error(mdaplotg(data, groupby = groupby, type = 'p', yticks = 1:10, yticklabels = yticklabels_p))
})


## using different colormaps and opacity colors

test_thatf('different colormaps can be used', {
   expect_silentf(mdaplotg(data, groupby = groupby_df, type = 'p', ))
   expect_silentf(mdaplotg(data, groupby = groupby_df, type = 'p', colmap = 'old'))
   expect_silentf(mdaplotg(data, groupby = groupby_df, type = 'p', colmap = 'gray'))
   expect_silentf(mdaplotg(data, groupby = groupby_df, type = 'p', colmap = 'jet'))
})

test_thatf('different colormaps can be used with single opacity value', {
   expect_silentf(mdaplotg(data, groupby = groupby_df, opacity = 0.5, type = 'p', ))
   expect_silentf(mdaplotg(data, groupby = groupby_df, opacity = 0.5, type = 'p', colmap = 'old'))
   expect_silentf(mdaplotg(data, groupby = groupby_df, opacity = 0.5, type = 'p', colmap = 'gray'))
   expect_silentf(mdaplotg(data, groupby = groupby_df, opacity = 0.5, type = 'p', colmap = 'jet'))
})

test_thatf('different colormaps can be used with individual opacity values', {
   expect_silentf(mdaplotg(data, groupby = groupby_df, opacity = c(0.2, 0.4, 0.6, 0.8), type = 'p', ))
   expect_silentf(mdaplotg(data, groupby = groupby_df, opacity = c(0.2, 0.4, 0.6, 0.8), type = 'p', colmap = 'old'))
   expect_silentf(mdaplotg(data, groupby = groupby_df, opacity = c(0.2, 0.4, 0.6, 0.8), type = 'p', colmap = 'gray'))
   expect_silentf(mdaplotg(data, groupby = groupby_df, opacity = c(0.2, 0.4, 0.6, 0.8), type = 'p', colmap = 'jet'))
})

## test grid parameters
test_thatf('grid can be shown', {
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'p', show.grid = TRUE))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'l', show.grid = TRUE))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'h', show.grid = TRUE))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'b', show.grid = TRUE))
})

test_thatf('grid can be hidden', {
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'p', show.grid = FALSE))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'l', show.grid = FALSE))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'h', show.grid = FALSE))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'b', show.grid = FALSE))
})

test_thatf('grid color and thickness can be changed', {
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'p', show.grid = TRUE, grid.col = 'red', grid.lwd = '2'))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'l', show.grid = TRUE, grid.col = 'red', grid.lwd = '2'))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'h', show.grid = TRUE, grid.col = 'red', grid.lwd = '2'))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'b', show.grid = TRUE, grid.col = 'red', grid.lwd = '2'))
})


## handling hidden data

data = mda.exclrows(data, data[, 'Beer'] > 300)
data = mda.exclcols(data, 'Height')

test_thatf('hidden data is not shown by default', {
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'p'))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'l'))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'b'))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'h'))
})

test_thatf('hidden data is shown as gray if needed and labels are produced correctly', {
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'p', show.excluded = TRUE, show.labels = T))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'l', show.excluded = TRUE, show.labels = T))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'b', show.excluded = TRUE, show.labels = T))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'h', show.excluded = TRUE, show.labels = T))
})

test_thatf('hidden data can be used with labels as values', {
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'p', show.excluded = TRUE, show.labels = T, labels= 'values'))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'l', show.excluded = TRUE, show.labels = T, labels= 'values'))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'b', show.excluded = TRUE, show.labels = T, labels= 'values'))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'h', show.excluded = TRUE, show.labels = T, labels= 'values'))
})

test_thatf('hidden data can be used with labels as names', {
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'p', show.excluded = TRUE, show.labels = T, labels = 'names'))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'l', show.excluded = TRUE, show.labels = T, labels = 'names'))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'b', show.excluded = TRUE, show.labels = T, labels = 'names'))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'h', show.excluded = TRUE, show.labels = T, labels = 'names'))
})

test_thatf('hidden data can be used with labels as indices', {
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'p', show.excluded = TRUE, show.labels = T, labels='indices'))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'l', show.excluded = TRUE, show.labels = T, labels='indices'))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'b', show.excluded = TRUE, show.labels = T, labels='indices'))
   expect_silentf(mdaplotg(data, groupby = groupby, type = 'h', show.excluded = TRUE, show.labels = T, labels='indices'))
})

context("mdaplotg: plots with list as data")
context("mdaplotg: plots with matrix as data without groupby")

