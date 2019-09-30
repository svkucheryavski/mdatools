context('mdaplotg: simple plots with groupby parameter')

data(people)
data = people
groupby = as.factor(people[, 'Sex'])
groupby_wrong = people[, 'Sex']

# any plot with groupby parameter
test_that('providing correct "groupy" parameter works for any plot', {
   expect_silent(mdaplotg(data, groupby = groupby, type = 'p'))
   expect_silent(mdaplotg(data, groupby = groupby, type = 'l'))
   expect_silent(mdaplotg(data, groupby = groupby, type = 'b'))
   expect_silent(mdaplotg(data, groupby = groupby, type = 'h'))
})

test_that('providing wrong "groupy" parameter raises error for any plot', {
   expect_error(mdaplotg(data, groupby = groupby_wrong, type = 'p'))
   expect_error(mdaplotg(data, groupby = groupby_wrong, type = 'l'))
   expect_error(mdaplotg(data, groupby = groupby_wrong, type = 'b'))
   expect_error(mdaplotg(data, groupby = groupby_wrong, type = 'h'))
})

test_that('"show.labels" is accepted', {
   expect_silent(mdaplotg(data, groupby = groupby, type = 'p', show.labels = T))
   expect_silent(mdaplotg(data, groupby = groupby, type = 'l', show.labels = T))
   expect_silent(mdaplotg(data, groupby = groupby, type = 'b', show.labels = T))
   expect_silent(mdaplotg(data, groupby = groupby, type = 'h', show.labels = T))
})

test_that('values can be used as labels', {
   expect_silent(mdaplotg(data, groupby = groupby, type = 'p', show.labels = T, labels = 'values'))
   expect_silent(mdaplotg(data, groupby = groupby, type = 'l', show.labels = T, labels = 'values'))
   expect_silent(mdaplotg(data, groupby = groupby, type = 'b', show.labels = T, labels = 'values'))
   expect_silent(mdaplotg(data, groupby = groupby, type = 'h', show.labels = T, labels = 'values'))
})

test_that('indices can be used as labels', {
   expect_silent(mdaplotg(data, groupby = groupby, type = 'p', show.labels = T, labels = 'indices'))
   expect_silent(mdaplotg(data, groupby = groupby, type = 'l', show.labels = T, labels = 'indices'))
   expect_silent(mdaplotg(data, groupby = groupby, type = 'b', show.labels = T, labels = 'indices'))
   expect_silent(mdaplotg(data, groupby = groupby, type = 'h', show.labels = T, labels = 'indices'))
})


test_that('legend can be hidden or have another position', {
   expect_silent(mdaplotg(data, groupby = groupby, type = 'p', show.legend = FALSE))
   expect_silent(mdaplotg(data, groupby = groupby, type = 'l', legend.position = 'topleft'))
   expect_silent(mdaplotg(data, groupby = groupby, type = 'b', legend.position = 'top'))
   expect_silent(mdaplotg(data, groupby = groupby, type = 'h', legend.position = 'bottomleft'))
})

test_that('legend values can be modified', {
   expect_silent(mdaplotg(data, groupby = groupby, type = 'p', legend = c('a', 'b')))
   expect_silent(mdaplotg(data, groupby = groupby, type = 'l', legend = c('a', 'b')))
   expect_silent(mdaplotg(data, groupby = groupby, type = 'b', legend = c('a', 'b')))
   expect_silent(mdaplotg(data, groupby = groupby, type = 'h', legend = c('a', 'b')))
})

test_that('if number of provided legend values is wrong - error raises', {
   expect_error(mdaplotg(data, groupby = groupby, type = 'p', legend = c('a', 'b', 'c', 'd')))
   expect_error(mdaplotg(data, groupby = groupby, type = 'l', legend = c('a', 'b', 'c', 'd')))
   expect_error(mdaplotg(data, groupby = groupby, type = 'b', legend = c('a', 'b', 'c', 'd')))
   expect_error(mdaplotg(data, groupby = groupby, type = 'h', legend = c('a', 'b', 'c', 'd')))
})

test_that('"cex" can be specified as a single value or individual for each group"', {
   expect_silent(mdaplotg(data, groupby = groupby, cex = 2, type = 'p'))
   expect_silent(mdaplotg(data, groupby = groupby, cex = 2, type = 'b'))
   expect_silent(mdaplotg(data, groupby = groupby, cex = 1:2, type = 'p'))
   expect_silent(mdaplotg(data, groupby = groupby, cex = 1:2, type = 'b'))
   expect_error( mdaplotg(data, groupby = groupby, cex = 1:3, type = 'b'))
})

test_that('"pch" can be specified as a single value or individual for each group"', {
   expect_silent(mdaplotg(data, groupby = groupby, pch = 2, type = 'p'))
   expect_silent(mdaplotg(data, groupby = groupby, pch = 2, type = 'b'))
   expect_silent(mdaplotg(data, groupby = groupby, pch = 1:2, type = 'p'))
   expect_silent(mdaplotg(data, groupby = groupby, pch = 1:2, type = 'b'))
   expect_error( mdaplotg(data, groupby = groupby, pch = 1:3, type = 'b'))
})

test_that('"lty" can be specified as a single value or individual for each group"', {
   expect_silent(mdaplotg(data, groupby = groupby, lty = 1, type = 'l'))
   expect_silent(mdaplotg(data, groupby = groupby, lty = 1, type = 'b'))
   expect_silent(mdaplotg(data, groupby = groupby, lty = 1:2, type = 'l'))
   expect_silent(mdaplotg(data, groupby = groupby, lty = 1:2, type = 'b'))
   expect_error( mdaplotg(data, groupby = groupby, lty = 1:3, type = 'l'))
})

test_that('"lwd" can be specified as a single value or individual for each group"', {
   expect_silent(mdaplotg(data, groupby = groupby, lwd = 1, type = 'l'))
   expect_silent(mdaplotg(data, groupby = groupby, lwd = 1, type = 'b'))
   expect_silent(mdaplotg(data, groupby = groupby, lwd = 1:2, type = 'l'))
   expect_silent(mdaplotg(data, groupby = groupby, lwd = 1:2, type = 'b'))
   expect_error( mdaplotg(data, groupby = groupby, lwd = 1:3, type = 'l'))
})

test_that('"col" can be specified as a single value or individual for each group"', {
   expect_silent(mdaplotg(data, groupby = groupby, col = c('red', 'green'), type = 'p'))
   expect_silent(mdaplotg(data, groupby = groupby, col = c('red', 'green'), type = 'h'))
   expect_silent(mdaplotg(data, groupby = groupby, col = c('red', 'green'), type = 'l'))
   expect_silent(mdaplotg(data, groupby = groupby, col = c('red', 'green'), type = 'b'))
   expect_silent( mdaplotg(data, groupby = groupby, col = 'red', type = 'p', pch = 1:2))
   expect_error( mdaplotg(data, groupby = groupby, col = c('red', 'green', 'blue'), type = 'p'))
})


## checking ticks and ticklabels

xticks = c(2, 4, 6, 8, 10, 12)
xticks_p = c(150, 170, 190)
xticklabels = colnames(people)[xticks]
xticklabels_p = c('Small', 'Normal', 'Large')

test_that('"xticks" can be user defined"', {
   expect_silent(mdaplotg(data, groupby = groupby, type = 'l', xticks = xticks))
   expect_silent(mdaplotg(data, groupby = groupby, type = 'b', xticks = xticks))
   expect_silent(mdaplotg(data, groupby = groupby, type = 'h', xticks = xticks))
   expect_silent(mdaplotg(data, groupby = groupby, type = 'p', xticks = xticks_p))
})

test_that('"xticklabels" can be user defined"', {
   expect_silent(mdaplotg(data, groupby = groupby, type = 'l', xticks = xticks, xticklabels = xticklabels))
   expect_silent(mdaplotg(data, groupby = groupby, type = 'b', xticks = xticks, xticklabels = xticklabels))
   expect_silent(mdaplotg(data, groupby = groupby, type = 'h', xticks = xticks, xticklabels = xticklabels))
   expect_silent(mdaplotg(data, groupby = groupby, type = 'p', xticks = xticks_p, xticklabels = xticklabels_p))
})

test_that('"xticklabels" can not be used without "xticks" or if they have different length"', {
   expect_error(mdaplotg(data, groupby = groupby, type = 'l', xticklabels = xticklabels))
   expect_error(mdaplotg(data, groupby = groupby, type = 'p', xticklabels = xticklabels_p))
   expect_error(mdaplotg(data, groupby = groupby, type = 'l', xticks = 1:10, xticklabels = xticklabels))
   expect_error(mdaplotg(data, groupby = groupby, type = 'p', xticks = 1:10, xticklabels = xticklabels_p))
})

yticks = c(10000, 30000)
yticks_p = c(50, 70, 90)
yticklabels = c('Mini', 'Mega')
yticklabels_p = c('Small', 'Normal', 'Large')

test_that('"yticks" can be user defined"', {
   expect_silent(mdaplotg(data, groupby = groupby, type = 'l', yticks = yticks))
   expect_silent(mdaplotg(data, groupby = groupby, type = 'b', yticks = yticks))
   expect_silent(mdaplotg(data, groupby = groupby, type = 'h', yticks = yticks))
   expect_silent(mdaplotg(data, groupby = groupby, type = 'p', yticks = yticks_p))
})

test_that('"yticklabels" can be user defined"', {
   expect_silent(mdaplotg(data, groupby = groupby, type = 'b', yticks = yticks,   yticklabels = yticklabels))
   expect_silent(mdaplotg(data, groupby = groupby, type = 'l', yticks = yticks,   yticklabels = yticklabels))
   expect_silent(mdaplotg(data, groupby = groupby, type = 'h', yticks = yticks,   yticklabels = yticklabels))
   expect_silent(mdaplotg(data, groupby = groupby, type = 'p', yticks = yticks_p, yticklabels = yticklabels_p))
})

test_that('"yticklabels" can not be used without "yticks" or if they have different length"', {
   expect_error(mdaplotg(data, groupby = groupby, type = 'l', yticklabels = yticklabels))
   expect_error(mdaplotg(data, groupby = groupby, type = 'p', yticklabels = yticklabels_p))
   expect_error(mdaplotg(data, groupby = groupby, type = 'l', yticks = 1:10, yticklabels = yticklabels))
   expect_error(mdaplotg(data, groupby = groupby, type = 'p', yticks = 1:10, yticklabels = yticklabels_p))
})