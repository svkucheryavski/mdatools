## Normal bar and group bar plots

data1 = matrix(c(1:3, 60, 30, 10, 55, 30, 15), ncol = 3)
colnames(data1) = c('x', 'cal', 'val')
data2 = matrix(c(1:20, 1:20, 80:61), ncol = 3)
colnames(data2) = c('x', 'cal', 'val')
data3 = matrix(c(1:3, 60, 30, 10, 55, 30, 15, 100, 60, 20), ncol = 4)
colnames(data3) = c('x', 'cal', 'val', 'cv')

par(mfrow = c(3, 2))
mdaplot(data1[, 1:2], type = 'h', col = 'red', show.labels = T)
mdaplotg(data1, type = 'h', legend = c('cal', 'val'), col = c('green', 'orange'))
mdaplot(data2[, 1:2], type = 'h', show.labels = T)
mdaplotg(data2, type = 'h')
mdaplot(data3[, 1:2], type = 'h')
mdaplotg(data3, type = 'h', show.labels = T, labels = data3[, -1])

## Bar plots with negative values

data1 = matrix(c(1:3, 60, -30, 10, 55, 30, -15), ncol = 3)
colnames(data1) = c('x', 'cal', 'val')
data2 = matrix(c(1:20, 1:20, 80:61), ncol = 3)
colnames(data2) = c('x', 'cal', 'val')
data3 = matrix(c(1:3, 60, 30, -10, 55, 30, -15, 100, -60, 20), ncol = 4)
colnames(data3) = c('x', 'cal', 'val', 'cv')

par(mfrow = c(2, 1))
mdaplot(data1[, 1:2], type = 'h', col = 'red', show.labels = T)
mdaplotg(data3, type = 'h', show.labels = T, labels = data3[, -1])

# Bar plots with data in lists
m1 = cbind(c(1:3), c(60, 30, -10))
rownames(m1) = c('A', 'B', 'C')
m2 = cbind(c(4:7), c(55, 30, -15, 11))
m3 = cbind(c(8:9), c(100, -60))

data3 = list(m1, m2, m3)
labels = list(data3[[1]][, 2], data3[[2]][, 2], data3[[3]][, 2])

par(mfrow = c(2, 1))
mdaplotg(data3, type = 'h', show.labels = T, labels = labels)
mdaplotg(data3, type = 'h', show.labels = T)

