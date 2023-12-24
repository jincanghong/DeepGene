
# library(kaos)
# data <- read.csv("inputs.csv", header = T, stringsAsFactors = F)
# print("数据维度：")
# print(dim(data))

# index <- read.csv("index.csv", header = T, stringsAsFactors = F)
# print("标签维度：")
# print(dim(index))
# index <- as.character(index[[1]])

# res <- 200
# for (i in seq(dim(data)[2])) {
#     print(sprintf("正在处理第%d个...", i))
#     cgr_list <- cgr(data[[i]], seq.base = c("A", "C", "T", "G"), res = res)
#     cgr_matrix <- cgr_list$matrix
#     write.csv(cgr_matrix, file = sprintf("FCGR_output/%s_%d.csv", index[i], res), quote = F, row.names = F)
# }



library(kaos)
cgr_list <- cgr(HIV, seq.base = c("a","g","t","c"), res = 100)
cgr_matrix <- cgr_list$matrix
# cgr.plot(cgr_list, mode = "matrix")
write.csv(HIV, file = "HIV.csv", quote = F, row.names = F)
write.csv(cgr_matrix, file = "test_output.csv", quote = F, row.names = F)