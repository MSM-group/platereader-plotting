install.packages("readr")
install.packages("writexl")

library(readr)
library(writexl)

txt_file <- "C:\\Users\\schweidi\\Documents\\GitHub\\platereader-plotting\\data\\hplctesttxt\\Mutation_LB_S146_alignment.txt"

data <- read_delim(txt_file, delim = "\t")

xlsx_file <- "C:\\Users\\schweidi\\Documents\\GitHub\\platereader-plotting\\data\\hplctest\\Mutation_LB_S170_alignment.xlsx"

write_xlsx(data, path = xlsx_file)
