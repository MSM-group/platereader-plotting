install.packages("readr")
install.packages("writexl")

library(readr)
library(writexl)

#insert file name

txt_file <- "C:\\Users\\schweidi\\Documents\\GitHub\\platereader-plotting\\data\\hplc_txt\\filename.txt"

data <- read_delim(txt_file, delim = "\t")

#insert file name

xlsx_file <- "C:\\Users\\schweidi\\Documents\\GitHub\\platereader-plotting\\data\\hplc_xlsx\\filename.xlsx"

write_xlsx(data, path = xlsx_file)
