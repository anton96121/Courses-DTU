library("rstudioapi")
path_file_dir = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path_file_dir)

