# Exploratory graphics for selection data
library(tidyverse)
library(data.table)

# Load in data
#list of files
file_list <- list.files(path = "./data/",
                        pattern = "*_selection_data.csv",
                        full.names = T)
selection_dat <- data_frame(bacteria = file_list) %>% # create a data frame
  # holding the file names
  mutate(file_contents = map(bacteria,          # read files into
                             ~ read_csv(.))) # a new data column
  )
selection_df <- unnest(selection_dat)

#change the file names to bacteria names
selection_df <- selection_df %>%
  mutate_if(is.character, str_replace_all, pattern = '^./data/bass_selection_data.csv', replacement = 'bass')
selection_df <- selection_df %>%
  mutate_if(is.character, str_replace_all, pattern = '^./data/ecoli_selection_data.csv', replacement = 'ecoli')
selection_df <- selection_df %>%
  mutate_if(is.character, str_replace_all, pattern = '^./data/strep_selection_data.csv', replacement = 'strep')
selection_df <- selection_df %>%
  mutate_if(is.character, str_replace_all, pattern = '^./data/sinoC_selection_data.csv', replacement = 'sinoC')
selection_df <- selection_df %>%
  mutate_if(is.character, str_replace_all, pattern = '^./data/pSymA_selection_data.csv', replacement = 'pSymA')
selection_df <- selection_df %>%
  mutate_if(is.character, str_replace_all, pattern = '^./data/pSymB_selection_data.csv', replacement = 'pSymB')
levels(selection_df$bacteria)
class(selection_df)
summarize(selection_df)
