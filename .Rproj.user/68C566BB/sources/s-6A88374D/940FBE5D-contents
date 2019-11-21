# Exploratory graphics for selection data
library(tidyverse)
library(data.table)

theme_set(theme_bw() + theme(strip.background =element_rect(fill="#e7e5e2")) +
            theme(strip.text = element_text(size =10)) +
            #theme(plot.title = element_text(hjust = 0.5), 
            theme(plot.title = element_text(), 
                  panel.background = element_rect(fill = "white", colour = NA), 
                  panel.grid.major = element_line(colour = "grey90", size = 0.2), 
                  panel.grid.minor = element_line(colour = "grey98", size = 0.5), 
                  panel.spacing = unit(0.25, "lines"), 
                  axis.text=element_text(size=10),
                  legend.position="top") 
)


# Load in data
#list of files
file_list <- list.files(path = "./data/",
                        pattern = "*_selection_data.csv",
                        full.names = T)
selection_dat <- tibble(bacteria = file_list) %>% # create a data frame
  # holding the file names
  mutate(file_contents = map(bacteria,          # read files into
                             ~ read_csv(.))) # a new data column
selection_df <- unnest(selection_dat, cols = c(file_contents))

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

#combine dN, dS, and omega data into one column for eventual plotting
selection_df <- selection_df %>% 
  gather(sel_class, val, dN:omega) %>%
  mutate(sel_class = gsub("sel_class", "", sel_class)) %>%
  arrange(bacteria, gene_name, sel_class, midpoint, tmp_pos)

#remove NAs to avoid warnings
selection_df_no_NA <- tidyr::drop_na(selection_df,val)

#reorder by median omega
selection_df_no_NA <- (selection_df_no_NA
             %>% group_by(bacteria)
             %>% mutate(omegaval=median(val[sel_class=="omega"]))
             %>% ungroup()
             %>% mutate(bacteria = reorder(bacteria,omegaval))
             %>% select(-omegaval)
)

#number of bacterial replicons
num_of_plots <- length(levels(selection_df_no_NA$bacteria))

#italic bacteria names
  levels(selection_df_no_NA$bacteria) <- c("ecoli" = expression(paste(italic("E.coli"), "")),
                          "bass" = expression(paste(italic("B. subtilis"), "")),
                          "sinoC" = expression(paste(italic("S.meliloti"), "")),
                          "pSymA" = expression(paste(italic("S.meliloti"), " pSymA")),
                          "strep" = expression(paste(italic("Streptomyces"), " Chromosome`")),
                          "pSymB" = expression(paste(italic("S.meliloti"), " pSymB")))

#choose colours
colours_arr <- rep(c("#CFE7C8","#D2E4DC","#A0747A"),num_of_plots)
#plot
vio_str<-(ggplot(selection_df_no_NA, aes(x=sel_class, y=val, fill=sel_class)) 
          + geom_violin() 
          + geom_jitter(position=position_jitter(0.2))
          + geom_boxplot(width=.1, outlier.shape=NA, fill = colours_arr) 
          + stat_boxplot(geom = "errorbar", width = 0.2) 
          + facet_wrap(~bacteria, labeller=label_parsed)
          + xlab("") 
          + ylab("Value") 
          + scale_color_manual(values=colours_arr)
          #make the omega a math symbol in both the legend and x-axis
          + scale_x_discrete(breaks = c("dN", "dS", "omega"),labels = c("dN","dS", expression(omega))) 
          #log scale
          #+ scale_y_continuous(trans='log10')
          + scale_fill_manual(values=c("#CFE7C8","#D2E4DC","#A0747A"), labels = c(" dN", " dS", expression(omega))) 
)

vio_str
