# Exploratory graphics for selection data
library(tidyverse)
library(data.table)
library(devtools)
library(cowplot)
library(ggpubr)
library(plyr)


options(scipen=10000)
theme_set(theme_bw() + theme(strip.background =element_rect(fill="#e7e5e2")) +
            theme(strip.text = element_text(size =10)) +
            #theme(plot.title = element_text(hjust = 0.5), 
            theme(plot.title = element_text(), 
                  panel.background = element_rect(fill = "white", colour = NA), 
                  panel.grid.major = element_line(colour = "grey90", size = 0.2), 
                  panel.grid.minor = element_line(colour = "grey98", size = 0.5), 
                  panel.spacing = unit(0.25, "lines"), 
                  axis.text=element_text(size=10),
                  legend.title = element_blank(),
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

###################################
#########################
#chunking the data into sections of the genome
#########################
#decreasing order
selection_df_medians <- arrange(selection_df_no_NA,tmp_pos)
nmax_pos <- max(selection_df_medians$tmp_pos)
nmin_pos <- min(selection_df_medians$tmp_pos)
nmin_pos
#lenght of section of genome
chunklen <- 100000
chunk_len_of_genome <- round_any(nmax_pos, chunklen, f=ceiling)
chunk_len_of_genome
chunks <- seq(0, chunk_len_of_genome, chunklen)
chunks

#expty vec to hold the rows that split the dat into X kb chunks
exp_rows_to_split_dat <- vector()
for (i in chunks) {
  exp_rows <-
    which(abs(selection_df_medians$tmp_pos-i)==min(abs(selection_df_medians$tmp_pos-i)))
  # finding the closest number to each 10kb without going over it
  exp_max_row <- max(exp_rows)
  #  actual_pos <- data_chrom_ordered$midpoint[max_row]
  #  rows_to_split_dat <- c(rows_to_split_dat, actual_pos)
  exp_rows_to_split_dat <- c(exp_rows_to_split_dat, exp_max_row)
}#for
#rows_to_split_dat <- rows_to_split_dat + 1
exp_rows_to_split_dat
#make new column with "groups" of these chunks
new_sections <- rep(chunks, times = exp_rows_to_split_dat)
head(new_sections)                    

selection_df_medians$new_sections <- new_sections



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
                          "sinoC" = expression(paste(italic("S.meliloti"), " Chromosome")),
                          "pSymA" = expression(paste(italic("S.meliloti"), " pSymA")),
                          "strep" = expression(paste(italic("Streptomyces"), "")),
                          "pSymB" = expression(paste(italic("S.meliloti"), " pSymB")))

#choose colours
colours_arr <- rep(c("#CFE7C8","#D2E4DC","#A0747A"),num_of_plots)
#plot
set.seed(1738)
vio_str_box <-(ggplot(selection_df_no_NA, aes(x=sel_class, y=val, fill = sel_class, colour = sel_class)) 
          + geom_jitter(position=position_jitter(0.2))
          + geom_violin(colour = "black", fill = NA) 
          #+ geom_boxplot(width=.1, outlier.shape=NA, fill = colours_arr) 
          + geom_boxplot(width=.1, outlier.shape=NA, colour = "black", fill = colours_arr) 
          + stat_boxplot(geom = "errorbar", width = 0.2, colour = "black") 
          #omega = 1 reference line
          +  geom_hline(yintercept=1, linetype="dashed", color = "#023C40")
          + facet_wrap(~bacteria, labeller=label_parsed)
          + xlab("") 
          + ylab("Value") 
          + scale_color_manual(values=c("#CFE7C8","#D2E4DC","#A0747A"),labels = c(" dN", " dS", expression(omega)))
          #+ scale_color_manual(values=c("#CFE7C8","#D2E4DC","#A0747A"))
          #make the omega a math symbol in x-axis
          + scale_x_discrete(breaks = c("dN", "dS", "omega"),labels = c("dN","dS", expression(omega))) 
          #log scale and removing trailing zeros from y-axis labels
          + scale_y_continuous(trans='log10',labels = function(x) ifelse(x == 0, "0", x))
          #+ scale_fill_manual(values=c("#CFE7C8","#D2E4DC","#A0747A"), labels = c(" dN", " dS", expression(omega))) 
)

vio_str_box


#hex plot
library(hexbin)
omega_data <- selection_df[which(selection_df$sel_class == "omega"),]
ecol_dat <- selection_df[which(selection_df$bacteria == "ecoli"),]
pa_dat <- omega_data[which(omega_data$bacteria == "ecoli"),]

hex_p <- (ggplot(omega_data, aes(x=tmp_pos, y=val)) 
          + facet_wrap(~bacteria, labeller=label_parsed)
          + geom_hex(bins=30) 
          #omega = 1 reference line
          +  geom_hline(yintercept=1, linetype="dashed", color = "red")
          + xlab("Genomic Position") 
          + ylab("Value")
          + scale_y_continuous(trans='log10',labels = function(x) ifelse(x == 0, "0", x))
)
hex_p

pa_hex <- (ggplot(pa_omeg, aes(x=tmp_pos, y=val)) 
          #+ geom_hex(bins=30) 
          + geom_point() 
          #omega = 1 reference line
          +  geom_hline(yintercept=1, linetype="dashed", color = "red")
          + xlab("Genomic Position") 
          + ylab("Value")
          + scale_y_continuous(trans='log10',labels = function(x) ifelse(x == 0, "0", x))
)
pa_hex


# Scatter plot colored by groups ("Species")
sp <- ggscatter(pa_dat, x = "tmp_pos", y = "val",
#sp <- geom_smooth(method = lm, val ~ tmp_pos, x = "tmp_pos", y = "val", data = pa_dat,
                color = "sel_class", palette = "jco",
                size = 3, alpha = 0.6)+
  geom_smooth(method = lm) +
          scale_y_continuous(trans='log10') +
  border()                                         
# Marginal density plot of x (top panel) and y (right panel)
xplot <- ggdensity(pa_dat, "tmp_pos", fill = "sel_class",
                   palette = "jco")
yplot <- ggdensity(pa_dat, "val", fill = "sel_class", 
                   palette = "jco")+
  #this transformation does not look right....
          scale_y_continuous(trans='log10') +
  rotate()
# Cleaning the plots
sp <- sp + rremove("legend")
yplot <- yplot + clean_theme() + rremove("legend") 
xplot <- xplot + clean_theme() + rremove("legend")
# Arranging the plot using cowplot
library(cowplot)
plot_grid(xplot, NULL, sp, yplot, ncol = 2, align = "hv", 
          rel_widths = c(2, 1), rel_heights = c(1, 2))

#########################
#chunking the data into sections of the genome
#########################
#adding a fake row for position 0
pa_dat <- add_row(pa_dat, tmp_pos = 0)
pa_dat <- arrange(pa_dat, tmp_pos)
nmax_pos <- max(pa_dat$tmp_pos)
nmin_pos <- min(pa_dat$tmp_pos)
nmin_pos
#lenght of section of genome
chunklen <- 100000
chunk_len_of_genome <- round_any(nmax_pos, chunklen, f=ceiling)
chunk_len_of_genome
chunks <- seq(nmin_pos, chunk_len_of_genome, chunklen)
chunks

#expty vec to hold the rows that split the dat into X kb chunks
exp_rows_to_split_dat <- vector()
for (i in chunks) {
  exp_rows <-
    which(abs(pa_dat$tmp_pos-i)==min(abs(pa_dat$tmp_pos-i)))
  # finding the closest number to each 10kb without going over it
  exp_max_row <- max(exp_rows)
  #  actual_pos <- data_chrom_ordered$midpoint[max_row]
  #  rows_to_split_dat <- c(rows_to_split_dat, actual_pos)
  exp_rows_to_split_dat <- c(exp_rows_to_split_dat, exp_max_row)
}#for
#rows_to_split_dat <- rows_to_split_dat + 1
exp_rows_to_split_dat
#split data by the above specified rows
#sometimes the closest number was technically in the next 10kb chunk of seq
#but its fine.
data_just_exp <- pa_dat$val
#list_dat_sets <- split(data_chrom_ordered, findInterval(1:nrow(data_chrom_ordered), rows_to_split_dat))
list_dat_sets_just_exp <- split(data_just_exp, findInterval(1:nrow(pa_dat), exp_rows_to_split_dat))
list_dat_sets_just_exp

#get median expression in each of the sections
list_total_exp_level <- lapply(list_dat_sets_just_exp, median)
print("list_total_exp_level")
list_total_exp_level
#as df
median_gene_exp_10kb <- as.data.frame(matrix(unlist(list_total_exp_level), byrow = F))
median_gene_exp_10kb <- cbind(median_gene_exp_10kb, new_pos)
write.csv(median_gene_exp_10kb)
#write.table(median_gene_exp_10kb, 'median_10kb_exp_data.csv', sep = "\t")

(ggplot(data = median_gene_exp_10kb, aes(x = new_pos, y = V1))
+ geom_point()
  )

###################################
#########################
#chunking the data into sections of the genome
#########################
#decreasing order
selection_df_medians <- arrange(selection_df_no_NA,tmp_pos)
nmax_pos <- max(selection_df_medians$tmp_pos)
nmin_pos <- min(selection_df_medians$tmp_pos)
nmin_pos
#lenght of section of genome
chunklen <- 100000
chunk_len_of_genome <- round_any(nmax_pos, chunklen, f=ceiling)
chunk_len_of_genome
chunks <- seq(0, chunk_len_of_genome, chunklen)
chunks

#expty vec to hold the rows that split the dat into X kb chunks
exp_rows_to_split_dat <- vector()
for (i in chunks) {
  exp_rows <-
    which(abs(selection_df_medians$tmp_pos-i)==min(abs(selection_df_medians$tmp_pos-i)))
  # finding the closest number to each 10kb without going over it
  exp_max_row <- max(exp_rows)
  #  actual_pos <- data_chrom_ordered$midpoint[max_row]
  #  rows_to_split_dat <- c(rows_to_split_dat, actual_pos)
  exp_rows_to_split_dat <- c(exp_rows_to_split_dat, exp_max_row)
}#for
#rows_to_split_dat <- rows_to_split_dat + 1
exp_rows_to_split_dat
#make new column with "groups" of these chunks
new_sections <- rep(chunks, times = exp_rows_to_split_dat)
new_sections                    

selection_df_medians$new_sections <- new_sections
