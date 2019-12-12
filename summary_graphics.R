# Exploratory graphics for selection data
library(tidyverse)
library(data.table)
library(devtools)
library(cowplot)
library(ggpubr)
library(plyr)
library(dplyr)
library(directlabels)
library(gridExtra)
library(grid)


options(scipen=10000)
theme_set(theme_bw() + theme(strip.background =element_rect(fill="#e7e5e2")) +
            theme(strip.text = element_text(size =12)) +
            #theme(plot.title = element_text(hjust = 0.5), 
            theme(plot.title = element_text(hjust = 0.5, size = 15), 
                  panel.background = element_rect(fill = "white", colour = NA), 
                  panel.grid.major = element_line(colour = "grey90", size = 0.2), 
                  panel.grid.minor = element_line(colour = "grey98", size = 0.5), 
                  panel.spacing = unit(0.25, "lines"), 
                  strip.background = element_rect(fill = "#E6E1EA"),
                  axis.text=element_text(size=15),
                  axis.text.y.right = element_text(size=15),
                  axis.title = element_text(size = 15),
                  legend.title = element_blank(),
                  legend.text = element_text(size = 15),
                  legend.key = element_blank(),
                  legend.background=element_blank(),
                  legend.position="top") 
)

strwrap_strip_text = function(p, pad=0.05) { 
  # get facet font attributes
  th = theme_get()
  if (length(p$theme) > 0L)
    th = th + p$theme
  
  require("grid")
  grobs <- ggplotGrob(p)
  
  # wrap strip x text
  if ((class(p$facet)[1] == "grid" && !is.null(names(p$facet$cols))) ||
      class(p$facet)[1] == "wrap")
  {
    ps = calc_element("strip.text.x", th)[["size"]]
    family = calc_element("strip.text.x", th)[["family"]]
    face = calc_element("strip.text.x", th)[["face"]]
    
    if (class(p$facet)[1] == "wrap") {
      nm = names(p$facet$facets)
    } else {
      nm = names(p$facet$cols)
    }
    
    # get number of facet columns
    levs = levels(factor(p$data[[nm]]))
    npanels = length(levs)
    if (class(p$facet)[1] == "wrap") {
      cols = n2mfrow(npanels)[1]
    } else {
      cols = npanels
    }
    
    # get plot width
    sum = sum(sapply(grobs$width, function(x) convertWidth(x, "in")))
    panels_width = par("din")[1] - sum  # inches
    # determine strwrap width
    panel_width = panels_width / cols
    mx_ind = which.max(nchar(levs))
    char_width = strwidth(levs[mx_ind], units="inches", cex=ps / par("ps"), 
                          family=family, font=gpar(fontface=face)$font) / 
      nchar(levs[mx_ind])
    width = floor((panel_width - pad)/ char_width)  # characters
    
    # wrap facet text
    p$data[[nm]] = unlist(lapply(strwrap(p$data[[nm]], width=width, 
                                         simplify=FALSE), paste, collapse="\n"))
  }
  
  if (class(p$facet)[1] == "grid" && !is.null(names(p$facet$rows))) {  
    ps = calc_element("strip.text.y", th)[["size"]]
    family = calc_element("strip.text.y", th)[["family"]]
    face = calc_element("strip.text.y", th)[["face"]]
    
    nm = names(p$facet$rows)
    
    # get number of facet columns
    levs = levels(factor(p$data[[nm]]))
    rows = length(levs)
    
    # get plot height
    sum = sum(sapply(grobs$height, function(x) convertWidth(x, "in")))
    panels_height = par("din")[2] - sum  # inches
    # determine strwrap width
    panels_height = panels_height / rows
    mx_ind = which.max(nchar(levs))
    char_height = strwidth(levs[mx_ind], units="inches", cex=ps / par("ps"), 
                           family=family, font=gpar(fontface=face)$font) / 
      nchar(levs[mx_ind])
    width = floor((panels_height - pad)/ char_height)  # characters
    
    # wrap facet text
    p$data[[nm]] = unlist(lapply(strwrap(p$data[[nm]], width=width, 
                                         simplify=FALSE), paste, collapse="\n"))
  }
  
  invisible(p)
}


#read in data'
sel_dat <- read.csv("./data/all_bac_selection_data.csv", header = FALSE)
#drop first column because it is just the row numbers
sel_dat <- sel_dat[,-1]
#col names
colnames(sel_dat) <- c("value", "position", "class", "bacteria")
#remove NAs
sel_dat <- sel_dat[!is.na(sel_dat$value),]
#number of bacterial replicons
num_of_plots <- length(levels(sel_dat$bacteria))

###################################
#########################
#chunking the data into sections of the genome
#########################
#decreasing order
sel_dat_ord <- arrange(sel_dat,position)
nmax_pos <- max(sel_dat_ord$position)
nmin_pos <- min(sel_dat_ord$position)
nmin_pos
#lenght of section of genome
chunklen <- 10000
chunk_len_of_genome <- round_any(nmax_pos, chunklen, f=ceiling)
chunk_len_of_genome
chunks <- seq(chunklen, chunk_len_of_genome, chunklen)
chunks

#expty vec to hold the rows that split the dat into X kb chunks
exp_rows_to_split_dat <- vector()
for (i in chunks) {
  exp_rows <-
    which(abs(sel_dat_ord$position-i)==min(abs(sel_dat$position-i)))
  # finding the closest number to each 10kb without going over it
  exp_max_row <- max(exp_rows)
  #  actual_pos <- data_chrom_ordered$midpoint[max_row]
  #  rows_to_split_dat <- c(rows_to_split_dat, actual_pos)
  exp_rows_to_split_dat <- c(exp_rows_to_split_dat, exp_max_row)
}#for
#rows_to_split_dat <- rows_to_split_dat + 1
exp_rows_to_split_dat
diff <- c(exp_rows_to_split_dat[-1],NA)
tmp_df <- as.data.frame(cbind(exp_rows_to_split_dat,diff))
tmp_df$rep_vec <- tmp_df$diff - tmp_df$exp_rows_to_split_dat
head(tmp_df)
tail(tmp_df)
#make rep vector
rep_vec <- c(exp_rows_to_split_dat[1],tmp_df$rep_vec[-length(tmp_df$rep_vec)])
rep_vec
#rep_vec <- rep_vec[which(rep_vec != 0)]
#rep_vec
#make new column with "groups" of these chunks
new_sections <- rep(chunks, times = rep_vec)
sel_dat_ord$new_sections <- new_sections
head(sel_dat_ord)

#mean for each genome chunk
median_sel_dat <- sel_dat_ord %>%
  dplyr::group_by(bacteria,class,new_sections) %>%
  dplyr::summarise(value.mean = mean(value), value.sd = sd(value))
median_sel_dat
median_sel_dat <- tidyr::drop_na(median_sel_dat)

#plot just ecoli
ecol_dat <- median_sel_dat[which(median_sel_dat$bacteria == "strep"),]
head(ecol_dat)
#scale the genomic position
ecol_dat$new_sections <- ecol_dat$new_sections / 1000000


## The errorbars overlapped, so use position_dodge to move them horizontally
#pd <- position_dodge(0.1) # move them .05 to the left and right
#
#ecol_dat$class <- factor(ecol_dat$class, levels = c("dS", "omega", "dN"))
#
#levels(ecol_dat$class) <- c("omega" = expression(omega),
#                            "dS" = " dS",
#                              "dN" = " dN")
#levels(ecol_dat$class)
##make fake variable so I can subset the data
#ecol_dat$fake_class <- factor(ifelse(ecol_dat$class == "omega", 'ome', 'rates'))
#
#
#
#head(ecol_dat[which(ecol_dat$class == "omega"),])
#head(ecol_dat[which(ecol_dat$class == " dN"),])
#head(ecol_dat[which(ecol_dat$class == " dS"),])
#levels(ecol_dat$fake_class)
##levels(ecol_dat$fake_class) <- factor(ecol_dat$fake_class, levels= c("omeg", "rates"))
ecol_omeg <- ecol_dat[which(ecol_dat$class == "omega"),]
head(ecol_omeg)
class(ecol_omeg)
#make new column to colour points with omega > 1
res <- ecol_omeg %>% mutate(category=cut(value.mean, breaks=c(-Inf, 1, Inf), labels=c("#C29979","#594637")))
res <- ecol_omeg %>% mutate(category=cut(value.mean, breaks=c(-Inf, 1, Inf), labels=c("#C29979","#6A5442")))

ecol_rates <- ecol_dat[which(ecol_dat$class == "dS" | ecol_dat$class == "dN"),]
head(ecol_rates)
#choose colours
#colours_arr <- c("#B0413E","#928CAB")
#colours_arr <- c("#B7524F","#A09ABC")
#colours_arr <- c("#B7524F","#928CAB")
colours_arr <- c("#6494AA","#A09ABC")
#graph with only dS and dN
rate_g <- (ggplot(ecol_rates, aes(x=new_sections, y=value.mean, colour=class))
          #  geom_errorbar(aes(ymin=value.mean-sd, ymax=value.mean+sd), width=.1, position=pd) +
          + geom_point(alpha = 0.75)
          + geom_smooth()
          #labels for the colours
          + annotate("text", x=5.4,y=0.015,label="dS", colour = "#A09ABC", size = 6)
          + annotate("text", x=5.4,y=0.0035,label="dN", colour = "#6494AA", size = 6)
          + scale_y_continuous(trans='log10',labels = function(x) ifelse(x == 0, "0", x), breaks=c(0.001,0.1, 1, 10))
          + scale_color_manual(values=colours_arr)
          #remove redundant xaxis
          + theme(axis.title.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.text.x = element_blank())
          #axis labels
          + xlab("Distance from the Origin of Replication (Mbp)")
          + ylab("Mean Expected #\nof Substitutions per 10Kbp\n")
          + ggtitle(expression(paste(italic("Streptomyces"), " Chromosome")))
          + theme(legend.position = "none",)
)

rate_g

colours_arr <- res$category
colours_arr <- res$category
omeg_g <- (ggplot(res, aes(x=new_sections, y=value.mean))
           + geom_point(alpha = 0.75, fill = colours_arr, colour = colours_arr)
           + geom_smooth(colour = "#C29979")
           #labels for the colours
           + annotate("text", x=5.4,y=0.4,label="omega", parse=TRUE, colour = "#C29979", size = 10)
           + scale_y_continuous(trans='log10',labels = function(x) ifelse(x == 0, "0", x), breaks=c(0.001,0.01,0.1, 1, 10))
           #+ scale_color_manual(values=colours_arr)
           #omega = 1 reference line
           +  geom_hline(yintercept=1, linetype="dashed", color = "black", size=1)
           #axis labels
           + ggtitle(expression(paste(italic("Streptomyces"), " Chromosome")))
           + ylab("Mean Ratio\n(dN/dS)\nper 10Kbp")
           #remove x axis labels
           + theme(axis.title.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.text.x = element_blank(),
                   #remove legend
                   legend.position = "none")
)
omeg_g
#colours_arr <- c("#C29979")
#colours_arr <- c("#B18C6E")
#omeg_g <- (ggplot(ecol_omeg, aes(x=new_sections, y=value.mean, colour=class))
omeg_g <- (ggplot(res, aes(x=new_sections, y=value.mean))
           #  geom_errorbar(aes(ymin=value.mean-sd, ymax=value.mean+sd), width=.1, position=pd) +
           + geom_point(alpha = 0.75, fill = colours_arr, colour = colours_arr)
           + geom_smooth(colour = "#C29979")
           #labels for the colours
           + annotate("text", x=5.4,y=0.4,label="omega", parse=TRUE, colour = "#B18C6E", size = 10)
           + scale_y_continuous(trans='log10',labels = function(x) ifelse(x == 0, "0", x), breaks=c(0.001,0.1, 1, 10))
           + scale_color_manual(values=colours_arr)
           #omega = 1 reference line
           +  geom_hline(yintercept=1, linetype="dashed", color = "black", size=1)
           #axis labels
           + xlab("Distance from the Origin of Replication (Mbp)")
           + ylab("Mean Ratio (dN/dS) per 10Kbp")
           #+ ggtitle(expression(paste(italic("Streptomyces"), " Chromosome")))
           + theme(legend.position = "none")
)

omeg_g

#arrange the graphs on one. since facet will not let you re-lable each axis in a facet
grid.newpage()
grid.draw(rbind(ggplotGrob(rate_g), ggplotGrob(omeg_g)))


###get separate datasets for the rates and omega
rates_dat <- sel_dat[which(sel_dat$class == "dN" | sel_dat$class == "dS"),]
omega_dat <- sel_dat[which(sel_dat$class == "omega"),]
### separate those again by the non-multirepliconic and multirepliconic bacteria
single_rates_dat <- rates_dat[which(rates_dat$bacteria == "ecoli" | rates_dat$bacteria == "bass" | rates_dat$bacteria == "strep"),]
multi_rates_dat <- rates_dat[which(rates_dat$bacteria == "sinoC" | rates_dat$bacteria == "pSymA" | rates_dat$bacteria == "pSymB"),]

single_omega_dat <- omega_dat[which(omega_dat$bacteria == "ecoli" | omega_dat$bacteria == "bass" | omega_dat$bacteria == "strep"),]
multi_omega_dat <- omega_dat[which(omega_dat$bacteria == "sinoC" | omega_dat$bacteria == "pSymA" | omega_dat$bacteria == "pSymB"),]


#get levels into order we want
levels(single_rates_dat$bacteria)
single_rates_dat$bacteria <- fct_relevel(single_rates_dat$bacteria,"ecoli","bass","strep","sinoC","pSymA","pSymB")
levels(single_rates_dat$bacteria)
#italic bacteria names
levels(single_rates_dat$bacteria) <- c("ecoli" = expression(paste(italic("E.coli"), "")),
                              "bass" = expression(paste(italic("B. subtilis"), "")),
                              "strep" = expression(paste(italic("Streptomyces"), "")),
                              "sinoC" = expression(paste(italic("S.meliloti"), " Chromosome")),
                              "pSymA" = expression(paste(italic("S.meliloti"), " pSymA")),
                              "pSymB" = expression(paste(italic("S.meliloti"), " pSymB")))

single_omega_dat$bacteria <- fct_relevel(single_omega_dat$bacteria,"ecoli","bass","strep","sinoC","pSymA","pSymB")
#italic bacteria names
levels(single_omega_dat$bacteria) <- c("ecoli" = expression(paste(italic("E.coli"), "")),
                                       "bass" = expression(paste(italic("B. subtilis"), "")),
                                       "strep" = expression(paste(italic("Streptomyces"), "")),
                                       "sinoC" = expression(paste(italic("S.meliloti"), " Chromosome")),
                                       "pSymA" = expression(paste(italic("S.meliloti"), " pSymA")),
                                       "pSymB" = expression(paste(italic("S.meliloti"), " pSymB")))


multi_rates_dat$bacteria <- fct_relevel(multi_rates_dat$bacteria,"ecoli","bass","strep","sinoC","pSymA","pSymB")
#italic bacteria names
levels(multi_rates_dat$bacteria) <- c("ecoli" = expression(paste(italic("E.coli"), "")),
                                       "bass" = expression(paste(italic("B. subtilis"), "")),
                                       "strep" = expression(paste(italic("Streptomyces"), "")),
                                       "sinoC" = expression(paste(italic("S.meliloti"), " Chromosome")),
                                       "pSymA" = expression(paste(italic("S.meliloti"), " pSymA")),
                                       "pSymB" = expression(paste(italic("S.meliloti"), " pSymB")))

multi_omega_dat$bacteria <- fct_relevel(multi_omega_dat$bacteria,"ecoli","bass","strep","sinoC","pSymA","pSymB")
#italic bacteria names
levels(multi_omega_dat$bacteria) <- c("ecoli" = expression(paste(italic("E.coli"), "")),
                                      "bass" = expression(paste(italic("B. subtilis"), "")),
                                      "strep" = expression(paste(italic("Streptomyces"), "")),
                                      "sinoC" = expression(paste(italic("S.meliloti"), " Chromosome")),
                                      "pSymA" = expression(paste(italic("S.meliloti"), " pSymA")),
                                      "pSymB" = expression(paste(italic("S.meliloti"), " pSymB")))


#### separate facet for rates for non-multi rep bacteria
colours_arr <- rep(c("#8BC1C1","#928CAB"), num_of_plots/2)
single_rates_summary <-(ggplot(single_rates_dat, aes(x=class, y=value, fill = class, colour = class)) 
               + geom_jitter(position=position_jitter(0.2))
               + geom_violin(colour = "black", fill = NA) 
               + geom_boxplot(width=.1, outlier.shape=NA, colour = "black", fill = colours_arr) 
               + stat_boxplot(geom = "errorbar", width = 0.2, colour = "black") 
               + facet_wrap(~bacteria, labeller=label_parsed)
               + xlab("") 
               + ylab("Expected Number\nof Substitutions\nper Site") 
               + scale_color_manual(values=c("#8BC1C1","#928CAB"),labels = c(" dN", " dS"))
               + scale_x_discrete(breaks = c("dN", "dS"),labels = c("dN","dS")) 
               #log scale and removing trailing zeros from y-axis labels
               + scale_y_continuous(trans='log10',labels = function(x) ifelse(x == 0, "0", x), breaks=c(0.001,0.1, 1, 10,1000))
                 #remove all legends
                 + theme(legend.position = "none",
                         #removing space at bottom of plot
                         plot.margin=unit(c(1,1,-0.5,1), "cm"),
                         #remove space between facest
                         panel.spacing.x=unit(0, "lines"))
)

single_rates_summary

### summary graph for omega values for non-multi-rep bacteria
colours_arr <- rep(c("#B18C6E"), num_of_plots/2)
single_omega_summary <-(ggplot(single_omega_dat, aes(x=class, y=value, fill = class, colour = class)) 
                 + geom_jitter(position=position_jitter(0.2))
                 + geom_violin(colour = "black", fill = NA) 
                 + geom_boxplot(width=.1, outlier.shape=NA, colour = "black", fill = colours_arr) 
                 + stat_boxplot(geom = "errorbar", width = 0.2, colour = "black") 
                 + facet_wrap(~bacteria, labeller=label_parsed)
                 #omega = 1 reference line only on the omega data
                 +  geom_hline(data = data.frame(yint=1,fake_class="ome"), aes(yintercept= yint), linetype="dashed", color = "black")
                 + xlab("") 
                 + ylab("Ratio of dN/dS") 
                 + scale_color_manual(values=c("#C29979"),labels = c(expression(omega)))
                 #make the omega a math symbol in x-axis
                 + scale_x_discrete(breaks = c("omega"),labels = c(expression(omega))) 
                 #log scale and removing trailing zeros from y-axis labels
                 ##remove the axis ticks
                 #+ theme(axis.ticks.y = element_blank(),
                 #        axis.title.y = element_blank(),
                 #        axis.text.y = element_blank() )
                 #and have label at omega = 1 for reference
                 + scale_y_continuous(trans='log10',labels = function(x) ifelse(x == 0, "0", x), breaks=c(0.001,0.1, 1, 10,1000))
                 #+ scale_y_continuous(trans='log10',labels = function(x) ifelse(x == 0, "0", x), breaks=c(0.001,0.1, 1, 10,1000),
                  #                  #adding a second axis for omega!
                  #                    sec.axis = sec_axis(trans = ~ . * 1,
                  #                                       name = 'Ratio of dN/dS \n',labels = function(x) ifelse(x == 0, "0", x), breaks=c(0.001,0.1, 1, 10,1000)))
                 #remove all legends and facet labels (redundant)
                 + theme(legend.position = "none",
                         strip.background = element_blank(),
                         strip.text.x = element_blank(),
                         #specify plot margins
                         plot.margin=unit(c(0.1,1,0.5,1), "cm"))
)

single_omega_summary

#arrange the graphs for non-multi rep bac
grid.newpage()
grid.draw(rbind(ggplotGrob(single_rates_summary), ggplotGrob(single_omega_summary)))


##colours_arr <- c("#CFE7C8","#D2E4DC","#A0747A")
##colours_arr <- c("#7A306C","#8E8DBE","#5EB26D")
##colours_arr <- c("#A0747A","#5EB26D","#8E8DBE")
#colours_arr <- c("#5EB26D","#A0747A","#8E8DBE")
#span = 0.3
#distg <- (ggplot(ecol_dat, aes(x=new_sections, y=value.mean, colour=class))
##  geom_errorbar(aes(ymin=value.mean-sd, ymax=value.mean+sd), width=.1, position=pd) +
#  + geom_point()
#  + geom_smooth(method = lm)
#  #labels for the colours
#  + annotate("text", x=5.4,y=0.4,label="omega", parse=TRUE, colour = "#A0747A", size = 10)
#  + annotate("text", x=5.4,y=0.015,label="dS", colour = "#5EB26D", size = 6)
#  + annotate("text", x=5.4,y=0.0035,label="dN", colour = "#8E8DBE", size = 6)
#  + scale_y_continuous(trans='log10',labels = function(x) ifelse(x == 0, "0", x), breaks=c(0.001,0.1, 1, 10))
#  #omega = 1 reference line
#  +  geom_hline(yintercept=1, linetype="dashed", color = "black", size=1)
#  + scale_color_manual(values=colours_arr)
#  #axis labels
#  + xlab("Distance from the Origin of Replication (Mbp)")
#  + ylab("Mean Expected Number of Substitutions per 10Kbp")
#  + ggtitle(expression(paste(italic("Streptomyces"), " Chromosome")))
#  + theme(legend.position = "none")
#)
#
#pdf("strep_selection.pdf")
#distg
#
#  
#dev.off()
#
#
#(distg + facet_grid(fake_class ~ .)
#  #(distg + facet_grid(fake_class ~ . ,scales = "free_x", space = "free_x",labeller=labeller(treatment = labels))
#  #removing the facet labels
#  + theme(
#    strip.background = element_blank(),
#    #    strip.text.y  = element_blank()
#  )
#)

strep_tmp <- median_sel_dat[which(median_sel_dat$bacteria == "strep"),]
levels(strep_tmp$class)
strep_tmp$fake_class <- factor(ifelse(strep_tmp$class == "omega", 'ome', 'rates'))
levels(strep_tmp$fake_class)
class(strep_tmp$fake_class)
strep_tmp$fake_class <- fct_relevel(strep_tmp$fake_class,"rates","ome")
levels(strep_tmp$fake_class)
head(strep_tmp)
strep_tmp$new_sections <- strep_tmp$new_sections / 1000000


colours_arr <- c("#6494AA","#A09ABC", "black")
#graph with only dS and dN
facet_g <- (ggplot(strep_tmp, aes(x=new_sections, y=value.mean, colour=class))
           #  geom_errorbar(aes(ymin=value.mean-sd, ymax=value.mean+sd), width=.1, position=pd) +
           + geom_point(alpha = 0.75)
           + geom_smooth()
           + facet_grid(fake_class ~ ., labeller=label_parsed, scales='free')
           #labels for the colours
           + annotate("text", x=5.4,y=0.015,label="dS", colour = "#A09ABC", size = 6)
           + annotate("text", x=5.4,y=0.0035,label="dN", colour = "#6494AA", size = 6)
           + scale_y_continuous(trans='log10',labels = function(x) ifelse(x == 0, "0", x), breaks=c(0.001,0.1, 1, 10))
           + scale_color_manual(values=colours_arr)
           #axis labels
           + xlab("Distance from the Origin of Replication (Mbp)")
           + ylab("Mean Expected Number of Substitutions per 10Kbp")
           + ggtitle(expression(paste(italic("Streptomyces"), " Chromosome")))
           + theme(legend.position = "none")
)

facet_g



#get levels into order we want
levels(sel_dat$bacteria)
sel_dat$bacteria <- fct_relevel(sel_dat$bacteria,"ecoli","bass","strep","sinoC","pSymA","pSymB")
levels(sel_dat$bacteria)
#italic bacteria names
levels(sel_dat$bacteria) <- c("ecoli" = expression(paste(italic("E.coli"), "")),
                              "bass" = expression(paste(italic("B. subtilis"), "")),
                              "strep" = expression(paste(italic("Streptomyces"), "")),
                              "sinoC" = expression(paste(italic("S.meliloti"), "\n Chromosome")),
                              "pSymA" = expression(paste(italic("S.meliloti"), " pSymA")),
                              "pSymB" = expression(paste(italic("S.meliloti"), " pSymB")))



class(sel_dat$bacteria)
levels(sel_dat$bacteria)

#make a fake factor to separate dS and dN from omega
#so it will facet the way I want
tmp_data <- sel_dat
levels(tmp_data$class)
tmp_data$fake_class <- factor(ifelse(tmp_data$class == "omega", 'ome', 'rates'))
tmp_data$fake_class <- fct_relevel(tmp_data$fake_class, "rates", "ome")
levels(tmp_data$fake_class)
class(tmp_data$fake_class)
#levels(tmp_data$fake_class) <- factor(tmp_data$fake_class, levels= c("omeg", "rates"))

#another fake factor
tmp_data$bac_class <- factor(paste(tmp_data$bacteria, tmp_data$fake_class, sep=""))
tmp_data$bac_class <- fct_relevel(tmp_data$bac_class, "ecolirates", "ecoliome",
                                  "bassrates", "bassome",
                                  "streprates", "strepome",
                                  "sinoCrates", "sinoCome",
                                  "pSymArates", "pSymAome",
                                  "pSymBrates", "pSymBome")
class(tmp_data$bac_class)
levels(tmp_data$bac_class)


#plot it
#choose colours
#colours_arr <- rep(c("#5EB26D","#8E8DBE","#A0747A"),num_of_plots)
#colours_arr <- rep(c("#8E8DBE","#89C794","#A0747A"),num_of_plots)
colours_arr <- rep(c( "#8BC1C1","#928CAB","#C29979"),num_of_plots)
#rates_arr <- rep(c( "#8BC1C1","#928CAB"),num_of_plots)
#omeg_arr <- rep(c("#C29979"), num_of_plots)
#colours_arr <- c(rates_arr, omeg_arr)
#plot
set.seed(1738);
vio_str_box <-(ggplot(tmp_data, aes(x=class, y=value, fill = class, colour = class)) 
               + geom_jitter(position=position_jitter(0.2))
               + geom_violin(colour = "black", fill = NA) 
               #+ geom_boxplot(width=.1, outlier.shape=NA, fill = colours_arr) 
               + geom_boxplot(width=.1, outlier.shape=NA, colour = "black", fill = colours_arr) 
               + stat_boxplot(geom = "errorbar", width = 0.2, colour = "black") 
               #+ facet_grid(fake_class ~bacteria, labeller=label_parsed)
               + facet_wrap(~bac_class, labeller=label_parsed, scales = "free_x", space = "free_x", ncol = 6)
               #omega = 1 reference line only on the omega data
               #+  geom_hline(yintercept=1, linetype="dashed", color = "black")
               +  geom_hline(data = data.frame(yint=1,fake_class="ome"), aes(yintercept= yint), linetype="dashed", color = "black")
               + xlab("") 
               + ylab("Expected Number of Substitutions per Site") 
               #+ scale_color_manual(values=c("#8E8DBE","#89C794","#A0747A"),labels = c(" dN", " dS", expression(omega)))
               #+ scale_color_manual(values=c("#8BC1C1","#928CAB","#BE6361"),labels = c(" dN", " dS", expression(omega)))
               + scale_color_manual(values=c("#8BC1C1","#928CAB","#C29979"),labels = c(" dN", " dS", expression(omega)))
               #make the omega a math symbol in x-axis
               + scale_x_discrete(breaks = c("dN", "dS", "omega"),labels = c("dN","dS", expression(omega))) 
               #log scale and removing trailing zeros from y-axis labels
               #and have label at omega = 1 for reference
               + scale_y_continuous(trans='log10',labels = function(x) ifelse(x == 0, "0", x), breaks=c(0.001,0.1, 1, 10,1000),
                                    #adding a second axis for omega!
                                    sec.axis = sec_axis(trans = ~ . * 1,
                                                        name = 'Ratio of dN/dS \n',labels = function(x) ifelse(x == 0, "0", x), breaks=c(0.001,0.1, 1, 10,1000)))
               #remove weird second legend
               + guides(fill=FALSE)
               #remove the y facet labels
               + theme(
                 strip.background.y = element_blank(),
                 strip.text.y = element_blank()
               #  strip.background = element_blank(),
               #  strip.text = element_blank()
               )
)

vio_str_box

strwrap_strip_text(vio_str_box)




#choose colours
#colours_arr <- rep(c("#5EB26D","#8E8DBE","#A0747A"),num_of_plots)
#colours_arr <- rep(c("#8E8DBE","#89C794","#A0747A"),num_of_plots)
colours_arr <- rep(c( "#8BC1C1","#928CAB","#C29979"),num_of_plots)
#plot
set.seed(1738);
vio_str_box <-(ggplot(sel_dat, aes(x=class, y=value, fill = class, colour = class)) 
               + geom_jitter(position=position_jitter(0.2))
               + geom_violin(colour = "black", fill = NA) 
               #+ geom_boxplot(width=.1, outlier.shape=NA, fill = colours_arr) 
               + geom_boxplot(width=.1, outlier.shape=NA, colour = "black", fill = colours_arr) 
               + stat_boxplot(geom = "errorbar", width = 0.2, colour = "black") 
               #omega = 1 reference line
               +  geom_hline(yintercept=1, linetype="dashed", color = "black")
               + facet_wrap(~bacteria, labeller=label_parsed)
               + xlab("") 
               + ylab("Expected Number of Substitutions per Site") 
               #+ scale_color_manual(values=c("#8E8DBE","#89C794","#A0747A"),labels = c(" dN", " dS", expression(omega)))
               #+ scale_color_manual(values=c("#8BC1C1","#928CAB","#BE6361"),labels = c(" dN", " dS", expression(omega)))
               + scale_color_manual(values=c("#8BC1C1","#928CAB","#C29979"),labels = c(" dN", " dS", expression(omega)))
               #make the omega a math symbol in x-axis
               + scale_x_discrete(breaks = c("dN", "dS", "omega"),labels = c("dN","dS", expression(omega))) 
               #log scale and removing trailing zeros from y-axis labels
               #and have label at omega = 1 for reference
               + scale_y_continuous(trans='log10',labels = function(x) ifelse(x == 0, "0", x), breaks=c(0.001,0.1, 1, 10,1000),
                                    #adding a second axis for omega!
                                    sec.axis = sec_axis(trans = ~ . * 1,
                                                        name = 'Ratio of dN/dS \n',labels = function(x) ifelse(x == 0, "0", x), breaks=c(0.001,0.1, 1, 10,1000)))
               #remove weird second legend
               + guides(fill=FALSE)
)
pdf("selection_vio_box.pdf")
vio_str_box

dev.off()


#for (i in (length(exp_rows_to_split_dat)-1)){
#  tmp <- exp_rows_to_split_dat[i+1] - exp_rows_to_split_dat[i]
#  rep_vec <- c(rep_vec,tmp)
#}
#rep_vec
##make new column with "groups" of these chunks
#new_sections <- rep(chunks, times = exp_rows_to_split_dat)
#head(new_sections)                    
#length(new_sections)
#sel_dat_ord$new_sections <- new_sections
#
#
## Load in data
##list of files
#file_list <- list.files(path = "./data/",
#                        pattern = "*_selection_data.csv",
#                        full.names = T)
#selection_dat <- tibble(bacteria = file_list) %>% # create a data frame
#  # holding the file names
#  mutate(file_contents = map(bacteria,          # read files into
#                             ~ read_csv(.))) # a new data column
#selection_df <- unnest(selection_dat, cols = c(file_contents))
#
##change the file names to bacteria names
#selection_df <- selection_df %>%
#  mutate_if(is.character, str_replace_all, pattern = '^./data/bass_selection_data.csv', replacement = 'bass')
#selection_df <- selection_df %>%
#  mutate_if(is.character, str_replace_all, pattern = '^./data/ecoli_selection_data.csv', replacement = 'ecoli')
#selection_df <- selection_df %>%
#  mutate_if(is.character, str_replace_all, pattern = '^./data/strep_selection_data.csv', replacement = 'strep')
#selection_df <- selection_df %>%
#  mutate_if(is.character, str_replace_all, pattern = '^./data/sinoC_selection_data.csv', replacement = 'sinoC')
#selection_df <- selection_df %>%
#  mutate_if(is.character, str_replace_all, pattern = '^./data/pSymA_selection_data.csv', replacement = 'pSymA')
#selection_df <- selection_df %>%
#  mutate_if(is.character, str_replace_all, pattern = '^./data/pSymB_selection_data.csv', replacement = 'pSymB')
#
##combine dN, dS, and omega data into one column for eventual plotting
#selection_df <- selection_df %>% 
#  gather(sel_class, val, dN:omega) %>%
#  mutate(sel_class = gsub("sel_class", "", sel_class)) %>%
#  arrange(bacteria, gene_name, sel_class, midpoint, tmp_pos)
#
##remove NAs to avoid warnings
#selection_df_no_NA <- tidyr::drop_na(selection_df,val)
#
####################################
##########################
##chunking the data into sections of the genome
##########################
##decreasing order
#selection_df_medians <- arrange(selection_df_no_NA,tmp_pos)
#nmax_pos <- max(selection_df_medians$tmp_pos)
#nmin_pos <- min(selection_df_medians$tmp_pos)
#nmin_pos
##lenght of section of genome
#chunklen <- 100000
#chunk_len_of_genome <- round_any(nmax_pos, chunklen, f=ceiling)
#chunk_len_of_genome
#chunks <- seq(0, chunk_len_of_genome, chunklen)
#chunks
#
##expty vec to hold the rows that split the dat into X kb chunks
#exp_rows_to_split_dat <- vector()
#for (i in chunks) {
#  exp_rows <-
#    which(abs(selection_df_medians$tmp_pos-i)==min(abs(selection_df_medians$tmp_pos-i)))
#  # finding the closest number to each 10kb without going over it
#  exp_max_row <- max(exp_rows)
#  #  actual_pos <- data_chrom_ordered$midpoint[max_row]
#  #  rows_to_split_dat <- c(rows_to_split_dat, actual_pos)
#  exp_rows_to_split_dat <- c(exp_rows_to_split_dat, exp_max_row)
#}#for
##rows_to_split_dat <- rows_to_split_dat + 1
#exp_rows_to_split_dat
##make new column with "groups" of these chunks
#new_sections <- rep(chunks, times = exp_rows_to_split_dat)
#head(new_sections)                    
#
#selection_df_medians$new_sections <- new_sections
#
#
#
##reorder by median omega
#selection_df_no_NA <- (selection_df_no_NA
#             %>% group_by(bacteria)
#             %>% mutate(omegaval=median(val[sel_class=="omega"]))
#             %>% ungroup()
#             %>% mutate(bacteria = reorder(bacteria,omegaval))
#             %>% select(-omegaval)
#)
#
##number of bacterial replicons
#num_of_plots <- length(levels(selection_df_no_NA$bacteria))
#
##italic bacteria names
#  levels(selection_df_no_NA$bacteria) <- c("ecoli" = expression(paste(italic("E.coli"), "")),
#                          "bass" = expression(paste(italic("B. subtilis"), "")),
#                          "sinoC" = expression(paste(italic("S.meliloti"), " Chromosome")),
#                          "pSymA" = expression(paste(italic("S.meliloti"), " pSymA")),
#                          "strep" = expression(paste(italic("Streptomyces"), "")),
#                          "pSymB" = expression(paste(italic("S.meliloti"), " pSymB")))
#
##choose colours
#colours_arr <- rep(c("#CFE7C8","#D2E4DC","#A0747A"),num_of_plots)
##plot
#set.seed(1738)
#vio_str_box <-(ggplot(selection_df_no_NA, aes(x=sel_class, y=val, fill = sel_class, colour = sel_class)) 
#          + geom_jitter(position=position_jitter(0.2))
#          + geom_violin(colour = "black", fill = NA) 
#          #+ geom_boxplot(width=.1, outlier.shape=NA, fill = colours_arr) 
#          + geom_boxplot(width=.1, outlier.shape=NA, colour = "black", fill = colours_arr) 
#          + stat_boxplot(geom = "errorbar", width = 0.2, colour = "black") 
#          #omega = 1 reference line
#          +  geom_hline(yintercept=1, linetype="dashed", color = "#023C40")
#          + facet_wrap(~bacteria, labeller=label_parsed)
#          + xlab("") 
#          + ylab("Value") 
#          + scale_color_manual(values=c("#CFE7C8","#D2E4DC","#A0747A"),labels = c(" dN", " dS", expression(omega)))
#          #+ scale_color_manual(values=c("#CFE7C8","#D2E4DC","#A0747A"))
#          #make the omega a math symbol in x-axis
#          + scale_x_discrete(breaks = c("dN", "dS", "omega"),labels = c("dN","dS", expression(omega))) 
#          #log scale and removing trailing zeros from y-axis labels
#          + scale_y_continuous(trans='log10',labels = function(x) ifelse(x == 0, "0", x))
#          #+ scale_fill_manual(values=c("#CFE7C8","#D2E4DC","#A0747A"), labels = c(" dN", " dS", expression(omega))) 
#)
#
#vio_str_box
#
#
##hex plot
#library(hexbin)
#omega_data <- selection_df[which(selection_df$sel_class == "omega"),]
#ecol_dat <- selection_df[which(selection_df$bacteria == "ecoli"),]
#pa_dat <- omega_data[which(omega_data$bacteria == "ecoli"),]
#
#hex_p <- (ggplot(omega_data, aes(x=tmp_pos, y=val)) 
#          + facet_wrap(~bacteria, labeller=label_parsed)
#          + geom_hex(bins=30) 
#          #omega = 1 reference line
#          +  geom_hline(yintercept=1, linetype="dashed", color = "red")
#          + xlab("Genomic Position") 
#          + ylab("Value")
#          + scale_y_continuous(trans='log10',labels = function(x) ifelse(x == 0, "0", x))
#)
#hex_p
#
#pa_hex <- (ggplot(pa_omeg, aes(x=tmp_pos, y=val)) 
#          #+ geom_hex(bins=30) 
#          + geom_point() 
#          #omega = 1 reference line
#          +  geom_hline(yintercept=1, linetype="dashed", color = "red")
#          + xlab("Genomic Position") 
#          + ylab("Value")
#          + scale_y_continuous(trans='log10',labels = function(x) ifelse(x == 0, "0", x))
#)
#pa_hex
#
#
## Scatter plot colored by groups ("Species")
#sp <- ggscatter(pa_dat, x = "tmp_pos", y = "val",
##sp <- geom_smooth(method = lm, val ~ tmp_pos, x = "tmp_pos", y = "val", data = pa_dat,
#                color = "sel_class", palette = "jco",
#                size = 3, alpha = 0.6)+
#  geom_smooth(method = lm) +
#          scale_y_continuous(trans='log10') +
#  border()                                         
## Marginal density plot of x (top panel) and y (right panel)
#xplot <- ggdensity(pa_dat, "tmp_pos", fill = "sel_class",
#                   palette = "jco")
#yplot <- ggdensity(pa_dat, "val", fill = "sel_class", 
#                   palette = "jco")+
#  #this transformation does not look right....
#          scale_y_continuous(trans='log10') +
#  rotate()
## Cleaning the plots
#sp <- sp + rremove("legend")
#yplot <- yplot + clean_theme() + rremove("legend") 
#xplot <- xplot + clean_theme() + rremove("legend")
## Arranging the plot using cowplot
#library(cowplot)
#plot_grid(xplot, NULL, sp, yplot, ncol = 2, align = "hv", 
#          rel_widths = c(2, 1), rel_heights = c(1, 2))
#
##########################
##chunking the data into sections of the genome
##########################
##adding a fake row for position 0
#pa_dat <- add_row(pa_dat, tmp_pos = 0)
#pa_dat <- arrange(pa_dat, tmp_pos)
#nmax_pos <- max(pa_dat$tmp_pos)
#nmin_pos <- min(pa_dat$tmp_pos)
#nmin_pos
##lenght of section of genome
#chunklen <- 100000
#chunk_len_of_genome <- round_any(nmax_pos, chunklen, f=ceiling)
#chunk_len_of_genome
#chunks <- seq(nmin_pos, chunk_len_of_genome, chunklen)
#chunks
#
##expty vec to hold the rows that split the dat into X kb chunks
#exp_rows_to_split_dat <- vector()
#for (i in chunks) {
#  exp_rows <-
#    which(abs(pa_dat$tmp_pos-i)==min(abs(pa_dat$tmp_pos-i)))
#  # finding the closest number to each 10kb without going over it
#  exp_max_row <- max(exp_rows)
#  #  actual_pos <- data_chrom_ordered$midpoint[max_row]
#  #  rows_to_split_dat <- c(rows_to_split_dat, actual_pos)
#  exp_rows_to_split_dat <- c(exp_rows_to_split_dat, exp_max_row)
#}#for
##rows_to_split_dat <- rows_to_split_dat + 1
#exp_rows_to_split_dat
##split data by the above specified rows
##sometimes the closest number was technically in the next 10kb chunk of seq
##but its fine.
#data_just_exp <- pa_dat$val
##list_dat_sets <- split(data_chrom_ordered, findInterval(1:nrow(data_chrom_ordered), rows_to_split_dat))
#list_dat_sets_just_exp <- split(data_just_exp, findInterval(1:nrow(pa_dat), exp_rows_to_split_dat))
#list_dat_sets_just_exp
#
##get median expression in each of the sections
#list_total_exp_level <- lapply(list_dat_sets_just_exp, median)
#print("list_total_exp_level")
#list_total_exp_level
##as df
#median_gene_exp_10kb <- as.data.frame(matrix(unlist(list_total_exp_level), byrow = F))
#median_gene_exp_10kb <- cbind(median_gene_exp_10kb, new_pos)
#write.csv(median_gene_exp_10kb)
##write.table(median_gene_exp_10kb, 'median_10kb_exp_data.csv', sep = "\t")
#
#(ggplot(data = median_gene_exp_10kb, aes(x = new_pos, y = V1))
#+ geom_point()
#  )
#
####################################
##########################
##chunking the data into sections of the genome
##########################
##decreasing order
#selection_df_medians <- arrange(selection_df_no_NA,tmp_pos)
#nmax_pos <- max(selection_df_medians$tmp_pos)
#nmin_pos <- min(selection_df_medians$tmp_pos)
#nmin_pos
##lenght of section of genome
#chunklen <- 100000
#chunk_len_of_genome <- round_any(nmax_pos, chunklen, f=ceiling)
#chunk_len_of_genome
#chunks <- seq(0, chunk_len_of_genome, chunklen)
#chunks
#
##expty vec to hold the rows that split the dat into X kb chunks
#exp_rows_to_split_dat <- vector()
#for (i in chunks) {
#  exp_rows <-
#    which(abs(selection_df_medians$tmp_pos-i)==min(abs(selection_df_medians$tmp_pos-i)))
#  # finding the closest number to each 10kb without going over it
#  exp_max_row <- max(exp_rows)
#  #  actual_pos <- data_chrom_ordered$midpoint[max_row]
#  #  rows_to_split_dat <- c(rows_to_split_dat, actual_pos)
#  exp_rows_to_split_dat <- c(exp_rows_to_split_dat, exp_max_row)
#}#for
##rows_to_split_dat <- rows_to_split_dat + 1
#exp_rows_to_split_dat
##make new column with "groups" of these chunks
#new_sections <- rep(chunks, times = exp_rows_to_split_dat)
#new_sections                    
#
#selection_df_medians$new_sections <- new_sections
#
#
#
#