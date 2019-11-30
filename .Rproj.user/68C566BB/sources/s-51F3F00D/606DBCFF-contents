##### fake graphs for spatial genomic trends####
library(ggplot2)
require(ggplot2)
require("ggplot2")
library(grid)
library(gridBase)

x <- seq(1,10,1)
y <- seq(1,10,1)
rev_y <- rev(y)
exp <- rep("exp", 10)
mut <- rep("mut", 10)
ess <- rep("ess", 10)
cons <- rep("cons", 10)
df_mut <- data.frame(cbind(x,y,mut))
df_exp <- data.frame(cbind(x,rev_y,exp))
df_ess <- data.frame(cbind(x,rev_y,ess))
df_cons <- data.frame(cbind(x,rev_y,cons))
df_pos <- data.frame(cbind(x,y))
df_neg <- data.frame(cbind(x,rev_y))
colnames(df_mut)[3] <- "trend"
colnames(df_exp)[3] <- "trend"
colnames(df_exp)[2] <- "y"
colnames(df_cons)[2] <- "y"
colnames(df_ess)[2] <- "y"
colnames(df_ess)[3] <- "trend"
colnames(df_cons)[3] <- "trend"

df_all <- data.frame(rbind(df_ess, df_exp, df_cons, df_mut))

#ggplot(data=df_all, aes(x=x, y=x, group=1)) +
#  geom_smooth(method = "lm", se =F, formula = y~x) +
#  facet_wrap(~trend, ncol = 2) +geom_line(color="#1A936F", size=4)+
#  geom_point(pch=".", color="#1A936F") +
#  geom_line(color="#CBE896", size=2)+
#  geom_point(pch=".", color="#CBE896")
#
#  geom_point()
#  

#define colours
txt_col <- "#1D3461"
outline_col <- "#926A6F"
main_col <- "#A0747A"


mut_rate <- ggplot(data=df_pos, aes(x=x, y=x, group=1)) +
  geom_line(color= outline_col, size=4)+
#  geom_line(color="#FE5F55", size=4)+
  geom_point(pch=".", color= main_col) +
#  geom_point(pch=".", color="#FE5F55") +
#  geom_line(color="#FE5F55", size=2)+
  geom_line(color=main_col, size=2)+
#  geom_point(pch=".", color="#CBE896") +
#  geom_point(pch=".", color="#FE5F55") +
  theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(size=30)) +
  ggtitle(expression(paste(underline("Mutation Rate")))) +
  labs(x = "", y = "Mutation Rate") +
#  labs(x = "Genomic Position", y = "Mutation Rate") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_text(color = txt_col),
        axis.title.y=element_text(color = txt_col),
        axis.ticks.y=element_blank(),
#        axis.ticks.y=element_blank()) +
        plot.title = element_text(color = txt_col))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "#074F57"))


conserv <- ggplot(data=df_neg, aes(x=x, y=rev_y, group=1)) +
#  geom_line(color="#FE5F55", size=4)+
  #geom_line(color="#1A936F", size=4)+
  geom_line(color= outline_col, size=4)+
  #geom_point(pch=".", color="#1A936F") +
  geom_point(pch=".", color= main_col) +
#  geom_point(pch=".", color="#FE5F55") +
  #geom_line(color="#CBE896", size=2)+
  geom_line(color=main_col, size=2)+
#  geom_line(color="#FE5F55", size=2)+
  #geom_point(pch=".", color="#CBE896") +
  geom_point(pch=".", color= main_col) +
#  geom_point(pch=".", color="#FE5F55") +
  theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(size=30)) +
ggtitle(expression(paste(underline("Gene Conservation")))) +
#  annotate("text", x = 7, y = 9.5, label = "Gene Conservation", color = "#074F57", size = 13) +
  labs(x = "Genomic Position", y = "Conservation") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_text(color = txt_col),
        axis.title.y=element_text(color = txt_col),
        axis.ticks.y=element_blank(),
#        axis.ticks.y=element_blank()) +
        plot.title = element_text(color = txt_col))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = txt_col))

expression <- ggplot(data=df_neg, aes(x=x, y=rev_y, group=1)) +
#  geom_line(color="#FE5F55", size=4)+
  geom_line(color=outline_col, size=4)+
#  geom_point(pch=".", color="#FE5F55") +
  geom_point(pch=".", color=outline_col) +
  geom_line(color=main_col, size=2)+
#  geom_line(color="#FE5F55", size=2)+
#  geom_point(pch=".", color="#FE5F55") +
  geom_point(pch=".", color=main_col) +
  theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(size=30)) +
  ggtitle(expression(paste(underline("Gene Expression")))) +
#  annotate("text", x = 7.5, y = 9.5, label = "Gene Expression", color = "#074F57", size = 13) +
  labs(x = "Genomic Position", y = "Expression") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_text(color =txt_col ),
        axis.title.y=element_text(color =txt_col ),
        axis.ticks.y=element_blank(),
#        axis.ticks.y=element_blank()) +
        plot.title = element_text(color = txt_col))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = txt_col))

dosage <- ggplot(data=df_neg, aes(x=x, y=rev_y, group=1)) +
  geom_line(color=outline_col, size=4)+
  #  geom_line(color="#FE5F55", size=4)+
  geom_point(pch=".", color=outline_col) +
  #  geom_point(pch=".", color="#FE5F55") +
  geom_line(color=main_col, size=2)+
  #  geom_line(color="#FE5F55", size=2)+
  geom_point(pch=".", color=main_col) +
  #  geom_point(pch=".", color="#FE5F55") +
  theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(size=30)) +
  ggtitle(expression(paste(underline("Gene Dosage")))) +
  labs(x = "", y = "Dosage") +
  #  labs(x = "Genomic Position", y = "Essentiality") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
                axis.title.x=element_text(color = txt_col),
                axis.title.y=element_text(color = txt_col),
                axis.ticks.y=element_blank(),
#        axis.ticks.y=element_blank()) +
          plot.title = element_text(color = txt_col))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = txt_col))








essential <- ggplot(data=df_neg, aes(x=x, y=rev_y, group=1)) +
  geom_line(color=outline_col, size=4)+
#  geom_line(color="#FE5F55", size=4)+
  geom_point(pch=".", color=outline_col) +
#  geom_point(pch=".", color="#FE5F55") +
  geom_line(color=main_col, size=2)+
#  geom_line(color="#FE5F55", size=2)+
  geom_point(pch=".", color=main_col) +
#  geom_point(pch=".", color="#FE5F55") +
  theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(size=30)) +
  ggtitle(expression(paste(underline("Gene Essentiality")))) +
  labs(x = "", y = "Essentiality") +
#  labs(x = "Genomic Position", y = "Essentiality") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_text(color = txt_col),
        axis.title.y=element_text(color = txt_col),
        axis.ticks.y=element_blank(),
#        axis.ticks.y=element_blank()) +
        plot.title = element_text(color = txt_col))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = txt_col))

exp_results <- ggplot(data=df_neg, aes(x=x, y=rev_y, group=1)) +
  geom_line(color=outline_col, size=4)+
  geom_point(pch=".", color=outline_col) +
  geom_line(color=main_col, size=2)+
  geom_point(pch=".", color=main_col) +
  theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(size=30)) +
  ggtitle("Gene Expression") +
  labs(x = "Genomic Position", y = "Expression") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_text(color = txt_col),
        axis.title.y=element_text(color = txt_col),
        axis.ticks.y=element_blank(),
        plot.title = element_text(color = txt_col))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = txt_col))
pdf("exp_results_graph.pdf")
expression
dev.off()

sub_results <- ggplot(data=df_neg, aes(x=x, y=rev_y, group=1)) +
  geom_line(color=outline_col, size=4)+
  geom_point(pch=".", color=outline_col) +
  geom_line(color=main_col, size=2)+
  geom_point(pch=".", color=main_col) +
  theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(size=30)) +
  ggtitle(expression(paste(underline("Substitutions")))) +
  labs(x = "", y = "Number of Substitutions") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_text(color = txt_col),
        axis.title.y=element_text(color =txt_col),
        axis.ticks.y=element_blank(),
        plot.title = element_text(color = txt_col))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = txt_col))
pdf("sub_results_graph.pdf")
sub_results
dev.off()

pdf("exp_graph.pdf")
expression
dev.off()


pdf("ess_graph.pdf")
essential
dev.off()


pdf("cons_graph.pdf")
conserv
dev.off()

pdf("mut_graph.pdf")
mut_rate
dev.off()


pdf("dosage_graph.pdf")
dosage
dev.off()
