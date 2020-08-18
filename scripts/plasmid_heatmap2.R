library(pheatmap)
library(ggplot2)
library(magrittr)
library(dplyr)
library(reshape2)
library(readr)
library(OneR)
library(microbenchmark)
library(tidytree)
library(ggtree)
library(magrittr)
library(dplyr)
library(readr)
library(reshape2)
library(ComplexHeatmap)
library(ggimage)

### IMPORTANT ###
### Check lines where there are hardcoded changes, such as line 21,
### if you plan to modify this script's inputs.

#Replace the variable below with the path to your SG17-135 repo
path_to_repo <-
        "/Users/maxcummins/Dropbox/Doctorate/Manuscripts/AVC171/AVC171"

#Changes working directory
setwd(path_to_repo)

path_to_tree <- "analysis/snippy/HC50_1106/fasttree/AVC171.clean.fullcore.tree"

path_to_abricate_data <- "analysis/abricate/"

#Read in the abricate genotype data sheet (small number of rows for colname reassignment)
abricate_hits <-
        read_delim(
                path_to_abricate_data,
                "\t",
                escape_double = FALSE,
                trim_ws = TRUE,
                n_max = 10
        )

#Colname reassignment
colnames(abricate_hits)[c(1, 10:11)] <-
        c("name", "perc_coverage", "perc_identity")
abricate_hits_colnames <- colnames(abricate_hits)

#Re-read in PAI abricate genotype data sheet
abricate_hits <-
        read_delim(
                path_to_abricate_data,
                "\t",
                escape_double = FALSE,
                trim_ws = TRUE,
                col_names = FALSE,
                skip = 1
        )

#Remove cases where there are multiple headers from concatenation of abricate reports
abricate_hits <- abricate_hits %>% filter(X2 != "SEQUENCE")

#Colname reassignment
colnames(abricate_hits) <- abricate_hits_colnames

#Convert percent coverage and identity to numeric type to allow filtering
abricate_hits$perc_coverage <- as.numeric(abricate_hits$perc_coverage)
abricate_hits$perc_identity <- as.numeric(abricate_hits$perc_identity)

#Filter to perc_identity > 95%
#abricate_hits <-
abricate_hits <- abricate_hits %>% filter(perc_identity > 90)

#Trim excess characters the assembly names and reassign this to rownames
abricate_hits$name <- gsub("\\..*", "", abricate_hits$name)

#Replace "SAL_HC4750AA_AS" with SG17-135
#abricate_hits$name <- gsub("SAL_HC4750AA_AS", "SG17-135", abricate_hits$name)

abricate_hits$newcov <- gsub("\\/.*","", abricate_hits$COVERAGE)
abricate_hits$length_gene <- gsub(".*\\/","", abricate_hits$COVERAGE)

abricate_hits$newcov <- gsub("-",":", abricate_hits$newcov)

new_df <- abricate_hits %>% group_by(name, GENE, length_gene) %>% filter(perc_coverage > 0.5) %>% summarise(start =paste(sort(unique(newcov)), collapse=","), end = paste(sort(unique(newcov)), collapse=",")) #%>% filter(grepl("SPI-1_", GENE))
#new_df <- abricate_hits %>% group_by(name, GENE, length_gene) %>% summarise(start =paste(sort(unique(newcov)), collapse=","), end = paste(sort(unique(newcov)), collapse=",")) #%>% filter(grepl("SPI-1_", GENE))

new_df$end <- gsub("[0-9]+:","", new_df$end)
new_df$start <- gsub(":[0-9]+","", new_df$start)

spl <-strsplit(as.character(new_df$start), ",")

hit_counts <- abricate_hits %>% filter(perc_coverage > 0.5) %>% group_by(name, GENE) %>% summarise(counts = n())

max_hits <- max(hit_counts$counts)

if(max_hits > 100){
        message("More than 100 hits detected for reference plasmid for a given strain - this script may not work because of this")
}

start_coord <- data.frame(name = new_df$name, gene = new_df$GENE,
                          length_gene = new_df$length_gene,
                          chunk1 = sapply(spl, "[", 1),
                          chunk2 = sapply(spl, "[", 2),
                          chunk3 = sapply(spl, "[", 3),
                          chunk4 = sapply(spl, "[", 4),
                          chunk5 = sapply(spl, "[", 5),
                          chunk6 = sapply(spl, "[", 6),
                          chunk7 = sapply(spl, "[", 7),
                          chunk8 = sapply(spl, "[", 8),
                          chunk9 = sapply(spl, "[", 9),
                          chunk10 = sapply(spl, "[", 10),
                          chunk11 = sapply(spl, "[", 11),
                          chunk12= sapply(spl, "[", 12),
                          chunk13 = sapply(spl, "[", 13),
                          chunk14 = sapply(spl, "[", 14),
                          chunk15 = sapply(spl, "[", 15),
                          chunk16 = sapply(spl, "[", 16),
                          chunk17 = sapply(spl, "[", 17),
                          chunk18 = sapply(spl, "[", 18),
                          chunk19 = sapply(spl, "[", 19),
                          chunk20= sapply(spl, "[", 20),
                          chunk21 = sapply(spl, "[", 21),
                          chunk22 = sapply(spl, "[", 22),
                          chunk23 = sapply(spl, "[", 23),
                          chunk24 = sapply(spl, "[", 24),
                          chunk25 = sapply(spl, "[", 25),
                          chunk26 = sapply(spl, "[", 26),
                          chunk27 = sapply(spl, "[", 27),
                          chunk28= sapply(spl, "[", 28),
                          chunk29 = sapply(spl, "[", 29),
                          chunk30 = sapply(spl, "[", 30),
                          chunk31 = sapply(spl, "[", 31),
                          chunk32 = sapply(spl, "[", 32),
                          chunk33 = sapply(spl, "[", 33),
                          chunk34 = sapply(spl, "[", 34),
                          chunk35 = sapply(spl, "[", 35),
                          chunk36= sapply(spl, "[", 36),
                          chunk37 = sapply(spl, "[", 37),
                          chunk38 = sapply(spl, "[", 38),
                          chunk39 = sapply(spl, "[", 39),
                          chunk40 = sapply(spl, "[", 40),
                          chunk41 = sapply(spl, "[", 41),
                          chunk42 = sapply(spl, "[", 42),
                          chunk43 = sapply(spl, "[", 43),
                          chunk44= sapply(spl, "[", 44),
                          chunk45 = sapply(spl, "[", 45),
                          chunk46 = sapply(spl, "[", 46),
                          chunk47 = sapply(spl, "[", 47),
                          chunk48 = sapply(spl, "[", 48),
                          chunk49 = sapply(spl, "[", 49),
                          chunk50= sapply(spl, "[", 50),
                          chunk51 = sapply(spl, "[", 51),
                          chunk52 = sapply(spl, "[", 52),
                          chunk53 = sapply(spl, "[", 53),
                          chunk54 = sapply(spl, "[", 54),
                          chunk55 = sapply(spl, "[", 55),
                          chunk56 = sapply(spl, "[", 56),
                          chunk57 = sapply(spl, "[", 57),
                          chunk58 = sapply(spl, "[", 58),
                          chunk59 = sapply(spl, "[", 59),
                          chunk60 = sapply(spl, "[", 60),
                          chunk61 = sapply(spl, "[", 61),
                          chunk62 = sapply(spl, "[", 62),
                          chunk53 = sapply(spl, "[", 53),
                          chunk64 = sapply(spl, "[", 64),
                          chunk65 = sapply(spl, "[", 65),
                          chunk66 = sapply(spl, "[", 66),
                          chunk67 = sapply(spl, "[", 67),
                          chunk68 = sapply(spl, "[", 68),
                          chunk69 = sapply(spl, "[", 69),
                          chunk70 = sapply(spl, "[", 70),
                          chunk71 = sapply(spl, "[", 71),
                          chunk72 = sapply(spl, "[", 72),
                          chunk73 = sapply(spl, "[", 73),
                          chunk74 = sapply(spl, "[", 74),
                          chunk75 = sapply(spl, "[", 75),
                          chunk76 = sapply(spl, "[", 76),
                          chunk77 = sapply(spl, "[", 77),
                          chunk78 = sapply(spl, "[", 78),
                          chunk79 = sapply(spl, "[", 79),
                          chunk80 = sapply(spl, "[", 80),
                          chunk81 = sapply(spl, "[", 81),
                          chunk82 = sapply(spl, "[", 82),
                          chunk83 = sapply(spl, "[", 83),
                          chunk84 = sapply(spl, "[", 84),
                          chunk85 = sapply(spl, "[", 85),
                          chunk86 = sapply(spl, "[", 86),
                          chunk87 = sapply(spl, "[", 87),
                          chunk88 = sapply(spl, "[", 88),
                          chunk89 = sapply(spl, "[", 89),
                          chunk90 = sapply(spl, "[", 90),
                          chunk91 = sapply(spl, "[", 91),
                          chunk92 = sapply(spl, "[", 92),
                          chunk93 = sapply(spl, "[", 93),
                          chunk94 = sapply(spl, "[", 94),
                          chunk95 = sapply(spl, "[", 95),
                          chunk96 = sapply(spl, "[", 96),
                          chunk97 = sapply(spl, "[", 97),
                          chunk98 = sapply(spl, "[", 98),
                          chunk99 = sapply(spl, "[", 99),
                          chunk100 = sapply(spl, "[", 100)
)

start_coord <- melt(start_coord, id=1:3, value.name = "start")

start_coord <- start_coord %>% select(-starts_with("variable")) 

spl <-strsplit(as.character(new_df$end), ",")
end_coord <- data.frame(name = new_df$name, gene = new_df$GENE,
                        length_gene = new_df$length_gene,
                        chunk1 = sapply(spl, "[", 1),
                        chunk2 = sapply(spl, "[", 2),
                        chunk3 = sapply(spl, "[", 3),
                        chunk4 = sapply(spl, "[", 4),
                        chunk5 = sapply(spl, "[", 5),
                        chunk6 = sapply(spl, "[", 6),
                        chunk7 = sapply(spl, "[", 7),
                        chunk8 = sapply(spl, "[", 8),
                        chunk9 = sapply(spl, "[", 9),
                        chunk10 = sapply(spl, "[", 10),
                        chunk11 = sapply(spl, "[", 11),
                        chunk12= sapply(spl, "[", 12),
                        chunk13 = sapply(spl, "[", 13),
                        chunk14 = sapply(spl, "[", 14),
                        chunk15 = sapply(spl, "[", 15),
                        chunk16 = sapply(spl, "[", 16),
                        chunk17 = sapply(spl, "[", 17),
                        chunk18 = sapply(spl, "[", 18),
                        chunk19 = sapply(spl, "[", 19),
                        chunk20= sapply(spl, "[", 20),
                        chunk21 = sapply(spl, "[", 21),
                        chunk22 = sapply(spl, "[", 22),
                        chunk23 = sapply(spl, "[", 23),
                        chunk24 = sapply(spl, "[", 24),
                        chunk25 = sapply(spl, "[", 25),
                        chunk26 = sapply(spl, "[", 26),
                        chunk27 = sapply(spl, "[", 27),
                        chunk28= sapply(spl, "[", 28),
                        chunk29 = sapply(spl, "[", 29),
                        chunk30 = sapply(spl, "[", 30),
                        chunk31 = sapply(spl, "[", 31),
                        chunk32 = sapply(spl, "[", 32),
                        chunk33 = sapply(spl, "[", 33),
                        chunk34 = sapply(spl, "[", 34),
                        chunk35 = sapply(spl, "[", 35),
                        chunk36= sapply(spl, "[", 36),
                        chunk37 = sapply(spl, "[", 37),
                        chunk38 = sapply(spl, "[", 38),
                        chunk39 = sapply(spl, "[", 39),
                        chunk40 = sapply(spl, "[", 40),
                        chunk41 = sapply(spl, "[", 41),
                        chunk42 = sapply(spl, "[", 42),
                        chunk43 = sapply(spl, "[", 43),
                        chunk44= sapply(spl, "[", 44),
                        chunk45 = sapply(spl, "[", 45),
                        chunk46 = sapply(spl, "[", 46),
                        chunk47 = sapply(spl, "[", 47),
                        chunk48 = sapply(spl, "[", 48),
                        chunk49 = sapply(spl, "[", 49),
                        chunk50= sapply(spl, "[", 50),
                        chunk51 = sapply(spl, "[", 51),
                        chunk52 = sapply(spl, "[", 52),
                        chunk53 = sapply(spl, "[", 53),
                        chunk54 = sapply(spl, "[", 54),
                        chunk55 = sapply(spl, "[", 55),
                        chunk56 = sapply(spl, "[", 56),
                        chunk57 = sapply(spl, "[", 57),
                        chunk58 = sapply(spl, "[", 58),
                        chunk59 = sapply(spl, "[", 59),
                        chunk60 = sapply(spl, "[", 60),
                        chunk61 = sapply(spl, "[", 61),
                        chunk62 = sapply(spl, "[", 62),
                        chunk53 = sapply(spl, "[", 53),
                        chunk64 = sapply(spl, "[", 64),
                        chunk65 = sapply(spl, "[", 65),
                        chunk66 = sapply(spl, "[", 66),
                        chunk67 = sapply(spl, "[", 67),
                        chunk68 = sapply(spl, "[", 68),
                        chunk69 = sapply(spl, "[", 69),
                        chunk70 = sapply(spl, "[", 70),
                        chunk71 = sapply(spl, "[", 71),
                        chunk72 = sapply(spl, "[", 72),
                        chunk73 = sapply(spl, "[", 73),
                        chunk74 = sapply(spl, "[", 74),
                        chunk75 = sapply(spl, "[", 75),
                        chunk76 = sapply(spl, "[", 76),
                        chunk77 = sapply(spl, "[", 77),
                        chunk78 = sapply(spl, "[", 78),
                        chunk79 = sapply(spl, "[", 79),
                        chunk80 = sapply(spl, "[", 80),
                        chunk81 = sapply(spl, "[", 81),
                        chunk82 = sapply(spl, "[", 82),
                        chunk83 = sapply(spl, "[", 83),
                        chunk84 = sapply(spl, "[", 84),
                        chunk85 = sapply(spl, "[", 85),
                        chunk86 = sapply(spl, "[", 86),
                        chunk87 = sapply(spl, "[", 87),
                        chunk88 = sapply(spl, "[", 88),
                        chunk89 = sapply(spl, "[", 89),
                        chunk90 = sapply(spl, "[", 90),
                        chunk91 = sapply(spl, "[", 91),
                        chunk92 = sapply(spl, "[", 92),
                        chunk93 = sapply(spl, "[", 93),
                        chunk94 = sapply(spl, "[", 94),
                        chunk95 = sapply(spl, "[", 95),
                        chunk96 = sapply(spl, "[", 96),
                        chunk97 = sapply(spl, "[", 97),
                        chunk98 = sapply(spl, "[", 98),
                        chunk99 = sapply(spl, "[", 99),
                        chunk100 = sapply(spl, "[", 100)
)

end_coord <- melt(end_coord, id=1:3, value.name = "end")

end_coord <- end_coord %>% select(-starts_with("variable")) 

coords <- start_coord

coords$end <- end_coord$end

coords <- coords[complete.cases(coords),]

unique(coords$length_gene)

coords$start <- as.numeric(coords$start)
coords$end <- as.numeric(coords$end)
coords$length_gene <- as.numeric(coords$length_gene)


coords$percentage <-  (((coords$end-coords$start)+1)/coords$length_gene)*100

test <- coords# %>% filter(name == "SAL_AB7542AA_AS", gene == "SPI-12_NC_006905_P4") %>% arrange(desc(end))

list_ <- NULL

for(sample in unique(test$name)){
        test2 <- test %>% filter(name == sample)
        for(gene_ in unique(test$gene)){
                test3 <- test2 %>% filter(gene == gene_)
                length_of_gene <- test3$length_gene[1]
                if(is.na(length_of_gene) == FALSE){
                        range_matrix <- rep(0, times = length_of_gene)
                        for(hit in 1:nrow(test3)){
                                start_ <- test3[hit,4]
                                end_ <- test3[hit,5]
                                range_matrix[start_:end_] <- 1
                                range_matrix[range_matrix > 1] <- 1
                        }
                }
                newline <- c(sample, range_matrix)
                list_ <- rbind(list_,newline)
        }
        
}

list_2 <- as.data.frame(list_, stringsAsFactors = FALSE)

base_ <- as.data.frame(list_2[,1:2])

cols <- c(2:ncol(list_2))

base_$V2 <- apply(list_2[ , cols ] , 1 , paste, collapse = "")

listy_ <- NULL

x <- 1

for(i in 1:54){
        d <- unlist(list_2[i,2:ncol(list_2)])
        bins <- split(d, ceiling(seq_along(d)/100))
        bins <- lapply(bins, as.numeric)
        binsums <- lapply(bins, sum)
        listy_ <- rbind(listy_, unlist(binsums))
        print(paste("another one:",x))
        x <- x + 1
        
}

listy_ <- rbind(listy_, rep(x = 100, times = ncol(listy_)))

rownames(listy_) <- c(list_2$V1,"AVC171")

pheatmap(listy_, cluster_cols = FALSE, fontsize_col = 2)


abc <- length(list_)/3

df <- data.frame(matrix(unlist(list_), nrow = length(unique(abricate_hits$name)), byrow=T), stringsAsFactors = F)

colnames(df) <- c("name","GENE","Coverage_percentage")

df$Coverage_percentage[is.na(df$Coverage_percentage)] <- 0

df$Coverage_percentage <- as.numeric(df$Coverage_percentage)

final_table <- dcast(df, name ~ GENE)

final_final_table <- final_table[1:nrow(final_table),2:ncol(final_table)]

final_final_table_2 <- final_final_table

final_final_table[final_final_table < 60] <- 0
final_final_table[final_final_table >= 60] <- 1

rownames(final_final_table) <- final_table$name

final_table <- final_final_table

#write.csv(final_table, "analysis/PAIs_present_absent.csv")

pheatmap(final_final_table, fontsize_row = 2)

