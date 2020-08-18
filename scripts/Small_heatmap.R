library(ggtree)

tree_path <- "analysis/snippy/AVC171_chromosome/fasttree/AVC171_chromosome.clean.fullcore.tree"
abricate_path <- "analysis/abricate/ST95_all/abricate.txt"      
pointfinder_path <- "analysis/pointfinder/all/ST95_all_pointfinder.txt"
ColV_data_path <- "analysis/abricate/ST95_all/colV_abricate.txt"
output_name <- "HC50_1106"
refname <- "AVC171"
cgMLST_path <- "metadata/curated_metadata_all.txt"

source("scripts/abricateR.R")

abricateR(
        file = abricate_path,
        output = output_name,
        identity = 90,
        length = 90,
        writecsv = TRUE,
        pointfinder_data = pointfinder_path,
        ColV_Liu_data = ColV_data_path
)

#Read in the tree file
tree <-
        read.tree(file = tree_path)

#trim the names of the assemblies in the tree tip labels
tree$tip.label <- gsub("\\..*", "", tree$tip.label)

#reassign strain SG17-135's assembly barcode to its strain name (tree tip label)
tree$tip.label <- gsub("Reference", refname, tree$tip.label)

#Read in cgMLST data
Metadata <- read_delim(cgMLST_path,
                       "\t",
                       escape_double = FALSE,
                       trim_ws = TRUE)

#Add new column working_name
Metadata$working_name <- Metadata$`Assembly barcode`

#replace the assembly barcode with strain name
Metadata$working_name <-
        gsub("ESC_SA8243AA_AS", "AVC171", Metadata$working_name)

#Filter metadata table to only contain strains from the tree
Metadata <- Metadata %>% filter(working_name %in% tree$tip.label)

#Change the column order to put working_name first
Metadata <- Metadata %>% select(working_name, everything())

#set rowname to be equal to working name
rownames(Metadata) <- Metadata$working_name

#Remove spaces from column names for metadata - this usually causes issues
colnames(Metadata) <- gsub(" ", "_", colnames(Metadata))

#Provides clade numbers to colour by.
#If you dont know what the node labels are use the line below under "get node labels"
#You may need to change the sizing..
#tree2 <- groupClade(tree, c(87, 86))

Metadata$Flag <- Metadata$Country
Metadata$Flag <- gsub("Australia", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Flags/Australia.jpg", Metadata$Flag)
Metadata$Flag <- gsub("Denmark", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Flags/Denmark.jpg", Metadata$Flag)
Metadata$Flag <- gsub("United States", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Flags/USA.png", Metadata$Flag)
Metadata$Flag <- gsub("^Ireland$", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Flags/Ireland.png", Metadata$Flag)
Metadata$Flag <- gsub("Vietnam", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Flags/Vietnam.png", Metadata$Flag)
Metadata$Flag <- gsub("United Kingdom", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Flags/United_Kindgom.png", Metadata$Flag)
Metadata$Flag <- gsub("Mexico", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Flags/Mexico.png", Metadata$Flag)
Metadata$Flag <- gsub("Ghana", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Flags/Ghana.png", Metadata$Flag)
Metadata$Flag <- gsub("United Arab Emirates", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Flags/UAE.png", Metadata$Flag)
Metadata$Flag <- gsub("Scotland", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Flags/Scotland.png", Metadata$Flag)
Metadata$Flag <- gsub("Germany", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Flags/Germany.png", Metadata$Flag)
Metadata$Flag <- gsub("Northern Ireland", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Flags/Northern_Ireland.jpg", Metadata$Flag)
Metadata$Flag <- gsub("Japan", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Flags/japan.png", Metadata$Flag)
Metadata$Flag <- gsub("Sweden", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Flags/sweden.png", Metadata$Flag)
Metadata$Flag <- gsub("Norway", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Flags/norway.png", Metadata$Flag)
Metadata$Flag <- gsub("China", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Flags/china.png", Metadata$Flag)
Metadata$Flag <- gsub("Nepal", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Flags/nepal.png", Metadata$Flag)
Metadata$Flag <- gsub("Saudi Arabia", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Flags/saudi-arabia.png", Metadata$Flag)
Metadata$Flag <- gsub("Netherlands", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Flags/netherlands.png", Metadata$Flag)
Metadata$Flag <- gsub("New Zealand", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Flags/new-zealand.png", Metadata$Flag)
Metadata$Flag <- gsub("Hungary", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Flags/hungary.png", Metadata$Flag)
Metadata$Flag <- gsub("France", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Flags/france.png", Metadata$Flag)
Metadata$Flag <- gsub("Singapore", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Flags/singapore.png", Metadata$Flag)
Metadata$Flag <- gsub("Croatia", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Flags/croatia.png", Metadata$Flag)
Metadata$Flag <- gsub("Sri Lanka", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Flags/sri-lanka.png", Metadata$Flag)
Metadata$Flag <- gsub("Finland", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Flags/finland.png", Metadata$Flag)
Metadata$Flag <- gsub("Canada", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Flags/canada.png", Metadata$Flag)
Metadata$Flag <- gsub("India", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Flags/india.png", Metadata$Flag)
Metadata$Flag <- gsub("Italy", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Flags/italy.png", Metadata$Flag)




Metadata$Flag[is.na(Metadata$Flag)] <- "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/question.png"

#
Metadata$Classification_img <- Metadata$Classification
Metadata$Classification_img <- gsub("SEPEC", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/blood.png", Metadata$Classification_img)
Metadata$Classification_img[is.na(Metadata$Classification_img)] <- "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/question.png"
Metadata$Classification_img <- gsub("Faecal", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/colon.png", Metadata$Classification_img)
Metadata$Classification_img <- gsub("RMAE", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/meat.png", Metadata$Classification_img)
Metadata$Classification_img <- gsub("APEC", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/poultry.png", Metadata$Classification_img)
Metadata$Classification_img <- gsub("UPEC", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/urine.png", Metadata$Classification_img)
Metadata$Classification_img <- gsub("ExPEC", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/ExPEC.png", Metadata$Classification_img)
Metadata$Classification_img <- gsub("Environmental", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/earth.png", Metadata$Classification_img)
#
Metadata$Pathogen <- Metadata$Classification
Metadata$Pathogen <- gsub("ExPEC", "Systemic", Metadata$Pathogen)
Metadata$Pathogen <- gsub("APEC", "Systemic", Metadata$Pathogen) 
Metadata$Pathogen <- gsub("UPEC", "Urine", Metadata$Pathogen) 
Metadata$Pathogen <- gsub("SEPEC", "Systemic", Metadata$Pathogen) 
Metadata$Pathogen <- gsub("RMAE", "Flora", Metadata$Pathogen) 
Metadata$Pathogen <- gsub("Faecal", "Flora", Metadata$Pathogen)
#
Metadata$Pathogen_img <- Metadata$Pathogen
Metadata$Pathogen_img <- gsub("Flora", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/colon.png", Metadata$Pathogen_img)
Metadata$Pathogen_img <- gsub("Systemic", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/blood.png", Metadata$Pathogen_img)
Metadata$Pathogen_img <- gsub("Urine", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/urine.png", Metadata$Pathogen_img)
Metadata$Pathogen_img[is.na(Metadata$Pathogen_img)] <- "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/question.png"
#
Metadata$Revised_Source_Niche_img <- Metadata$Revised_Source_Niche
Metadata$Revised_Source_Niche_img <- gsub("Canine", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/dog.png", Metadata$Revised_Source_Niche_img)
Metadata$Revised_Source_Niche_img <- gsub("Poultry", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/poultry.png", Metadata$Revised_Source_Niche_img)
Metadata$Revised_Source_Niche_img <- gsub("Human", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/human.png", Metadata$Revised_Source_Niche_img)
Metadata$Revised_Source_Niche_img[is.na(Metadata$Revised_Source_Niche_img)] <- "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/question.png"

##Designate strains as HC5-4181 or Other
Metadata$HC50_or_other <- gsub("^1106$*", "HC50-1106", Metadata$HC50)
Metadata$HC50_or_other <-
        gsub("^[0-9].*", "Other", Metadata$HC50_or_other)

df <- HC50_1106_simple_summary_N90L90

df <- df %>% select(ColV, everything(), -name)

df$ColV <- as.character(df$ColV)

df[df > 1] <- 1

colsum <- cbind(colnames(df),colSums(df)) %>% as.data.frame()

colsum$V2 <- as.numeric(colsum$V2)

df <- data.frame(lapply(df, as.character), stringsAsFactors=FALSE)

#Remove unwanted AMR genes from card
df <- df %>% select(-starts_with("card_acr"),
              -starts_with("card_bac"),
              -starts_with("card_bae"),
              -starts_with("card_cpxA"),
              -starts_with("card_CRP"),
              -starts_with("card_emr"),
              -starts_with("card_eptA"),
              -starts_with("card_Escherichia_coli_acrA"),
              -starts_with("card_Escherichia_coli_amp"),
              -starts_with("card_Escherichia_coli_emrE"),
              -starts_with("card_Escherichia_coli_mdfA"),
              -starts_with("card_evg"),
              -starts_with("card_gad"),
              -starts_with("card_H-NS"),
              -starts_with("card_kdpE"),
              -starts_with("card_marA"),
              -starts_with("card_mdt"),
              -starts_with("card_msbA"),
              -starts_with("card_pmrF"),
              -starts_with("card_tolC"),
              -starts_with("card_ugd"),
              -starts_with("card_yoj"))

df <- df %>% select(-starts_with("EC_custom_malX"),
              -starts_with("EC_custom_VGI"),
              -starts_with("EC_custom_malX"),
              -starts_with("EC_custom_yeeT"))

df <- df %>% select(-starts_with("ISfinder_Feb_2020"),
                    ISfinder_Feb_2020_IS26.X00011,
                    ISfinder_Feb_2020_IS15DII.M12900,
                    starts_with
              
              

df$ColV[df$ColV > 0] <- 3

rownames(df) <- HC50_1106_simple_summary_N90L90$name

p <- ggtree(tree) %<+%
        Metadata +
        geom_tiplab(size = 2,
                    align = TRUE,
                    offset = 0.007,
                    aes(color = Revised_Source_Niche)) #+
#geom_tippoint(size = 3, aes = (color = HC50_or_other))

colval <- c("white", #0
            "red", #1
            "purple", #3
            #"brown", # Bovine
            "gold", # Canine
          # "green", # Environment
            "blue", # Human
           # "black", #Other
            "red" # Poultry
            # "grey" #?
)

gheatmap(
        p = p,
        data = df,
        colnames_offset_y = -0.1,
        font.size = 0.7,
        hjust = 0,
        colnames_position = "top",
        #colnames = FALSE,
        colnames_angle = 90,
        offset = 0.180,
        width = 3,
        color = 'grey'
)+ ggplot2::ylim(NA, 60) +
        theme(legend.position = "bottom",
              legend.title = element_blank(),
              legend.text = element_text(size=14),
              legend.box = "horizontal") +
        scale_fill_manual(
                aesthetics = c("colour", "fill"),
                values = colval,
                na.value = 'grey') +
        geom_tiplab(aes(image = Flag), geom="image", size = 0.01, align = TRUE, linetype = NULL, offset = 0.120
        ) +
        geom_tiplab(aes(image = Revised_Source_Niche_img), geom="image", size = 0.01, align = TRUE, linetype = NULL, offset = 0.1325
            ) +
        geom_tiplab(aes(image = Pathogen_img), geom="image", size = 0.01, align = TRUE, linetype = NULL, offset = 0.1450
                    ) 

#