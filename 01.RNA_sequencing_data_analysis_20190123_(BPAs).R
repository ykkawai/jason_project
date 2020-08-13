## 20190123 (BPC2 only) ##
##### Package and data loading #####

library(dplyr)
library(limma)
library(VennDiagram)

data_190123 <- readxl::read_excel("C:/Users/User/Desktop/RNA-sequencing/02.BPAs RNA-seq/Count (BPAs).xlsx", sheet = 2)
data_190123 ## DMSO and BPC2

# This database have two seperated exposure data.
# Therefore, I need to seperate them for data analysis.

##### Relative gene expression level (190123) #####

# BPC2_0.1 and BPC2_1, upregulation (> 5)

data_190123_BPC2_0.1_u <- select(data_190123, Geneid, Genbank, DMSO, BPC2_0.1)
data_190123_BPC2_0.1_u <- mutate(data_190123_BPC2_0.1_u, fold_change_BPC2_0.1 = BPC2_0.1/DMSO)
data_190123_BPC2_0.1_u <- na.omit(data_190123_BPC2_0.1_u)
data_190123_BPC2_0.1_u <- data_190123_BPC2_0.1_u[!is.infinite(data_190123_BPC2_0.1_u$fold_change_BPC2_0.1), ]
data_190123_BPC2_0.1_u <- filter(data_190123_BPC2_0.1_u, fold_change_BPC2_0.1 > 5)
data_190123_BPC2_0.1_u <- data_190123_BPC2_0.1_u %>%  arrange(desc(fold_change_BPC2_0.1))
data_190123_BPC2_0.1_u
data_190123_BPC2_0.1_u <- select(data_190123_BPC2_0.1_u, Geneid, Genbank, fold_change_BPC2_0.1)
data_190123_BPC2_0.1_u

data_190123_BPC2_1_u <- select(data_190123, Geneid, Genbank, DMSO, BPC2_1)
data_190123_BPC2_1_u <- mutate(data_190123_BPC2_1_u, fold_change_BPC2_1 = BPC2_1/DMSO)
data_190123_BPC2_1_u <- na.omit(data_190123_BPC2_1_u)
data_190123_BPC2_1_u <- data_190123_BPC2_1_u[!is.infinite(data_190123_BPC2_1_u$fold_change_BPC2_1), ]
data_190123_BPC2_1_u <- filter(data_190123_BPC2_1_u, fold_change_BPC2_1 > 5)
data_190123_BPC2_1_u <- data_190123_BPC2_1_u %>%  arrange(desc(fold_change_BPC2_1))
data_190123_BPC2_1_u
data_190123_BPC2_1_u <- select(data_190123_BPC2_1_u, Geneid, Genbank, fold_change_BPC2_1)
data_190123_BPC2_1_u


##### Merge all the upregulated data #####

db_totalu_190123 <- merge(data_190123_BPC2_0.1_u, data_190123_BPC2_1_u, by = c("Geneid", "Genbank"), all = T)
db_totalu_190123


##### Venndiagram with the prepared database (VennDiagram package) #####

### Group by 2 (chemicals)

# BPC2_0.1 and BPC2_1 - Upregulation

db_totalu_190123[is.na(db_totalu_190123)] <- -1 # Change all the NA values to -1 for calculation
db_totalu_190123

db_totalu_190123_BPC2_0.1 <- filter(db_totalu_190123, fold_change_BPC2_0.1 > 5)
db_totalu_190123_BPC2_0.1 <- arrange(db_totalu_190123_BPC2_0.1, Geneid)
db_totalu_190123_BPC2_0.1 <- db_totalu_190123_BPC2_0.1 %>% select(Geneid)
nrow(db_totalu_190123_BPC2_0.1) # 571

db_totalu_190123_BPC2_1 <- filter(db_totalu_190123, fold_change_BPC2_1 > 5)
db_totalu_190123_BPC2_1 <- arrange(db_totalu_190123_BPC2_1, Geneid)
db_totalu_190123_BPC2_1 <- db_totalu_190123_BPC2_1 %>% select(Geneid)
nrow(db_totalu_190123_BPC2_1) # 540

BPC2_list_u = list(db_totalu_190123_BPC2_0.1[, 1], db_totalu_190123_BPC2_1[, 1]) # Check the number of data
length(BPC2_list_u[[1]])
length(BPC2_list_u[[2]])

colors = c("blue", "green")

venn.diagram(x = list(db_totalu_190123_BPC2_0.1[, 1], db_totalu_190123_BPC2_1[, 1]),
             category.names = c("BPC2_0.1", "BPC2_1"),
             filename = "file_BPC2_u.png", 
             scaled = F,
             col = "black",
             fill = colors,
             cat.col = "black",
             cat.cex = 1,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-5, 5),
             cat.dist = c(0.055, 0.055),
             margin = 0.15)


######################################################################################################

##### Relative gene expression level (190123) #####

# BPC2_0.1 and BPC2_1, Downregulation (< 0.2)

data_190123

data_190123_BPC2_0.1_d <- select(data_190123, Geneid, Genbank, DMSO, BPC2_0.1)
data_190123_BPC2_0.1_d <- mutate(data_190123_BPC2_0.1_d, fold_change_BPC2_0.1 = BPC2_0.1/DMSO)
data_190123_BPC2_0.1_d <- na.omit(data_190123_BPC2_0.1_d)
data_190123_BPC2_0.1_d <- data_190123_BPC2_0.1_d[!is.infinite(data_190123_BPC2_0.1_d$fold_change_BPC2_0.1), ]
data_190123_BPC2_0.1_d <- filter(data_190123_BPC2_0.1_d, fold_change_BPC2_0.1 < 0.2)
data_190123_BPC2_0.1_d <- data_190123_BPC2_0.1_d %>%  arrange(desc(fold_change_BPC2_0.1))
data_190123_BPC2_0.1_d
data_190123_BPC2_0.1_d <- select(data_190123_BPC2_0.1_d, Geneid, Genbank, fold_change_BPC2_0.1)
data_190123_BPC2_0.1_d

data_190123_BPC2_1_d <- select(data_190123, Geneid, Genbank, DMSO, BPC2_1)
data_190123_BPC2_1_d <- mutate(data_190123_BPC2_1_d, fold_change_BPC2_1 = BPC2_1/DMSO)
data_190123_BPC2_1_d <- na.omit(data_190123_BPC2_1_d)
data_190123_BPC2_1_d <- data_190123_BPC2_1_d[!is.infinite(data_190123_BPC2_1_d$fold_change_BPC2_1), ]
data_190123_BPC2_1_d <- filter(data_190123_BPC2_1_d, fold_change_BPC2_1 < 0.2)
data_190123_BPC2_1_d <- data_190123_BPC2_1_d %>%  arrange(desc(fold_change_BPC2_1))
data_190123_BPC2_1_d
data_190123_BPC2_1_d <- select(data_190123_BPC2_1_d, Geneid, Genbank, fold_change_BPC2_1)
data_190123_BPC2_1_d


##### Merge all the upregulated data #####

db_totald_190123 <- merge(data_190123_BPC2_0.1_d, data_190123_BPC2_1_d, by = c("Geneid", "Genbank"), all = T)
db_totald_190123


##### Venndiagram with the prepared database (VennDiagram package) #####

### Group by 2 (chemicals)

# BPC2_0.1 and BPC2_1 - Downregulation

db_totald_190123[is.na(db_totald_190123)] <- 9999999 # Change all the NA values to 9999999 for calculation
db_totald_190123

db_totald_190123_BPC2_0.1 <- filter(db_totald_190123, fold_change_BPC2_0.1 < 0.2)
db_totald_190123_BPC2_0.1 <- arrange(db_totald_190123_BPC2_0.1, Geneid)
db_totald_190123_BPC2_0.1 <- db_totald_190123_BPC2_0.1 %>% select(Geneid)
nrow(db_totald_190123_BPC2_0.1) # 3315

db_totald_190123_BPC2_1 <- filter(db_totald_190123, fold_change_BPC2_1 < 0.2)
db_totald_190123_BPC2_1 <- arrange(db_totald_190123_BPC2_1, Geneid)
db_totald_190123_BPC2_1 <- db_totald_190123_BPC2_1 %>% select(Geneid)
nrow(db_totald_190123_BPC2_1) # 3147

BPC2_list_d = list(db_totald_190123_BPC2_0.1[, 1], db_totald_190123_BPC2_1[, 1]) # Check the number of data
length(BPC2_list_d[[1]])
length(BPC2_list_d[[2]])

colors = c("blue", "green")

venn.diagram(x = list(db_totald_190123_BPC2_0.1[, 1], db_totald_190123_BPC2_1[, 1]),
             category.names = c("BPC2_0.1", "BPC2_1"),
             filename = "file_BPC2_d.png", 
             scaled = F,
             col = "black",
             fill = colors,
             cat.col = "black",
             cat.cex = 1,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-5, 5),
             cat.dist = c(0.055, 0.055),
             margin = 0.15)
