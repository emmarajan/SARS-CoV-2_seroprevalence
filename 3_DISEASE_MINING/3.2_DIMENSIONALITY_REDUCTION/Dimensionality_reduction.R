#Load necessary packages

library(ggplot2)
library(umap)
library(factoextra)
library(ggpubr)

################################################################################ Prepare the script

# load data
# a data frame containing patient IDs, posterior probability, sex, age, seropositivity, and arrayed ICD-10 codes
data <- read.csv("ICD10_sex_age_matrix.csv")
rnames <- data[,2]

# a data frame containing patient IDs, posterior probability, sex, age, seropositivity, and arrayed ICD-10 codes only for patients with ICD-10 entries
data_selection <- read.csv("ICD10_sex_age_matrix_onlyifICD10.csv") 
rnames_selection <- data_selection[,2]


################################################################################ Fig. S8B
### UMAP COSINE WITH AGE AND SEX

# prepare for plotting with ggplot2

groups_01 <- as.factor(data$Serostatus_numerical[1:36431])
label_01 <- unique(groups_01)
print(groups_01)
print(label_01)

custom.config = umap.defaults
custom.config$metric = 'cosine'
custom.config

set.seed(42)
umap_01 <- umap::umap(data.matrix(data[,5:202]), config=custom.config)

df_01 <- data.frame(x = umap_01$layout[,1],
                 y = umap_01$layout[,2],
                 color = groups_01)

# plot with ggplot2

gg01 <- ggplot(df_01, aes(x, y, color = groups_01)) +
  geom_point(size = 2, stroke = 0) +
  theme_bw() + 
  ylab("UMAP Dimension 2") + 
  xlab("UMAP Dimension 1") + 
  scale_color_manual(values=c("#6A8ED2", "red")) +
  labs(colour="CoV2 seropositivity")

################################################################################ Fig. S8C
### PCA WITH AGE AND SEX

# Indicate range of data that should be used for analysis

mat_data_02 <- data.matrix(data[,5:202])
print(mat_data_02)
rownames(mat_data_02) <- rnames          
print(rownames(mat_data_02))

# Use column "class" to colour label different groups

groups_02 <- as.factor(data$Serostatus_numerical[1:36431])
print(groups_02)

# Compute PCA with prcomp

set.seed(42) # reproducibility
mat_data_02.pca <- prcomp(mat_data_02, center = TRUE)
print(mat_data_02.pca)

# visualize eigenvalues

fviz_eig(mat_data_02.pca)

# summary method
summary(mat_data_02.pca)

# Illustrate using factoextra

gg02 <- fviz_pca_ind(mat_data_02.pca,
             col.ind = groups_02,
             palette = c("#6A8ED2",  "red"), 
             addEllipses = FALSE,
             legend.title = "CoV2 seropositivity", 
             repel = TRUE, 
             geom.ind = c("point")
)

################################################################################ Fig. S8D
### UMAP COSINE WITHOUT AGE/SEX

# prepare for plotting with ggplot2

groups_03 <- as.factor(data$Serostatus_numerical[1:36431])
label_03 <- unique(groups_03)
print(groups_03)
print(label_03)

custom.config = umap.defaults
custom.config$metric = 'cosine'
custom.config

set.seed(42)
umap_03 <- umap::umap(data.matrix(data[,6:202]), config=custom.config)

df_03 <- data.frame(x = umap_03$layout[,1],
                 y = umap_03$layout[,2],
                 color = groups_03)

# plot with ggplot2

gg03 <- ggplot(df_03, aes(x, y, color = groups_03)) +
  geom_point(size = 2, stroke = 0) +
  theme_bw() + 
  ylab("UMAP Dimension 2") + 
  xlab("UMAP Dimension 1") + 
  scale_color_manual(values=c("#6A8ED2", "red")) +
  labs(colour="CoV2 seropositivity")

################################################################################ Fig. S8E
### PCA WITHOUT AGE/SEX

# Indicate range of data that should be used for analysis

mat_data_04 <- data.matrix(data[,6:202])
print(mat_data_04)
rownames(mat_data_04) <- rnames                  
print(rownames(mat_data_04))

# Use column "class" to colour label different groups

groups_04 <- as.factor(data$Serostatus_numerical[1:36431])
print(groups_04)

# Compute PCA with prcomp

set.seed(42)
mat_data_04.pca <- prcomp(mat_data_04, center = TRUE)
print(mat_data_04.pca)

# visualize eigenvalues

fviz_eig(mat_data_04.pca)

# summary method

summary(mat_data_04.pca)

# Illustrate using factoextra

gg04 <- fviz_pca_ind(mat_data_04.pca,
             col.ind = groups_04,
             palette = c("#6A8ED2",  "red"),
             addEllipses = FALSE,
             legend.title = "CoV2 seropositivity", 
             repel = TRUE, 
             geom.ind = c("point")
)

################################################################################ Fig. S8F
### UMAP COSINE AFTER EXCLUSION OF PATIENTS WITHOUT ICD-10 code entries

# prepare for plotting with ggplot2

groups_05 <- as.factor(data_selection$Serostatus_numeric[1:12730])
label_05 <- unique(groups_05)
print(groups_05)
print(label_05)

custom.config = umap.defaults
custom.config$metric = 'cosine'
custom.config

set.seed(42)
umap_05 <- umap::umap(data.matrix(data_selection[,5:202]), config=custom.config)

df_05 <- data.frame(x = umap_05$layout[,1],
                    y = umap_05$layout[,2],
                    color = groups_05)

# plot with ggplot2

gg05 <- ggplot(df_05, aes(x, y, color = groups_05)) +
  geom_point(size = 2, stroke = 0) +
  theme_bw() + 
  ylab("UMAP Dimension 2") + 
  xlab("UMAP Dimension 1") + 
  scale_color_manual(values=c("#6A8ED2", "red")) +
  labs(colour="CoV2 seropositivity")

################################################################################ Fig. S8G
### UMAP EUCLIDEAN AFTER EXCLUSION OF PATIENTS WITHOUT ICD-10 code entries

# prepare for plotting with ggplot2

groups_06 <- as.factor(data_selection$Serostatus_numeric[1:12730])
label_06 <- unique(groups_06)
print(groups_06)
print(label_06)

set.seed(42)
umap_06 <- umap::umap(data.matrix(data_selection[,5:202]))

df_06 <- data.frame(x = umap_06$layout[,1],
                    y = umap_06$layout[,2],
                    color = groups_06)

# plot with ggplot2

gg06 <- ggplot(df_06, aes(x, y, color = groups_06)) +
  geom_point(size = 2, stroke = 0) +
  theme_bw() + 
  ylab("UMAP Dimension 2") + 
  xlab("UMAP Dimension 1") + 
  scale_color_manual(values=c("#6A8ED2", "red")) +
  labs(colour="CoV2 seropositivity")


################################################################################
### Bring all the plots onto one page

ggarrange(gg01, gg02, gg03, gg04, gg05, gg06, 
          labels = c("B", "C", 'D', 'E', 'F', 'G'),
          ncol = 2, nrow = 3,
          common.legend = TRUE)
