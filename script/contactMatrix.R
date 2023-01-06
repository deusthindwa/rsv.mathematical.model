
## read in matrix of effective contacts relevant to airborne infectious diseases
## https://www.nature.com/articles/s41467-020-20544-y#Sec2
## R may not show everything but has 85 rows and 85 columns
USAcontactnew <- read.csv("https://raw.githubusercontent.com/mobs-lab/mixing-patterns/main/data/contact_matrices/United_States_country_level_M_overall_contact_matrix_85.csv",header = F)

# visualize contact matrix
n.cols = 100
nice.cols <-  colorRampPalette(brewer.pal(9, "YlOrRd"))(n.cols)
heatmap(as.matrix(USAcontactnew)/sum(diag(as.matrix(USAcontactnew))), Rowv = NA, Colv = NA, scale = 'none', col = nice.cols)

## create a symmetric contact matrix
symmetricUSA <- 0.5 * (USAcontactnew + t(USAcontactnew)) #add the transposed matrix to original matrix to create symmetry

## expand the contact under 1 year old
expendUSA_1 <- matrix(data = symmetricUSA[1,1], nrow = 5,ncol = 5) #age group <1y only, repeat 5 times to create rows, hence a new matrix
expendUSA_2 <- matrix(data = as.numeric(rep(symmetricUSA[1,], 5)), byrow = T, nrow = 5, ncol = 85) #all age groups, repeat 5 times to create rows, hence a matrix
expendUSA <- cbind(expendUSA_1, expendUSA_2)
expendUSA <- rbind(expendUSA, cbind(t(expendUSA_2), as.matrix(symmetricUSA))) #overall goal is to create 5 rows and 5 columns (age groups <1y) to symmetricUSA dataset and copy upwards row 1 (<1y)

#rename the columns noting that row 6 upwards will have same values for now as has been copied upwards (see code above)
colnames(expendUSA) = c("<2m", "2-3m", "4-5m", "6-7m", "8-9m", "10-11m", "1Y",
                        rep("2-4Y", 3), rep("5-9Y", 5),
                        rep("10-19Y",10), rep("20-39Y", 20),
                        rep("40-59Y", 20), rep("60Y+", 25))

# aggregate ages into 13 age groups; average the number of contacts within each age group
ave_mat_x = aggregate( expendUSA, by = list(colnames(expendUSA)), FUN = 'sum' )
ave_mat = aggregate( t(ave_mat_x[,-1]), by = list(colnames(expendUSA)), FUN = 'mean') #except the first row which is the collomn names
rownames(ave_mat) <- ave_mat[,1] #reassign collumn names
colnames(ave_mat)[2:length(ave_mat)] <- rownames(ave_mat)

#arrange the column names
ageorder <- c("Group.1", "<2m", "2-3m", "4-5m", "6-7m", "8-9m", "10-11m", "1Y", "2-4Y", "5-9Y", "10-19Y", "20-39Y", "40-59Y", "60Y+")
contactUSA <- select(ave_mat, ageorder)

#arrange the row names
contactUSA <- 
  contactUSA %>% 
  mutate(Group.1 = factor(Group.1, levels = c("<2m","2-3m","4-5m","6-7m","8-9m","10-11m","1Y","2-4Y","5-9Y","10-19Y","20-39Y","40-59Y","60Y+"))) %>%
  arrange(Group.1)

#assign as matrix and remove the first column seeming as duplicated
contactUSA <- as.matrix(contactUSA[,-1])

## create a symmetric contact matrix (Don't understand the logic here)
contactUSA[lower.tri(contactUSA, diag = FALSE)] <- 0
contactUSA <- contactUSA + t(contactUSA)
diag(contactUSA) <- 0.5*diag(contactUSA)
rownames(contactUSA) <-  c("<2m", "2-3m", "4-5m", "6-7m", "8-9m", "10-11m", "1Y", "2-4Y", "5-9Y", "10-19Y", "20-39Y", "40-59Y", "60Y+")

# visualize contact matrix
heatmapUSA <- as.matrix(contactUSA)
heatmap(heatmapUSA/sum(diag(heatmapUSA)), Rowv=NA, Colv=NA, scale='none', col=nice.cols)

# scale it to represent per capita contact per day
heatmapUSA <- heatmapUSA*rowSums(expendUSA)[1]/rowSums(heatmapUSA)[1]
