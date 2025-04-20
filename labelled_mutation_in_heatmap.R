# libraries
library(ggplot2)
library(reshape2)
library(dplyr)
#install.packages("gplots")
library(gplots)

# conversion of three letter amino acid to one letter 
aa_conversion <- c(
  "ALA" = "A", "CYS" = "C", "ASP" = "D", "GLU" = "E", "PHE" = "F", 
  "GLY" = "G", "HIS" = "H", "ILE" = "I", "LYS" = "K", "LEU" = "L", 
  "MET" = "M", "ASN" = "N", "PRO" = "P", "GLN" = "Q", "ARG" = "R", 
  "SER" = "S", "THR" = "T", "VAL" = "V", "TRP" = "W", "TYR" = "Y"
)

# file path
main_dir <- "~/Desktop/stj/last-cfg-mgl-onlyRBD/scanning-output-RBM-from-onlyRBD/"

# .txt list
file_list <- list.files(main_dir, pattern = "\\.txt$", full.names = TRUE)

# sequence and first position
sequence <- "SNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQ"
start_position <- 438
sequence_positions <- paste0(substring(sequence, 1:nchar(sequence), 1:nchar(sequence)), 
                             start_position:(start_position + nchar(sequence) - 1))

# Picking the unique mutation for Y axis 
mutations <- c()
# empty data.frame
all_data <- data.frame()

# run for whole file
for (file in file_list) {
  # Take the position and amino acids name from the file name 
  file_name <- basename(file)
  file_amino_acid <- substr(basename(file), 4, 4)   #  4th character of file names gives the amino acid in the sequence 
  file_position <- as.numeric(substr(file_name, 6, 8))  # 505, 6-8 character position
  
  # reading file
  data <- read.table(file, header = FALSE, sep = "\t", col.names = c("Mutation", "Energy"))
  
  # Seperate original amino acid and mutation 
  data$Original <- substr(data$Mutation, 1, 3)  # First 3 original aminoacid 
  data$Mutated_to <- substr(data$Mutation, 5, nchar(data$Mutation))  # 505R, after 4th character
  data$Pozisyon <- substr(data$Mutation, 5, 7) # only position 505
  data$Mutasyon <- substr(data$Mutation, 8, 8) # mutation R
  data$Mutations <- substr(data$Mutation, 1, 8)
  
  # Long version of data
  df_long <- data.frame(
    Mutated_to = data$Mutated_to, 
    Pozisyon_mt = data$Pozisyon,
    Position = file_amino_acid,
    harf_mut = data$Mutasyon, 
    Energy = data$Energy
  )
  
  # Delete the 2nd self mutation 
  filtered <- df_long %>%
    group_by(Pozisyon_mt, harf_mut) %>%
    filter(!(Pozisyon_mt == harf_mut & row_number() > 1)) 
  
  mutations <- unique(c(mutations, filtered$harf_mut))
  
  # Merge the all data 
  all_data <- rbind(all_data, filtered)
  
}

# Creating amino acid sequence for X axis 
deneme <-all_data$Pozisyon_mt <- factor(all_data$Pozisyon_mt, levels = as.character(start_position:(start_position + 68)), labels = sequence_positions)

# Limit the energy value in between -10 and +10 
all_data$Energy <- ifelse(all_data$Energy < -10, -10, ifelse(all_data$Energy > 10, 10, all_data$Energy))

# convert it to wide format for heatmap
heatmap_data <- dcast(all_data, harf_mut ~ Pozisyon_mt, value.var = "Energy", fun.aggregate = mean)


# Select specific rows (o, e and H row)
selected_rows <- heatmap_data %>%
  filter(harf_mut %in% c("o", "e", "H"))

# Sum the all row 
summed_row <- colSums(selected_rows[,-1], na.rm = TRUE)
# delete all row (o, e ve H)
heatmap_data <- heatmap_data %>%
  filter(!harf_mut %in% c("o", "e", "H"))
# add new row
heatmap_data <- rbind(heatmap_data, c(harf_mut = "H", as.list(summed_row)))


# fill the gaps
heatmap_data[is.na(heatmap_data)] <- 0

# Max and min energy 
min_energy <- min(all_data$Energy, na.rm = TRUE)
max_energy <- max(all_data$Energy, na.rm = TRUE)

# Put the mutation name in Y axis 
rownames(heatmap_data) <- heatmap_data$harf_mut
heatmap_data <- heatmap_data[,-1]  

heatmap_data <- as.matrix(heatmap_data)


labs <- matrix("", nrow = 69, ncol = 69)

#  1. Y (mutation aa) - 2. X (orginal aa) 


labs[19, 33] <- "*"
labs[5, 63] <- "*"
labs[8, 63] <- "*"
labs[9, 63] <- "*"
labs[10, 63] <- "*"
labs[18, 63] <- "*"
labs[19, 63] <- "*"
labs[5, 16] <- "*"
labs[10, 8] <- "*"
labs[12, 8] <- "*"
labs[18, 8] <- "*"
labs[12, 66] <- "*"
labs[8, 69] <- "*"
labs[8, 67] <- "*"
labs[11, 67] <- "*"
labs[12, 67] <- "*"
labs[12, 21] <- "*"
labs[9, 57] <- "*"
labs[14, 57] <- "*"
labs[8, 2] <- "*"
labs[10, 2] <- "*"
labs[14, 2] <- "*"
labs[5, 56] <- "*"
labs[9, 56] <- "*"
labs[9, 61] <- "*"
labs[10, 61] <- "*"


# Define colors for specific x-axis labels
highlight_posit <- c("G446", "Y449", "Y453", "L455", "F456", "Y473", "A475","S477", 
                     "G476", "F486", "N487", "Y489", "Q493", "G496", "Q498", 
                     "T500", "N501", "G502", "Y505")

# Create a vector for x-axis label colors, with red for highlighted positions and black elsewhere
cols <- rep("black", length(sequence_positions))
cols[sequence_positions %in% highlight_posit] <- "purple"  # Set red for highlighted positions

# Draw the heatmap with the colored x-axis labels
heatmap.2(as.matrix(heatmap_data), 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none",  
          col = bluered(25),  
          trace = "none",  
          main = "Energy Differences of RBM from only RBD", 
          xlab = "Residues within 5 Angstroms",  
          ylab = "Mutated Amino Acid-interesting mutation as stars (0/-1 range of corplot)",  
          key = TRUE,  
          cexCol = 0.7,  
          srtCol = 45,
          keysize = 1.5,  
          margins = c(5, 8),
          cellnote = labs, 
          notecol = "black",
          notecex = 2,
          labCol = sequence_positions,  # Custom x-axis labels
          colCol = cols)               # Color the x-axis labels
