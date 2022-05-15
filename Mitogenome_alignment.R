dir.create("data_output")
dir.create("data_output_nt")
dir.create("data_output_aa")
library(seqinr)
library(dplyr)
library(spgs) # use for getting the reverse-complementary sequences of ND6

# get access to the original mitogenome alignment that has all the elements
Mitogenome_aligned <- read.fasta(file = "mitogenome_alignment.txt",
                           seqtype = "DNA", as.string = TRUE,forceDNAtolower = FALSE,set.attributes = FALSE)
seq_align <- seq(1:length(Mitogenome_aligned))
seq_names <- names(Mitogenome_aligned)

# Read in position range of each mitogenome element (compiled by hand using one annotated mitogenome as the reference)
Mito_elem_range <- read.csv("Mito_elem_range.csv")
Elem_num<- nrow(Mito_elem_range)
Elem_len <- Mito_elem_range$Elem_end-Mito_elem_range$Elem_start+1 # to get the length of each element

#retrieve and save the individual alignment of all genes except ND6
Final_seq <-vector(mode = "list", length = length(Mitogenome_aligned))
for (i  in 1:(Elem_num-1)) {
  Element_seq <- lapply(seq_align, function(x) {substr(Mitogenome_aligned[x],Mito_elem_range[i,]$Elem_start,Mito_elem_range[i,]$Elem_end)})
  write.fasta(sequences = Element_seq, names = seq_names, file.out = paste0("data_output_nt/",Mito_elem_range$Elem_name[i],"_alignment",".fasta"))
  Element_aa <- lapply(seq_align, function(x) {translate(s2c(as.character(Element_seq[x])), frame = 0, sens = "F", numcode = 2, NAstring = "X", ambiguous = FALSE)})
  write.fasta(sequences = Element_aa, names = seq_names, file.out = paste0("data_output_aa/",Mito_elem_range$Elem_name[i],"_aa_alignment",".fasta"))
  Final_seq <- mapply(c, Final_seq, Element_seq,SIMPLIFY=F)
}

Concat_seq <- lapply(seq_align, function(x) paste(Final_seq[[x]][1:(Elem_num-1)],collapse='')) # concatenate all genes but ND6, which is the last one

ND6_seq <- lapply(seq_align, function(x) {substr(Mitogenome_aligned[x],Mito_elem_range[Elem_num,]$Elem_start,Mito_elem_range[Elem_num,]$Elem_end)})
ND6_seq_RC <- sapply(seq_align, function(x) reverseComplement(ND6_seq[x], case = "as is")) # to get the reverse-complementary seq of ND6
write.fasta(sequences = ND6_seq_RC, names = seq_names, file.out = "data_output_nt/ND6r_alignment.fasta")
ND6_aa <- lapply(seq_align, function(x) {translate(s2c(as.character(ND6_seq_RC[x])), frame = 0, sens = "F", numcode = 2, NAstring = "X", ambiguous = FALSE)})
write.fasta(sequences = ND6_aa, names = seq_names, file.out = "data_output_aa/ND6r_aa_alignment.fasta")

All_seq <- mapply(c, Concat_seq, ND6_seq_RC,SIMPLIFY=F)
All_concat_seq <- lapply(seq_align, function(x) paste(All_seq[[x]][1:2],collapse='')) # all genes with ND6 at the end
write.fasta(sequences = All_concat_seq, names = seq_names, file.out = "data_output/Selected_gene_alignment.fasta")


##### codes below are for generating partition files
Par_len <- 0
Elem_tab <-vector(mode = "list", length = 0)
for (i  in 1:(Elem_num)) {
  Par_len <- Par_len+Elem_len[i]
  Elem_tab <- append(Elem_tab, Par_len)
  Elem_tab_plus1 <- as.numeric(Elem_tab)+1 # to calcuate the begin position of the next gene (all but the first gene, which begins at 1)
}

begin_pos_all_except_first <- Elem_tab_plus1[1:Elem_num-1]
begin_pos1 <- c("1") # the first gene begins at 1
begin_pos <- append(begin_pos1,begin_pos_all_except_first)
Elem_tab_pos <- bind_cols(Mito_elem_range$Elem_name,data.frame(as.numeric(begin_pos)), data.frame(as.numeric(Elem_tab)))
colnames(Elem_tab_pos) <- c("Element","Begin_pos","End_pos")

Element_begin_end <- paste0(Elem_tab_pos$Element," = ",Elem_tab_pos$Begin_pos, "-", Elem_tab_pos$End_pos,";")
write.table(Element_begin_end,file = "data_output/Element_begin_end.txt", append = FALSE, quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)

Elem_tab_pos2 <- Elem_tab_pos[rep(1:Elem_num, each = 3), ] # duplicate each line twice
Elem_tab_row <- nrow(Elem_tab_pos2) # now have three times of rows

#Begin position
Elem_tab_pos2[((1:Elem_tab_row)+1) %% 3 == 0, ]$Begin_pos <- as.numeric(Elem_tab_pos2[((1:Elem_tab_row)+1) %% 3 == 0, ]$Begin_pos)+1 # numbers on each of the 1st duplicated+1
Elem_tab_pos2[1:Elem_tab_row %% 3 == 0, ]$Begin_pos <- as.numeric(Elem_tab_pos2[1:Elem_tab_row %% 3 == 0, ]$Begin_pos)+2 # numbers on each of the 2nd duplicated+2
#End position
Elem_tab_pos2[((1:Elem_tab_row)+2) %% 3 == 0, ]$End_pos <- as.numeric(Elem_tab_pos2[((1:Elem_tab_row)+2) %% 3 == 0, ]$End_pos)-2 # each original end -2
Elem_tab_pos2[((1:Elem_tab_row)+1) %% 3 == 0, ]$End_pos <- as.numeric(Elem_tab_pos2[((1:Elem_tab_row)+1) %% 3 == 0, ]$End_pos)-1 # numbers on each of the 1st duplicated+1-2
Elem_tab_pos2[1:Elem_tab_row %% 3 == 0, ]$End_pos <- as.numeric(Elem_tab_pos2[1:Elem_tab_row %% 3 == 0, ]$End_pos) # numbers on each of the 2nd duplicated+2-2
#Element names
Elem_tab_pos2$Element <- as.character(Elem_tab_pos2$Element)
Elem_tab_pos2[((1:Elem_tab_row)+2) %% 3 == 0, ]$Element <- paste0(Elem_tab_pos2[((1:Elem_tab_row)+2) %% 3 == 0, ]$Element,"_codon1")
Elem_tab_pos2[((1:Elem_tab_row)+1) %% 3 == 0, ]$Element <- paste0(Elem_tab_pos2[((1:Elem_tab_row)+1) %% 3 == 0, ]$Element,"_codon2")
Elem_tab_pos2[1:Elem_tab_row %% 3 == 0, ]$Element <- paste0(Elem_tab_pos2[1:Elem_tab_row %% 3 == 0, ]$Element,"_codon3")

Final_partition_file <- paste0(Elem_tab_pos2$Element," = ",Elem_tab_pos2$Begin_pos, "-", Elem_tab_pos2$End_pos, "\\3;")
write.table(Final_partition_file,file = "data_output/Partition_file.txt", append = FALSE, quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)


##### codes below are for checking unexpected stop codons in protein-coding genes
Codon_check <- Elem_len%%3==0 # check if the length of protein genes can be divided by 3
Codon_check_results <- bind_cols(Mito_elem_range$Elem_name, data.frame(Codon_check))
colnames(Codon_check_results) <- c("Element","Divisible by 3")
write.table(Codon_check_results,file = "data_output/Codon_check_results.csv", append = FALSE, quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)

All_concat_aa <- lapply(seq_align, function(x) {translate(s2c(as.character(All_concat_seq[x])), frame = 0, sens = "F", numcode = 2, NAstring = "X", ambiguous = FALSE)})

unexpect_stop_loc <- character(0)
for (i  in seq_align) {
  Asterisk_list <- grepl("\\*", All_concat_aa[i][[1]][1:(length(All_concat_aa[i][[1]]))]) # find all * in aa sequences
  Asterisk_loc <- as.integer(which(Asterisk_list==TRUE)*3) %in% as.integer(Elem_tab_pos$End_pos)
  unexpect_stop <- which(Asterisk_loc==FALSE) # unexpected stop codons are not found at the end of a gene
  unexpect_stop_loc<- append(unexpect_stop_loc, paste0("#", i, "_", seq_names[i]," site ", as.character(which(Asterisk_list==TRUE)*3)[unexpect_stop]))
}

seq_numb_only <- gsub("site.*","", unexpect_stop_loc)
Unexpec_stop_only <- as.numeric(gsub(".*site","", unexpect_stop_loc))
Start_row_results <- lapply(Unexpec_stop_only, function(x) {Elem_tab_pos$Begin_pos<=x & Elem_tab_pos$End_pos>=x})
Start_row_number <- as.numeric(lapply(Start_row_results,function(x){which(x==T)}))
Element_unexpect_stop <- sapply(Start_row_number,function(x) Elem_tab_pos$Element[x])
Stop_elem_only <- as.character(Element_unexpect_stop)

unexpect_stop_compile <- mapply(c, seq_numb_only, Unexpec_stop_only,Stop_elem_only, SIMPLIFY=F)
dat <- as.data.frame(unexpect_stop_compile)
write.table(t(dat),file = "data_output/unexpected_stop_codon_loc.csv", append = FALSE, na = "", quote = FALSE, sep = ",",
            col.names = c("Sequence#", "Location","Gene"), row.names = FALSE)



