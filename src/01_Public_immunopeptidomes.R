### ---------------------------------------------- Public immunopeptidomes annotation  ----------------------------------------------
# description:  Curate public MHC-I immunopeptidome data into a small dataset with user criteria
# 
# input:        - IEDB, MHCIatlas, and SysMHC datasets of T-cell epitopes
# 
# output:       - Filtered human self-peptides detected by mass-spectrometry, affinity predicted with NetMHCPan
#               - Affinity distributions for length-haplotype combinations
# author:       Yehor Horokhovskyi

library(Biostrings)
library(data.table)
library(dplyr)
library(dtplyr)
library(eulerr)
library(ggplot2)
library(ggthemes)
library(grid)
library(purrr)
library(RColorBrewer)
library(stringr)
library(tidyr)
library(readr)

# Multi-threading
library(doParallel)
n_cores = 8
cl <- registerDoParallel(makeCluster(n_cores))

wd <- "/home/yhorokh/data/Results_reports/wd/AffinitySPI_cutoffs/"
dir.create(paste0(wd,"results"))
dir.create(paste0(wd,"results/Public_IP_MHC_affinity"))
dir.create(paste0(wd,"results/Public_IP_MHC_affinity/kmers"))
dir.create(paste0(wd,"results/Public_IP_MHC_affinity/kmers_predicted"))

# MHC-I affinity prediction
netMHCpan <- "/home/yhorokh/data/Results_reports/Mouse_lymphoma_target_11.2020/bin/netMHCpan-4.1/netMHCpan"
chunk_length <- 10^5
netMHCpan_alleles <- fread("data/alleles.txt", col.names = c("V1", "V2"), nThread = n_cores) %>%
  .[, V1 := str_squish(V1)] %>%
  .[, V2 := str_squish(V2)]

### Alleles to predict. If alleles = T, use all available
alleles <- T
# alleles <- c("HLA-A0101", "HLA-A0201", "HLA-A0301", "HLA-A1101", "HLA-A2301", "HLA-A2402", "HLA-A3101","HLA-A6801", 
#              "HLA-B0801", "HLA-B1401", "HLA-B1501", "HLA-B2705", "HLA-B3501", "HLA-B4001", "HLA-B4402", "HLA-B4403",
#              "HLA-B5101", "HLA-C0802", "HLA-C1506", "HLA-A0319", "HLA-A6802", "HLA-B1503", "HLA-C1203")
### Which data to use
Use_strictly_monoallelic <- TRUE
min_peptides_per_allele.length <- 100
min_length = 8
max_length = 15

### IP databases
IP_DB <- list()

# IEDB
IP_DB$IEDB <- fread(paste0(wd, "/data/Public_Immunopeptidomes/mhc_ligand_full_03.02.2023.csv.gz"), sep = ",", nThread = n_cores, check.names = T)

if (Use_strictly_monoallelic == FALSE) {
  # MHCIatlas
  IP_DB$MHCIatlas <- list.files(paste0(wd, "data/Public_Immunopeptidomes/MHCIatlas/"), full.names = T) %>% 
    as.list() %>%
    lapply(fread, sep = ",", nThread = n_cores) %>%
    rbindlist()
  
  # SysMHC dataset
  IP_DB$SysMHC <- list.files(paste0(wd, "data/Public_Immunopeptidomes/SYSTEMHCATLAS/"), full.names = T) %>%
    as.list() %>%
    lapply(fread, colClasses=cols(
      dir0 = col_character(),
      SysteMHC_ID = col_character(),
      SampleID = col_character(),
      protein_id = col_character(),
      search_hit = col_character(),
      prob = col_double(),
      annotation_score = col_double(),
      numNeighbors = col_double(),
      spectral_counts = col_double(),
      top_allele = col_character(),
      MHCClass = col_character(),
      Organism = col_character()
    ), sep = ",", nThread = n_cores) %>%
    rbindlist()
}

# For plots
gg <- list()

# -------------------------------------------- Pre-processing --------------------------------------------
### Create a common immune peptide sequence database. Keep human MHC-I self-peptides
### IEDB
IP_DB$IEDB <- IP_DB$IEDB %>%
  lazy_dt() %>%
  dplyr::filter(str_detect(Host, "Homo sapiens")) %>%
  dplyr::filter(MHC.3 == "I") %>%
  dplyr::filter(grepl("Homo sapiens", Epitope.10)) %>%
  dplyr::filter(!is.na(Epitope.2)) %>%
  dplyr::filter(Assay.1 %in% c("mass spectrometry", 
                               "secreted MHC/mass spectrometry", 
                               "cellular MHC/mass spectrometry")) %>%
  dplyr::select(Reference.7, `Antigen.Processing.Cells`, MHC, Epitope.2) %>%
  dplyr::rename(Experiment = Reference.7,
                Material = `Antigen.Processing.Cells`,
                MHC = MHC,
                Peptide = Epitope.2) %>%
  as.data.table()

if (!Use_strictly_monoallelic) {
  ### SysMHC
  IP_DB$SysMHC <- IP_DB$SysMHC %>%
    .[Organism == "Human",] %>%
    .[MHCClass == "I",] %>%
    .[,Material := str_split_fixed(dir0, "_", Inf)[,6]] %>%
    lazy_dt() %>%
    dplyr::filter(!grepl("CONT_", protein_id, fixed = T))  %>%
    dplyr::select(dir0, Material, top_allele, search_hit) %>%
    dplyr::rename(Experiment = dir0,
                  Material = Material,
                  MHC = top_allele,
                  Peptide = search_hit) %>%
    as.data.table()
  
  ### IP mouse
  IP_DB$MHCIatlas <- IP_DB$MHCIatlas %>%
    lazy_dt() %>%
    rename(Accession = accession) %>%
    dplyr::filter(grepl("HUMAN", Accession, fixed = T)) %>%
    dplyr::filter(!grepl("#CONTAM#", Accession, fixed = T)) %>%
    dplyr::filter(!grepl("CONT_", Accession, fixed = T)) %>%
    dplyr::filter(!grepl("iRT", Accession, fixed = T)) %>%
    dplyr::filter(!is.na(Accession)) %>%
    dplyr::select(File_name, Tissue, Best_Patient.Allele, sequence) %>%
    dplyr::rename(Experiment = File_name,
                  Material = Tissue,
                  MHC = Best_Patient.Allele,
                  Peptide = sequence) %>%
    as.data.table()
}

IP_DB <- rbindlist(IP_DB, idcol = "Database") %>%
  .[,Peptide := str_split_fixed(Peptide, pattern = fixed(" + "), Inf)[,1]] %>%
  .[,Peptide := str_replace_all(Peptide, pattern = "I", replacement = "L")] %>%
  .[,Peptide := str_squish(Peptide)] %>%
  .[,length := nchar(Peptide)] %>%
  .[length %in% min_length:max_length,]  %>%
  .[!str_detect(Peptide, "X|U|\\*"),]  %>%
  unique() 

### Before filtering
IP_DB$MHC %>% table() %>% sort() %>% print()
length_IP_DB_unfiltered <- n_distinct(IP_DB$Peptide)

# Keep only alleles of interest
IP_DB <- IP_DB %>%
  .[,MHC := str_remove_all(MHC, "\\*|\\:|\\)|\\(")] %>%
  .[,MHC := fifelse(!str_starts(MHC, "HLA-"), paste0("HLA-", MHC), MHC)] %>%
  .[,MHC := str_replace_all(MHC, "HLA-HLA-", "HLA-")] %>%
  .[,MHC := str_squish(MHC)] %>%
  .[MHC %chin% netMHCpan_alleles[["V1"]] | MHC %chin% netMHCpan_alleles[["V2"]]]

### after filtering by NetMHCPan alleles
allele_freq <- IP_DB %>%
  lazy_dt() %>%
  group_by(MHC, length) %>%
  summarise(unique_peptides = n_distinct(Peptide)) %>%
  as.data.table()
allele_freq %>% 
  print()

keep_alleles <- allele_freq[unique_peptides >= min_peptides_per_allele.length, .(MHC, length)]

### Keep alleles with minimal number of peptides
IP_DB <- keep_alleles[IP_DB, on =  .(MHC, length), nomatch = NULL]
IP_DB

if (length(alleles) == 1) {
  if (alleles == T) {
    print("All available alleles will be used")
  }
} else {
  IP_DB <- IP_DB %>%
    .[MHC %chin% alleles,]
}

# EDA
nchar(IP_DB$Peptide) %>% table() %>% plot(main = "Peptide length distribution")
IP_DB$Material %>% table()
IP_DB$MHC %>% table()
IP_DB$Database %>% table()

# -------------------------------------------- Predict MHC-I affinity --------------------------------------------
### Select relevant alleles
kmers <- IP_DB %>%
  .[, .(Peptide, length, MHC)] %>%
  unique() %>%
  split(by = c("length", "MHC"))

# Export kmers for search
for (k in seq_along(kmers)) {
  seq_in <- kmers[[k]]$Peptide
  
  n_chunks = ceiling(ceiling(length(seq_in)) / chunk_length)
  for (ch in 1:n_chunks) {
    first = 1 + ((-1 + ch) * chunk_length)
    last = chunk_length + ((-1 + ch) * chunk_length)
    
    if (last > length(seq_in)) {
      last <- length(seq_in)
    }
    fwrite(as.data.table(seq_in[first:last]), col.names = F, append = F, 
           file = paste0(wd, "results/Public_IP_MHC_affinity/kmers/",names(kmers)[k],"_ch_",ch,".tsv"))
  }
}

# Generate commands
files <- list.files(paste0(wd, "/results/Public_IP_MHC_affinity/kmers/"), ".tsv", 
                    full.names = T, recursive = F)

files_names <- list.files(paste0(wd, "/results/Public_IP_MHC_affinity/kmers/"), ".tsv", 
                          full.names = F, recursive = F) %>%
  str_remove_all(pattern = ".tsv") %>%
  str_remove_all(pattern = "_mers") %>%
  str_remove_all(pattern = "_ch") %>%
  str_split_fixed(pattern = "_", n = Inf) %>%
  as.data.frame() %>%
  rename(pep_length=V1,
         chunk=V2) %>%
  separate(pep_length, into = c("pep_length", "allele"), sep = "\\.") %>%
  cbind(file=files) %>%
  arrange(-file.size(file)) %>%
  as_tibble()

files_names$cmds <- paste(netMHCpan,
                          "-BA", "-inptype 1",
                          "-a", files_names$allele,
                          "-l", files_names$pep_length, 
                          "-p -f",
                          files_names$file,
                          ">", paste0(wd, "/results/Public_IP_MHC_affinity/kmers_predicted/", 
                                      files_names$pep_length, "_mers_", files_names$allele, "_ch_", files_names$chunk, ".txt"),
                          "-v")
files_names$output <- str_split_fixed(files_names$cmds, "/kmers_predicted/", Inf)[,2] %>%
  str_remove_all(pattern = " -v")
files_names
files_names$cmds[1]

# Remove processed files
files_ready <- list.files(paste0(wd, "/results/Public_IP_MHC_affinity/kmers_predicted/"), ".txt", 
                          full.names = F, recursive = F)
files_names <- filter(files_names, !output %in% files_ready)
files_names

# Execute the commands in parallel
foreach(i = 1:length(files_names$cmds)) %dopar% system(paste0(files_names$cmds[i]))

# Read in the files 
affinity_pred <- list.files(paste0(wd, "/results/Public_IP_MHC_affinity/kmers_predicted/"), ".txt", 
           full.names = T, recursive = F) %>%
  lapply(FUN = function(x){
    print(x)
  
  if (file.exists(x)) {
    # Read and filter netMHCpan output
    binders_df <- try(fread(x, nThread = n_cores, sep = "^",
                            blank.lines.skip = TRUE, col.names = "V1",
                            skip = " Pos"))
    if (!"try-error" %in% class(binders_df)) {
      binders_df <- binders_df$V1 %>%
        str_replace_all(pattern = "[[:space:]]+", " ") %>%
        str_replace_all(pattern = "\\*", "") %>%
        str_split_fixed(pattern = " ", n = Inf) %>%
        as.data.table() %>%
        lazy_dt() %>%
        rename(MHC = V2, peptide = V3, `Aff(nM)` = V16) %>%
        select(MHC, peptide, `Aff(nM)`) %>%
        mutate(`Aff(nM)` = as.numeric(`Aff(nM)`)) %>%
        filter((!is.na(`Aff(nM)`))) %>%
        mutate(peptide = str_remove_all(peptide, "X")) %>%
        as.data.table()
      # print(head(binders_df))
    } else {
      binders_df = data.table(MHC = character(),
                              peptide = character(),
                              `Aff(nM)` = numeric())
      print("Warning: empty input!")
      print("Make sure that allele is present in netMHCpan-4.1/data/allelenames")
    }
    return(binders_df)
  } 
}) %>%
  rbindlist() %>%
  .[, MHC := str_remove_all(MHC, "\\*")] %>%
  .[!MHC %in% c("PEPLIST","PEPLIST."),] %>%
  .[!is.na(MHC),] %>%
  .[!is.na(`Aff(nM)`),] %>%
  setnames(old = "peptide", new = "Peptide") %>%
  setkey(Peptide)
affinity_pred

### ------------------------------------------- Plots -------------------------------------------
# Plot affinity distribution for the corresponding alleles
IP_DB_affinity <- IP_DB %>%
  setkey(Peptide) %>%
  merge(y = affinity_pred, by = c("Peptide", "MHC"), all.x = T, all.y = T, allow.cartesian = T) %>%
  na.omit() %>%
  unique() %>%
  as_tibble() %>%
  relocate(Database, Material, Peptide, length, MHC, `Aff(nM)`, Experiment)

IP_DB_affinity$length <- factor(IP_DB_affinity$length, levels = min_length:max_length)
IP_DB_affinity$MHC <- factor(IP_DB_affinity$MHC, levels = sort(unique(IP_DB_affinity$MHC)))

mycolors = c(rev(colorRampPalette(brewer.pal(9, "Blues"))(length(which(str_starts(levels(IP_DB_affinity$MHC), "HLA-A"))) + 3)), 
             rev(colorRampPalette(brewer.pal(9, "Oranges"))(length(which(str_starts(levels(IP_DB_affinity$MHC), "HLA-B")))  + 3)), 
             rev(colorRampPalette(brewer.pal(9, "Purples"))(length(which(str_starts(levels(IP_DB_affinity$MHC), "HLA-C")))  + 3)))

gg$IP_DB_peptide_affinity <- IP_DB_affinity %>%
  mutate(MHC_type = str_extract(MHC, "HLA-A|HLA-B|HLA-C")) %>%
  ggplot(aes(x=`Aff(nM)`)) + 
  geom_density(aes(y=..density.., color=MHC, linetype=Database), stat="density", position="identity", alpha=0.5) + 
  facet_grid(length ~ MHC_type) + 
  theme_bw() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5), 
        legend.position = "none") + 
  xlab("Aff(nM)") + 
  ylab("Density") +
  # scale_x_log10(limits = c(NA, 500)) +
  ggtitle(label = "IC50 predicted with NetMHCpan 4.1")  +
  scale_color_manual(values = mycolors)
gg$IP_DB_peptide_affinity

gg$IP_DB_peptide_length <-  IP_DB_affinity %>%
  select(MHC, Peptide, length) %>%
  unique() %>%
  mutate(length = as.numeric(as.character(length))) %>%
  ggplot(aes(x=length, fill=MHC)) + 
  geom_histogram(stat="bin", position="dodge", alpha=1, bins=1) + 
  geom_hline(yintercept = 1000, color = "grey50", show.legend = T) +
  theme_bw() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5), 
        legend.position = "none") + 
  xlab("length") + 
  ylab("count") +
  ggtitle(label = "# Unique peptides per HLA-I haplotype") +
  scale_fill_manual(values = mycolors) 
gg$IP_DB_peptide_length

affinity_thresholds <- IP_DB_affinity %>%
  select(MHC, Peptide, `Aff(nM)`, length) %>%
  unique() %>%
  filter(!is.na(`Aff(nM)`))%>%
  mutate(length = as.numeric(as.character(length))) %>%
  group_by(MHC, length) %>%
  mutate(`Unique peptides` = n_distinct(Peptide),
         percentile_50 = quantile(`Aff(nM)`, probs = c(0.5)) %>% round(),
         percentile_60 = quantile(`Aff(nM)`, probs = c(0.6)) %>% round(),
         percentile_75 = quantile(`Aff(nM)`, probs = c(0.75)) %>% round(),
         percentile_95 = quantile(`Aff(nM)`, probs = c(0.95))  %>% round(),
         percentile_99 = quantile(`Aff(nM)`, probs = c(0.99))  %>% round()) %>%
  select(MHC, length, `Unique peptides`, contains("percentile")) %>%
  unique() %>%
  arrange(MHC, length) 
print.data.frame(affinity_thresholds)
print(paste0(round(100 * n_distinct(IP_DB_affinity$Peptide) / length_IP_DB_unfiltered), "% Epitopes annotated from initial data"))

# -------------------------------------------- Export tables --------------------------------------------
write_csv(IP_DB, paste0(wd, "results/Public_IP_databases.csv"))
write_csv(IP_DB_affinity, paste0(wd, "results/Public_IP_databases_affinity.csv"))
write_csv(affinity_thresholds, paste0(wd, "results/Public_IP_affinity_thresholds.csv"))

# -------------------------------------------- EDA --------------------------------------------
pdf(file = paste0(wd,"/results/IP_DB_overlap.pdf"),
    width = 12, height = 8)

nchar(IP_DB$Peptide) %>% table() %>% plot(main = "Peptide length distribution")
IP_DB$Material %>% table() %>% plot()
IP_DB$MHC %>% table() %>% t() %>% plot()
IP_DB$Database %>% table() %>% plot()

### Explore overlaps - By database
# peptide
set_list <- IP_DB %>%
  as_tibble() %>%
  select(Peptide, Database) %>%
  split(f = as.factor(.$Database)) %>%
  map(~discard(flatten_chr(.), . == "")[-1]) %>%
  lapply(unique)
plot(euler(set_list, shape = "ellipse"), quantities = TRUE, main = "Unique pepties by Database")

# ### Explore overlaps - By Material
# peptide
set_list <- IP_DB %>%
  as_tibble() %>%
  select(Peptide, Material) %>%
  split(f = as.factor(.$Material)) %>%
  map(~discard(flatten_chr(.), . == "")[-1]) %>%
  lapply(unique)
UpSet.set_list <- set_list %>%
  UpSetR::fromList() %>%
  UpSetR::upset(nsets = length(set_list),
                keep.order = T)
UpSet.set_list
grid.text("Unique Peptide sequence overlaps",x = 0.65, y=0.95, gp=gpar(fontsize=14))

for (i in seq_along(gg)) {
  plot(gg[[i]])
}

dev.off()

# -------------------------------------------- Export ggplots --------------------------------------------
for (i in seq_along(gg)) {
  g <- gg[[i]]
  ggsave(g, filename = paste0(wd, "/results/",names(gg)[[i]],".png"), width = 20, height = 10, pointsize = 15, dpi = "retina")
  ggsave(g, filename = paste0(wd, "/results/",names(gg)[[i]],".svg"), width = 20, height = 10)
  ggsave(g, filename = paste0(wd, "/results/",names(gg)[[i]],".pdf"), width = 20, height = 10)
}

