# INS CDS 进化树（仅输出分支长度图）

suppressPackageStartupMessages({
  library(Biostrings)
  library(muscle)
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(grid)
})

options(ignore.negative.edge = TRUE)

input_fasta <- "INS_cds.fasta"
output_png <- "INS_tree.png"

dna_seqs <- Biostrings::readDNAStringSet(input_fasta)
alignment <- muscle::muscle(dna_seqs, verbose = FALSE)
dna_bin <- as.DNAbin(alignment)

dist_matrix <- ape::dist.dna(dna_bin, model = "TN93", pairwise.deletion = TRUE)
tree <- ape::fastme.bal(dist_matrix)

# 简单定根：取树上最远的 tip 作为 outgroup
D <- stats::cophenetic(tree)
root_tip <- rownames(D)[which(D == max(D), arr.ind = TRUE)[1, 1]]
tree <- ape::root(tree, outgroup = root_tip, resolve.root = TRUE)

species_names <- c(
  "human_tx1" = "Homo sapiens (人-1)",
  "human_tx2" = "Homo sapiens (人-2)",
  "human_tx3" = "Homo sapiens (人-3)",
  "human_tx4" = "Homo sapiens (人-4)",
  "mouse_tx1" = "Mus musculus (小鼠-1)",
  "mouse_tx2" = "Mus musculus (小鼠-2)",
  "mouse_tx3" = "Mus musculus (小鼠-3)",
  "rat_tx1" = "Rattus norvegicus (大鼠-1)",
  "rat_tx2" = "Rattus norvegicus (大鼠-2)",
  "cow_tx1" = "Bos taurus (牛-1)",
  "cow_tx2" = "Bos taurus (牛-2)",
  "pig_tx1" = "Sus scrofa (猪-1)",
  "pig_tx2" = "Sus scrofa (猪-2)",
  "dog_tx1" = "Canis lupus familiaris (狗)",
  "rabbit_tx1" = "Oryctolagus cuniculus (兔)",
  "chimpanzee_tx1" = "Pan troglodytes (黑猩猩-1)",
  "chimpanzee_tx2" = "Pan troglodytes (黑猩猩-2)",
  "cat_tx1" = "Felis catus (猫-1)",
  "cat_tx2" = "Felis catus (猫-2)",
  "macaque_tx1" = "Macaca mulatta (猕猴)"
)

groups <- c(
  "human_tx1" = "Primates (灵长类)", "human_tx2" = "Primates (灵长类)", "human_tx3" = "Primates (灵长类)", "human_tx4" = "Primates (灵长类)",
  "chimpanzee_tx1" = "Primates (灵长类)", "chimpanzee_tx2" = "Primates (灵长类)",
  "macaque_tx1" = "Primates (灵长类)",
  "mouse_tx1" = "Rodentia (啮齿类)", "mouse_tx2" = "Rodentia (啮齿类)", "mouse_tx3" = "Rodentia (啮齿类)",
  "rat_tx1" = "Rodentia (啮齿类)", "rat_tx2" = "Rodentia (啮齿类)",
  "cow_tx1" = "Artiodactyla (偶蹄类)", "cow_tx2" = "Artiodactyla (偶蹄类)",
  "pig_tx1" = "Artiodactyla (偶蹄类)", "pig_tx2" = "Artiodactyla (偶蹄类)",
  "dog_tx1" = "Carnivora (食肉类)", "cat_tx1" = "Carnivora (食肉类)", "cat_tx2" = "Carnivora (食肉类)",
  "rabbit_tx1" = "Lagomorpha (兔类)"
)

original_labels <- tree$tip.label
full_names <- unname(species_names[original_labels])
label_map <- data.frame(
  label = full_names,
  group = unname(groups[original_labels]),
  stringsAsFactors = FALSE
)

tree$tip.label <- full_names

group_colors <- c(
  "Primates (灵长类)" = "#E41A1C",
  "Rodentia (啮齿类)" = "#377EB8",
  "Artiodactyla (偶蹄类)" = "#4DAF4A",
  "Carnivora (食肉类)" = "#FF7F00",
  "Lagomorpha (兔类)" = "#A65628"
)

# 压缩横向比例尺：越大越“短”
max_x <- max(ape::node.depth.edgelength(tree), na.rm = TRUE)
compress_factor <- 1.35
x_limits <- c(0, max_x * compress_factor)
x_breaks <- seq(0, ceiling(max_x / 0.05) * 0.05, by = 0.05)

p <- ggtree(tree, branch.length = "branch.length") %<+% label_map +
  geom_tiplab(
    aes(color = group),
    size = 7,
    offset = max(tree$edge.length, na.rm = TRUE) * 0.02,
    fontface = "italic"
  ) +
  scale_color_manual(values = group_colors, name = "Taxonomic Group") +
  theme_tree2() +
  theme(text = element_text(size = 14)) +
  xlab("Genetic Distance") +
  scale_x_continuous(limits = x_limits, breaks = x_breaks) +
  coord_cartesian(clip = "off") +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c("left", "top"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20),
    legend.key.size = unit(1.3, "lines"),
    legend.spacing.y = unit(1.6, "lines"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.subtitle = element_blank(),
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    plot.margin = margin(2, 2, 2, 2, unit = "cm")
  ) +
  labs(title = "Phylogenetic Tree of INS Gene (CDS) with Branch Lengths") +
  guides(color = guide_legend(override.aes = list(size = 5)))

ggsave(output_png, plot = p, width = 18, height = 12, dpi = 300, bg = "white")
message("Saved: ", output_png)
