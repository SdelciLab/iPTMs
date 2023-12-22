# protein-coverage proteomehd2
library(data.table);library(stringr);library(ggplot2)

isoform_fasta <- seqinr::read.fasta('/home/v1skourt/proteomehd2/fasta/micro_prot_decoy.fasta',
                                    seqtype = "AA",as.string = TRUE)
isoform_fasta = data.table(Annotation  = unlist(seqinr::getAnnot(isoform_fasta)) |>  str_remove(' [:print:]*$') |> str_remove('>'),
                           Seq = seqinr::getSequence(isoform_fasta, as.string = T) |>  unlist())
pgFDR_pGroups =fread(here::here('in','datasets',
                              'picked_protein_group_no_remap',
                              'proteinGroups_picked_protein_group_no_remap.txt'))
PG_all_cont = pgFDR_pGroups$`Protein IDs` |> strsplit(';') |> unlist() |> stringr::str_subset('orf|ribo|ORF|cogn',negate = T) |> stringr::str_remove('[:print:]*_') |> 
  table() |> names()
PG_all_cont = paste0('_',PG_all_cont[stringr::str_detect(PG_all_cont,'HUMAN',negate = T)])
pgFDR_pGroups = pgFDR_pGroups[,`Potential contaminant`:= str_detect(`Protein IDs`,paste(PG_all_cont,collapse = '|'))]
pgFDR_pGroups = pgFDR_pGroups[`Q-value`<0.01]
pgFDR_pGroups = pgFDR_pGroups[Reverse != '+' ]
pgFDR_pGroups = pgFDR_pGroups[`Potential contaminant`==F]

pgFDR_pGroups = str_split(pgFDR_pGroups[,`Protein IDs`] ,';')

length_PG_per_group = purrr::map(.x = pgFDR_pGroups,~
                                   purrr::map_dbl(.x= .x,~isoform_fasta[Annotation==.x,Seq] |> str_length()))
longest_PG = purrr::map_dbl(.x = length_PG_per_group,~head(which(.x ==max(.x)),1))
longest_PGs = purrr::map2_chr(.x = pgFDR_pGroups,.y = longest_PG,~.x[.y])

Protein_coverage = isoform_fasta[Annotation %in% longest_PGs]

psm_percolator = fread(here::here('in','datasets','pep_to_prot_mapping.txt'),header = F)
names(psm_percolator) = c('sequence','Uniprot')
psm_percolator[,sequence:= str_remove(sequence,'^n')]
count_coverage <- function(Annotation,sequence,psms_percolator){
  
  # Annotation = Protein_coverage[2,1] |> unlist()
  Annotation = Annotation |> str_remove('^[:print:]*?\\|') |> str_remove('\\|[:print:]*') 
  # sequence  = Protein_coverage[2,2] |> unlist()
  psm_percolator_tmp = psm_percolator[str_detect(Uniprot,glue::glue('(^|;){Annotation}')),sequence]
  positions_tmp = c()
  # print(Annotation)
  # print(sequence)
  for(i in psm_percolator_tmp){
    # sometimes it's missing for leucine and isoleucine
    # i = psm_percolator_tmp[1]
    location = str_locate(sequence,i)
    if(!any(is.na(location))){
      positions_tmp = c(positions_tmp,location[1]:location[2])}
  }
  positions_tmp = unique(positions_tmp)
  perc_coverage = length(positions_tmp)/str_length(sequence)
  perc_coverage
}
Protein_coverage[,perc_coverage:= count_coverage(Annotation,Seq,psms_percolator),by = Annotation]
fwrite(Protein_coverage,here::here('out','datasets','pgFDR_prot_coverage.tsv' ))

ggplot(Protein_coverage,aes(x = perc_coverage))+
  geom_histogram(bins = 50)+
  ggtitle('protein coverage prohd2')
