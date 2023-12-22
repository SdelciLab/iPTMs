# HUPO missing proteins
library(stringr)
library(data.table)
library(ggplot2)
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10081950/
Missing_prots = openxlsx::read.xlsx('/home/v1skourt/Documents/NIHMS1867254-supplement-Suppl__Table_1.xlsx', startRow = 2)
Missing_uniprots = Missing_prots$acc..code |> str_remove_all('NX_')
isoform_fasta <- seqinr::read.fasta('/home/v1skourt/proteomehd2/fasta/micro_prot_decoy.fasta',
                                    seqtype = "AA",as.string = TRUE)
isoform_fasta = data.table(Annotation  = unlist(seqinr::getAnnot(isoform_fasta)) |>  str_remove(' [:print:]*$') |> str_remove('>'),
                           Seq = seqinr::getSequence(isoform_fasta, as.string = T) |>  unlist())
proteomehd2 =fread(here::here('in','datasets',
                              'picked_protein_group_no_remap',
                              'proteinGroups_picked_protein_group_no_remap.txt'))
proteomehd2_proteins = proteomehd2[Reverse != '+'][`Q-value`<0.01,`Protein IDs`] |> str_subset(';',negate = T)
proteomehd2_proteins = proteomehd2_proteins |> str_match('\\|([:print:]*)\\|')
proteomehd2_proteins = proteomehd2_proteins[,2]
proteomehd2_proteins=  proteomehd2_proteins[!is.na(proteomehd2_proteins)]
found_missing = Missing_uniprots[Missing_uniprots%in% proteomehd2_proteins ]

proteomehd2_pep_mapping =fread(here::here('in','datasets','pep_to_prot_mapping.txt'), header = F)
proteomehd2_pep_mapping = proteomehd2_pep_mapping[str_detect(V2,';', negate = T )][str_detect(V2,'REV', negate = T )]
proteomehd2_pep_mapping_missing = proteomehd2_pep_mapping[V2 %in% found_missing][,stripped_pept:= str_remove(V1,'^n')]

detect_overlapping = function(list_of_peptides,peptide_of_interest){
  list_of_peptides = list_of_peptides[list_of_peptides != peptide_of_interest]
  any(str_detect(list_of_peptides,peptide_of_interest))
}

proteomehd2_pep_mapping_missing[,overlapping:=detect_overlapping(proteomehd2_pep_mapping_missing$stripped_pept,
                                                                 stripped_pept) , by = stripped_pept]
proteomehd2_pep_mapping_missing[overlapping==F][,number:=1] |> 
  ggplot(aes(x = V2,y= number))+
  geom_col()+theme_bw()+
  geom_hline(yintercept = 2, linetype = 'dotted' ,colour = 'red')+
  coord_flip()+
  ggtitle('HUPO requires two non-overlapping unique peptides for these missing uniprots')

plot_coverage_prot <- function(Uniprot_tmp,mapping, isoform_fasta){
  # Uniprot_tmp = 'Q8NET4'
  # mapping = proteomehd2_pep_mapping
  mapping_tmp = mapping[V2 == Uniprot_tmp] |> unique()
  fasta = isoform_fasta[str_detect(Annotation,'rev',negate = T)]
  fasta[,Uniprot:=str_remove(Annotation,'^[:print:]*?\\|') |> str_remove('\\|[:print:]*$')]
  fasta_tmp = fasta[Uniprot==Uniprot_tmp,Seq]
  coverage_plot = data.table(peptide = 'Whole_prot',
                             position = 1:(stringr::str_length(fasta_tmp)))
  for(i in 1:nrow(mapping_tmp)){
    tmp_peptide =mapping_tmp[i,V1] |> str_remove('^n')
    positions = str_locate(fasta_tmp,tmp_peptide)
    coverage_plot = rbind(coverage_plot,
                          data.table(
                            peptide = tmp_peptide,
                            position = positions[1,1]:positions[1,2]
                          ))
    
  }
  levels = coverage_plot[,.(min_position = min(position)),by = peptide][order(min_position),peptide]
  coverage_plot[,`:=`( peptide =factor(peptide,levels = levels),
                       position = as.factor(position))]
  coverage_plot
  
}

dt = plot_coverage_prot('single_ribo1318_HUMAN',proteomehd2_pep_mapping,isoform_fasta)


ggplot(dt,aes(y =peptide ,x= position))+
  geom_point()+theme(axis.text.x = element_blank(),
                     axis.ticks = element_blank())+
  theme_bw()+  ggtitle('single_ribo1318_HUMAN')
