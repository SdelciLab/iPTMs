# Analysis of HRAS Mut data
library(data.table);library(stringr);library(ggplot2)
list_figures <- list()
HRAS_psms = fread(here::here('Datasets','Raw','HRAS_all_psms.csv'))
HRAS_psms[,sample:= stringr::str_remove_all(sample,'^[:print:]*\\/\\/') |> stringr::str_remove_all('\\/psm.tsv')]
HRAS_psms[Peptide =='LVVVGAGGVGK'] |> View()
Samples =openxlsx::read.xlsx(here::here('Datasets','Raw','SampleAnnotations.xlsx')) |> as.data.table()
Samples[,`:=`(sample_name = paste(Biological.Replicate,Cell.Type, sep= '_') |> stringr::str_replace_all('-','_'),
              `Raw.File.Name` =`Raw.File.Name` |> stringr::str_replace_all('-','_'))]
setnames(Samples,'Raw.File.Name','sample')

HRAS_psms = Samples[,.(sample,sample_name)][HRAS_psms, on = 'sample']
Psm_per_prot = HRAS_psms[Gene != ''][,.(N_psms= .N), by = .(Gene,sample_name)]
Psm_per_prot[, bait := fifelse(Gene =='HRAS','bait','prey')]
Psm_per_prot[,cell_line := stringr::str_remove_all(sample_name,'[:print:]*_._')]
Psm_per_prot[,HRAS_state := fifelse(stringr::str_detect(sample_name,'wt'), 'WT','G12D')]
HRAS_enrichment <- Psm_per_prot |> ggplot(aes(x = log10(N_psms), y= sample_name, fill = HRAS_state ))+
    ggridges::geom_density_ridges(alpha = 0.4,  aes(colour = HRAS_state),rel_min_height=0.001)+theme_bw()+
    geom_point(data = Psm_per_prot[bait == 'bait'], aes(colour = HRAS_state),size = 2.5)+
    facet_wrap('cell_line',nrow = 3,scales = 'free_y')+
    scale_fill_manual(values = c('#EDA4B2','#5BCFFA'))+
    scale_colour_manual(values = c('#EDA4B2','#5BCFFA'))+
    ggtitle('Bait HRAS is enriched with numerous Psms in IP-MS samples')
ggsave(here::here('Output','HRAS_enrichment.pdf'),HRAS_enrichment)
list_figures$HRAS_enrichment = HRAS_enrichment
isoform_fasta <- seqinr::read.fasta(here::here('Datasets','Raw','micro_prot_decoy.fasta'),
                                    seqtype = "AA",as.string = TRUE)
isoform_fasta = data.table(Annotation  = unlist(seqinr::getAnnot(isoform_fasta)) |>  str_remove(' [:print:]*$') |> str_remove('>'),
                           Seq = seqinr::getSequence(isoform_fasta, as.string = T) |>  unlist())
isoform_fasta= isoform_fasta[str_detect(Annotation,'rev',negate = T)]
isoform_fasta[,Annotation:= stringr::str_remove_all(Annotation,'^[:print:]*?\\|') |> stringr::str_remove_all('\\|[:print:]*?$')]
Uniprots_FDR =HRAS_psms$`Protein ID` |> unique()

Protein_coverage = isoform_fasta[Annotation %in% Uniprots_FDR]

psm_percolator = HRAS_psms[,.(Peptide,`Protein ID`)] |> unique()
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
HRAS_coverage = Protein_coverage[Annotation == 'P01112',perc_coverage]*100
HRAS_coverage = Protein_coverage |> 
    ggplot(aes(x= perc_coverage*100))+
    theme_bw()+
    geom_histogram(binwidth =1, alpha  = 0.5)+
    geom_vline(xintercept = HRAS_coverage, colour = 'black', linetype="longdash" )+
    ggtitle('HRAS bait has very high protein coverage')
ggsave(here::here('Output','HRAS_relative_coverage.pdf'),HRAS_coverage)
list_figures$HRAS_coverage = HRAS_coverage
plot_coverage_prot <- function(Uniprot_tmp,mapping, isoform_fasta){
    # Uniprot_tmp = 'Q8NET4'
    # mapping = proteomehd2_pep_mapping
    mapping_tmp = mapping[Uniprot == Uniprot_tmp] |> unique()
    fasta = isoform_fasta[str_detect(Annotation,'rev',negate = T)]
    fasta[,Uniprot:=str_remove(Annotation,'^[:print:]*?\\|') |> str_remove('\\|[:print:]*$')]
    fasta_tmp = fasta[Uniprot==Uniprot_tmp,Seq]
    coverage_plot = data.table(peptide = 'Whole_prot',
                               position = 1:(stringr::str_length(fasta_tmp)))
    for(i in 1:nrow(mapping_tmp)){
        tmp_peptide =mapping_tmp[i,sequence] |> str_remove('^n')
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
isoforms = HRAS_psms[Gene=='HRAS',Protein] |> 
    unlist() |> stringr::str_extract('\\|[:print:]*?\\|') |> stringr::str_remove_all('\\|') |> unique()
dt = plot_coverage_prot('P01112',psm_percolator,isoform_fasta)
HRAS_peptidoforms = HRAS_psms[`Protein ID` %in% isoforms][,.(Peptide,`Delta Mass`)]
HRAS_peptidoforms[, mass_shift :=round(`Delta Mass`,digits = 1)]
HRAS_peptidoforms = HRAS_peptidoforms[!between(mass_shift,-0.5,2)][,.(N_pepts=.N), by = .(Peptide,mass_shift)
                                      ][N_pepts>1][,.(N_variation = .N), by = Peptide]
fwrite(HRAS_peptidoforms, here::here('Datasets','Processed','HRAS_peptidoforms.csv'))
HRAS_pepts = HRAS_psms[`Protein ID` %in% isoforms][,.(Peptide,`Delta Mass`,sample_name)]
HRAS_pepts[,`Delta Mass`:=round(`Delta Mass`,digits = 1)]
HRAS_pepts[,Mod_pept:= paste(Peptide,`Delta Mass`,sep = '_')]
HRAS_pepts = HRAS_pepts[,.(N_pepts=.N), by = .(Mod_pept,sample_name)]
HRAS_pepts[,total_pept_sample:=sum(N_pepts), by = sample_name]
HRAS_pepts[,denom := min(total_pept_sample)]
HRAS_pepts[,norm_pepts := round(N_pepts*denom/total_pept_sample), by = Mod_pept]
abundant = HRAS_pepts[,.(total_psms = sum(norm_pepts)), by = Mod_pept]
HRAS_pepts[,.(sample_name,total_pept_sample)] |> unique()
HRAS_pepts_mod_state = HRAS_pepts[,.(Mod_pept)
][,`:=`(state = str_remove_all(Mod_pept,'^[:print:]*_') |>
            as.numeric(),
        Pept = str_remove_all(Mod_pept,'_[:print:]*$'))] |> 
    unique()
HRAS_pepts_mod_state[,modified := !between(state,-0.5,2)]
HRAS_pepts_mod_state = HRAS_pepts_mod_state[,.(Pept,modified)] |> unique()
dt[,detected_state := fcase(
    peptide %in%   HRAS_pepts_mod_state[modified==T,Pept ] & 
        peptide %in% HRAS_pepts_mod_state[modified==F,Pept ], 'modified_and_unmodified',
    peptide %in%   HRAS_pepts_mod_state[modified==T,Pept ],'only_modified',
    peptide %in% HRAS_pepts_mod_state[modified==F,Pept ], 'only_unmodified',
    peptide == 'Whole_prot',NA_character_
), by = position]


HRAS_coverage_mod = ggplot(dt,aes(y =peptide ,x= position, colour = detected_state))+theme_bw()+ 
    geom_point()+theme(axis.text.x = element_blank(),
                       axis.ticks = element_blank())+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))+
     ggtitle('HRAS protein coverage at 98%')+
    scale_colour_manual(values = c('#A37261','#543B2D','#E3C1B5'))
ggsave(here::here('Output','HRAS_peptide_coverage.pdf'),HRAS_coverage_mod)
list_figures$HRAS_coverage_mod = HRAS_coverage_mod
# isoforms = c(isoforms,stringr::str_remove(isoforms,'-.')) |> na.omit() |> unique()
# HRAS_psms[`Protein ID` %in% isoforms][,.(Peptide,`Delta Mass`,`Assigned Modifications`,`Observed Modifications`, `MSFragger Localization`)] |>
#     unique() |> View()
# HRAS_psms[`Protein ID` %in% isoforms][stringr::str_detect(Peptide,'MPDEPEEPVVAVSSPAVPPPTK')][,.(sample,Peptide,`Delta Mass`,`Assigned Modifications`,`Observed Modifications`, `MSFragger Localization`,`Mapped Genes`)] |>
#     unique() |> View()
# isoform_relative_abundance = HRAS_psms[`Protein ID` %in% isoforms][,.(sample,`Mapped Proteins`)
# ][,.(N_iso = .N), by = .(sample,`Mapped Proteins`)
# ][,total_pepts:= sum(N_iso), by = sample
# ][, perc_isoform:= N_iso/total_pepts
# ][,.(sample,perc_isoform,`Mapped Proteins`)] |>
#     dcast(sample~`Mapped Proteins`, value.var = 'perc_isoform')
# isoform_relative_abundance = isoform_relative_abundance |> tibble::column_to_rownames('sample') |> as.matrix() |> t()
# isoform_relative_abundance[is.na(isoform_relative_abundance)]=0
# pheatmap::pheatmap(isoform_relative_abundance, scale = 'row')

# finding_similar_samples = HRAS_pepts |> dcast(Mod_pept ~ sample_name, value.var = 'N_pepts')
# finding_similar_samples = finding_similar_samples[,-1] |> as.matrix() |> is.finite() |> apply(2,as.numeric)
# pheatmap::pheatmap(finding_similar_samples)
coln_annot = data.table(samples = HRAS_pepts$sample_name |> unique())
coln_annot[,`:=`(Cell_line = stringr::str_remove_all(samples,'^[:print:]*_._'),
                 G12_state = fifelse(stringr::str_detect(samples,'wt'),'G12_wt','G12D_mut'))]
ann_colors = list(
    G12_state = c(G12_wt="grey50", G12D_mut="white"),
    Cell_line = c(CAL_33 = '#58C8F2',SCC_25   = 'white',HET_1A = '#EDA4B2'))
coln_annot= coln_annot |> tibble::column_to_rownames('samples')
finding_similar_samples = HRAS_pepts[Mod_pept %in% abundant[total_psms>30,Mod_pept]] |>
    dcast(Mod_pept ~ sample_name , value.var = 'norm_pepts')
# finding_similar_samples = finding_similar_samples[,-1] |> as.matrix() |> is.finite() |> apply(2,as.numeric)
heatmap_to_plot = finding_similar_samples[,-1] |> as.matrix() |> log2()
rownames(heatmap_to_plot) = finding_similar_samples[,1] |> unlist()
heatmap_to_plot[is.na(heatmap_to_plot)| is.infinite(heatmap_to_plot)] <- 0
pdf(here::here('Output','HRAS_ptms_heatmap.pdf'), width = 10, height = 8)
HRAS_heatmap = pheatmap::pheatmap(heatmap_to_plot, cluster_rows = F,scale = 'none', show_colnames = F, 
                   color = c('white',rev(monochromeR::generate_palette("darkred", modification = "go_darker", blend_colour = 'white',
                                                             n_colours = 40, view_palette = F))), 
                   annotation_col = coln_annot,annotation_colors = ann_colors,
                   main = '                         Modified peptides cluster the samples by cell line and HRAS mutation state')
dev.off()
list_figures$HRAS_heatmap <- HRAS_heatmap
fwrite(HRAS_pepts[norm_pepts>0],here::here('Datasets','Raw','HRAS_normalised_PSMs.tsv'))
location = 'C:\\Users\\skourtis\\OneDrive - CRG - Centre de Regulacio Genomica\\Bioinformatics Projects\\Proteoform_IPMS\\Datasets\\Processed\\'
saints_location  = 'C:\\Users\\skourtis\\Downloads\\SAINTexpress_Windows_executables_v3.4__2014_Aug_19\\SAINTexpress-spc'

for(i in unique(coln_annot$Cell_line)){
    
    cell_line_tmp = i
    # prey file
    Prey = fread(here::here('Datasets','Raw','HRAS_normalised_PSMs.tsv'))
    Prey = Prey[stringr::str_detect(sample_name,cell_line_tmp)]
    Prey[,length:= stringr::str_remove(Mod_pept,'_[:print:]*$') |> stringr::str_length(), by  = Mod_pept]
    Prey[,name:= Mod_pept]
    Prey = Prey[,.(Mod_pept,length,name)] |> unique()
    fwrite(Prey,here::here('Datasets','Processed',glue::glue('{cell_line_tmp}_prey.tsv')), sep = '\t', col.names = F)
    # interaction file
    Prey = fread(here::here('Datasets','Raw','HRAS_normalised_PSMs.tsv'))
    Prey = Prey[stringr::str_detect(sample_name,cell_line_tmp)]
    Prey[,sample_name :=stringr::str_remove(sample_name,glue::glue('_{cell_line_tmp}')) ]
    Prey[,IP_name := sample_name]
    Prey[,bait_name := stringr::str_remove(sample_name,'_.$')]
    Prey[,preyname := Mod_pept]
    Prey[,counts := norm_pepts]
    Prey = Prey[,.(IP_name,bait_name,preyname,counts)] 
    fwrite(Prey,here::here('Datasets','Processed',glue::glue('{cell_line_tmp}_interaction.tsv')), sep = '\t', col.names = F)
    
    # bait file
    for(type in c('forw','rev')){
        conditions= c('HRAS_g12d','HRAS_wt')
        if(type == 'rev'){
            conditions= rev(conditions)
        }
        Bait= data.table(Sample = paste(rep(conditions, each = 3), rep(1:3,times = 2),sep = '_'),
                         Condition = rep(c('HRAS_g12d','HRAS_wt'), each = 3),
                         Test = rep(c('T','C'),each = 3))
        fwrite(Bait,here::here('Datasets','Processed',glue::glue('{cell_line_tmp}_bait_{type}.tsv')), sep = '\t', col.names = F)
        command = glue::glue('{saints_location} -L4 "{location}{i}_interaction.tsv" "{location}{i}_prey.tsv" "{location}{i}_bait_{type}.tsv"')
        # file.copy(from = here::here('list.txt'),glue::glue('{location}\\{i}_{type}.tsv'))
        
        system(command)
        file.rename(here::here('list.txt'), glue::glue('{location}\\{i}_saints_results_{type}.tsv'))
    }
 
       
}

result_files =  list.files(location, pattern = 'saints_results') |>
    stringr::str_subset('forw')  |> 
    str_subset('HET|CAL|SCC')
significant_peptides = data.table()
for(file_n in result_files){
    cell_line = str_remove(file_n,'_[:print:]*$')
    cell_line_position = glue::glue('{cell_line}_position')
    proteoforms = rbind(fread(here::here('Datasets','Processed',file_n)),
                        fread(here::here('Datasets','Processed',stringr::str_replace(file_n,'forw','rev'))))
    proteoforms = proteoforms[order(BFDR )][,head(.SD,1), by = Prey]
    proteoforms[,LogFC:= log2(FoldChange)]
    baits = proteoforms$Bait |> unique() 
    first_bait = baits[1]
    second_bait = baits[2]
    proteoforms[,mut_pept:= stringr::str_detect(Prey,'LVVVGAGG')]
    proteoforms[,LogFC:= fifelse(Bait !=first_bait,LogFC*(-1),LogFC)]
    significant_peptides = rbind(significant_peptides,proteoforms[BFDR <0.01][,comparison:=cell_line])
    max_show = proteoforms$LogFC |> abs() |> sort() |> tail(8) |> min()
    peptides_mod = proteoforms[,.(LogFC,BFDR,Prey)] |> unique()
    peptides_mod[,peptide:=str_remove(Prey,'_[:print:]*$')]
    dt_comparison = dt[peptide %in% c(peptides_mod$peptide,'Whole_prot')]
    peptides_mod = peptides_mod[dt_comparison, on = 'peptide',allow.cartesian=TRUE]
    peptides_mod[,min_location:= min(as.numeric(position)), by = Prey]
    peptides_mod = peptides_mod[order(min_location,Prey)]
    peptides_mod[,Prey:= factor(Prey,levels = unique(peptides_mod$Prey) )]
    q = ggplot(peptides_mod,aes(y =Prey ,x= position, colour = LogFC, size = -log10(BFDR+0.001)))+theme_bw()+ 
        geom_point()+theme(axis.text.x = element_blank(),
                           axis.ticks = element_blank())+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.y = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))+
        ggtitle('KRAS protein coverage at 88%')+
        scale_colour_gradient2( mid = 'white',high = 'darkred',low = 'darkblue' )+
        scale_size_continuous(range = c(0.01,5))+
        ggtitle('Position of modified peptides relevant for isoform switching',
                subtitle = glue::glue('IP-MS between {first_bait} vs {second_bait}'))
    ggsave(here::here('Output',glue::glue('HRAS_G12D_{cell_line}_sequence_plot.pdf')),q)
    list_figures[[cell_line_position]] = q
    p = ggplot(proteoforms,aes(x = LogFC, y= -log10(BFDR+0.001), label = Prey, colour = mut_pept))+
        geom_point()+
        ggrepel::geom_label_repel(data = proteoforms[BFDR<0.05], alpha = 1)+
        theme_bw()+
        scale_colour_manual(values =c('grey20','darkred'))+
        ggtitle('IP-MS between HRAS G12D and WT',
                subtitle = cell_line)
    ggsave(here::here('Output',glue::glue('HRAS_G12D_{cell_line}_proteoform_volcano.pdf')),p)
    list_figures[[cell_line]] = p
}

fwrite(significant_peptides,here::here('Datasets','Processed','HRAS_study_significant_SAINTSexpress.csv'))
