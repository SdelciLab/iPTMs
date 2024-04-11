# Analysis of HRAS Mut data
library(data.table);library(stringr)

BRD4_psms = fread(here::here('Datasets','Raw','BRD4_MTHFD1_targets_psms.tsv'))
BRD4_psms[str_detect(Peptide,'SEPFSPSL'),.(sample,`Delta Mass`,`Observed Modifications`,`MSFragger Localization`)] |> View() 
BRD4_psms[,sample:= stringr::str_remove_all(sample,'^[:print:]*\\/\\/') |> stringr::str_remove_all('\\/psm.tsv')]

Samples = data.table(`Raw.File.Name` = unique(BRD4_psms$sample))
Samples[,`:=`(sample_name = `Raw.File.Name` |> 
                  stringr::str_remove_all('^[:print:]*poolLys[:print:]*?_') |> 
                  stringr::str_remove_all('rep.') )]
setnames(Samples,'Raw.File.Name','sample')
Samples =Samples[order(sample_name)]
Samples[,sample_name := paste(rep(c('MEG','MOLM','MV4','IgG','K562'),each =3),rep(1:3,times = 5),sep = '_')]
BRD4_psms = Samples[,.(sample,sample_name)][BRD4_psms, on = 'sample']
Psm_per_prot = BRD4_psms[Gene != ''][,.(N_psms= .N), by = .(Gene,sample_name)]
Psm_per_prot[, bait := fifelse(Gene =='BRD4','bait','prey')]
# Psm_per_prot[,inhibitor := as.factor(fifelse(stringr::str_detect(sample_name,'i'), 'inh','control'))]
Psm_per_prot[,BRD4_state := as.factor(stringr::str_remove(sample_name,'_.'))]
BRD4_enrichment =  Psm_per_prot |> ggplot(aes(x = log10(N_psms), y= sample_name, fill = BRD4_state ))+
    ggridges::geom_density_ridges(alpha = 0.4,  aes(colour = BRD4_state),rel_min_height=0.01)+theme_bw()+
    geom_point(data = Psm_per_prot[bait == 'bait'], aes(colour = BRD4_state),size = 2.5)+
    facet_wrap('BRD4_state',nrow = 5,scales = 'free_y')+
    scale_fill_manual(values = c('#C39A96','#A95C26','#D18547','#C6A886','#7F522B'))+
    scale_colour_manual(values = c('#C39A96','#A95C26','#D18547','#C6A886','#7F522B'))+
    ggtitle('Bait BRD4 is enriched with numerous Psms in IP-MS samples')
ggsave(here::here('Output','BRD4_psm_enrichment.pdf'),BRD4_enrichment)
list_figures$BRD4_enrichment =BRD4_enrichment
isoform_fasta <- seqinr::read.fasta(here::here('Datasets','Raw','micro_prot_decoy.fasta'),
                                    seqtype = "AA",as.string = TRUE)
isoform_fasta = data.table(Annotation  = unlist(seqinr::getAnnot(isoform_fasta)) |>  str_remove(' [:print:]*$') |> str_remove('>'),
                           Seq = seqinr::getSequence(isoform_fasta, as.string = T) |>  unlist())
isoform_fasta= isoform_fasta[str_detect(Annotation,'rev',negate = T)]
isoform_fasta[,Annotation:= stringr::str_remove_all(Annotation,'^[:print:]*?\\|') |> stringr::str_remove_all('\\|[:print:]*?$')]
Uniprots_FDR =BRD4_psms$`Protein ID` |> unique()

Protein_coverage = isoform_fasta[Annotation %in% Uniprots_FDR]

psm_percolator = BRD4_psms[,.(Peptide,`Protein ID`)] |> unique()
names(psm_percolator) = c('sequence','Uniprot')
psm_percolator[,sequence:= str_remove(sequence,'^n')]
count_coverage <- function(Annotation,sequence){
    
    # Annotation = 'O60885'
    Annotation = Annotation |> str_remove('^[:print:]*?\\|') |> str_remove('\\|[:print:]*') 
    # sequence  = Protein_coverage[Annotation =='O60885',Seq] |> unlist()
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
Protein_coverage[,perc_coverage:= count_coverage(Annotation,Seq),by = Annotation]
BRD4_coverage = Protein_coverage[Annotation == 'O60885',perc_coverage]*100
BRD4_coerage = Protein_coverage |> 
    ggplot(aes(x= perc_coverage*100))+
    theme_bw()+
    geom_histogram(binwidth =1, alpha  = 0.5)+
    geom_vline(xintercept = BRD4_coverage, colour = 'black', linetype="longdash" )+
    ggtitle('BRD4 bait has a high protein coverage')
ggsave(here::here('Output','BRD4i_relative_perc_coverage.pdf'),BRD4_coerage)
list_figures$BRD4_coverage = BRD4_coerage
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
isoforms = BRD4_psms[Gene=='BRD4',Protein] |> 
    unlist() |> stringr::str_extract('\\|[:print:]*?\\|') |> stringr::str_remove_all('\\|') |> unique()

dt = plot_coverage_prot('O60885',psm_percolator,isoform_fasta)
BRD4_peptidoforms = BRD4_psms[`Protein ID` %in% isoforms][,.(Peptide,`Delta Mass`)]
BRD4_peptidoforms[, mass_shift :=round(`Delta Mass`,digits = 1)]
BRD4_peptidoforms = BRD4_peptidoforms[!between(mass_shift,-0.5,2)][,.(N_pepts=.N), by = .(Peptide,mass_shift)
][N_pepts>1][,.(N_variation = .N), by = Peptide]
fwrite(BRD4_peptidoforms, here::here('Datasets','Processed','BRD4_peptidoforms'))


BRD4_pepts = BRD4_psms[`Protein ID` %in% isoforms][,.(Peptide,`Delta Mass`,sample_name)]
BRD4_pepts[,`Delta Mass`:=round(`Delta Mass`,digits = 1)]
BRD4_pepts[,Mod_pept:= paste(Peptide,`Delta Mass`,sep = '_')]
BRD4_pepts = BRD4_pepts[,.(N_pepts=.N), by = .(Mod_pept,sample_name)]
# BRD4_pepts[,mut:= str_detect(sample_name,'C12')]
BRD4_pepts[,total_pept_sample:=sum(N_pepts), by = .(sample_name)]
BRD4_pepts = BRD4_pepts[str_detect(sample_name, 'IgG', negate = T)]
BRD4_pepts[,denom := min(total_pept_sample)]
BRD4_pepts[,norm_pepts := round(N_pepts*denom/total_pept_sample,digits = 0), by = Mod_pept]
abundant = BRD4_pepts[,.(total_psms = sum(norm_pepts)), by = Mod_pept]
BRD4_pepts[,.(sample_name,total_pept_sample)] |> unique()
BRD4_pepts[,.(total_psms = sum(norm_pepts)), by = sample_name]
BRD4_pepts_mod_state = BRD4_pepts[,.(Mod_pept)
                                  ][,`:=`(state = str_remove_all(Mod_pept,'^[:print:]*_') |>
                                              as.numeric(),
                                          Pept = str_remove_all(Mod_pept,'_[:print:]*$'))] |> 
    unique()
BRD4_pepts_mod_state[,modified := !between(state,-0.5,2)]
BRD4_pepts_mod_state = BRD4_pepts_mod_state[,.(Pept,modified)] |> unique()
dt[,detected_state := fcase(
  peptide %in%   BRD4_pepts_mod_state[modified==T,Pept ] & 
      peptide %in% BRD4_pepts_mod_state[modified==F,Pept ], 'modified_and_unmodified',
  peptide %in%   BRD4_pepts_mod_state[modified==T,Pept ],'only_modified',
  peptide %in% BRD4_pepts_mod_state[modified==F,Pept ], 'only_unmodified',
  peptide == 'Whole_prot',NA_character_
), by = position]

BRD4_mod_pept = ggplot(dt,aes(y =peptide ,x= position, colour = detected_state))+theme_bw()+ 
    geom_point()+theme(axis.text.x = element_blank(),
                       axis.ticks = element_blank())+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))+
     ggtitle(glue::glue('BRD4 protein coverage at {round(BRD4_coverage)}%'))+
    scale_colour_manual(values = c('#A37261','#543B2D','#E3C1B5'))
ggsave(here::here('Output','BRD4_prot_coverage_modification_state.pdf'),BRD4_mod_pept)
list_figures$BRD4_mod_pept = BRD4_mod_pept
# isoforms = c(isoforms,stringr::str_remove(isoforms,'-.')) |> na.omit() |> unique()
# BRD4_psms[`Protein ID` %in% isoforms][,.(Peptide,`Delta Mass`,`Assigned Modifications`,`Observed Modifications`, `MSFragger Localization`)] |>
#     unique() |> View()
# BRD4_psms[`Protein ID` %in% isoforms][stringr::str_detect(Peptide,'MPDEPEEPVVAVSSPAVPPPTK')][,.(sample,Peptide,`Delta Mass`,`Assigned Modifications`,`Observed Modifications`, `MSFragger Localization`,`Mapped Genes`)] |>
#     unique() |> View()
# isoform_relative_abundance = BRD4_psms[`Protein ID` %in% isoforms][,.(sample,`Mapped Proteins`)
# ][,.(N_iso = .N), by = .(sample,`Mapped Proteins`)
# ][,total_pepts:= sum(N_iso), by = sample
# ][, perc_isoform:= N_iso/total_pepts
# ][,.(sample,perc_isoform,`Mapped Proteins`)] |>
#     dcast(sample~`Mapped Proteins`, value.var = 'perc_isoform')
# isoform_relative_abundance = isoform_relative_abundance |> tibble::column_to_rownames('sample') |> as.matrix() |> t()
# isoform_relative_abundance[is.na(isoform_relative_abundance)]=0
# pheatmap::pheatmap(isoform_relative_abundance, scale = 'row')

# finding_similar_samples = BRD4_pepts |> dcast(Mod_pept ~ sample_name, value.var = 'N_pepts')
# finding_similar_samples = finding_similar_samples[,-1] |> as.matrix() |> is.finite() |> apply(2,as.numeric)
# pheatmap::pheatmap(finding_similar_samples)
coln_annot = data.table(samples = BRD4_pepts[,sample_name] |> unique())
coln_annot[,`:=`(cell_line = str_remove(samples,'_.$'))]
ann_colors = list(
    cell_line = c(IgG = '#C39A96', MV4= '#A95C26', MOLM = '#D18547', K562 = '#C6A886',MEG = '#7F522B'))
coln_annot= coln_annot |> tibble::column_to_rownames('samples')
finding_similar_samples = BRD4_pepts[Mod_pept %in% abundant[total_psms>14,Mod_pept]][str_detect(sample_name,'IgG',negate = T)] |>
    dcast(Mod_pept ~ sample_name , value.var = 'norm_pepts')
# finding_similar_samples = finding_similar_samples[,-1] |> as.matrix() |> is.finite() |> apply(2,as.numeric)
heatmap_to_plot = finding_similar_samples[,-1] |> as.matrix() |> log2()
rownames(heatmap_to_plot) = finding_similar_samples[,1] |> unlist()
heatmap_to_plot[is.na(heatmap_to_plot) |is.infinite(heatmap_to_plot) ] <- 0
pdf(here::here('Output','BRD4_ptms_heatmap.pdf'), width = 10, height = 8)
BRD4_heatmap = pheatmap::pheatmap(heatmap_to_plot, cluster_rows = F,scale = 'none', show_colnames = F, 
                   color = c('white',rev(monochromeR::generate_palette("darkred", modification = "go_darker", blend_colour = 'white',
                                                                       n_colours = 40, view_palette = F))), 
                   annotation_col = coln_annot,annotation_colors = ann_colors,
                   main = '                         Modified peptides cluster the samples by cell_line mutation state')
dev.off()
list_figures$BRD4_heatmap= BRD4_heatmap
fwrite(BRD4_pepts[norm_pepts>0],here::here('Datasets','Raw','BRD4_normalised_PSMs.tsv'))
location = 'C:\\Users\\skourtis\\OneDrive - CRG - Centre de Regulacio Genomica\\Bioinformatics Projects\\Proteoform_IPMS\\Datasets\\Processed\\'
saints_location  = 'C:\\Users\\skourtis\\Downloads\\SAINTexpress_Windows_executables_v3.4__2014_Aug_19\\SAINTexpress-spc'


cell_line_tmp = 'BRD4s'
Samples[,bait_name:=str_remove(sample_name,'_.$')]
conditions = Samples$bait_name |> unique()
conditions = str_subset(conditions,'IgG',negate = T)
conditions_comparison = expand.grid(conditions,conditions) |> as.data.table()
conditions_comparison = conditions_comparison[as.numeric(Var1) > as.numeric(Var2)]
for(row_n in 1:nrow(conditions_comparison)){
    conditions= conditions_comparison[row_n] |> as.matrix() |> as.vector()

    # prey file
    # Prey = Prey[stringr::str_detect(sample,'03|06',negate = T)]
    Prey = fread(here::here('Datasets','Raw','BRD4_normalised_PSMs.tsv'))
    
    Prey[,length:= stringr::str_remove(Mod_pept,'_[:print:]*$') |> stringr::str_length(), by  = Mod_pept]
    Prey[,name:= Mod_pept]
    Prey = Prey[,.(Mod_pept,length,name)] |> unique()
    fwrite(Prey,here::here('Datasets','Processed',glue::glue('{cell_line_tmp}_prey.tsv')), sep = '\t', col.names = F)
    # interaction file
    Prey = fread(here::here('Datasets','Raw','BRD4_normalised_PSMs.tsv'))
    Prey[,`:=`(bait_name = sample_name |> str_remove('_.$'),
               IP_name = sample_name)]
    Prey[,preyname := Mod_pept]
    Prey[,counts := norm_pepts]
    Prey = Prey[,.(IP_name,bait_name,preyname,counts)] 
    Prey= Prey[bait_name %in% conditions]
    Prey = Prey[,.(counts = sum(counts)), by = .(IP_name,bait_name,preyname)]
    fwrite(Prey,here::here('Datasets','Processed',glue::glue('{cell_line_tmp}_{row_n}_interaction.tsv')), sep = '\t', col.names = F)
    
    # bait file
    
    for(type in c('forw','rev')){
        # conditions= c('MDAMB231','LM2')
        if(type == 'rev'){
            conditions_tmp= rev(conditions)
        }else{
            conditions_tmp = conditions
        }
        Bait= data.table(Sample = paste(rep(conditions_tmp, each = 3), rep(1:3,times = 2),sep = '_'),
                         Condition = rep(conditions, each = 3),
                         Test = rep(c('T','C'),each = 3))
        # Bait = Bait[stringr::str_detect(Sample,'_3',negate = T)]
        fwrite(Bait,here::here('Datasets','Processed',glue::glue('{cell_line_tmp}_bait_{row_n}_{type}.tsv')), sep = '\t', col.names = F)
        command = glue::glue('{saints_location} -L4 "{location}{cell_line_tmp}_{row_n}_interaction.tsv" "{location}{cell_line_tmp}_prey.tsv" "{location}{cell_line_tmp}_bait_{row_n}_{type}.tsv"')
        # file.copy(from = here::here('list.txt'),glue::glue('{location}\\{i}_{type}.tsv'))
        
        system(command)
        file.rename(here::here('list.txt'), glue::glue('{location}\\{cell_line_tmp}_{row_n}_saints_results_{type}.tsv'))
    }
}


result_files =  list.files(location, pattern = 'saints_results') |>
    stringr::str_subset('forw') |> 
    str_subset(cell_line_tmp)
for(file_n in result_files){
    cell_line = str_remove(file_n,'_saints_results_forw.tsv')
    cell_line_pos = glue::glue('{cell_line}_pos')
    proteoforms = rbind(fread(here::here('Datasets','Processed',file_n)),
                        fread(here::here('Datasets','Processed',stringr::str_replace(file_n,'forw','rev'))))
    proteoforms = proteoforms[order(BFDR )][,head(.SD,1), by = Prey]
    proteoforms[,LogFC:= log2(FoldChange)]
    baits = proteoforms$Bait |> unique() 
    first_bait = baits[1]
    second_bait = baits[2]
    # if(file_n %in% result_files[c(1,6)]){
    proteoforms[,mut_pept:= stringr::str_detect(Prey,'_80') & BFDR<0.01 ]
    proteoforms[,LogFC:= fifelse(Bait !=first_bait,LogFC*(-1),LogFC)]
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
    ggsave(here::here('Output',glue::glue('{first_bait}_{second_bait}_sequence_plot.pdf')),q)
    list_figures[[cell_line_pos]] = q
    p = ggplot(proteoforms,aes(x = LogFC, y= -log10(BFDR+0.001), label = Prey, colour = mut_pept))+
        geom_jitter()+
        ggrepel::geom_label_repel(data = proteoforms[BFDR<0.05 &mut_pept == T ], alpha = 1, max.overlaps = 10)+
        theme_bw()+
        scale_colour_manual(values =c('grey50','darkred'))+
        ggtitle(glue::glue('IP-MS between  {first_bait}vs{second_bait} '))
    ggsave(here::here('Output',glue::glue('{first_bait}_{second_bait}_proteoform_volcano.pdf')),p)
    list_figures[[cell_line]] = p
}


