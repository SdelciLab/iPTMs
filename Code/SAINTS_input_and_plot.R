# create SAINTSexpress files
library(data.table);library(ggplot2)
# https://ebi.ac.uk/pride/archive/projects/PXD019469
Samples =openxlsx::read.xlsx(here::here('Datasets','Raw','SampleAnnotations.xlsx')) |> as.data.table()
location = 'C:\\Users\\skourtis\\OneDrive - CRG - Centre de Regulacio Genomica\\Bioinformatics Projects\\Proteoform_IPMS\\Datasets\\Processed\\'
saints_location  = 'C:\\Users\\skourtis\\Downloads\\SAINTexpress_Windows_executables_v3.4__2014_Aug_19\\SAINTexpress-spc'
for(i in unique(Samples$Cell.Type)){
  
cell_line_tmp = i
Samples_cell_lines = Samples[`Cell.Type` == cell_line_tmp]
setnames(Samples_cell_lines,'Raw.File.Name','sample')
Samples_cell_lines[,sample := stringr::str_replace(sample,'-','_')]

# prey file
Prey = fread(here::here('Datasets','Raw','proteoform_pept_PSMs.tsv'))
Prey = Samples_cell_lines[Prey, on = 'sample', nomatch = NULL] 
Prey[,length:= stringr::str_remove(Mod_pept,'_[:print:]*$') |> stringr::str_length(), by  = Mod_pept]
Prey[,name:= Mod_pept]
Prey = Prey[,.(Mod_pept,length,name)] |> unique()
fwrite(Prey,here::here('Datasets','Processed',glue::glue('{cell_line_tmp}_prey.tsv')), sep = '\t', col.names = F)
# interaction file
Prey = fread(here::here('Datasets','Raw','proteoform_pept_PSMs.tsv'))
Prey = Samples_cell_lines[Prey, on = 'sample', nomatch = NULL] 
Prey[,IP_name := stringr::str_replace(Biological.Replicate,'-','_')]
Prey[,bait_name := Bait]
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
    file.rename(here::here('list.txt'), glue::glue('{location}\\{i}_{type}.tsv'))
}

}
    for(cell_line in unique(Samples$Cell.Type)){
proteoforms = rbind(fread(here::here('Datasets','Processed',glue::glue('{cell_line}_forw.tsv'))),
                    fread(here::here('Datasets','Processed',glue::glue('{cell_line}_rev.tsv'))))
proteoforms = proteoforms[order(BFDR )][,head(.SD,1), by = Prey]
proteoforms[,LogFC:= log2(FoldChange)]
proteoforms[,mut_pept:= stringr::str_detect(Prey,'LVVVGAGG')]
proteoforms[,LogFC:= fifelse(Bait =='HRAS_wt',LogFC*(-1),LogFC)]
max_show = proteoforms$LogFC |> abs() |> sort() |> tail(8) |> min()
p = ggplot(proteoforms,aes(x = LogFC, y= -log10(BFDR+0.001), label = Prey, colour = mut_pept))+
    geom_point()+
    ggrepel::geom_label_repel(data = proteoforms[abs(LogFC)>max_show & BFDR<0.01], alpha = 1)+
    theme_bw()+
    scale_colour_manual(values =c('grey20','darkred'))+
    ggtitle('IP-MS between HRAS G12D and WT',
            subtitle = cell_line)
ggsave(here::here('Output',glue::glue('{cell_line}_proteoform_volcano.pdf')),p)
}

 
      
      
      