# create SAINTSexpress files
library(data.table);library(ggplot2)
# https://ebi.ac.uk/pride/archive/projects/PXD019469
Samples =fread(here::here('Datasets','Raw','README.csv')) |> as.data.table()
location = 'C:\\Users\\skourtis\\OneDrive - CRG - Centre de Regulacio Genomica\\Bioinformatics Projects\\Proteoform_IPMS\\Datasets\\Processed\\'
saints_location  = 'C:\\Users\\skourtis\\Downloads\\SAINTexpress_Windows_executables_v3.4__2014_Aug_19\\SAINTexpress-spc'
Samples[,sample:= stringr::str_remove(FileName,'.raw$') |> stringr::str_remove('_0._10pto$')]
cell_line_tmp = 'A549'

# prey file
Prey = fread(here::here('Datasets','Raw','KRAS_proteoform_pept_PSMs.tsv'))
Prey[,sample:= stringr::str_remove(sample,'_0._10pto$')]
Prey = Samples[Prey, on = 'sample', nomatch = NULL]
Prey[,length:= stringr::str_remove(Mod_pept,'_[:print:]*$') |> stringr::str_length(), by  = Mod_pept]
Prey[,name:= Mod_pept]
Prey = Prey[,.(Mod_pept,length,name)] |> unique()
fwrite(Prey,here::here('Datasets','Processed',glue::glue('{cell_line_tmp}_prey.tsv')), sep = '\t', col.names = F)
# interaction file
Prey = fread(here::here('Datasets','Raw','KRAS_proteoform_pept_PSMs.tsv'))
Prey[,sample:= stringr::str_remove(sample,'_0._10pto$')]
Prey = Samples[Prey, on = 'sample', nomatch = NULL] 
Prey[,Replicate := fcase(
    Replicate %in% c(1,4),1,
    Replicate %in% c(2, 5),2,
    Replicate %in% c(3,6),3
)]
Prey[,IP_name := paste(Condition,Replicate,sep = '_'),by=FileName]
Prey[,bait_name := Condition]
Prey[,preyname := Mod_pept]
Prey[,counts := norm_pepts]
Prey = Prey[,.(IP_name,bait_name,preyname,counts)] 
fwrite(Prey,here::here('Datasets','Processed',glue::glue('{cell_line_tmp}_interaction.tsv')), sep = '\t', col.names = F)

# bait file
for(type in c('forw','rev')){
    conditions= c('KRAS_cytosol','KRAS_chromatin')
    if(type == 'rev'){
        conditions= rev(conditions)
    }
    Bait= data.table(Sample = paste(rep(conditions, each = 3), rep(1:3,times = 2),sep = '_'),
                     Condition = rep(c('KRAS_cytosol','KRAS_chromatin'), each = 3),
                     Test = rep(c('T','C'),each = 3))
    fwrite(Bait,here::here('Datasets','Processed',glue::glue('{cell_line_tmp}_bait_{type}.tsv')), sep = '\t', col.names = F)
    command = glue::glue('{saints_location} -L4 "{location}{cell_line_tmp}_interaction.tsv" "{location}{cell_line_tmp}_prey.tsv" "{location}{cell_line_tmp}_bait_{type}.tsv"')
    system(command)
    file.rename(here::here('list.txt'), glue::glue('{location}\\{cell_line_tmp}_{type}.tsv'))
}


cell_line = 'A549'
proteoforms = rbind(fread(here::here('Datasets','Processed',glue::glue('{cell_line}_forw.tsv'))),
                    fread(here::here('Datasets','Processed',glue::glue('{cell_line}_rev.tsv'))))
proteoforms = proteoforms[order(BFDR )][,head(.SD,1), by = Prey]
proteoforms[,LogFC:= log2(FoldChange)]
proteoforms[,LogFC:= fifelse(Bait =='KRAS_cytosol',LogFC*(-1),LogFC)]
max_show = proteoforms$LogFC |> abs() |> sort() |> tail(8) |> min()
p = ggplot(proteoforms,aes(x = LogFC, y= -log10(BFDR+0.001), label = Prey))+
    geom_point()+
    ggrepel::geom_label_repel(data = proteoforms[abs(LogFC)>max_show])+
    theme_bw()+
    ggtitle('IP-MS between HRAS G12D and WT',
            subtitle = cell_line)
ggsave(here::here('Output',glue::glue('{cell_line}_proteoform_volcano.pdf')),p)

 
      
      
      