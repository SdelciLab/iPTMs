# create SAINTSexpress files
library(data.table);library(ggplot2)
# https://ebi.ac.uk/pride/archive/projects/PXD019469
Samples =openxlsx::read.xlsx(here::here('Datasets','Raw','SampleAnnotations.xlsx')) |> as.data.table()
location = 'C:\\Users\\skourtis\\OneDrive - CRG - Centre de Regulacio Genomica\\Bioinformatics Projects\\Proteoform_IPMS\\Datasets\\Processed\\'
saints_location  = 'C:\\Users\\skourtis\\Downloads\\SAINTexpress_Windows_executables_v3.4__2014_Aug_19\\SAINTexpress-spc'

cell_line_tmp = 'MDAMB231'

# prey file
Prey = fread(here::here('Datasets','Raw','HNRNPC_proteoform_pept_PSMs.tsv'))
# Prey = Prey[stringr::str_detect(sample,'03|06',negate = T)]
Prey[,length:= stringr::str_remove(Mod_pept,'_[:print:]*$') |> stringr::str_length(), by  = Mod_pept]
Prey[,name:= Mod_pept]
Prey = Prey[,.(Mod_pept,length,name)] |> unique()
fwrite(Prey,here::here('Datasets','Processed',glue::glue('{cell_line_tmp}_prey.tsv')), sep = '\t', col.names = F)
# interaction file
Prey = fread(here::here('Datasets','Raw','HNRNPC_proteoform_pept_PSMs.tsv'))
# Prey = Prey[stringr::str_detect(sample,'03|06',negate = T)]
Prey[,IP_name := fcase(
    stringr::str_detect(sample,'1$'),'MDAMB231_1',
    stringr::str_detect(sample,'2$'),'MDAMB231_2',
    stringr::str_detect(sample,'3$'),'MDAMB231_3',
    stringr::str_detect(sample,'4$'), 'LM2_1',
    stringr::str_detect(sample,'5$'),'LM2_2',
    stringr::str_detect(sample,'6$'),'LM2_3'
    
)]
Prey[,bait_name := stringr::str_remove(IP_name,'_.')]
Prey[,preyname := Mod_pept]
Prey[,counts := norm_pepts]
Prey = Prey[,.(IP_name,bait_name,preyname,counts)] 

fwrite(Prey,here::here('Datasets','Processed',glue::glue('{cell_line_tmp}_interaction.tsv')), sep = '\t', col.names = F)

# bait file
for(type in c('forw','rev')){
    conditions= c('MDAMB231','LM2')
    if(type == 'rev'){
        conditions= rev(conditions)
    }
    Bait= data.table(Sample = paste(rep(conditions, each = 3), rep(1:3,times = 2),sep = '_'),
                     Condition = rep(c('MDAMB231','LM2'), each = 3),
                     Test = rep(c('T','C'),each = 3))
    # Bait = Bait[stringr::str_detect(Sample,'_3',negate = T)]
    fwrite(Bait,here::here('Datasets','Processed',glue::glue('{cell_line_tmp}_bait_{type}.tsv')), sep = '\t', col.names = F)
    command = glue::glue('{saints_location} -L4 "{location}{cell_line_tmp}_interaction.tsv" "{location}{cell_line_tmp}_prey.tsv" "{location}{cell_line_tmp}_bait_{type}.tsv"')
    # file.copy(from = here::here('list.txt'),glue::glue('{location}\\{i}_{type}.tsv'))
    
    system(command)
    file.rename(here::here('list.txt'), glue::glue('{location}\\{cell_line_tmp}_{type}.tsv'))
}
cell_line = cell_line_tmp
proteoforms = rbind(fread(here::here('Datasets','Processed',glue::glue('{cell_line}_forw.tsv'))),
                    fread(here::here('Datasets','Processed',glue::glue('{cell_line}_rev.tsv'))))
proteoforms = proteoforms[order(BFDR )][,head(.SD,1), by = Prey]
proteoforms[,LogFC:= log2(FoldChange)]
proteoforms[,LogFC:= fifelse(Bait =='MDAMB231',LogFC*(-1),LogFC)]
# max_show = proteoforms$LogFC |> abs() |> sort() |> tail(8) |> min()
p = ggplot(proteoforms,aes(x = LogFC, y= -log10(BFDR+0.001), label = Prey))+
    geom_point()+
    ggrepel::geom_label_repel(data = proteoforms[BFDR <0.01], max.overlaps = 100)+
    theme_bw()+
    ggtitle('IP-MS between HNRNPC MDAMB231 and LM2',
            subtitle = cell_line)
ggsave(here::here('Output',glue::glue('{cell_line}_proteoform_volcano.pdf')),p, width = 10,height = 10)


 
      
      
      