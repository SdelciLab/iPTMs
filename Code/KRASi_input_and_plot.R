# create SAINTSexpress files
library(data.table);library(ggplot2)
# https://ebi.ac.uk/pride/archive/projects/PXD019469
location = 'C:\\Users\\skourtis\\OneDrive - CRG - Centre de Regulacio Genomica\\Bioinformatics Projects\\Proteoform_IPMS\\Datasets\\Processed\\'
saints_location  = 'C:\\Users\\skourtis\\Downloads\\SAINTexpress_Windows_executables_v3.4__2014_Aug_19\\SAINTexpress-spc'

Prey = fread(here::here('Datasets','Raw','KRASi_C12_proteoform_pept_PSMs.tsv'))
Prey = Prey[stringr::str_detect(sample,'_IP_')]
Samples =data.table(sample = unique(Prey$sample))
Samples[,bait_name:= sample |> stringr::str_remove('^[:print:]*?_ANolan_') |> stringr::str_remove('._S1_[:print:]*$') |> 
            stringr::str_remove('_$')]
Samples[,replicate:= stringr::str_extract(sample,'._S1') |> stringr::str_remove('_S1')]
Samples[,IP_name := paste(bait_name,replicate,sep = '_') ]
Samples = Samples[stringr::str_detect(IP_name,'EV',negate = T)]
cell_line_tmp = 'HEK_kras_C12'
BRD = c('')
# for(BRD in BRDs){
    conditions = Samples$bait_name |> stringr::str_subset(BRD) |> unique()
    conditions_comparison = expand.grid(conditions,conditions) |> as.data.table()
    conditions_comparison = conditions_comparison[as.numeric(Var1) > as.numeric(Var2)]
    for(row_n in 1:nrow(conditions_comparison)){
        conditions= conditions_comparison[row_n] |> as.matrix() |> as.vector()
# prey file
# Prey = Prey[stringr::str_detect(sample,'03|06',negate = T)]
        Prey = fread(here::here('Datasets','Raw','KRASi_C12_proteoform_pept_PSMs.tsv'))
        
Prey[,length:= stringr::str_remove(Mod_pept,'_[:print:]*$') |> stringr::str_length(), by  = Mod_pept]
Prey[,name:= Mod_pept]
Prey = Prey[,.(Mod_pept,length,name)] |> unique()
fwrite(Prey,here::here('Datasets','Processed',glue::glue('{cell_line_tmp}_prey.tsv')), sep = '\t', col.names = F)
# interaction file
Prey = fread(here::here('Datasets','Raw','KRASi_C12_proteoform_pept_PSMs.tsv'))
Prey = Prey[Samples,on = 'sample']
Prey[,preyname := Mod_pept]
Prey[,counts := norm_pepts]
Prey = Prey[,.(IP_name,bait_name,preyname,counts)] 
Prey= Prey[bait_name %in% conditions]
Prey = Prey[,.(counts = sum(counts)), by = .(IP_name,bait_name,preyname)]
fwrite(Prey,here::here('Datasets','Processed',glue::glue('{cell_line_tmp}_{BRD}_{row_n}_interaction.tsv')), sep = '\t', col.names = F)

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
    fwrite(Bait,here::here('Datasets','Processed',glue::glue('{cell_line_tmp}_bait_{BRD}_{row_n}_{type}.tsv')), sep = '\t', col.names = F)
    command = glue::glue('{saints_location} -L4 "{location}{cell_line_tmp}_{BRD}_{row_n}_interaction.tsv" "{location}{cell_line_tmp}_prey.tsv" "{location}{cell_line_tmp}_bait_{BRD}_{row_n}_{type}.tsv"')
    # file.copy(from = here::here('list.txt'),glue::glue('{location}\\{i}_{type}.tsv'))
    
    system(command)
    file.rename(here::here('list.txt'), glue::glue('{location}\\{cell_line_tmp}_{BRD}_{row_n}_{type}_results.tsv'))
}
    }


    result_files =  list.files(location, pattern = 'results') |>
        stringr::str_subset('forw') |> 
        stringr::str_subset('HEK_kras')
PTMS = data.table()
for(file_n in result_files){

proteoforms = rbind(fread(here::here('Datasets','Processed',file_n)),
                    fread(here::here('Datasets','Processed',file_n |> stringr::str_replace('forw','rev'))))
proteoforms = proteoforms[order(BFDR )][,head(.SD,1), by = Prey]
proteoforms[,LogFC:= log2(FoldChange)]
baits = proteoforms$Bait |> unique() 
first_bait = baits[1]
second_bait = baits[2]
proteoforms[,LogFC:= fifelse(Bait ==first_bait,LogFC*(-1),LogFC)]
# max_show = proteoforms$LogFC |> abs() |> sort() |> tail(8) |> min()
PTMS_tmp =proteoforms[BFDR<0.02,Prey] |> stringr::str_remove_all('_[:print:]*$') |> table() |> 
    tibble::enframe(name = 'Prey', value = 'Number_significant') |> as.data.frame() |>  as.data.table()
PTMS_tmp[,comparison := file_n]
print(PTMS_tmp |> colnames())
PTMS = rbind(PTMS,PTMS_tmp)
p = ggplot(proteoforms,aes(x = LogFC, y= -log10(BFDR+0.001), label = Prey))+
    geom_jitter()+
    ggrepel::geom_label_repel(data = proteoforms[BFDR <0.01], max.overlaps = 10)+
    ggrepel::geom_label_repel(data = proteoforms[ stringr::str_detect(Prey,'LVVVGAGGVGK_(0|46)')], max.overlaps = 100, colour = 'red')+
    theme_bw()+
    ggtitle(glue::glue('IP-MS between {second_bait} vs {first_bait} '))
ggsave(here::here('Output',glue::glue('{second_bait}_{first_bait}_proteoform_volcano.pdf')),p)
}
fwrite(PTMS, here::here('Datasets','Processed',glue::glue('hits_{file_n}')))

 
      
      
      