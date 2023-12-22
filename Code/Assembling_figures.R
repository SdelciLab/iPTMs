# Assembling figures
library("grid")
library("gridExtra")
library("pheatmap")
library("ggplot2")
theme_set(theme_bw(base_size = 10))
Fig_1 = ggpubr::ggarrange(plotlist  = list(
    ggplot(),
    list_figures$HRAS_enrichment+labs(title = NULL, subtitle = NULL)+
        xlab('Number of Psms per Protein (log10)')+ ylab('Sample replicates per cell line')+
                                               theme(legend.position="none",
                                                     axis.text.y=element_blank(),
                                                     axis.ticks.y=element_blank()),
                               list_figures$HRAS_heatmap[[4]],
                               list_figures$CAL+ ylab('-log10(BFDR)')+ 
        labs(title = NULL, subtitle = NULL)+
                                   theme(legend.position="none"),
                               list_figures$HET+ylab('-log10(BFDR)')+
        labs(title = NULL, subtitle = NULL)+
                                   theme(legend.position="none"),
                               list_figures$SCC+ylab('-log10(BFDR)')+
        labs(title = NULL, subtitle = NULL)+
                                   theme(legend.position="none")), nrow = 3, ncol = 2)
ggsave(here::here('Output','Fig_1.pdf'), Fig_1, width = 210, height = 297, units = 'mm')
Fig_S1 = ggpubr::ggarrange(plotlist  = 
                               list(list_figures$HRAS_coverage+
                                        labs(title = NULL, subtitle = NULL)+
                                        xlab('Protein Percentage coverage'),ggplot(),ggplot(),
                           list_figures$HRAS_coverage_mod+
                               theme(legend.position="bottom",
                                     axis.text.y=element_blank(),
                                     axis.ticks.y=element_blank())+
                               labs(title = NULL, subtitle = NULL)+
                               xlab('Peptide relative position'),ggplot(),ggplot(),
                           list_figures$HET_position+
                               theme(legend.position="bottom")+
                               labs(title = NULL, subtitle = NULL)+
                               xlab('Peptide relative position'),ggplot(),ggplot()), nrow = 3, ncol = 3,align = 'none')
ggsave(here::here('Output','Fig_S1.pdf'), Fig_S1, width = 210, height = 297, units = 'mm')

Fig_2 = ggpubr::ggarrange(plotlist  = list(
    list_figures$KRAS_enrichment+labs(title = NULL, subtitle = NULL)+
        xlab('Number of Psms per Protein (log10)')+ ylab('Sample replicates per cell line')+
        theme(legend.position="none",
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank()),
    list_figures$KRASi_1+ ylab('-log10(BFDR)')+ 
        labs(title = NULL, subtitle = NULL)+
        theme(legend.position="none"),
    ggplot(),
    list_figures$BRD4_enrichment+labs(title = NULL, subtitle = NULL)+
        xlab('Number of Psms per Protein (log10)')+ ylab('Sample replicates per cell line')+
        theme(legend.position="none",
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank()),
    list_figures$BRD4s_5+ylab('-log10(BFDR)')+
        labs(title = NULL, subtitle = NULL)+
        theme(legend.position="none"), ggplot(), ggplot(), ggplot(), ggplot()), nrow = 4, ncol = 3)
ggsave(here::here('Output','Fig_2.pdf'), Fig_2, width = 210, height = 297, units = 'mm')
Fig_S2 = ggpubr::ggarrange(plotlist  = 
                               list(
                                   list_figures$KRAS_cov+
                                        labs(title = NULL, subtitle = NULL)+
                                        xlab('Protein Percentage coverage'),    
                                    # list_figures$KRAS_heatmap[[4]],
                                   ggplot(),
                                   ggplot(),
                                    list_figures$KRAS_mod+
                                        theme(legend.position="none",
                                              axis.text.y=element_blank(),
                                              axis.ticks.y=element_blank())+
                                        labs(title = NULL, subtitle = NULL)+
                                        xlab('Peptide relative position'),
                                    list_figures$KRASi_6+ylab('-log10(BFDR)')+
                                        labs(title = NULL, subtitle = NULL)+
                                        theme(legend.position="none"),
                                    # list_figures$KRASi_1_position+
                                    #     theme(legend.position="bottom")+
                                    #     labs(title = NULL, subtitle = NULL)+
                                    #     xlab('Peptide relative position'),
                                    list_figures$BRD4_coverage+ labs(title = NULL, subtitle = NULL),
                                    # list_figures$BRD4_heatmap[[4]],
                                   ggplot(),
                                    list_figures$BRD4s_5_pos+
                                        theme(legend.position="none")+
                                        labs(title = NULL, subtitle = NULL)+
                                        xlab('Peptide relative position'),
                                    list_figures$BRD4s_3+ylab('-log10(BFDR)')+
                                        labs(title = NULL, subtitle = NULL)+
                                        theme(legend.position="none")),ncol =3, nrow = 4, align = 'none')
ggsave(here::here('Output','Fig_S2.pdf'), Fig_S2, width = 210, height = 297, units = 'mm')
                           