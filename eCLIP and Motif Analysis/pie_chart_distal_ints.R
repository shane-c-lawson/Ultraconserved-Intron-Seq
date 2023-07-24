slices <- c((2364-473-89-160-227-2-3), 473, (89+227+2), (160+3))
lbls <- c("Other", "Start/Stop Exon", "Cons_NMD_In", "Cons_NMD_Ex")
pct <- round(slices/sum(slices)*100, 1)
#lbls <- paste(lbls, pct) # add percents to labels
#lbls <- paste(lbls,"%",sep="") # ad % to labels
cols = c('lightblue', 'orange', "indianred", "seagreen")
pie(slices,labels = lbls, col=cols,#rainbow(length(lbls)),
    main="Exons with distant intron conservation", cex.main = 1)#, font.family='Arial')+